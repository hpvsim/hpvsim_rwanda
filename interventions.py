"""
Interventions for the Rwanda HPVsim model (v3 / Starsim API).

MIGRATION NOTES (v2 -> v3)
--------------------------
* Products: `hpv.default_dx('hpv')` -> `hpv.dx(name='hpv')`; a custom-CSV product
  -> `hpv.dx(df=pd.read_csv(...), hierarchy=[...])` (hierarchy order matters:
  most-severe result first, because hpv.dx uses hierarchy-min semantics).
* The v2 "therapeutic vaccine" (`default_tx('txvx1')` + `txvx_pars_*.csv`) is a
  per-state EFFICACY product -- i.e. an `hpv.tx` state-flip treatment, NOT the
  immunity-based `hpv.txvx`. We therefore model TxV faithfully as an `hpv.tx`
  product delivered by a treatment intervention (lesion regression by efficacy).
  (v2's product also conferred a small imm_init immunity; that secondary effect
  is not reproduced -- flagged in the migration report.)
* Intervention identity: `label=` -> `name=`; `sim.get_intervention(x)` ->
  `sim.interventions[x]`; a treatment's `name` must differ from its product name.
* Re-screen eligibility: `sim.people.date_screened` is gone -> use the screening
  module's `screened` / `ti_screened` state.
* `sim.yearvec[sim.t]` -> `sim.t.now('year')`.
* Coverage/eligibility outcomes: `screening.outcomes['positive']`,
  `triage.outcomes['ablation']`, `treatment.outcomes['unsuccessful']` all exist
  on the v3 modules (populated during the same timestep, so ordering in the
  intervention list matters: screen -> triage -> treat).
"""
import numpy as np
import pandas as pd
import starsim as ss
import hpvsim as hpv


# Hierarchies for the custom diagnostic CSVs (most-severe first).
_TX_ASSIGN_HIER = ['radiation', 'excision', 'ablation', 'none']
_HPV_DX_HIER = ['positive', 'inadequate', 'negative']


def _named(product, name):
    """Give a custom df-based product a module name.

    hpv.dx/hpv.tx built from a `df=` leave `name=None`, but Starsim registers
    every product as a People module keyed by `name`, so a name is required
    (and must be unique within a sim / distinct from the owning intervention).
    """
    product.name = name
    product.label = name
    return product


def make_vx(end_year=2100):
    """Routine-style bivalent vaccination of girls 11-12, scaling up from 2011.

    v2 used `campaign_vx(sex=0)` (girls only); v3 `campaign_vx` defaults to both
    sexes, so we pass `sex='f'` to match.
    """
    vx_years = np.arange(2011, end_year + 1)
    scaleup = [.2, .4, .6, .8, .9]
    final_cov = 0.9
    vx_cov = np.concatenate([scaleup + [final_cov] * (len(vx_years) - len(scaleup))])
    routine_vx = hpv.campaign_vx(
        name='routine_vx', product=_named(hpv.vx(name='bivalent'), 'routine_vx_prod'), sex='f',
        age_range=[11, 12], prob=vx_cov, interpolate=False, years=vx_years,
    )
    return routine_vx


def make_hpv_test():
    """Custom HPV DNA diagnostic from hpvdna.csv."""
    return _named(hpv.dx(df=pd.read_csv('hpvdna.csv'), hierarchy=_HPV_DX_HIER), 'hpvdna')


def make_st(primary='hpv', prev_screen_cov=0.1, future_screen_cov=0.18,
            screen_change_year=2025, age_range=[30, 50], start_year=2020,
            end_year=2100, future_treat_cov=0.75, txv_pars=None, txv=False,
            tx_assigner_csv='tx_assigner'):
    """Screening + triage + treatment (+ optional TxV) cascade."""
    # Determine screening years
    screen_years = np.arange(start_year, end_year + 1)
    final_prev_year = min(screen_change_year, end_year)
    prev_years = np.arange(start_year, final_prev_year + 1)
    future_years = np.arange(screen_change_year + 1, end_year + 1)
    n_prev_years = len(prev_years)
    n_future_years = len(future_years)

    # Define screening coverage (adjusted for age range width)
    screen_cov = np.array([prev_screen_cov] * n_prev_years + [future_screen_cov] * n_future_years)
    len_age_range = (age_range[1] - age_range[0]) / 2
    model_annual_screen_prob = 1 - (1 - screen_cov) ** (1 / len_age_range)

    # Routine screening -- re-screen if never screened or >10y since last screen.
    def screen_eligible(sim):
        scr = sim.interventions['screening']
        due = scr.screened & ((sim.ti - scr.ti_screened) > (10.0 / sim.t.dt_year))
        return ~scr.screened | due

    screening = hpv.routine_screening(
        name='screening', product=primary, prob=model_annual_screen_prob,
        eligibility=screen_eligible, years=screen_years, age_range=age_range,
        sex='f',
    )
    st_intvs = [screening]

    # Triage: assign treatment from a screen-positive.
    triage_end_year = 2030 if txv_pars == 'cin' else end_year
    triage_years = np.arange(start_year, triage_end_year + 1)
    triage_prob = np.array([0.9] * len(triage_years))

    tx_assigner = _named(hpv.dx(df=pd.read_csv(f'{tx_assigner_csv}.csv'),
                                hierarchy=_TX_ASSIGN_HIER), 'tx_assigner')
    assign_treatment = hpv.routine_triage(
        name='triage', product=tx_assigner, years=triage_years, prob=triage_prob,
        annual_prob=False, sex='f',
        eligibility=lambda sim: sim.interventions['screening'].outcomes['positive'],
    )

    # Ablation treatment. NB: intervention `name` must differ from the product
    # module name ('ablation'), so treatments use an `_rx` suffix.
    ablation = hpv.treat_num(
        name='ablation_rx', prob=future_treat_cov, product='ablation',
        eligibility=lambda sim: sim.interventions['triage'].outcomes['ablation'],
    )
    # Excision treatment (triage-excision OR failed-ablation)
    def excision_eligible(sim):
        return ss.uids(np.union1d(
            np.asarray(sim.interventions['triage'].outcomes['excision'], dtype=int),
            np.asarray(sim.interventions['ablation_rx'].outcomes['unsuccessful'], dtype=int),
        ))
    excision = hpv.treat_num(
        name='excision_rx', prob=future_treat_cov, product='excision',
        eligibility=excision_eligible,
    )
    # Radiation treatment (cancer)
    radiation = hpv.treat_num(
        name='radiation_rx', prob=1 / 4, product=hpv.radiation(),
        eligibility=lambda sim: sim.interventions['triage'].outcomes['radiation'],
    )

    st_intvs += [assign_treatment, ablation, excision, radiation]

    # Optional therapeutic vaccine, delivered as an hpv.tx state-flip treatment.
    if txv:
        txv_prod = _named(hpv.tx(df=pd.read_csv(f'txvx_pars_{txv_pars}.csv')), 'txvprod')

        def txv_eligible(sim):
            if sim.t.now('year') >= 2030:
                return sim.interventions['screening'].outcomes['positive']
            return np.array([], dtype=int)

        txv_intv = hpv.treat_num(
            name='txv', prob=0.9, product=txv_prod, eligibility=txv_eligible,
        )
        st_intvs += [txv_intv]

    return st_intvs


def make_mv_intvs(campaign_coverage=None, txv_pars=None, intro_year=2030,
                  campaign_age=[20, 50]):
    """One-time mass TxV campaign (delivered as an hpv.tx state-flip treatment)."""
    campaign_years = [intro_year]
    txv_prod = _named(hpv.tx(df=pd.read_csv(f'txvx_pars_{txv_pars}.csv')), 'txvprod')

    campaign_txvx = hpv.campaign_txvx(
        name='campaign txvx', prob=campaign_coverage, interpolate=False,
        years=campaign_years, age_range=campaign_age, product=txv_prod,
        eligibility=lambda sim: (sim.interventions['campaign txvx'].txvx_doses == 0),
    )

    # Add historical screening and treatment
    hist_intvs = make_st(screen_change_year=2026)
    return hist_intvs + [campaign_txvx]


def make_st_older(start_year=2027, screen_cov=0.4, treat_cov=1, age_range=[20, 50]):
    """HPV-Faster: one-time screen + treat + vaccinate campaign for older women."""
    primary = make_hpv_test()
    screening = hpv.campaign_screening(
        name='screening_older', prob=screen_cov, interpolate=False,
        years=[start_year], product=primary, age_range=age_range, sex='f',
    )

    tx_assigner = _named(hpv.dx(df=pd.read_csv('tx_assigner_no_triage.csv'),
                                hierarchy=_TX_ASSIGN_HIER), 'tx_assigner_older')
    assign_treatment = hpv.campaign_triage(
        name='triage_older', years=[start_year], prob=1, interpolate=False,
        product=tx_assigner, sex='f',
        eligibility=lambda sim: sim.interventions['screening_older'].outcomes['positive'],
    )

    # Explicit, distinctly-named products: this cascade is combined with the
    # normal make_st() cascade (below), which already registers 'ablation' /
    # 'excision' product modules; product module names must be unique per sim.
    ablation = hpv.treat_num(
        name='ablation_older', prob=treat_cov,
        product=_named(hpv.tx(name='ablation'), 'ablation_older_prod'),
        eligibility=lambda sim: sim.interventions['triage_older'].outcomes['ablation'],
    )
    def excision_older_eligible(sim):
        return ss.uids(np.union1d(
            np.asarray(sim.interventions['triage_older'].outcomes['excision'], dtype=int),
            np.asarray(sim.interventions['ablation_older'].outcomes['unsuccessful'], dtype=int),
        ))
    excision = hpv.treat_num(
        name='excision_older', prob=treat_cov,
        product=_named(hpv.tx(name='excision'), 'excision_older_prod'),
        eligibility=excision_older_eligible,
    )
    radiation = hpv.treat_num(
        name='radiation_older', prob=treat_cov, product=hpv.radiation(name='radiation_older_prod'),
        eligibility=lambda sim: sim.interventions['triage_older'].outcomes['radiation'],
    )

    # Vaccinate women screened this campaign (nonavalent).
    def mass_eligible(sim):
        scr = sim.interventions['screening_older']
        return scr.screened & (scr.ti_screened == sim.ti)
    mass_vx = hpv.campaign_vx(
        name='mass_vax', product=_named(hpv.vx(name='nonavalent'), 'mass_vax_prod'), sex='f',
        eligibility=mass_eligible, interpolate=False, age_range=age_range,
        prob=screen_cov, years=[start_year],
    )

    normal_intvs = make_st(screen_change_year=2026)
    return normal_intvs + [screening, assign_treatment, ablation, excision,
                           radiation, mass_vx]
