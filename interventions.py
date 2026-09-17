"""
Interventions for the Rwanda HPVsim model (v3.2)
"""

import numpy as np
import pandas as pd
import starsim as ss
import hpvsim as hpv


# Skip anyone re-screened within this many years.
SCREEN_GAP_YEARS = 10


def make_vx(end_year=2100):
    """Routine bivalent vaccination at ages 11-12: 20-90% scale-up 2011-15, 90% thereafter."""
    vx_years = np.arange(2011, end_year + 1)
    scaleup = [0.2, 0.4, 0.6, 0.8, 0.9]
    final_cov = 0.9
    vx_cov = np.concatenate([scaleup, [final_cov] * (len(vx_years) - len(scaleup))])
    return hpv.campaign_vx(
        name='routine_vx',
        product='bivalent',
        age_range=[11, 12],
        prob=vx_cov,
        years=vx_years,
    )


def make_hpv_test(name='hpv_test'):
    """HPV DNA diagnostic with Rwanda-specific per-genotype probabilities."""
    return hpv.dx(
        name=name,
        df=pd.read_csv('hpvdna.csv'),
        hierarchy=['positive', 'inadequate', 'negative'],
    )


def make_via_test(name='via_test'):
    """VIA (visual inspection with acetic acid) diagnostic.

    Sensitivity by stage: precin 15%, cin 55%, cancerous 80%,
    specificity ~95%. Used for the R2.6 primary-test sensitivity sweep.
    """
    return hpv.dx(
        name=name,
        df=pd.read_csv('via.csv'),
        hierarchy=['positive', 'inadequate', 'negative'],
    )


def make_st(future_screen_cov=0.18, coverage_change_year=2028, age_range=[30, 50],
            start_year=2020, end_year=2100, future_treat_cov=0.75,
            txv_pars=None, txv=False, tx_assigner_csv='tx_assigner',
            txv_start_year=2030, treat_capacity=None,
            txv_efficacy_mult=1.0):
    """
    Every scenario runs the status-quo screening program (S&T&T at 18%
    coverage, HPV DNA + VIA triage, 75% ablation attendance) from
    start_year..coverage_change_year-1. This is the shared 2020-2027
    "before the switch" period; all scenarios have identical dynamics
    in it.

    From coverage_change_year (default 2028) onwards, the variant-specific
    intervention era runs:

      tx_assigner_csv='tx_assigner' (default, VIA triage):
        chain per lesion = 0.9 (triage attends) x VIA_sens x
        future_treat_cov x 0.936 (ablation efficacy). VIA sens is
        30% for precin, 60% for CIN per tx_assigner.csv.

      tx_assigner_csv='tx_assigner_no_triage' (S&T direct):
        no separate VIA/ablation LTFU split - ablation given at the
        same visit as receiving results. future_treat_cov is force-
        overridden to 1.0 internally so LTFU isn't double-counted.
        Chain per lesion = 0.9 (single-visit attends) x 1.0 x 1.0 x
        0.936 = 84%.

    TxV (linked_txvx at prob=0.9) fires from txv_start_year onwards on
    screen positives. For txv_pars='cin' (lesion-regressing profile) the
    intv-era ablate/excise path stops at txv_start_year so TxV is the
    sole treatment path 2030+. For txv_pars='precin' the two run side-
    by-side.

    treat_capacity caps ablation+excision agents per timestep (None = no
    cap). txv_efficacy_mult scales the loaded txvx_pars CSV in-memory
    (clipped to [0,1]).
    """
    len_age_range = (age_range[1] - age_range[0]) / 2

    # Screening: one intervention across both eras with time-varying prob.
    # Status quo (SQ) years use 18% lifetime coverage; intv-era years use
    # future_screen_cov.
    sq_years = np.arange(start_year, coverage_change_year)          # 2020..2027
    intv_years = np.arange(coverage_change_year, end_year + 1)      # 2028..2100
    screen_years = np.concatenate([sq_years, intv_years])
    screen_cov = np.concatenate([
        np.full(len(sq_years), 0.18),
        np.full(len(intv_years), future_screen_cov),
    ])
    annual_screen_prob = 1 - (1 - screen_cov) ** (1 / len_age_range)

    def screen_eligible(sim):
        # Never screened, or last screened > gap ago.
        ti_s = sim.interventions['screening'].ti_screened
        gap_ti = SCREEN_GAP_YEARS / sim.t.dt_year
        stale = sim.ti > (ti_s + gap_ti)
        return (ti_s.isnan | stale).uids

    primary = make_hpv_test(name='hpv_test_routine')
    screening = hpv.routine_screening(
        name='screening',
        prob=annual_screen_prob,
        eligibility=screen_eligible,
        years=screen_years,
        product=primary,
        age_range=age_range,
    )
    st_intvs = [screening]

    screen_positive = lambda sim: sim.interventions['screening'].outcomes['positive']

    # === Status quo era: S&T&T at 18% with VIA + 75% treat_num ===
    if len(sq_years) > 0:
        tx_assigner_sq = hpv.dx(
            name='tx_assigner_sq_product',
            df=pd.read_csv('tx_assigner.csv'),
            hierarchy=['radiation', 'excision', 'ablation', 'none'],
        )
        assign_sq = hpv.routine_triage(
            name='tx_assigner_sq_intv',
            years=sq_years,
            prob=np.full(len(sq_years), 0.9),
            annual_prob=False,
            product=tx_assigner_sq,
            eligibility=screen_positive,
        )
        # SQ-era treatment products carry unique module_names to avoid
        # colliding with the intv-era treat_num instances that use the
        # shipped default products.
        ablation_sq_eligible = lambda sim: sim.interventions['tx_assigner_sq_intv'].outcomes['ablation']
        ablation_sq = hpv.treat_num(
            name='ablation_sq_intv',
            prob=0.75,
            product=hpv.tx(name='ablation', module_name='ablation_sq_prod'),
            eligibility=ablation_sq_eligible,
        )
        def excision_sq_eligible(sim):
            triage_out = sim.interventions['tx_assigner_sq_intv'].outcomes['excision']
            abl_fail = sim.interventions['ablation_sq_intv'].outcomes['unsuccessful']
            return triage_out | abl_fail
        excision_sq = hpv.treat_num(
            name='excision_sq_intv',
            prob=0.75,
            product=hpv.tx(name='excision', module_name='excision_sq_prod'),
            eligibility=excision_sq_eligible,
        )
        radiation_sq_eligible = lambda sim: sim.interventions['tx_assigner_sq_intv'].outcomes['radiation']
        radiation_sq = hpv.treat_num(
            name='radiation_sq_intv',
            prob=1/4,
            product=hpv.radiation(name='radiation_sq_prod'),
            eligibility=radiation_sq_eligible,
        )
        st_intvs += [assign_sq, ablation_sq, excision_sq, radiation_sq]

    # === Intervention era: variant-specific ===
    # For txv_pars='cin', triage/ablate stops at txv_start_year and TxV
    # becomes the sole treatment path.
    triage_end_year = min(txv_start_year, end_year) if txv_pars == 'cin' else end_year
    intv_active_years = np.arange(coverage_change_year, triage_end_year + 1)

    if len(intv_active_years) > 0:
        # tx_assigner_no_triage: same-visit ablation, no additional LTFU
        # beyond the 0.9 attendance step.
        if tx_assigner_csv == 'tx_assigner_no_triage':
            treat_prob = 1.0
        else:
            treat_prob = future_treat_cov

        tx_assigner = hpv.dx(
            name='tx_assigner_product',
            df=pd.read_csv(f'{tx_assigner_csv}.csv'),
            hierarchy=['radiation', 'excision', 'ablation', 'none'],
        )
        assign_treatment = hpv.routine_triage(
            name='tx_assigner_intv',
            years=intv_active_years,
            prob=np.full(len(intv_active_years), 0.9),
            annual_prob=False,
            product=tx_assigner,
            eligibility=screen_positive,
        )
        ablation_eligible = lambda sim: sim.interventions['tx_assigner_intv'].outcomes['ablation']
        ablation = hpv.treat_num(
            name='ablation_intv',
            prob=treat_prob,
            product='ablation',
            eligibility=ablation_eligible,
            max_capacity=treat_capacity,
        )
        def excision_eligible(sim):
            triage_out = sim.interventions['tx_assigner_intv'].outcomes['excision']
            abl_fail = sim.interventions['ablation_intv'].outcomes['unsuccessful']
            return triage_out | abl_fail
        excision = hpv.treat_num(
            name='excision_intv',
            prob=treat_prob,
            product='excision',
            eligibility=excision_eligible,
            max_capacity=treat_capacity,
        )
        radiation_eligible = lambda sim: sim.interventions['tx_assigner_intv'].outcomes['radiation']
        radiation = hpv.treat_num(
            name='radiation_intv',
            prob=1/4,  # extra dropoff for cancer treatment
            product=hpv.radiation(),
            eligibility=radiation_eligible,
        )
        st_intvs += [assign_treatment, ablation, excision, radiation]

    if txv:
        txv_df = pd.read_csv(f'txvx_pars_{txv_pars}.csv')
        if txv_efficacy_mult != 1.0:
            txv_df = txv_df.copy()
            txv_df['efficacy'] = (txv_df['efficacy'] * txv_efficacy_mult).clip(0, 1)
        txv_prod = hpv.txvx(
            df=txv_df,
            imm_init=ss.uniform(low=0.49, high=0.51),
        )
        def txv_eligible(sim):
            if sim.now.years >= txv_start_year:
                return sim.interventions['screening'].outcomes['positive']
            return ss.uids()
        st_intvs.append(hpv.linked_txvx(
            name='txv',
            prob=0.9,
            product=txv_prod,
            eligibility=txv_eligible,
        ))

    return st_intvs


def make_mv_intvs(campaign_coverage=None, txv_pars=None, intro_year=2030,
                  campaign_age=[20, 50], end_year=2100, st_kwargs=None,
                  txv_efficacy_mult=1.0):
    """Mass therapeutic vaccination campaign, layered on top of baseline S&T.

    st_kwargs is forwarded to the nested make_st() so the normalized-start
    scenario set can strip the pre-2030 screening history.
    """
    txv_df = pd.read_csv(f'txvx_pars_{txv_pars}.csv')
    if txv_efficacy_mult != 1.0:
        txv_df = txv_df.copy()
        txv_df['efficacy'] = (txv_df['efficacy'] * txv_efficacy_mult).clip(0, 1)
    txv_prod = hpv.txvx(
        df=txv_df,
        imm_init=ss.uniform(low=0.49, high=0.51),
    )

    def mv_eligible(sim):
        return ~sim.interventions.campaign_txvx.tx_vaccinated

    campaign_txvx = hpv.campaign_txvx(
        name='campaign_txvx',
        prob=campaign_coverage,
        years=[intro_year],
        age_range=campaign_age,
        product=txv_prod,
        eligibility=mv_eligible,
    )
    hist_intvs = make_st(end_year=end_year, **(st_kwargs or {}))
    return hist_intvs + [campaign_txvx]


def make_st_older(start_year=2028, screen_cov=0.4, treat_cov=1,
                  age_range=[20, 50], end_year=2100):
    """One-off screen-and-vax campaign for 20-50yo, layered on baseline S&T.

    campaign_triage prob=0.9 = 10% LTFU at the single results+treatment
    visit (matches the S&T same-visit attendance model). treat_cov=1 by
    default (no additional LTFU beyond the 0.9 attendance step)."""
    primary = make_hpv_test(name='hpv_test_older')
    screening = hpv.campaign_screening(
        name='screening_older',
        prob=screen_cov,
        years=[start_year],
        product=primary,
        age_range=age_range,
    )

    tx_assigner = hpv.dx(
        name='tx_assigner_older_product',
        df=pd.read_csv('tx_assigner_no_triage.csv'),
        hierarchy=['radiation', 'excision', 'ablation', 'none'],
    )
    screen_positive = lambda sim: sim.interventions['screening_older'].outcomes['positive']
    assign_treatment = hpv.campaign_triage(
        name='tx_assigner_older',
        years=[start_year],
        prob=0.9,  # 10% LTFU at same-visit results + treatment
        product=tx_assigner,
        eligibility=screen_positive,
    )

    # Older-cohort treatment products carry `module_name=` so they don't
    # collide with make_st's identically-shipped ablation/excision/radiation.
    ablation_eligible = lambda sim: sim.interventions['tx_assigner_older'].outcomes['ablation']
    ablation = hpv.treat_num(
        name='ablation_older',
        prob=treat_cov,
        product=hpv.tx(name='ablation', module_name='ablation_older_prod'),
        eligibility=ablation_eligible,
    )

    def excision_eligible(sim):
        triage_out = sim.interventions['tx_assigner_older'].outcomes['excision']
        abl_fail = sim.interventions['ablation_older'].outcomes['unsuccessful']
        return triage_out | abl_fail
    excision = hpv.treat_num(
        name='excision_older',
        prob=treat_cov,
        product=hpv.tx(name='excision', module_name='excision_older_prod'),
        eligibility=excision_eligible,
    )

    radiation_eligible = lambda sim: sim.interventions['tx_assigner_older'].outcomes['radiation']
    radiation = hpv.treat_num(
        name='radiation_older',
        prob=treat_cov/4,
        product=hpv.radiation(name='radiation_older_prod'),
        eligibility=radiation_eligible,
    )

    def mass_eligible(sim):
        # Just-screened by the older-cohort screen AND never vaccinated by any
        # prophylactic intervention.
        scr = sim.interventions['screening_older']
        just_screened_uids = (scr.ti_screened == sim.ti).uids
        for intv in sim.interventions.values():
            v = getattr(intv, 'vaccinated', None)
            if v is not None:
                just_screened_uids = just_screened_uids.remove(v.uids)
        return just_screened_uids

    mass_vx = hpv.campaign_vx(
        name='mass_vax',
        product='nonavalent',
        eligibility=mass_eligible,
        age_range=age_range,
        prob=screen_cov,
        years=[start_year],
    )

    normal_intvs = make_st(end_year=end_year)
    return normal_intvs + [
        screening, assign_treatment, ablation, excision, radiation, mass_vx,
    ]
