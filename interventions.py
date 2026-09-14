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


def make_st(primary=None, prev_screen_cov=0.1, future_screen_cov=0.18,
            screen_change_year=2027, age_range=[30, 50],
            start_year=2020, end_year=2100, future_treat_cov=0.75,
            txv_pars=None, txv=False, tx_assigner_csv='tx_assigner'):
    """
    Make screening and treatment interventions.
    """
    # Per-year prob split across pre / post the coverage change year.
    screen_years = np.arange(start_year, end_year + 1)
    final_prev_year = min(screen_change_year, end_year)
    prev_years = np.arange(start_year, final_prev_year + 1)
    future_years = np.arange(screen_change_year + 1, end_year + 1)
    screen_cov = np.array(
        [prev_screen_cov] * len(prev_years) + [future_screen_cov] * len(future_years)
    )
    # Convert lifetime coverage across the age window to an annual probability.
    len_age_range = (age_range[1] - age_range[0]) / 2
    model_annual_screen_prob = 1 - (1 - screen_cov) ** (1 / len_age_range)

    def screen_eligible(sim):
        # Never screened, or last screened > gap ago.
        ti_s = sim.interventions['screening'].ti_screened
        gap_ti = SCREEN_GAP_YEARS / sim.t.dt_year
        stale = sim.ti > ti_s + gap_ti
        return (ti_s.isnan | stale).uids

    if primary is None:
        primary = make_hpv_test(name='hpv_test_routine')
    screening = hpv.routine_screening(
        name='screening',
        prob=model_annual_screen_prob,
        eligibility=screen_eligible,
        years=screen_years,
        product=primary,
        age_range=age_range,
    )
    st_intvs = [screening]

    # If lesion-regressing therapeutic vaccine is on, stop the ablate/excise path
    # once it kicks in (2030).
    triage_end_year = min(2030, end_year) if txv_pars == 'cin' else end_year
    triage_years = np.arange(start_year, triage_end_year + 1)
    triage_prob = np.full(len(triage_years), 0.9)

    tx_assigner = hpv.dx(
        name='tx_assigner_product',
        df=pd.read_csv(f'{tx_assigner_csv}.csv'),
        hierarchy=['radiation', 'excision', 'ablation', 'none'],
    )
    screen_positive = lambda sim: sim.interventions['screening'].outcomes['positive']
    assign_treatment = hpv.routine_triage(
        name='tx_assigner_intv',
        years=triage_years,
        prob=triage_prob,
        annual_prob=False,
        product=tx_assigner,
        eligibility=screen_positive,
    )

    # Intervention names get an `_intv` suffix because v3 reserves the bare
    # product names ('ablation', 'excision', 'radiation', 'tx_assigner').
    ablation_eligible = lambda sim: sim.interventions['tx_assigner_intv'].outcomes['ablation']
    ablation = hpv.treat_num(
        name='ablation_intv',
        prob=future_treat_cov,
        product='ablation',
        eligibility=ablation_eligible,
    )

    def excision_eligible(sim):
        triage_out = sim.interventions['tx_assigner_intv'].outcomes['excision']
        abl_fail = sim.interventions['ablation_intv'].outcomes['unsuccessful']
        return triage_out | abl_fail  # TODO check this
    excision = hpv.treat_num(
        name='excision_intv',
        prob=future_treat_cov,
        product='excision',
        eligibility=excision_eligible,
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
        txv_prod = hpv.txvx(
            df=pd.read_csv(f'txvx_pars_{txv_pars}.csv'),
            imm_init=ss.uniform(low=0.49, high=0.51),
        )
        def txv_eligible(sim):
            if sim.now.years >= 2030:
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
                  campaign_age=[20, 50], end_year=2100):
    """Mass therapeutic vaccination campaign, layered on top of baseline S&T."""
    txv_prod = hpv.txvx(
        df=pd.read_csv(f'txvx_pars_{txv_pars}.csv'),
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
    hist_intvs = make_st(end_year=end_year)
    return hist_intvs + [campaign_txvx]


def make_st_older(start_year=2027, screen_cov=0.4, treat_cov=1,
                  age_range=[20, 50], end_year=2100):
    """One-off screen-and-vax campaign for 20-50yo, layered on baseline S&T."""
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
        prob=1,  # no LTFU by assumption
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
