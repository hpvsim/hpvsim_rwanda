"""
Interventions for the Rwanda HPVsim model
"""

import os
import hpvsim as hpv
import numpy as np
import pandas as pd


def make_vx(end_year=2100):
    vx_years = np.arange(2011, end_year + 1)
    scaleup = [.2, .4, .6, .8, .9]
    final_cov = 0.9
    vx_cov = np.concatenate([scaleup+[final_cov]*(len(vx_years)-len(scaleup))])
    routine_vx = hpv.routine_vx(product='bivalent', age_range=[11, 12], prob=vx_cov, years=vx_years)
    return routine_vx


def make_st(primary='hpv', prev_screen_cov=0.1, future_screen_cov=0.4, screen_change_year=2026,
            start_year=2020, end_year=2100, prev_treat_cov=0.3, txv_pars=None, future_treat_cov=0.7, treat_change_year=2030):
    """
    Make screening and treatment interventions
    """
    age_range = [30, 50]

    # Determine screening years
    screen_years = np.arange(start_year, end_year + 1)
    final_prev_year = min(screen_change_year, end_year)
    prev_years = np.arange(start_year, final_prev_year + 1)
    future_years = np.arange(screen_change_year + 1, end_year + 1)
    n_prev_years = len(prev_years)
    n_future_years = len(future_years)

    # Define screening coverage
    screen_cov = np.array([prev_screen_cov]*n_prev_years+[future_screen_cov]*n_future_years)
    len_age_range = (age_range[1]-age_range[0])/2
    model_annual_screen_prob = 1 - (1 - screen_cov)**(1/len_age_range)  # Adjusted for age range

    # Routine screening
    screen_eligible = lambda sim: np.isnan(sim.people.date_screened) | \
                                  (sim.t > (sim.people.date_screened + 5 / sim['dt']))
    screening = hpv.routine_screening(
        prob=model_annual_screen_prob,
        eligibility=screen_eligible,
        years=screen_years,
        product=primary,
        age_range=age_range,
        label='screening'
    )
   
    # Assign treatment - historical and status quo
    screen_positive = lambda sim: sim.get_intervention('screening').outcomes['positive']
    assign_treatment = hpv.routine_triage(
        start_year=start_year,
        end_year=final_prev_year,
        prob=1,
        annual_prob=False,
        product='tx_assigner',
        eligibility=screen_positive,
        label='tx assigner'
    )

    # Treatment options - previous
    # Ablation treatment
    ablation_eligible = lambda sim: sim.get_intervention('tx assigner').outcomes['ablation']
    ablation = hpv.treat_num(
        prob=prev_treat_cov,
        product='ablation',
        eligibility=ablation_eligible,
        label='ablation'
    )
    # Excision treatment
    excision_eligible = lambda sim: list(set(sim.get_intervention('tx assigner').outcomes['excision'].tolist()
                                            + sim.get_intervention('ablation').outcomes['unsuccessful'].tolist()))
    excision = hpv.treat_num(
        prob=prev_treat_cov,
        product='excision',
        eligibility=excision_eligible,
        label='excision'
    )
    # Radiation treatment
    radiation_eligible = lambda sim: sim.get_intervention('tx assigner').outcomes['radiation']
    radiation = hpv.treat_num(
        prob=prev_treat_cov/4,  # assume an additional dropoff in CaTx coverage
        product=hpv.radiation(),
        eligibility=radiation_eligible,
        label='radiation'
    )

    st_intvs = [
        screening,
        assign_treatment, ablation, excision, radiation,
    ]

    # Assign treatment - future
    if n_future_years > 0:
        txv_assigner = hpv.default_dx('txvx_assigner')
        txv_assigner.df = pd.read_csv('txv_assigner.csv')
        txv_assigner.hierarchy = ['radiation', 'txv', 'none']
        assign_treatment2 = hpv.routine_triage(
            start_year=screen_change_year,
            prob=1.0,
            annual_prob=False,
            product=txv_assigner,
            eligibility=screen_positive,
            label='txv_assigner'
        )

        # TxV treatment
        txv_eligible = lambda sim: sim.get_intervention('txv_assigner').outcomes['txv']
        txv_prod = hpv.default_tx('txvx2')
        if txv_pars is not None:
            if isinstance(txv_pars, str):
                txv_prod.df = pd.read_csv(txv_pars)
            elif isinstance(txv_pars, pd.DataFrame):
                txv_prod.df = txv_pars
        txv = hpv.linked_txvx(
            prob=future_treat_cov,
            product=txv_prod,
            eligibility=txv_eligible,
            label='txv'
        )

        # Radiation
        religible = lambda sim: sim.get_intervention('txv_assigner').outcomes['radiation']
        radiation2 = hpv.treat_num(
            prob=prev_treat_cov/4,  # assume an additional dropoff in CaTx coverage
            product=hpv.radiation(),
            eligibility=religible,
            label='radiation'
        )

        st_intvs += [assign_treatment2, txv, radiation2]

    return st_intvs

