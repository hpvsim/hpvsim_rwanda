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
    print("vx_years:", len(vx_years), vx_years)
    print("vx_cov:", len(vx_cov), vx_cov)
    print("Routine Vaccination")
    print(f"Years: {len(vx_years)}, Probabilities: {len(vx_cov)}")
   
    assert len(vx_years) == len(vx_cov)

    routine_vx = hpv.routine_vx(product='bivalent', age_range=[11, 12], prob=vx_cov, years=vx_years)
    print("Routine Vaccination")
    print(f"Years: {len(vx_years)}, Probabilities: {len(vx_cov)}")
    return routine_vx


def make_st(primary='hpv', prev_screen_cov=0.1, future_screen_cov=0.4, screen_change_year=2025,
            start_year=2020, end_year=2100,prev_treat_cov=0.3, txv_pars=None, future_treat_cov=0.7, treat_change_year=2100, dt=0.25):
    """
    Make screening and treatment interventions
    """
    print("end_year:", end_year)
    
    age_range = [30, 50]
    len_age_range = (age_range[1]-age_range[0])/2

    # Determine screening years
    screen_years = np.arange(start_year, end_year + 1)
    n_prev_years = screen_change_year - start_year + 1
    n_future_years = end_year - screen_change_year  
    print("start_year:", start_year)
    print("screen_change_year:", screen_change_year)
    print("end_year:", end_year)
    print("n_prev_years:", n_prev_years)
    print("n_future_years:", n_future_years)
    print("screen_years:", len(screen_years), screen_years)
    
    # Define screening coverage
    screen_cov = np.array([prev_screen_cov]*n_prev_years+[future_screen_cov]*n_future_years)
    # Debugging: Print lengths of screen_years and screen_cov
    print(f"Length of screen_years: {len(screen_years)}")
    print(f"Length of screen_cov: {len(screen_cov)}")

 # Ensure lengths match
    if len(screen_cov) > len(screen_years):  # Truncate screen_cov if it's longer
        screen_cov = screen_cov[:len(screen_years)]
    elif len(screen_years) > len(screen_cov):  # Extend screen_cov if it's shorter
        screen_cov = np.append(screen_cov, [screen_cov[-1]] * (len(screen_years) - len(screen_cov)))
    
    print(f"Adjusted Length of screen_cov: {len(screen_cov)}")

    print("screen_years:", len(screen_years), "screen_cov:", len(screen_cov))
    print("screen_cov:", len(screen_cov), screen_cov)
    assert len(screen_years) == len(screen_cov)

    # Calculate annual screening probability
    model_annual_screen_prob = 1 - (1 - screen_cov)**(1/1) # 1 year of screening
    # Interpolate prob to match the simulation's time step
    #timepoints = np.arange(start_year, end_year + dt, dt)  # Generate timepoints based on dt
    #model_annual_screen_prob = np.interp(timepoints, screen_years, model_annual_screen_prob)
    print(f"Length of model_annual_screen_prob: {len(model_annual_screen_prob)}")

    # Determine treatment coverage
    treat_years = np.arange(start_year, end_year + 1)  # Define treat_years here
    treat_change_year = min(treat_change_year, end_year)
    n_prev_years_treat = treat_change_year - start_year
    n_future_years_treat = end_year - treat_change_year 
    #treat_coverage = np.array([prev_treat_cov]*n_prev_years_treat+[future_treat_cov]*n_future_years_treat)
    treat_coverage = np.array(
    [prev_treat_cov] * n_prev_years_treat + [future_treat_cov] * (n_future_years_treat + 1)
 )

    # Ensure lengths match for treatment coverage
    if len(treat_coverage) > len(treat_years):  # Truncate treat_coverage if it's longer
       treat_coverage = treat_coverage[:len(treat_years)]
    elif len(treat_years) > len(treat_coverage):  # Extend treat_coverage if it's shorter
       treat_coverage = np.append(treat_coverage, [treat_coverage[-1]] * (len(treat_years) - len(treat_coverage)))
    print(f"Adjusted Length of treat_coverage: {len(treat_coverage)}")

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
        end_year=treat_change_year,
        prob=1,
        annual_prob=False,
        product='tx_assigner',
        eligibility=screen_positive,
        label='tx assigner'
    )

    # Assign treatment - future
    screen_positive = lambda sim: sim.get_intervention('screening').outcomes['positive']
    txv_assigner = hpv.default_dx('txvx_assigner')
    txv_assigner.df = pd.read_csv('txv_assigner.csv')
    txv_assigner.hierarchy = ['radiation', 'txv', 'none']
    assign_treatment2 = hpv.routine_triage(
        start_year=treat_change_year,
        prob=1.0,
        annual_prob=False,
        product=txv_assigner,
        eligibility=screen_positive,
        label='txv_assigner'
    )

    # Treatment options

    # Ablation treatment
    print("Ablation Treatment")
    print(f"Years: {len(treat_years)}, Probabilities: {len(treat_coverage)}")
    ablation_eligible = lambda sim: sim.get_intervention('tx assigner').outcomes['ablation']
    ablation = hpv.treat_num(
        prob=treat_coverage,
        product='ablation',
        eligibility=ablation_eligible,
        label='ablation'
    )
    # Excision treatment
    print("Excision Treatment")
    print(f"Years: {len(treat_years)}, Probabilities: {len(treat_coverage)}")
    excision_eligible = lambda sim: list(set(sim.get_intervention('tx assigner').outcomes['excision'].tolist()
                                            + sim.get_intervention('ablation').outcomes['unsuccessful'].tolist()))
    excision = hpv.treat_num(
        prob=treat_coverage,
        product='excision',
        eligibility=excision_eligible,
        label='excision'
    )
    # Radiation treatment
    print("Radiation Treatment")
    print(f"Years: {len(treat_years)}, Probabilities: {len(treat_coverage / 4)}")
    radiation_eligible = lambda sim: list(set(sim.get_intervention('tx assigner').outcomes['radiation'].tolist()
                                            + sim.get_intervention('txv_assigner').outcomes['radiation'].tolist()))
    radiation = hpv.treat_num(
        prob=treat_coverage/4,  # assume an additional dropoff in CaTx coverage
        product=hpv.radiation(),
        eligibility=radiation_eligible,
        label='radiation'
    )
    # TxV treatment
    txv_years = np.arange(start_year, end_year + 1)
    txv_prob = [future_treat_cov] * len(txv_years)  # Ensure lengths match
    print("TxV Treatment")
    print(f"Years: {len(txv_years)}, Probabilities: {len(txv_prob)}")
    txv_eligible = lambda sim: sim.get_intervention('txv_assigner').outcomes['txv']
    txv_prod = hpv.default_tx('txvx2')
    if txv_pars is not None:
        if isinstance(txv_pars, str):
            txv_prod.df = pd.read_csv(txv_pars)
        elif isinstance(txv_pars, pd.DataFrame):
            txv_prod.df = txv_pars
    txv = hpv.linked_txvx(
        prob=[future_treat_cov] * (end_year - start_year + 1),
        product=txv_prod,
        eligibility=txv_eligible,
        label='txv'
    )

    st_intvs = [
        screening,
        assign_treatment, ablation, excision, radiation,
        assign_treatment2, txv
    ]

    return st_intvs

