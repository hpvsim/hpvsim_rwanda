"""
Run scenarios
"""


# %% General settings

import os

os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

# Standard imports
import numpy as np
import sciris as sc
import hpvsim as hpv

# Imports from this repository
import run_sim as rs
from interventions import make_st, make_st_older, make_st_hiv, make_male_vx, make_mv_intvs

# Settings - used here and imported elsewhere
debug = 0
n_seeds = [10, 1][debug]  # How many seeds to run per cluster


# %% Create interventions


def make_txv_scenarios():
    scendict = dict()

    # Baseline
    scendict['No interventions'] = []
    scendict['Baseline'] = make_st(screen_change_year=2100)

    for cov in [0.18, 0.35, 0.7]:
        # Virus-clearing TxV within S&T. Triaged such that women with precin are given TxV
        # and women with CIN are given ablation
        scendict[f'TxV 90/50 within screening, {int(cov*100)}%'] = make_st(
            future_screen_cov=cov,
            txv_pars='precin',
            txv=True,
            screen_change_year=2027,)

        # Lesion-regressing TxV within S&T. Triaged such that women with both precin and CIN are given TxV
        scendict[f'TxV 50/90 within screening, {int(cov*100)}%'] = make_st(
            future_screen_cov=cov,
            txv_pars='cin',
            txv=True,
            screen_change_year=2027,)

        # Virus-clearing mass TxV
        scendict[f'Mass TxV 90/50, {int(cov*100)}%'] = make_mv_intvs(
            txv_pars='precin',
            campaign_coverage=cov,
            routine_coverage=cov,
        )

        # Virus-clearing mass TxV
        scendict[f'Mass TxV 50/90, {int(cov*100)}%'] = make_mv_intvs(
            txv_pars='cin',
            campaign_coverage=cov,
            routine_coverage=cov,
        )

    return scendict


def make_vx_scenarios():
    """
    Compare vaccination strategies:
        1. Test & mass vaccinate women aged 20-25
        2. Test & vaccinate WLHIV
    """
    scendict = dict()

    start_year = 2027
    mass_vx_age_range = [20, 50]
    cov_array = [.18, .35, .70]
    for cov_val in cov_array:

        # Screen, treat, & vaccinate older women
        mass_intvs = make_st_older(screen_cov=cov_val, age_range=mass_vx_age_range, start_year=start_year)
        scendict[f'HPV-Faster {cov_val*100:.0f}%'] = mass_intvs

        # Screen, treat, & vaccinate WLHIV
        hiv_intvs = make_st_hiv(screen_cov=cov_val, start_year=start_year)
        scendict[f'HIV+ vx {cov_val*100:.0f}%'] = hiv_intvs

        # Excluded strategies: scale up S&T
        # st_intvs = make_st(screen_change_year=start_year, future_screen_cov=cov_val)
        # scendict[f'S&T {cov_val*100:.0f}%'] = st_intvs
        # Excluded strategies: vaccinate males
        # intvs = make_male_vx(prob=cov_val)
        # scendict[f'Male vx {cov_val*100:.0f}%'] = intvs

    return scendict


def make_sims(scenarios=None, end=2100):
    """ Set up scenarios """

    all_msims = sc.autolist()
    for name, interventions in scenarios.items():
        sims = sc.autolist()
        for seed in range(n_seeds):
            if name == 'No interventions':
                add_vax = False
            else:
                add_vax = True
            sim = rs.make_sim(debug=debug, add_st=False, add_vax=add_vax, interventions=interventions, end=end, seed=seed)
            sim.label = name
            sims += sim
        all_msims += hpv.MultiSim(sims)

    msim = hpv.MultiSim.merge(all_msims, base=False)

    return msim


def run_sims(scenarios=None, end=2100, verbose=-1):
    """ Run the simulations """
    msim = make_sims(scenarios=scenarios, end=end)
    # for sim in msim.sims:
    #     sim.run()
    msim.run(verbose=verbose)
    return msim


# %% Run as a script
if __name__ == '__main__':

    T = sc.timer()
    do_run = True
    do_save = False
    do_process = True
    end = 2100

    # Run scenarios (usually on VMs, runs n_seeds in parallel over M scenarios)
    if do_run:
        scenarios = sc.mergedicts(make_txv_scenarios(), make_vx_scenarios())
        # scenarios = {k: v for k, v in scenarios.items() if k in ['S&T 18%', 'Mass vx 18%']}
        msim = run_sims(scenarios=scenarios, end=end)

        if do_process:

            metrics = ['year', 'asr_cancer_incidence', 'cancers', 'cancer_deaths', 'cancers_with_hiv', 'cancers_no_hiv', 'cancer_incidence_no_hiv', 'cancer_incidence_with_hiv']

            # Process results
            scen_labels = list(scenarios.keys())
            mlist = msim.split(chunks=len(scen_labels))

            msim_dict = sc.objdict()
            for si, scen_label in enumerate(scen_labels):
                reduced_sim = mlist[si].reduce(output=True)
                mres = sc.objdict({metric: reduced_sim.results[metric] for metric in metrics})

                for ii, intv in enumerate(reduced_sim['interventions']):
                    intv_label = intv.label
                    mres[intv_label] = reduced_sim['interventions'][ii].n_products_used
                msim_dict[scen_label] = mres

            sc.saveobj('results/st_scens.obj', msim_dict)

    print('Done.')
