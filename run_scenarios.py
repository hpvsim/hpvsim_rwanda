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
from interventions import make_st, make_st_20_25, make_st_hiv

# Settings - used here and imported elsewhere
debug = 0
n_seeds = [10, 1][debug]  # How many seeds to run per cluster


# %% Create interventions


def make_st_scenarios():
    scendict = dict()

    # Baseline
    scendict['No interventions'] = []
    scendict['Baseline'] = make_st(screen_change_year=2100)

    for future_screen_cov in [0.1, 0.5, 0.9]:
        for indication, effpars in {'cin': '50/90', 'precin': '90/50'}.items():
            scendict[f'TxV {effpars} with {int(future_screen_cov*100)}% screening'] = make_st(
                future_screen_cov=future_screen_cov,
                txv_pars=f'txvx_pars_{indication}.csv',
                txv=True,
                screen_change_year=2027,)

    return scendict


def make_vx_scenarios():
    """
    Compare vaccination strategies:
        1. Test & vaccinate WLHIV
        2. Test & mass vaccinate women aged 20-25
    """
    scendict = dict()

    start_year = 2027
    mass_vx_age_range = [20, 25]
    cov_array = [.1, .5, .9]
    for cov_val in cov_array:

        # Screen, treat, & vaccinate 20-25yos
        mass_intvs = make_st_20_25(screen_cov=cov_val, age_range=mass_vx_age_range, start_year=start_year)
        scendict[f'Mass vx {cov_val*100:.0f}%'] = mass_intvs
        hiv_intvs = make_st_hiv(screen_cov=cov_val, start_year=start_year)
        scendict[f'HIV+ vx {cov_val*100:.0f}%'] = hiv_intvs
        st_intvs = make_st(screen_change_year=start_year, future_screen_cov=cov_val, st=True)
        scendict[f'S&T {cov_val*100:.0f}%'] = st_intvs

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
        scenarios = sc.mergedicts(make_st_scenarios(), make_vx_scenarios())
        msim = run_sims(scenarios=scenarios, end=end)

        if do_process:

            metrics = ['year', 'asr_cancer_incidence', 'cancers', 'cancer_deaths']

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
