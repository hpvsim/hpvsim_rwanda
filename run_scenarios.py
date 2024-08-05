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
from interventions import make_st

# Settings - used here and imported elsewhere
debug = 0
n_seeds = [20, 1][debug]  # How many seeds to run per cluster


# %% Create interventions


def make_st_scenarios():
    scendict = dict()

    # Baseline
    scendict['Baseline'] = make_st()

    for future_screen_cov in [0.1, 0.5, 0.9]:
        for indication, effpars in {'cin': '50/90', 'precin': '90/50'}.items():
            scendict[f'TxV {effpars} with {int(future_screen_cov*100)}% screening'] = make_st(
                future_screen_cov=future_screen_cov,
                tvx_pars=f'txvx_pars_{indication}.csv',
                treat_change_year=2030
            )

    return scendict


def make_sims(st_scenarios=None, end=2100):
    """ Set up scenarios """

    all_msims = sc.autolist()
    for name, st_intv in st_scenarios.items():
        sims = sc.autolist()
        for seed in range(n_seeds):
            interventions = st_intv
            sim = rs.make_sim(debug=debug, add_st=False, interventions=interventions, end=end, seed=seed)
            sim.label = name
            sims += sim
        all_msims += hpv.MultiSim(sims)

    msim = hpv.MultiSim.merge(all_msims, base=False)

    return msim


def run_sims(st_scenarios=None, end=2100, verbose=0.2):
    """ Run the simulations """
    msim = make_sims(st_scenarios=st_scenarios, end=end)
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
        st_scenarios = make_st_scenarios()
        msim = run_sims(st_scenarios=st_scenarios, end=end)

        if do_process:

            metrics = ['year', 'asr_cancer_incidence', 'cancers', 'cancer_deaths']

            # Process results
            scen_labels = list(st_scenarios.keys())
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
