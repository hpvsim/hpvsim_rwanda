"""
This file is used to run calibrations for TxV in Rwanda.
Instructions: Go to the CONFIGURATIONS section on lines 29-36 to set up the script before running it.
"""

import os

# # Set the working directory
# os.chdir('D:\TxV_modeling\hpvsim_rwanda-main')
# # To check the working directory in Python
# print(os.getcwd())


# Additions to handle numpy multithreading

os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

# Standard imports
import numpy as np  # Add this at the top if not already imported
import sciris as sc
import hpvsim as hpv

# Imports from this repository
import run_sim as rs
import utils as ut

# CONFIGURATIONS TO BE SET BY USERS BEFORE RUNNING
to_run = [
    'run_calibration',  # Make sure this is uncommented if you want to _run_ the calibrations (usually on VMs)
    # 'plot_calibration',  # Make sure this is uncommented if you want to _plot_ the calibrations (usually locally)
]
debug = False  # If True, this will do smaller runs that can be run locally for debugging
do_save = True

# Run settings for calibration (dependent on debug)
n_trials = [1000, 10][debug]  # How many trials to run for calibration
n_workers = 100   # How many cores to use
n_to_save = 100  # How many results to save in the calibration
# storage = ["mysql://hpvsim_user@localhost/hpvsim_db", None][debug]  # Storage for calibrations
storage = None  # Use local storage for now

########################################################################
# Run calibration
########################################################################
def run_calib(n_trials=None, n_workers=None, do_plot=False, do_save=True, filestem=''):

    sim = rs.make_sim(calib=True, use_calib=False)

    dataloc = 'data/rwanda'  # Location of data files
    datafiles = [
        f'{dataloc}_cancer_cases.csv',  # Globocan
        f'{dataloc}_cancer_incidence_by_age_no_hiv.csv',
        f'{dataloc}_cancer_incidence_by_age_with_hiv.csv',
        f'{dataloc}_asr_cancer_incidence.csv',
        f'{dataloc}_cin_types.csv',
        f'{dataloc}_cancer_types.csv',
    ]
    
    # Define the calibration parameters
    calib_pars = dict(
        beta=[0.05, 0.02, 0.5, 0.02],
    )
    sexual_behavior_pars = dict(
        m_cross_layer=[0.3, 0.1, 0.7, 0.05],
        m_partners=dict(
            c=dict(par1=[0.5, 0.1, 0.6, 0.05])
        ),
        f_cross_layer=[0.4, 0.05, 0.7, 0.05],
        f_partners=dict(
            c=dict(par1=[0.2, 0.1, 0.6, 0.05])
        )
    )
    calib_pars = sc.mergedicts(calib_pars, sexual_behavior_pars)

    genotype_pars = dict(
        hi5=dict(
            cin_fn=dict(k=[.15, .1, .25, 0.01]),
        ),
        ohr=dict(
            cin_fn=dict(k=[.15, .1, .25, 0.01]),
        ),
    )

    hiv_pars = dict(
        rel_sus=dict(
            lt200=[2.25, 2, 5, 0.25],
            gt200=[2.25, 2, 4, 0.25]
        ),
        rel_sev=dict(
            lt200=[1.5, 2, 5, 0.25],
            gt200=[1.5, 1.25, 3, 0.25]
        ),
    )

    # Save some extra sim results
    extra_sim_result_keys = ['cancers', 'cancers_with_hiv', 'cancers_no_hiv',
                             'cancers_by_age_with_hiv', 'cancers_by_age_no_hiv',
                             'asr_cancer_incidence', 'cancer_incidence_by_age_with_hiv',
                             'cancer_incidence_by_age_no_hiv']

    calib = hpv.Calibration(sim, 
                            calib_pars=calib_pars, 
                            genotype_pars=genotype_pars,
                            hiv_pars=hiv_pars,
                            name=f'rwanda_calib',
                            datafiles=datafiles,
                            extra_sim_result_keys=extra_sim_result_keys,
                            total_trials=n_trials, n_workers=n_workers,
                            storage=storage
                            )
    calib.calibrate()
    filename = f'rwanda_calib{filestem}'
    if do_plot:
        calib.plot(do_save=True, fig_path=f'figures/{filename}.png')
    if do_save:
        calib_pars = calib.trial_pars_to_sim_pars(which_pars=0)
        sc.save(f'results/rwanda_pars{filestem}.obj', calib_pars)
        n_to_save = min(len(calib.sim_results), n_to_save)
        cal = ut.shrink_calib(calib, n_results=n_to_save)
        sc.saveobj(f'results/{filename}.obj', cal)

    print(f'Best pars are {calib.best_pars}')

    return sim, calib


########################################################################
# Load pre-run calibration
########################################################################
def load_calib(do_plot=True, filestem=''):
    filename = f'rwanda_calib{filestem}'
    calib = sc.load(f'results/rwanda_calib.obj')
    if do_plot:
        fig = calib.plot(res_to_plot=n_to_save, plot_type='sns.boxplot', do_save=False)
        fig.suptitle(f'Calibration results')
        fig.tight_layout()
        fig.savefig(f'figures/{filename}.png')

    return calib


# %% Run as a script
if __name__ == '__main__':

    T = sc.timer()

    # Run calibration - usually on VMs
    if 'run_calibration' in to_run:
        sim, calib = run_calib(n_trials=n_trials, n_workers=n_workers, do_save=do_save)

    # Load the calibration, plot it, and save the best parameters -- usually locally
    if 'plot_calibration' in to_run:
        calib = load_calib(do_plot=True)

    T.toc('Done')
