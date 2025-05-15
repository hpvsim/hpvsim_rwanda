"""
This file is used to run calibrations for TxV in Rwanda.

Instructions: Go to the CONFIGURATIONS section on lines 29-36 to set up the script before running it.
"""
#Set the working directory
import os
os.chdir('D:\TxV_modeling\hpvsim_rwanda-main')

#To check the working directory in Python
print(os.getcwd())


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
import run_sim_hiv as rs

# CONFIGURATIONS TO BE SET BY USERS BEFORE RUNNING
to_run = [
    'run_calibration',  # Make sure this is uncommented if you want to _run_ the calibrations (usually on VMs)
    'plot_calibration',  # Make sure this is uncommented if you want to _plot_ the calibrations (usually locally)
]
debug = True  # If True, this will do smaller runs that can be run locally for debugging
do_save = True

# Run settings for calibration (dependent on debug)
n_trials = [20000, 10][debug]  # How many trials to run for calibration
n_workers = 1   # How many cores to use
storage = "sqlite:///D:/TxV_modeling/rwanda_calib.db"  # Storage for calibrations 

########################################################################
# Run calibration
########################################################################
def make_priors():
                                                      
    return

def run_calib(location=None, n_trials=None, n_workers=None,
              do_plot=False, do_save=True, filestem=''):
    dflocation=location.replace(" ", "_")
    if location == 'rwanda':
        hiv_datafile = ['data/hiv_incidence_rwanda.csv',
                        'data/rwanda_female_hiv_mortality.csv',
                        'data/rwanda_male_hiv_mortality.csv']
        art_datafile = ['data/rwanda_art_coverage_by_age_males.csv',
                        'data/rwanda_art_coverage_by_age_females.csv']

    hiv_pars = {
        'cd4_trajectory': lambda y: np.full_like(np.atleast_1d(y), 0.5, dtype=float)
    }
                                                         ##### HIV/ART Data: Defined before rs.make_sim() because they are simulation inputs, whias Cancer Data are Defined after rs.make_sim() because they are calibration targets (e.g., the model adjusts parameters to match cancer incidence).
                                                    ##### Cancer Data are Calibration targets: Used to compare model outputs to real-world data during calibration (not as direct inputs to the simulation).
    sim = rs.make_sim(location, art_datafile=art_datafile, hiv_datafile=hiv_datafile, calib=True,add_vax=False, analyzers=[],hiv_pars=hiv_pars
    )
    #cd4_traj = sim['hiv_pars'].get('cd4_trajectory', None)
    #if not callable(cd4_traj):
       # val = cd4_traj if isinstance(cd4_traj, (float, int)) else 0.5
        #sim['hiv_pars']['cd4_trajectory'] = lambda y: np.full_like(np.atleast_1d(y), float(val), dtype=float)

    #sim['hiv_pars']['time_to_hiv_death_scale'] = lambda age: np.ones_like(age) * 1.0
    #sim['hiv_pars']['cd4_reconstitution'] = lambda months_on_ART: np.ones_like(months_on_ART) * 1.0  
    #sim.run()

    for intervention in sim['interventions']:
        print(f"Intervention: {type(intervention).__name__}")
    if hasattr(intervention, 'years'):
        print(f"Years: {intervention.years}")
    if hasattr(intervention, 'prob'):
        print(f"Probabilities: {intervention.prob}")
    
    datafiles = [
        f'data/{dflocation}_cancer_cases.csv', #Globocan
        f'data/{dflocation}_cancer_incidence_by_age_no_hiv.csv', 
        f'data/{dflocation}_cancer_incidence_by_age_with_hiv.csv', 
        f'data/{dflocation}_asr_cancer_incidence.csv',
        f'data/{dflocation}_cin_types.csv',
        f'data/{dflocation}_cancer_types.csv',
    ]
    
    # Define the calibration parameters
    calib_pars = dict(
        beta=[0.05, 0.02, 0.5, 0.02],
        own_imm_hr=[0.5, 0.25, 1, 0.05],
        age_risk=dict(risk=[3.2, 1, 4, 0.1],
                      age=[38, 30, 45, 1]),
        # sev_dist=dict(par1=[1, 1, 2, 0.1]),
        cell_imm_init=dict(par1=[0.2, 0.2, 0.8, 0.05]),
        # hpv_control_prob=[0,0,1, 0.25],
        # hpv_reactivation=[0.025, 0, 0.1, 0.025]
    )

    if location is None:
        sexual_behavior_pars = dict(
            m_cross_layer=[0.9, 0.5, 0.95, 0.05],
            m_partners=dict(
                c=dict(par1=[10, 5, 12, 1])
            ),
            f_cross_layer=[0.1, 0.05, 0.5, 0.05],
            f_partners=dict(
                c=dict(par1=[1, .5, 2, .1], par2=[.2, .1, 1, .05])
            )
        )
    else:
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
        hpv16=dict(
            dur_precin=dict(par1=[3, 1, 10, 0.5], par2=[9, 5, 15, 0.5]),
            cancer_fn=dict(transform_prob=[2e-3, 1e-3, 3e-3, 2e-4]),
            cin_fn=dict(k=[.35, .2, .4, 0.01]),
            dur_cin=dict(par1=[5, 4, 6, 0.5], par2=[20, 16, 24, 0.5]),
        ),
        hpv18=dict(
            dur_precin=dict(par1=[2.5, 1, 10, 0.5], par2=[9, 5, 15, 0.5]),
            cancer_fn=dict(transform_prob=[2e-3, 1e-3, 3e-3, 2e-4]),
            cin_fn=dict(k=[.4, .15, .35, 0.01]),
            dur_cin=dict(par1=[5, 4, 6, 0.5], par2=[20, 16, 24, 0.5]),
        ),
        hi5=dict(
            dur_precin=dict(par1=[2.5, 1, 10, 0.5], par2=[9, 5, 15, 0.5]),
            cancer_fn=dict(transform_prob=[1.5e-3, 0.5e-3, 2.5e-3, 2e-4]),
            cin_fn=dict(k=[.15, .1, .25, 0.01]),
            dur_cin=dict(par1=[4.5, 3.5, 5.5, 0.5], par2=[20, 16, 24, 0.5]),
        ),
        ohr=dict(
            dur_precin=dict(par1=[2.5, 1, 10, 0.5], par2=[9, 5, 15, 0.5]),
            cancer_fn=dict(transform_prob=[1.5e-3, 0.5e-3, 2.5e-3, 2e-4]),
            cin_fn=dict(k=[.15, .1, .25, 0.01]),
            dur_cin=dict(par1=[4.5, 3.5, 5.5, 0.5], par2=[20, 16, 24, 0.5]),
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
        rel_reactivation_prob=[3, 2, 5, 0.5],
        #time_to_hiv_death_scale=lambda age: np.ones_like(age) * 1.0  # or whatever scale is appropriate
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
                            name=f'{location}_calib',
                            datafiles=datafiles,
                            extra_sim_result_keys=extra_sim_result_keys,
                            total_trials=n_trials, n_workers=n_workers,
                            storage=storage
                            )
    calib.calibrate()
    filename = f'{location}_calib{filestem}'
    if do_plot:
        calib.plot(do_save=True, fig_path=f'figures/{filename}.png')
    if do_save:
        sc.saveobj(f'results/{filename}.obj', calib)

    print(f'Best pars are {calib.best_pars}')

    return sim, calib

########################################################################
# Load pre-run calibration
########################################################################
def load_calib(location=None, do_plot=True, which_pars=0, save_pars=True, filestem=''):
    fnlocation = location.replace(' ', '_')
    filename = f'{fnlocation}_calib{filestem}'
    calib = sc.load(f'results/rwanda_calib.obj')
    if do_plot:
        sc.fonts(add=sc.thisdir(aspath=True) / 'Libertinus Sans')
        sc.options(font='Libertinus Sans')
        fig = calib.plot(res_to_plot=200, plot_type='sns.boxplot', do_save=False)
        fig.suptitle(f'Calibration results, {location.capitalize()}')
        fig.tight_layout()
        fig.savefig(f'figures/{filename}.png')

    if save_pars:
       pars_path = f'results/{filename}_pars.obj'
       if os.path.exists(pars_path):
        calib_pars = sc.load(pars_path)
    else:
        print(f"Calibration parameter file {pars_path} not found. Run calibration first.")
        calib_pars = None 
    # If you want to save new parameters after calibration, keep this:
    calib_pars = calib.trial_pars_to_sim_pars(which_pars=which_pars)
    sc.save(pars_path, calib_pars)

    return calib


# %% Run as a script
if __name__ == '__main__':

    T = sc.timer()
    location = 'rwanda'

    # Run calibration - usually on VMs
    if 'run_calibration' in to_run:
        filestem = ''
        sim, calib = run_calib(location=location, n_trials=n_trials, n_workers=n_workers,
                                   do_save=do_save, do_plot=False, filestem=filestem)

    # Load the calibration, plot it, and save the best parameters -- usually locally
    if 'plot_calibration' in to_run:
        filestem = ''

        calib = load_calib(location=location, do_plot=True, save_pars=True, filestem=filestem)

        best_par_ind = calib.df.index[0]
        extra_sim_results = calib.extra_sim_results[best_par_ind]
        years = calib.sim.results['year']
        year_ind = sc.findinds(years, 1985)[0]

        import matplotlib.pylab as pl

        fig, axes = pl.subplots(3, 1)
        axes[0].plot(years[year_ind:], extra_sim_results['cancers_with_hiv'][year_ind:], label='HIV+')
        axes[0].plot(years[year_ind:], extra_sim_results['cancers_no_hiv'][year_ind:], label='HIV-')
        axes[0].plot(years[year_ind:], extra_sim_results['cancers'][year_ind:], label='Total')
        axes[0].set_title(f'Cancers over time')
        axes[0].legend()
        axes[1].plot(calib.sim.pars['age_bin_edges'][:-1],
                         extra_sim_results['cancer_incidence_by_age_with_hiv'][:, -2], label='HIV+')
        axes[1].plot(calib.sim.pars['age_bin_edges'][:-1],
                         extra_sim_results['cancer_incidence_by_age_no_hiv'][:, -2],
                         label='HIV-')
        axes[1].legend()

        axes[2].plot(years[year_ind:], extra_sim_results['asr_cancer_incidence'][year_ind:])

        fig.show()
        fig.savefig(f'cancers_hiv_calib.png')

    T.toc('Done')
