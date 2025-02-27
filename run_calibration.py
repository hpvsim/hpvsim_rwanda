import os
import sciris as sc
import hpvsim as hpv
import pylab as pl
import run_sim_hiv as rs  # Ensure this script is adapted for Rwanda

# Set environment variables to handle multithreading
os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

# User configurations
location = 'rwanda'
to_run = ['run_calibration']  # Change to 'plot_calibration' if plotting only
debug = False
do_save = True

# Calibration settings
n_trials = [50, 5][debug]
n_workers = [25, 1][debug]
storage = "sqlite:///hpvsim_rwanda.db"

def run_calib(location=None, n_trials=None, n_workers=None, do_plot=False, do_save=True, filestem=''):
    dflocation = location.replace(" ", "_")
    hiv_datafile = [
        'data/rwanda_hiv_incidence.csv',
        'data/rwanda_female_hiv_mortality.csv',
        'data/rwanda_male_hiv_mortality.csv'
    ]
    art_datafile = [
        'data/rwanda_art_coverage_by_age_males.csv',
        'data/rwanda_art_coverage_by_age_females.csv'
    ]
    
    sim = rs.make_sim(hiv_datafile=hiv_datafile, art_datafile=art_datafile, calib=True)
    
    datafiles = [
        'data/rwanda_cancer_cases.csv',
        'data/rwanda_asr_cancer_incidence.csv',
    ]
    
    # Update calibration parameters based on Rwanda data
    calib_pars = dict(
        beta=[0.04, 0.02, 0.5, 0.02],
        own_imm_hr=[0.6, 0.3, 1, 0.05],
        age_risk=dict(risk=[3.5, 1, 4, 0.1], age=[39, 30, 45, 1]),
        hpv_control_prob=[0, 0, 1, 0.25],
        hpv_reactivation=[0.03, 0, 0.1, 0.025]
    )
    
    sexual_behavior_pars = dict(
        m_cross_layer=[0.35, 0.1, 0.7, 0.05],
        m_partners=dict(c=dict(par1=[0.55, 0.1, 0.6, 0.05])),
        f_cross_layer=[0.45, 0.05, 0.7, 0.05],
        f_partners=dict(c=dict(par1=[0.25, 0.1, 0.6, 0.05]))
    )
    
    calib_pars = sc.mergedicts(calib_pars, sexual_behavior_pars)
    
    hiv_pars = dict(
        rel_sus=dict(lt200=[2.3, 2, 5, 0.25], gt200=[2.3, 2, 4, 0.25]),
        rel_sev=dict(lt200=[1.6, 1.25, 5, 0.25], gt200=[1.6, 1.25, 3, 0.25]),
        rel_reactivation_prob=[3.2, 2, 5, 0.5]
    )
    
    extra_sim_result_keys = ['cancers', 'cancers_with_hiv', 'cancers_no_hiv',
                             'cancers_by_age_with_hiv', 'cancers_by_age_no_hiv',
                             'asr_cancer_incidence', 'cancer_incidence_by_age_with_hiv',
                             'cancer_incidence_by_age_no_hiv']
    
    calib = hpv.Calibration(
        sim,
        calib_pars=calib_pars,
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
        sc.saveobj(f'results/{filename}_hiv.obj', calib)
    print(f'Best parameters: {calib.best_pars}')
    return sim, calib

if __name__ == '__main__':
    T = sc.timer()
    if 'run_calibration' in to_run:
        filestem = ''
        sim, calib = run_calib(location=location, n_trials=n_trials, n_workers=n_workers,
                               do_save=do_save, do_plot=False, filestem=filestem)
    T.toc('Calibration for Rwanda complete')