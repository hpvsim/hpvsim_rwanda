"""
Run calibration for HPVsim Rwanda.

Three modes:
  python run_calibration.py --run-sim    # run calibration (heavy, VM-side)
  python run_calibration.py              # extract plot-ready CSVs from existing rwanda_calib.obj
  python run_calibration.py --plot       # run local calibration plot via hpv.Calibration
"""
import argparse
import os

os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

import numpy as np
import pandas as pd
import sciris as sc
import hpvsim as hpv

import run_sim as rs
import utils as ut


# Run settings
debug = False
n_trials = [3000, 10][debug]
n_workers = 75
n_to_save = 500
storage = None

# Years over which extra_sim_results are extracted (matches hpv.Calibration conventions)
EXTRA_YEARS = np.arange(1960, 2026)
EXTRA_SI = int(np.where(EXTRA_YEARS == 2000)[0][0])


########################################################################
# Run calibration
########################################################################
def run_calib(n_trials=None, n_workers=None, do_plot=False, do_save=True,
              n_to_save=None, filestem=''):

    sim = rs.make_sim(calib=True, use_calib=False)

    dataloc = 'data/rwanda'
    datafiles = [
        f'{dataloc}_cancer_cases.csv',
        f'{dataloc}_cancer_incidence_by_age_no_hiv.csv',
        f'{dataloc}_cancer_incidence_by_age_with_hiv.csv',
        f'{dataloc}_asr_cancer_incidence.csv',
        f'{dataloc}_precin_types.csv',
        f'{dataloc}_cancer_types.csv',
    ]

    calib_pars = dict(
        beta=[0.05, 0.02, 0.5, 0.02],
        sev_dist=dict(par1=[1, 0.5, 1.5, 0.01]),
    )
    sexual_behavior_pars = dict(
        m_cross_layer=[0.3, 0.1, 0.7, 0.05],
        m_partners=dict(c=dict(par1=[0.5, 0.1, 0.6, 0.05])),
        f_cross_layer=[0.4, 0.05, 0.7, 0.05],
        f_partners=dict(c=dict(par1=[0.2, 0.1, 0.6, 0.05])),
    )
    calib_pars = sc.mergedicts(calib_pars, sexual_behavior_pars)

    genotype_pars = dict(
        hi5=dict(cin_fn=dict(k=[.15, .1, .25, 0.01])),
        ohr=dict(cin_fn=dict(k=[.15, .1, .25, 0.01])),
    )

    hiv_pars = dict(
        rel_sus=dict(
            lt200=[2.25, 2, 5, 0.25],
            gt200=[2.25, 2, 4, 0.25],
        ),
        rel_sev=dict(
            lt200=[2.25, 1.5, 5, 0.25],
            gt200=[2.25, 1.5, 5, 0.25],
        ),
        art_failure_prob=[0.1, 0.05, 0.3, 0.01],
    )

    extra_sim_result_keys = [
        'cancers', 'cancers_with_hiv', 'cancers_no_hiv',
        'cancers_by_age_with_hiv', 'cancers_by_age_no_hiv',
        'asr_cancer_incidence', 'cancer_incidence_with_hiv', 'cancer_incidence_no_hiv',
        'art_coverage', 'female_hiv_prevalence', 'male_hiv_prevalence',
        'n_females_with_hiv_alive', 'n_males_with_hiv_alive',
        'hiv_infections', 'hiv_deaths',
    ]

    calib = hpv.Calibration(
        sim, calib_pars=calib_pars, genotype_pars=genotype_pars, hiv_pars=hiv_pars,
        name='rwanda_calib', datafiles=datafiles,
        extra_sim_result_keys=extra_sim_result_keys,
        total_trials=n_trials, n_workers=n_workers, storage=storage,
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


def load_calib(do_plot=True, filestem=''):
    filename = f'rwanda_calib{filestem}'
    calib = sc.load(f'results/rwanda_calib.obj')
    if do_plot:
        fig = calib.plot(res_to_plot=None, plot_type='sns.boxplot', do_save=False)
        fig.suptitle('Calibration results')
        fig.tight_layout()
        fig.savefig(f'figures/{filename}.png')
    return calib


########################################################################
# Extract plot-ready CSVs from a calib object
########################################################################
def _bxp_stats(arr):
    q1, med, q3 = np.percentile(arr, [25, 50, 75])
    iqr = q3 - q1
    lo = float(arr[arr >= q1 - 1.5 * iqr].min())
    hi = float(arr[arr <= q3 + 1.5 * iqr].max())
    return dict(q1=float(q1), med=float(med), q3=float(q3),
                whislo=lo, whishi=hi)


def save_figS2_csvs(calib, resfolder='results'):
    os.makedirs(resfolder, exist_ok=True)
    analyzer_results = calib.analyzer_results
    sim_results = calib.sim_results
    extra_sim_results = calib.extra_sim_results

    # ---- Age-binned cancer metrics (boxplot stats per bin across runs) ----
    age_rkey_to_year = {
        'cancers': 2020,
        'cancer_incidence_no_hiv': 2017,
        'cancer_incidence_with_hiv': 2017,
    }
    for rkey, year in age_rkey_to_year.items():
        stacked = np.array([run[rkey][year] for run in analyzer_results])
        rows = [{'bin': bi, **_bxp_stats(stacked[:, bi])}
                for bi in range(stacked.shape[1])]
        pd.DataFrame(rows).to_csv(f'{resfolder}/figS2_{rkey}.csv', index=False)

    # ---- Time series: med + pi95 across runs ----
    ts_rkeys = ['asr_cancer_incidence', 'cancer_incidence_with_hiv',
                'cancer_incidence_no_hiv', 'art_coverage',
                'female_hiv_prevalence', 'male_hiv_prevalence',
                'hiv_infections', 'hiv_deaths']
    rows = []
    for rkey in ts_rkeys:
        stacked = np.array([run[rkey] for run in extra_sim_results])
        med = np.median(stacked, axis=0)
        lo = np.percentile(stacked, 2.5, axis=0)
        hi = np.percentile(stacked, 97.5, axis=0)
        for yi, yr in enumerate(EXTRA_YEARS):
            if yi < EXTRA_SI:
                continue
            rows.append({'year': int(yr), 'metric': rkey,
                         'med': float(med[yi]),
                         'pi95_low': float(lo[yi]),
                         'pi95_high': float(hi[yi])})
    pd.DataFrame(rows).to_csv(f'{resfolder}/figS2_timeseries.csv', index=False)

    # ---- Genotype distributions (boxplot stats per bin) ----
    for rkey in ['precin_genotype_dist', 'cancerous_genotype_dist']:
        stacked = np.array([run[rkey] for run in sim_results])
        rows = [{'bin': bi, **_bxp_stats(stacked[:, bi])}
                for bi in range(stacked.shape[1])]
        pd.DataFrame(rows).to_csv(f'{resfolder}/figS2_{rkey}.csv', index=False)

    # ---- Target (observed) data ----
    target_names = ['cancers', 'cancer_incidence_no_hiv', 'cancer_incidence_with_hiv',
                    'asr_cancer_incidence', 'precin_genotype_dist', 'cancerous_genotype_dist']
    for ti, name in enumerate(target_names):
        if ti < len(calib.target_data):
            calib.target_data[ti].to_csv(f'{resfolder}/figS2_target_{name}.csv', index=False)


# %% Run as a script
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-sim', action='store_true',
                        help='Run calibration (heavy, VM-side)')
    parser.add_argument('--plot', action='store_true',
                        help='Run local calibration diagnostic plot')
    parser.add_argument('--resfolder', default='results')
    args = parser.parse_args()

    T = sc.timer()
    if args.run_sim:
        sim, calib = run_calib(n_trials=n_trials, n_workers=n_workers,
                               n_to_save=n_to_save, do_save=True)
    else:
        calib = sc.load(f'{args.resfolder}/rwanda_calib.obj')

    if args.plot:
        fig = calib.plot(res_to_plot=None, plot_type='sns.boxplot', do_save=False)
        fig.suptitle('Calibration results')
        fig.tight_layout()
        os.makedirs('figures', exist_ok=True)
        fig.savefig('figures/rwanda_calib.png')

    save_figS2_csvs(calib, resfolder=args.resfolder)
    print(f'Saved figS2_*.csv to {args.resfolder}/')
    T.toc('Done')
