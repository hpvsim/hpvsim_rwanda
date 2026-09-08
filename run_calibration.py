"""
Run calibration for HPVsim Rwanda (v3.2 port).

Three modes:
  python run_calibration.py --run-sim    # run calibration (heavy, VM-side)
  python run_calibration.py --plot       # produce hpv.plot_calibration figure
  python run_calibration.py              # (extract CSVs; stub until step 5)
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


# Run settings
debug = False
n_trials = [3000, 10][debug]
n_workers = 75
n_to_save = 500


# File-level weights spread across the target columns each file produces.
# The default eval_fn multiplies per-column gof by weight; a file whose 5
# splits across 16 age bins contributes 5/16 per bin, matching the per-file
# intent regardless of how many bins the loader produces.
_FILE_WEIGHTS = {
    'cancer_cases':              (0.5,  lambda c: c.startswith('all_hpv.cancers.')),
    'asr_cancer_incidence':      (10.0, lambda c: c == 'all_hpv.asr_cancer_incidence'),
    'cancer_incidence_with_hiv': (0.5,  lambda c: c == 'all_hpv.cancer_incidence_with_hiv'),
    'cancer_types':              (0.5,  lambda c: c.startswith('by_genotype.cancerous_genotype_dist.')),
    'precin_types':              (0.5,  lambda c: c.startswith('by_genotype.precin_genotype_dist.')),
}


def _build_weights(datafiles):
    from hpvsim.data.loaders import load_calib_data
    columns = list(load_calib_data(datafiles).columns)
    weights = {}
    for file_w, matches in _FILE_WEIGHTS.values():
        cols = [c for c in columns if matches(c)]
        if cols:
            per_col = file_w / len(cols)
            for c in cols:
                weights[c] = per_col
    return weights


########################################################################
# Run calibration
########################################################################
def run_calib(n_trials=None, n_workers=None, do_plot=False, do_save=True,
              n_to_save=None, filestem=''):
    """Fit Rwanda parameters on cancer + genotype + HIV targets."""

    sim = rs.make_sim(calib=True, use_calib=False)

    dataloc = 'data/rwanda'
    datafiles = [
        f'{dataloc}_cancer_cases.csv',
        f'{dataloc}_cancer_incidence_with_hiv.csv',
        f'{dataloc}_asr_cancer_incidence.csv',
        f'{dataloc}_precin_types.csv',
        f'{dataloc}_cancer_types.csv',
    ]

    # v3 nested calib_pars: top-level keys are scopes (broadcast, per-genotype,
    # network, hiv); leaves are ``[best, low, high, step]`` lists.
    # cin_fn.k priors sit tight (+/-0.05) around the shipped defaults per
    # genotype -- natural history doesn't vary much between countries.
    calib_pars = dict(
        beta=[0.05, 0.02, 0.5, 0.02],

        hpv16=dict(cin_fn=dict(k=[0.30, 0.25, 0.35, 0.01])),  # default 0.30
        hpv18=dict(cin_fn=dict(k=[0.25, 0.20, 0.30, 0.01])),  # default 0.25
        hi5=dict(cin_fn=dict(k=[0.20, 0.15, 0.25, 0.01])),    # default 0.20
        ohr=dict(cin_fn=dict(k=[0.20, 0.15, 0.25, 0.01])),    # default 0.20

        # Sexual network. Cross-layer probabilities are ANNUAL since
        # hpvsim 2.3.0; v2.2.6 values were per-timestep at dt=0.25, so
        # converted via ``1 - (1 - p_step)**(1/dt)`` = ``1 - (1 - p)**4``:
        #   m_cross_layer: v2 [0.3, 0.1, 0.7]  -> v3 [0.76, 0.35, 0.95]
        #   f_cross_layer: v2 [0.4, 0.05, 0.7] -> v3 [0.87, 0.20, 0.95]
        # Bounds snapped to 0.05-divisible ranges so Optuna doesn't clamp
        # them silently. ``*_partners_casual`` are Poisson-lambda;
        # ``Pars.update`` routes a scalar into ``Dist.set(lam=...)``.
        network=dict(
            m_cross_layer=[0.76, 0.35, 0.95, 0.05],
            f_cross_layer=[0.87, 0.20, 0.95, 0.05],
            m_partners_casual=[0.5, 0.1, 0.6, 0.05],
            f_partners_casual=[0.2, 0.1, 0.6, 0.05],
        ),

        # HIV. v2's ``rel_sus.lt200`` / ``.gt200`` map to v3's ``rel_sus_lo``
        # (CD4 < cd4_threshold=200) and ``rel_sus_hi`` (CD4 >= 200). v2's
        # ``art_failure_prob`` becomes v3's ``p_effective_art = 1 - failure``.
        hiv=dict(
            rel_sus_lo=[2.25, 2, 5, 0.25],
            rel_sus_hi=[2.25, 2, 4, 0.25],
            rel_sev_lo=[2.25, 1.5, 5, 0.25],
            rel_sev_hi=[2.25, 1.5, 5, 0.25],
            p_effective_art=[0.9, 0.7, 0.95, 0.01],
        ),
    )

    # reseed=True makes rand_seed a searched par; the best trial saves both
    # (pars, rand_seed) so the fit is exactly reproducible. Otherwise, all
    # trials use one shared seed and the fit is a single stochastic draw.
    calib = hpv.Calibration(
        sim, calib_pars=calib_pars, data=datafiles,
        weights=_build_weights(datafiles),
        total_trials=n_trials, n_workers=n_workers,
        reseed=True,
        label='rwanda_calib',
    )
    calib.calibrate()

    filename = f'rwanda_calib{filestem}'
    if do_plot:
        os.makedirs('figures', exist_ok=True)
        hpv.plot_calibration(calib, fig_path=f'figures/{filename}.png')
    if do_save:
        # Two-tier artifact: raw_results/ is the full untracked calib
        # (needed for hpv.make_calib_sims re-runs); results/ is the shrunk,
        # tracked version.
        os.makedirs('raw_results', exist_ok=True)
        os.makedirs('results', exist_ok=True)
        sc.saveobj(f'raw_results/{filename}.obj', calib)
        shrunk = calib.shrink(n_results=n_to_save or 500)
        sc.saveobj(f'results/{filename}.obj', shrunk)
        sc.saveobj(f'results/rwanda_pars{filestem}.obj', calib.best_pars)

    print(f'Best pars are {calib.best_pars}')
    return sim, calib


def load_calib(filestem=''):
    """Load the shrunk calibration artifact committed under ``results/``."""
    return sc.load(f'results/rwanda_calib{filestem}.obj')


########################################################################
# Extract plot-ready CSVs from a calib object
########################################################################

# v3 result names -> v2 metric labels used in figS2_timeseries.csv and by
# plot_figS2_calib.py.
_HIV_TS_MAP = {
    'art_coverage':          'p_on_art',
    'female_hiv_prevalence': 'prevalence_f',
    'male_hiv_prevalence':   'prevalence_m',
    'hiv_infections':        'new_infections',
    'hiv_deaths':            'new_deaths',
}
_HPV_TS_KEYS = ('asr_cancer_incidence', 'cancer_incidence_with_hiv',
                'cancer_incidence_no_hiv')
_GENOTYPES = ('hpv16', 'hpv18', 'hi5', 'ohr')


def _snapshot(sim):
    ah = sim.results.all_hpv
    hiv = sim.results.hiv
    out = {k: ah[k].annualize().values for k in _HPV_TS_KEYS}
    for v2_key, v3_key in _HIV_TS_MAP.items():
        out[v2_key] = hiv[v3_key].annualize().values
    out['years'] = ah[_HPV_TS_KEYS[0]].annualize().timevec.years

    i2020 = sc.findnearest(sim.t.yearvec, 2020)
    precin_df = hpv.results_by_genotype(sim, 'n_precin', normalize=True)
    cancer_df = hpv.results_by_genotype(sim, 'cum_cancers', normalize=True)
    out['precin_genotype_dist'] = np.array([precin_df.iloc[i2020].get(g, 0.0) for g in _GENOTYPES])
    out['cancerous_genotype_dist'] = np.array([cancer_df.iloc[i2020].get(g, 0.0) for g in _GENOTYPES])

    by_age = sim.analyzers.get('all_hpv_by_age')
    if by_age is not None:
        df = by_age.to_dataframe('cancers')
        out['cancers_by_age'] = df.loc[2020].values
    return out


def _bxp_stats(arr):
    q1, med, q3 = np.percentile(arr, [25, 50, 75])
    iqr = q3 - q1
    lo = float(arr[arr >= q1 - 1.5 * iqr].min())
    hi = float(arr[arr <= q3 + 1.5 * iqr].max())
    return dict(q1=float(q1), med=float(med), q3=float(q3), whislo=lo, whishi=hi)


# Reruns extend one year past the calibration horizon so annualize()'s final
# bin isn't a partial year; the extra year is dropped when writing the CSVs.
_PLOT_YEAR_MAX = 2025


def save_calib_results(calib, resfolder='results', n=50, n_workers=None):
    """Rerun top-n trials via hpv.make_calib_sims, aggregate, write CSVs.

    HIV-strat age-binned panels (v2 cancers_by_age_{with,no}_hiv) are not
    produced: v3's hpv.by_age has no HIV strata. Deferred to step 5.
    """
    os.makedirs(resfolder, exist_ok=True)
    print(f'Rerunning top-{n} trials...')
    outs = hpv.make_calib_sims(
        calib, n=n, extract_fn=_snapshot, n_workers=n_workers,
        sim_kwargs=dict(stop=_PLOT_YEAR_MAX + 1),
    )

    years = outs[0]['years']
    keep = years <= _PLOT_YEAR_MAX

    ts_rows = []
    for rkey in list(_HPV_TS_KEYS) + list(_HIV_TS_MAP):
        stack = np.array([o[rkey] for o in outs])[:, keep]
        med, lo, hi = np.nanmedian(stack, axis=0), np.nanpercentile(stack, 2.5, axis=0), np.nanpercentile(stack, 97.5, axis=0)
        for yi, yr in enumerate(years[keep]):
            ts_rows.append(dict(year=float(yr), metric=rkey, med=float(med[yi]),
                                pi95_low=float(lo[yi]), pi95_high=float(hi[yi])))
    pd.DataFrame(ts_rows).to_csv(f'{resfolder}/figS2_timeseries.csv', index=False)

    for rkey in ('precin_genotype_dist', 'cancerous_genotype_dist'):
        stack = np.array([o[rkey] for o in outs])
        pd.DataFrame([dict(bin=bi, **_bxp_stats(stack[:, bi])) for bi in range(stack.shape[1])]) \
          .to_csv(f'{resfolder}/figS2_{rkey}.csv', index=False)

    if all('cancers_by_age' in o for o in outs):
        stack = np.array([o['cancers_by_age'] for o in outs])
        pd.DataFrame([dict(bin=bi, **_bxp_stats(stack[:, bi])) for bi in range(stack.shape[1])]) \
          .to_csv(f'{resfolder}/figS2_cancers.csv', index=False)

    _write_target_csvs(resfolder)
    print(f'Wrote figS2_*.csv to {resfolder}/')


def _write_target_csvs(resfolder):
    """Reshape source datafiles into figS2_target_<key>.csv (value[, year])."""
    src = 'data/rwanda'
    # Age-binned cancers 2020 (target for `cancers` panel).
    cases = pd.read_csv(f'{src}_cancer_cases.csv').sort_values('age')
    cases[['value']].to_csv(f'{resfolder}/figS2_target_cancers.csv', index=False)
    # HIV-stratified age panels: write targets so the plot can render them
    # if/when the matching model CSVs land.
    for stratum, out_key in [('no', 'cancer_incidence_no_hiv'),
                              ('with', 'cancer_incidence_with_hiv')]:
        df = pd.read_csv(f'{src}_cancer_incidence_by_age_{stratum}_hiv.csv').sort_values('age')
        df[['value']].to_csv(f'{resfolder}/figS2_target_{out_key}.csv', index=False)
    # ASR (single-year scalar target).
    asr = pd.read_csv(f'{src}_asr_cancer_incidence.csv')[['year', 'value']]
    asr.to_csv(f'{resfolder}/figS2_target_asr_cancer_incidence.csv', index=False)
    # Genotype distributions.
    for src_key, out_key in [('precin_types', 'precin_genotype_dist'),
                              ('cancer_types', 'cancerous_genotype_dist')]:
        df = pd.read_csv(f'{src}_{src_key}.csv')
        # match genotype order used in the model extract
        gorder = {g: i for i, g in enumerate(_GENOTYPES)}
        df['_ord'] = df['genotype'].map(lambda g: gorder.get(_normalize_geno(g), 99))
        df = df.sort_values('_ord')
        df[['value']].to_csv(f'{resfolder}/figS2_target_{out_key}.csv', index=False)


def _normalize_geno(g):
    """CSV genotype label -> the module name used on sim.diseases (hpv16, hi5, ...)."""
    m = {'16': 'hpv16', '18': 'hpv18', 'Hi5': 'hi5', 'hi5': 'hi5', 'OHR': 'ohr', 'ohr': 'ohr'}
    return m.get(str(g), str(g))


# %% Run as a script
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-sim', action='store_true',
                        help='Run calibration (heavy, VM-side)')
    parser.add_argument('--plot', action='store_true',
                        help='Run local calibration diagnostic plot')
    parser.add_argument('--extract-csvs', action='store_true',
                        help='Extract plot-ready CSVs (stubbed pending step 5)')
    parser.add_argument('--resfolder', default='results')
    args = parser.parse_args()

    T = sc.timer()
    if args.run_sim:
        sim, calib = run_calib(n_trials=n_trials, n_workers=n_workers,
                               n_to_save=n_to_save, do_save=True)
    else:
        calib = sc.load(f'{args.resfolder}/rwanda_calib.obj')

    if args.plot:
        os.makedirs('figures', exist_ok=True)
        hpv.plot_calibration(calib, fig_path='figures/rwanda_calib.png')

    if args.extract_csvs:
        save_calib_results(calib, resfolder=args.resfolder)
    T.toc('Done')
