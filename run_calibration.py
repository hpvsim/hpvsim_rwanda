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

import sciris as sc
import hpvsim as hpv

import run_sim as rs


# Run settings
debug = False
n_trials = [3000, 10][debug]
n_workers = 75
n_to_save = 500


########################################################################
# Run calibration
########################################################################
def run_calib(n_trials=None, n_workers=None, do_plot=False, do_save=True,
              n_to_save=None, filestem=''):
    """Fit Rwanda parameters on cancer + genotype + HIV targets."""

    sim = rs.make_sim(calib=True, use_calib=False)

    dataloc = 'data/rwanda'
    # v3 hpv.by_age ships `cancers` but not the HIV-stratified age keys v2
    # fit via AgeResults (`cancers_by_age_with_hiv` / `_no_hiv`). Fit the
    # scalar HIV-strat incidence here; the age x HIV Fig S2 panels come
    # back post-calibration via a custom analyzer at step 5 of the plan.
    datafiles = [
        f'{dataloc}_cancer_cases.csv',
        f'{dataloc}_cancer_incidence_with_hiv.csv',
        f'{dataloc}_asr_cancer_incidence.csv',
        f'{dataloc}_precin_types.csv',
        f'{dataloc}_cancer_types.csv',
    ]

    # v3 nested calib_pars: top-level keys are scopes (broadcast, hi5, ohr,
    # network, hiv); leaves are ``[best, low, high, step]`` lists.
    calib_pars = dict(
        # Broadcast to every HPV genotype
        beta=[0.05, 0.02, 0.5, 0.02],

        # Per-genotype natural-history slope for the two pooled genotypes.
        # v2 also carried a broadcast ``sev_dist.par1``; v3's per-genotype
        # ``cin_fn.k`` covers the same effect, and Rwanda only ever
        # calibrated hi5 / ohr, so the broadcast has no replacement here.
        hi5=dict(cin_fn=dict(k=[0.15, 0.1, 0.25, 0.01])),
        ohr=dict(cin_fn=dict(k=[0.15, 0.1, 0.25, 0.01])),

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

    calib = hpv.Calibration(
        sim, calib_pars=calib_pars, data=datafiles,
        total_trials=n_trials, n_workers=n_workers,
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
def save_figS2_csvs(calib, resfolder='results'):
    """Extract plot-ready CSVs from a calibration artifact.

    v2 read directly from ``calib.analyzer_results`` / ``sim_results`` /
    ``extra_sim_results``; v3 hpv.Calibration stores only per-trial
    mismatch and eval-column values. Age-binned cancers, per-run genotype
    distributions, and time-series like ``asr_cancer_incidence`` have to
    be rebuilt by rerunning the top-N trials via
    ``hpv.make_calib_sims(calib, n=..., extract_fn=...)``.

    Deferred to step 5 of docs/revision_plan.md (port plot scripts).
    """
    raise NotImplementedError(
        'save_figS2_csvs pending v3.2 rewrite (docs/revision_plan.md step 5). '
        'Use hpv.make_calib_sims(calib, n=n_to_save, extract_fn=...) to pull '
        'per-trial time series, age-binned cancers, and genotype distributions.'
    )


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
        save_figS2_csvs(calib, resfolder=args.resfolder)
        print(f'Saved figS2_*.csv to {args.resfolder}/')
    T.toc('Done')
