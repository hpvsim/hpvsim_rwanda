"""
Run scenarios (v3.2 port).

Two modes:
  python run_scenarios.py --run-sim   # run msim + save plot-ready CSVs (VM)
  python run_scenarios.py             # re-extract CSVs from existing st_scens.obj
"""


# %% General settings

import argparse
import os

os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

# Standard imports
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss

# Imports from this repository
import run_sim as rs
from interventions import make_st, make_st_older, make_mv_intvs


# Time-series cancer metrics (v3 module scope: sim.results.all_hpv.<key>)
TS_METRICS = ['asr_cancer_incidence', 'cancer_incidence_with_hiv', 'cancer_incidence_no_hiv']

# Per-step count metrics with low/high bounds. Column names kept v2-style so
# the plot scripts don't need to change; v3 sim.results.all_hpv exposes them
# under new_cancers / new_cancer_deaths, mapped via _V3_ALIAS below.
CUM_METRICS_BOUNDED = ['cancers', 'cancers_with_hiv', 'cancers_no_hiv', 'cancer_deaths']
_V3_ALIAS = {'cancers': 'new_cancers', 'cancer_deaths': 'new_cancer_deaths'}

# Program-level counters extracted from per-intervention results below.
CUM_METRICS_UNBOUNDED = ['ablations', 'txvs', 'vaccinations', 'screens', 'excisions',
                         'leeps', 'cancer_treatments']

CUM_START_YEAR = 2025

# v2 read `intv.n_products_used` uniformly; v3 splits per intervention type:
#   BaseVaccination  -> new_doses
#   BaseScreening    -> n_screened
#   BaseTreatment    -> new_cin_treated (CIN) OR new_cancer_treated (radiation)
#   BaseTxVx         -> new_txvx_doses
# The intervention `name`s here match those set in interventions.py.
INTV_TO_METRIC = {
    # baseline S&T (from make_st)
    'screening':       ('screens',           'n_screened'),
    'ablation_intv':   ('ablations',         'new_cin_treated'),
    'excision_intv':   ('leeps',             'new_cin_treated'),
    'radiation_intv':  ('cancer_treatments', 'new_cancer_treated'),
    'txv':             ('txvs',              'new_txvx_doses'),
    # mass therapeutic-vax campaign (from make_mv_intvs)
    'campaign_txvx':   ('txvs',              'new_txvx_doses'),
    # older-cohort screen-and-vax (from make_st_older)
    'screening_older': ('screens',           'n_screened'),
    'ablation_older':  ('ablations',         'new_cin_treated'),
    'excision_older':  ('excisions',         'new_cin_treated'),
    'radiation_older': ('cancer_treatments', 'new_cancer_treated'),
    'mass_vax':        ('vaccinations',      'new_doses'),
    # routine childhood HPV vax (from make_vx) — attached to every non-baseline scenario
    'routine_vx':      ('vaccinations',      'new_doses'),
}


# Settings - used here and imported elsewhere
debug = 0
n_seeds = [10, 1][debug]  # How many seeds to run per cluster


# %% Create interventions
def make_baselines(end_year=2100):
    """
    Baseline scenarios:
        1. No interventions
        2. Status quo screening + treatment
    """
    scendict = dict()
    scendict['No interventions'] = []
    scendict['Baseline'] = make_st(future_screen_cov=0.18, screen_change_year=2025, end_year=end_year)
    return scendict


def make_campaign_scenarios(end_year=2100):
    """
    Scenarios for mass one-time campaigns:
        1. Mass delivery of virus-clearing TxV
        2. Mass delivery of lesion-regressing TxV
        3. "HPV-faster", with screening+treatment+vaccination
    """
    scendict = dict()
    age_range = [20, 50]

    for cov in [0.18, 0.35, 0.7]:
        scendict[f'Mass TxV 90/0, {int(cov*100)}%'] = make_mv_intvs(
            txv_pars='precin', campaign_coverage=cov, end_year=end_year,
        )
        scendict[f'Mass TxV 50/90, {int(cov*100)}%'] = make_mv_intvs(
            txv_pars='cin', campaign_coverage=cov, end_year=end_year,
        )
        mass_intvs = make_st_older(screen_cov=cov, age_range=age_range,
                                   start_year=2026, end_year=end_year)
        scendict[f'HPV-Faster {cov*100:.0f}%'] = mass_intvs

    return scendict


def make_st_scenarios(end_year=2100):
    """
    Compare screen-and-treat vaccination strategies at three coverage levels.
    """
    scendict = dict()

    start_year = 2026
    for cov_val in [.18, .35, .70]:
        st_intvs = make_st(screen_change_year=start_year, future_screen_cov=cov_val,
                           end_year=end_year)
        scendict[f'S&T&T {cov_val*100:.0f}%'] = st_intvs

        st_intvs = make_st(screen_change_year=start_year, future_screen_cov=cov_val,
                           tx_assigner_csv='tx_assigner_no_triage', end_year=end_year)
        scendict[f'S&T {cov_val*100:.0f}%'] = st_intvs

        st_intvs = make_st(screen_change_year=start_year, future_screen_cov=cov_val,
                           txv_pars='precin', txv=True, end_year=end_year)
        scendict[f'S&TxV&T&T {cov_val*100:.0f}%'] = st_intvs

        st_intvs = make_st(screen_change_year=start_year, future_screen_cov=cov_val,
                           txv_pars='cin', txv=True, end_year=end_year)
        scendict[f'S&TxV {cov_val*100:.0f}%'] = st_intvs

    return scendict


def make_sims(scenarios=None, end=2100):
    """One flat list of sims (n_scenarios * n_seeds) labelled by scenario.

    v2 built one MultiSim per scenario and then `hpv.MultiSim.merge`d them
    into a super-MultiSim; v3 `ss.MultiSim` has neither `merge` nor
    `split`, so we flatten into one MultiSim and slice back at reduce time
    using scenario boundaries.
    """
    sims = sc.autolist()
    for name, interventions in scenarios.items():
        add_vax = name != 'No interventions'
        for seed in range(n_seeds):
            sim = rs.make_sim(
                debug=debug,
                add_st=False,
                add_vax=add_vax,
                interventions=interventions,
                stop=end,
                seed=seed,
            )
            sim.label = name
            sims += sim
    return ss.MultiSim(sims)


def run_sims(scenarios=None, end=2100, verbose=-1):
    """ Run the simulations """
    msim = make_sims(scenarios=scenarios, end=end)
    msim.run(verbose=verbose)
    return msim


def process_msim(msim, scenarios):
    """Reduce per-scenario slices → dict of year + metric arrays with low/high."""
    scen_labels = list(scenarios.keys())
    msim_dict = sc.objdict()
    for si, scen_label in enumerate(scen_labels):
        scen_sims = list(msim.sims[si * n_seeds : (si + 1) * n_seeds])
        scen_msim = ss.MultiSim(scen_sims)
        reduced_sim = scen_msim.reduce(output=True)

        # sim.results.timevec is a DateArray of ss.date; .t.yearvec is the
        # matching float-year ndarray downstream code (save_csvs, plots) expects.
        year = np.asarray(reduced_sim.t.yearvec)
        mres = sc.objdict(year=year)
        for metric in TS_METRICS + CUM_METRICS_BOUNDED:
            mres[metric] = reduced_sim.results.all_hpv[_V3_ALIAS.get(metric, metric)]

        # Zero-init program buckets, then sum per-intervention counters in.
        for bucket in set(m for (m, _) in INTV_TO_METRIC.values()):
            mres[bucket] = np.zeros_like(year, dtype=float)
        for intv_name, (bucket, counter) in INTV_TO_METRIC.items():
            intv = reduced_sim.interventions.get(intv_name)
            if intv is None:
                continue
            mres[bucket] = mres[bucket] + np.asarray(intv.results[counter])

        msim_dict[scen_label] = mres

    return msim_dict


def save_csvs(msim_dict, resfolder='results'):
    """Extract two plot-ready CSVs from an msim_dict.

    scens_timeseries.csv — year, scenario, metric, value, low, high
                           (for asr + cancer_incidence_with_hiv + cancer_incidence_no_hiv)
    scens_cumulative.csv — scenario, metric, value[, low, high]
                           (sums 2025-2100 for cancers, cancers_with_hiv, ablations, txvs, vaccinations, ...)
    """
    os.makedirs(resfolder, exist_ok=True)

    # Time series (only plotted metrics, full year range)
    ts_rows = []
    for scen_label, mres in msim_dict.items():
        years = np.asarray(mres.year)
        for metric in TS_METRICS:
            r = mres[metric]
            for yi, yr in enumerate(years):
                ts_rows.append({
                    'scenario': scen_label, 'year': float(yr), 'metric': metric,
                    'value': float(r[yi]),
                    'low': float(r.low[yi]),
                    'high': float(r.high[yi]),
                })
    pd.DataFrame(ts_rows).to_csv(f'{resfolder}/scens_timeseries.csv', index=False)

    # Cumulative sums from CUM_START_YEAR → end
    cum_rows = []
    for scen_label, mres in msim_dict.items():
        years = np.asarray(mres.year)
        fi = int(np.where(years == CUM_START_YEAR)[0][0])
        for metric in CUM_METRICS_BOUNDED:
            r = mres[metric]
            cum_rows.append({
                'scenario': scen_label, 'metric': metric,
                'value': float(np.sum(r.values[fi:])),
                'low': float(np.sum(r.low[fi:])),
                'high': float(np.sum(r.high[fi:])),
            })
        for metric in CUM_METRICS_UNBOUNDED:
            r = np.asarray(mres[metric])
            cum_rows.append({
                'scenario': scen_label, 'metric': metric,
                'value': float(np.sum(r[fi:])),
                'low': np.nan, 'high': np.nan,
            })
    pd.DataFrame(cum_rows).to_csv(f'{resfolder}/scens_cumulative.csv', index=False)


# %% Run as a script
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--run-sim', action='store_true',
                        help='Run scenarios on the VM (heavy); otherwise only re-extract CSVs')
    parser.add_argument('--end', type=int, default=2100)
    parser.add_argument('--resfolder', default='results')
    args = parser.parse_args()

    T = sc.timer()
    scenarios = sc.mergedicts(make_baselines(args.end),
                              make_st_scenarios(args.end),
                              make_campaign_scenarios(args.end))

    if args.run_sim:
        msim = run_sims(scenarios=scenarios, end=args.end)
        msim_dict = process_msim(msim, scenarios)
        sc.saveobj(f'{args.resfolder}/st_scens.obj', msim_dict)
    else:
        msim_dict = sc.loadobj(f'{args.resfolder}/st_scens.obj')

    save_csvs(msim_dict, resfolder=args.resfolder)
    print(f'Saved scens_timeseries.csv + scens_cumulative.csv to {args.resfolder}/')
    T.toc('Done')
