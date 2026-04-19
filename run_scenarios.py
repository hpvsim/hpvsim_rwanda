"""
Run scenarios.

Two modes:
  python run_scenarios.py --run-sim   # run msim + save plot-ready CSVs (VM)
  python run_scenarios.py             # re-extract CSVs from an existing st_scens.obj
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
import hpvsim as hpv

# Imports from this repository
import run_sim as rs
from interventions import make_st, make_st_older, make_mv_intvs


# Metrics with low/high bounds (hpv.Result objects)
TS_METRICS = ['asr_cancer_incidence', 'cancer_incidence_with_hiv', 'cancer_incidence_no_hiv']
CUM_METRICS_BOUNDED = ['cancers', 'cancers_with_hiv', 'cancers_no_hiv', 'cancer_deaths']
CUM_METRICS_UNBOUNDED = ['ablations', 'txvs', 'vaccinations', 'screens', 'excisions',
                         'leeps', 'cancer_treatments']
CUM_START_YEAR = 2025


# Settings - used here and imported elsewhere
debug = 0
n_seeds = [10, 1][debug]  # How many seeds to run per cluster


# %% Create interventions
def make_baselines():
    """
    Baseline scenarios:
        1. No interventions
        2. Status quo screening + treatment
    """
    scendict = dict()
    scendict['No interventions'] = []
    scendict['Baseline'] = make_st(future_screen_cov=0.18, screen_change_year=2025)
    return scendict


def make_campaign_scenarios():
    """
    Scenarios for mass one-time campaigns:
        1. Mass delivery of virus-clearing TxV
        2. Mass delivery of lesion-regressing TxV
        3. "HPV-faster", with screening+treatment+vaccination
    """
    scendict = dict()
    age_range = [20, 50]

    for cov in [0.18, 0.35, 0.7]:

        # Virus-clearing mass TxV
        scendict[f'Mass TxV 90/0, {int(cov*100)}%'] = make_mv_intvs(
            txv_pars='precin',
            campaign_coverage=cov,
        )

        # Virus-clearing mass TxV
        scendict[f'Mass TxV 50/90, {int(cov*100)}%'] = make_mv_intvs(
            txv_pars='cin',
            campaign_coverage=cov,
        )

        # Screen, treat, & vaccinate older women
        mass_intvs = make_st_older(screen_cov=cov, age_range=age_range, start_year=2026)
        scendict[f'HPV-Faster {cov*100:.0f}%'] = mass_intvs

    return scendict


def make_st_scenarios():
    """
    Compare vaccination strategies:
    """
    scendict = dict()

    start_year = 2026
    cov_array = [.18, .35, .70]
    for cov_val in cov_array:

        # Scale up S&T&T - default algorithm, screening + triage + treatment
        st_intvs = make_st(screen_change_year=start_year, future_screen_cov=cov_val)
        scendict[f'S&T&T {cov_val*100:.0f}%'] = st_intvs

        # Scale up S&T - streamlined algorithm, screening + treatment only
        st_intvs = make_st(screen_change_year=start_year, future_screen_cov=cov_val, tx_assigner_csv='tx_assigner_no_triage')
        scendict[f'S&T {cov_val*100:.0f}%'] = st_intvs

        # Scale up S&TxV&T&T - enhanced algorithm with virus-clearing TxV used in addition to triage and treatment
        st_intvs = make_st(
            screen_change_year=start_year,
            future_screen_cov=cov_val,
            txv_pars='precin',
            txv=True)
        scendict[f'S&TxV&T&T {cov_val*100:.0f}%'] = st_intvs

        # Scale up S&TxV - enhanced and streamlined algorithm with virus-clearing TxV used instead of treatment
        st_intvs = make_st(
            screen_change_year=start_year,
            future_screen_cov=cov_val,
            txv_pars='cin',
            txv=True)
        scendict[f'S&TxV {cov_val*100:.0f}%'] = st_intvs

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
            sim = rs.make_sim(
                debug=debug,
                add_st=False,
                add_vax=add_vax,
                interventions=interventions,
                end=end,
                seed=seed,
            )
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


def process_msim(msim, scenarios):
    """Reduce msim → per-scenario dict of year + metric arrays (with low/high)."""
    metrics = ['year'] + TS_METRICS + CUM_METRICS_BOUNDED

    scen_labels = list(scenarios.keys())
    mlist = msim.split(chunks=len(scen_labels))

    msim_dict = sc.objdict()
    for si, scen_label in enumerate(scen_labels):
        reduced_sim = mlist[si].reduce(output=True)
        mres = sc.objdict({metric: reduced_sim.results[metric] for metric in metrics})

        # Sum intervention product counts across matching interventions
        programs = {
            'mass_vax': 'vaccinations',
            'screening': 'screens',
            'ablation': 'ablations',
            'excision': 'leeps',
            'radiation': 'cancer_treatments',
            'txv': 'txvs',
            'campaign txvx': 'txvs',
            'ablation_older': 'ablations',
            'excision_older': 'excisions',
            'radiation_older': 'cancer_treatments',
        }
        for intv_name in set(programs.values()):
            mres[intv_name] = np.zeros_like(mres.year)
        for intv_name, df_key in programs.items():
            if reduced_sim.get_intervention(intv_name, die=False) is not None:
                mres[df_key] += reduced_sim.get_intervention(intv_name).n_products_used.values

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
    scenarios = sc.mergedicts(make_baselines(), make_st_scenarios(), make_campaign_scenarios())

    if args.run_sim:
        msim = run_sims(scenarios=scenarios, end=args.end)
        msim_dict = process_msim(msim, scenarios)
        sc.saveobj(f'{args.resfolder}/st_scens.obj', msim_dict)
    else:
        msim_dict = sc.loadobj(f'{args.resfolder}/st_scens.obj')

    save_csvs(msim_dict, resfolder=args.resfolder)
    print(f'Saved scens_timeseries.csv + scens_cumulative.csv to {args.resfolder}/')
    T.toc('Done')
