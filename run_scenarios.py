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
import pandas as pd
import sciris as sc
import starsim as ss

# Imports from this repository
import run_sim as rs
from interventions import make_st, make_st_older, make_mv_intvs


# Time-series cancer metrics (v3 module scope: sim.results.all_hpv.<key>)
TS_METRICS = ['asr_cancer_incidence', 'cancer_incidence_with_hiv', 'cancer_incidence_no_hiv']

# Cancer flow metrics with plot-visible bounds. Column names kept v2-style
# (cancers / cancer_deaths); v3 exposes them as new_cancers / new_cancer_deaths
# on sim.results.all_hpv, mapped via _V3_ALIAS.
CUM_METRICS_BOUNDED = ['cancers', 'cancers_with_hiv', 'cancers_no_hiv', 'cancer_deaths']
_V3_ALIAS = {'cancers': 'new_cancers', 'cancer_deaths': 'new_cancer_deaths'}

CUM_START_YEAR = 2025

# v3 per-intervention flow counter for each program bucket. Names match
# the intervention `name=` in interventions.py.
INTV_TO_METRIC = {
    # baseline S&T (from make_st)
    'screening':       ('screens',           'new_screens'),
    'ablation_intv':   ('ablations',         'new_cin_treated'),
    'excision_intv':   ('leeps',             'new_cin_treated'),
    'radiation_intv':  ('cancer_treatments', 'new_cancer_treated'),
    'txv':             ('txvs',              'new_txvx_doses'),
    # mass therapeutic-vax campaign (from make_mv_intvs)
    'campaign_txvx':   ('txvs',              'new_txvx_doses'),
    # older-cohort screen-and-vax (from make_st_older)
    'screening_older': ('screens',           'new_screens'),
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


def _result_rows(result, metric):
    """Long-format annual rows for one ss.Result: (year, metric, value)."""
    df = result.annualize().to_df()
    return pd.DataFrame({
        'year': pd.to_datetime(df['timevec']).dt.year,
        'metric': metric,
        'value': df['value'].astype(float),
    })


def process_msim(msim, scenarios):
    """Long-format ensemble: rows are (scenario, sim, year, metric, value)."""
    frames = []
    for si, scen_label in enumerate(scenarios):
        for sim_idx, sim in enumerate(msim.sims[si * n_seeds : (si + 1) * n_seeds]):
            for metric in TS_METRICS + CUM_METRICS_BOUNDED:
                frames.append(_result_rows(
                    sim.results.all_hpv[_V3_ALIAS.get(metric, metric)], metric,
                ).assign(scenario=scen_label, sim=sim_idx))
            for intv_name, (bucket, counter) in INTV_TO_METRIC.items():
                if intv_name not in sim.interventions:
                    continue
                frames.append(_result_rows(
                    sim.interventions[intv_name].results[counter], bucket,
                ).assign(scenario=scen_label, sim=sim_idx))
    return pd.concat(frames, ignore_index=True)


def save_csvs(long, resfolder='results'):
    """Two plot-ready CSVs from the long-format ensemble DataFrame."""
    os.makedirs(resfolder, exist_ok=True)

    q = {'value': 'median',
         'low':   lambda s: s.quantile(0.10),
         'high':  lambda s: s.quantile(0.90)}

    ts = (long[long.metric.isin(TS_METRICS)]
          .groupby(['scenario', 'metric', 'year'])['value'].agg(**q)
          .reset_index())
    ts.to_csv(f'{resfolder}/scens_timeseries.csv', index=False)

    per_sim = (long[long.year >= CUM_START_YEAR]
               .groupby(['scenario', 'sim', 'metric'])['value'].sum()
               .reset_index())
    cum = (per_sim.groupby(['scenario', 'metric'])['value'].agg(**q)
                  .reset_index())
    # Fill 0 for scenarios that didn't run a given program (e.g. txvs for
    # non-TxV scenarios); plot scripts assume every (scenario, metric) exists.
    scenarios = long['scenario'].unique()
    metrics = long['metric'].unique()
    grid = pd.MultiIndex.from_product([scenarios, metrics], names=['scenario', 'metric'])
    cum = cum.set_index(['scenario', 'metric']).reindex(grid, fill_value=0).reset_index()
    cum.to_csv(f'{resfolder}/scens_cumulative.csv', index=False)


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
