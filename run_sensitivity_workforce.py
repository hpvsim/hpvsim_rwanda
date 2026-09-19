"""
Workforce capacity sensitivity for R2.5.

Reruns the S&T-family scenarios (S&T&T, S&T, S&TxV&T&T, S&TxV at
18/35/70% coverage) under three workforce caps on treat_num.max_capacity:
10, 15, 30 agents/timestep. Also includes the no-cap S&T&T 18% baseline
(named exactly 'S&T&T 18%' so the paired-diff pipeline finds it) and
'No interventions' for context.

Baseline reference (S&T&T 18% at 2028, no cap): ~34,700 ablations/year
in population-scale terms. The three caps roughly correspond to 1.5x,
2.5x, 5x that baseline.

Campaign scenarios (HPV-Faster, Mass TxV) are excluded: those assume the
workforce constraint is relaxed during the one-off campaign year.

Memory strategy: each parallel worker runs its sim, extracts long-format
rows, and returns only the DataFrame; the parent process never holds the
sim objects. This keeps peak parent-side memory bounded by n_sims x
per-sim-rows (KB) rather than n_sims x per-sim-state (MB) as the
MultiSim path in run_scenarios.run_sims does.
"""

import argparse
import os

os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

import pandas as pd
import sciris as sc

import run_sim as rs
from interventions import make_st
from run_scenarios import (
    _top_pars, _result_rows, save_csvs, CUM_START_YEAR,
    TS_METRICS, CUM_METRICS_BOUNDED, INTV_TO_METRIC, _V3_ALIAS, n_reps,
)


CAPS = [10, 15, 30]  # agents/ti; approx 1.5x, 2.5x, 5x baseline 34.7K/yr


def make_workforce_scenarios(end_year=2100, caps=CAPS):
    """S&T-family scenarios at three coverage levels under three workforce caps."""
    scendict = dict()
    scendict['No interventions'] = []
    scendict['S&T&T 18%'] = make_st(future_screen_cov=0.18, end_year=end_year)

    for cap in caps:
        cap_lbl = f'cap{cap}'
        for cov in [0.18, 0.35, 0.70]:
            cov_pct = int(cov * 100)
            scendict[f'S&T&T {cov_pct}% ({cap_lbl})'] = make_st(
                future_screen_cov=cov, treat_capacity=cap, end_year=end_year,
            )
            scendict[f'S&T {cov_pct}% ({cap_lbl})'] = make_st(
                future_screen_cov=cov, tx_assigner_csv='tx_assigner_no_triage',
                treat_capacity=cap, end_year=end_year,
            )
            scendict[f'S&TxV&T&T {cov_pct}% ({cap_lbl})'] = make_st(
                future_screen_cov=cov, txv_pars='precin', txv=True,
                treat_capacity=cap, end_year=end_year,
            )
            scendict[f'S&TxV {cov_pct}% ({cap_lbl})'] = make_st(
                future_screen_cov=cov, txv_pars='cin', txv=True,
                treat_capacity=cap, end_year=end_year,
            )
    return scendict


def _run_and_extract(scenario_name, add_vax, interventions, calib_pars,
                     end, sim_idx):
    """Run a single sim, extract per-year rows, return DataFrame only."""
    sim = rs.make_sim(
        add_st=False,
        add_vax=add_vax,
        interventions=interventions,
        stop=end,
        calib_pars=dict(calib_pars),
        use_calib=False,
    )
    sim.label = scenario_name
    sim.run(verbose=-1)

    frames = []
    for metric in TS_METRICS + CUM_METRICS_BOUNDED:
        frames.append(_result_rows(
            sim.results.all_hpv[_V3_ALIAS.get(metric, metric)], metric,
        ).assign(scenario=scenario_name, sim=sim_idx))
    for intv_name, (bucket, counter) in INTV_TO_METRIC.items():
        if intv_name not in sim.interventions:
            continue
        frames.append(_result_rows(
            sim.interventions[intv_name].results[counter], bucket,
        ).assign(scenario=scenario_name, sim=sim_idx))
    return pd.concat(frames, ignore_index=True)


def run_and_extract_scenarios(scenarios, end=2100, top_pars=None):
    """Run all (scenario, rep) combinations in parallel, return long DataFrame."""
    if top_pars is None:
        top_pars = _top_pars(n_reps)
    iterkwargs = []
    for name, interventions in scenarios.items():
        add_vax = name != 'No interventions'
        for sim_idx, trial_pars in enumerate(top_pars):
            iterkwargs.append(dict(
                scenario_name=name,
                add_vax=add_vax,
                interventions=interventions,
                calib_pars=trial_pars,
                end=end,
                sim_idx=sim_idx,
            ))
    frames = sc.parallelize(_run_and_extract, iterkwargs=iterkwargs)
    return pd.concat(frames, ignore_index=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-sim', action='store_true',
                        help='Run scenarios on the VM (heavy); otherwise only re-extract CSVs')
    parser.add_argument('--end', type=int, default=2100)
    parser.add_argument('--resfolder', default='results/sens_workforce')
    args = parser.parse_args()

    T = sc.timer()
    scenarios = make_workforce_scenarios(end_year=args.end)
    os.makedirs(args.resfolder, exist_ok=True)
    obj_path = f'{args.resfolder}/sens_workforce.obj'

    if args.run_sim:
        long_df = run_and_extract_scenarios(scenarios=scenarios, end=args.end)
        sc.saveobj(obj_path, long_df)
    else:
        long_df = sc.loadobj(obj_path)

    save_csvs(long_df, resfolder=args.resfolder, cum_start_year=CUM_START_YEAR)
    print(f'Saved scens_*.csv to {args.resfolder}/ (cum accounting from {CUM_START_YEAR})')
    T.toc('Done')
