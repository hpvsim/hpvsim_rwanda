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

from interventions import make_st
from run_scenarios import (
    run_sims, process_msim, save_csvs, CUM_START_YEAR,
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
        msim = run_sims(scenarios=scenarios, end=args.end)
        msim_dict = process_msim(msim, scenarios)
        sc.saveobj(obj_path, msim_dict)
    else:
        msim_dict = sc.loadobj(obj_path)

    save_csvs(msim_dict, resfolder=args.resfolder, cum_start_year=CUM_START_YEAR)
    print(f'Saved scens_*.csv to {args.resfolder}/ (cum accounting from {CUM_START_YEAR})')
    T.toc('Done')
