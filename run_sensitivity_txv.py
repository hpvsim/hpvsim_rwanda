"""
Sensitivity sweeps for the reviewer response (paper revision).

Currently defines Sweep A (R1.3 + R2.1): vary the TxV introduction year
across [2030, 2035, 2040, 2045, 2050] for all four TxV-carrying arms
(S&TxV, S&TxV&T&T, Mass TxV 50/90, Mass TxV 90/0) plus the S&T&T 18%
baseline required by the paired-diff pipeline.

Two modes, mirroring run_scenarios.py:
  python run_sensitivity.py --run-sim   # run msim + save plot-ready CSVs (VM)
  python run_sensitivity.py             # re-extract CSVs from existing .obj
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

from interventions import make_st, make_mv_intvs
from run_scenarios import (
    run_sims, process_msim, save_csvs, CUM_START_YEAR,
)


TXV_INTRO_YEARS = [2030, 2035, 2040, 2045, 2050]
TXV_COV = 0.70


def make_txv_year_scenarios(end_year=2100, txv_years=TXV_INTRO_YEARS, cov=TXV_COV):
    """Sweep A: TxV intro year across the four TxV-carrying arms.

    Includes S&T&T 18% so the paired-diff pipeline in save_csvs can find
    its BASELINE_SCEN row.
    """
    scendict = dict()
    scendict['S&T&T 18%'] = make_st(future_screen_cov=0.18, end_year=end_year)

    for yr in txv_years:
        scendict[f'S&TxV&T&T {cov*100:.0f}% (TxV {yr})'] = make_st(
            future_screen_cov=cov, txv_pars='precin', txv=True,
            txv_start_year=yr, end_year=end_year,
        )
        scendict[f'S&TxV {cov*100:.0f}% (TxV {yr})'] = make_st(
            future_screen_cov=cov, txv_pars='cin', txv=True,
            txv_start_year=yr, end_year=end_year,
        )
        scendict[f'Mass TxV 90/0 {cov*100:.0f}% (TxV {yr})'] = make_mv_intvs(
            txv_pars='precin', campaign_coverage=cov,
            intro_year=yr, end_year=end_year,
        )
        scendict[f'Mass TxV 50/90 {cov*100:.0f}% (TxV {yr})'] = make_mv_intvs(
            txv_pars='cin', campaign_coverage=cov,
            intro_year=yr, end_year=end_year,
        )

    return scendict


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--run-sim', action='store_true',
                        help='Run scenarios on the VM (heavy); otherwise only re-extract CSVs')
    parser.add_argument('--end', type=int, default=2100)
    parser.add_argument('--resfolder', default='results/sens_txv_year')
    args = parser.parse_args()

    T = sc.timer()
    scenarios = make_txv_year_scenarios(end_year=args.end)

    os.makedirs(args.resfolder, exist_ok=True)
    obj_path = f'{args.resfolder}/sens_txv_year.obj'

    if args.run_sim:
        msim = run_sims(scenarios=scenarios, end=args.end)
        msim_dict = process_msim(msim, scenarios)
        sc.saveobj(obj_path, msim_dict)
    else:
        msim_dict = sc.loadobj(obj_path)

    save_csvs(msim_dict, resfolder=args.resfolder, cum_start_year=CUM_START_YEAR)
    print(f'Saved scens_*.csv to {args.resfolder}/ (cum accounting from {CUM_START_YEAR})')
    T.toc('Done')
