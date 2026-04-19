"""
Cross-version comparison for hpvsim_rwanda scenario outputs.

Compares cumulative cancers and ASR elimination year across baseline folders
(e.g. v2.2.6_baseline vs a future v2.3.0_baseline).

Usage:
  python compare_baselines.py --baselines v2.2.6_baseline v2.3.0_baseline
  python compare_baselines.py --baselines v2.2.6_baseline v2.3.0_baseline \\
                              --scenarios "S&T&T 18%" "S&TxV 70%" "HPV-Faster 70%"
"""
import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sciris as sc

import utils as ut


DEFAULT_SCENARIOS = [
    'Baseline', 'S&T&T 18%', 'S&T 70%', 'S&TxV 70%',
    'HPV-Faster 70%', 'Mass TxV 50/90, 70%',
]


def elim_year(ts_df, scenario, metric='asr_cancer_incidence',
              threshold=4, search_from=2020):
    years, values, _, _ = ut.get_ts(ts_df, scenario, metric)
    smoothed = np.convolve(values, np.ones(5), 'valid') / 5
    years_s = years[4:]
    mask = years_s >= search_from
    below = np.where((smoothed < threshold) & mask)[0]
    if len(below) == 0:
        return None
    return float(years_s[below[0]])


def compare(baselines, scenarios, resroot='results',
            outpath='figures/compare_baselines.png'):
    cum_rows = []
    ts_rows = []
    elim_rows = []

    for base in baselines:
        resfolder = f'{resroot}/{base}'
        ts_df, cum_df = ut.load_scens(resfolder)
        for scen in scenarios:
            v, lb, ub = ut.get_cum(cum_df, scen, 'cancers')
            cum_rows.append({'baseline': base, 'scenario': scen,
                             'cancers': v, 'low': lb, 'high': ub,
                             'elim_year': elim_year(ts_df, scen)})
        for scen in scenarios:
            yrs, vals, lo, hi = ut.get_ts(ts_df, scen, 'asr_cancer_incidence')
            for yi, yr in enumerate(yrs):
                ts_rows.append({'baseline': base, 'scenario': scen,
                                'year': float(yr), 'asr': float(vals[yi]),
                                'low': float(lo[yi]), 'high': float(hi[yi])})

    cum_out = pd.DataFrame(cum_rows)
    ts_out = pd.DataFrame(ts_rows)

    # --- Plot: one panel per scenario, ASR over time overlaid across baselines ---
    ut.set_font(14)
    n = len(scenarios)
    ncols = min(3, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows), squeeze=False)

    base_colors = sc.gridcolors(len(baselines))
    for si, scen in enumerate(scenarios):
        ax = axes[si // ncols, si % ncols]
        for bi, base in enumerate(baselines):
            sub = ts_out[(ts_out.baseline == base) & (ts_out.scenario == scen)].sort_values('year')
            sub = sub[(sub.year >= 2016) & (sub.year <= 2100)]
            smooth = np.convolve(sub.asr.values, np.ones(5), 'valid') / 5
            ax.plot(sub.year.values[4:], smooth, color=base_colors[bi], label=base)
        ax.axhline(4, color='k', ls='--', lw=0.5)
        ax.set_title(scen)
        ax.set_ylim(bottom=0)
        if si == 0:
            ax.legend(frameon=False)

    for si in range(n, nrows * ncols):
        axes[si // ncols, si % ncols].set_visible(False)

    fig.suptitle('ASR cancer incidence by baseline version')
    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.savefig(outpath, dpi=100)

    print('\n=== Cumulative cancers 2025-2100 (cancers | elim year) ===')
    pivot = cum_out.pivot_table(index='scenario', columns='baseline',
                                values=['cancers', 'elim_year'])
    print(pivot.to_string())

    return cum_out, ts_out


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--baselines', nargs='+', required=True,
                        help='Baseline folder names under results/ (e.g. v2.2.6_baseline)')
    parser.add_argument('--scenarios', nargs='+', default=DEFAULT_SCENARIOS)
    parser.add_argument('--resroot', default='results')
    parser.add_argument('--outpath', default='figures/compare_baselines.png')
    args = parser.parse_args()

    compare(args.baselines, args.scenarios, resroot=args.resroot, outpath=args.outpath)
    print('\nDone.')
