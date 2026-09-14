"""
Fig 5: cumulative cancers and cancers averted across all seven strategies.

Plots from plot-ready CSVs produced by `run_scenarios.py --run-sim`.
"""
import argparse
import os

import numpy as np
import pylab as pl
import sciris as sc

import utils as ut


ALL_STRATEGIES = [
    ('S&T', 'Screen+\ntreat'),
    ('S&T&T', 'Screen+\ntriage+\ntreat'),
    ('S&TxV&T&T', 'Screen+\nTxV 90/0+\ntriage+treat'),
    ('S&TxV', 'Screen+\nTxV 50/90'),
    ('HPV-Faster', 'HPV-Faster'),
    ('Mass TxV 90/0,', 'Mass TxV\n90/0'),
    ('Mass TxV 50/90,', 'Mass TxV\n50/90'),
]


def plot_fig5(resfolder='results', outpath='figures/fig5_comparison.png', poster=False):
    fs = 20 if not poster else 24
    ut.set_font(fs)
    _, cum_df, paired_df, _ = ut.load_scens(resfolder)

    figsize = (14, 12)
    fig = pl.figure(layout="tight", figsize=figsize)
    gs = fig.add_gridspec(2, 1)
    text_font = 14 if not poster else 20

    coverage_levels = ['18%', '35%', '70%']
    coverage_colors = sc.vectocolor(len(coverage_levels)).tolist()

    # ---- Top: cumulative cancers ----
    ax = fig.add_subplot(gs[0])
    bar_width = 0.28
    x_base = np.arange(len(ALL_STRATEGIES))
    offsets = [-bar_width, 0, bar_width]

    all_cum = []
    for cov_idx, cov in enumerate(coverage_levels):
        keys = [f'{strat_key} {cov}' for strat_key, _ in ALL_STRATEGIES]
        stats = [ut.get_cum(cum_df, k, 'cancers') for k in keys]
        med, lo, hi = zip(*stats)
        med, lo, hi = list(med), list(lo), list(hi)
        for k, v in zip(keys, med):
            print(f'{k}: {v:.0f} cancers')
        all_cum.extend(hi)
        bars = ax.bar(x_base + offsets[cov_idx], med, width=bar_width,
                      color=coverage_colors[cov_idx], label=cov,
                      yerr=ut.yerr(med, lo, hi), capsize=3, ecolor='0.4')
        for bar, h_bar in zip(bars, hi):
            ax.text(bar.get_x() + bar.get_width() / 2., h_bar,
                    f'{sc.sigfig(bar.get_height() / 1e3, 3)}',
                    ha='center', va='bottom', fontsize=text_font)

    ax.set_xticks(x_base)
    ax.set_xticklabels([label for _, label in ALL_STRATEGIES])
    ax.set_title('Cumulative cancers 2025-2100'); sc.SIticks()
    ax.set_ylim(bottom=0, top=max(all_cum) * 1.15)
    ax.legend(title='Coverage', loc='upper right', frameon=False, fontsize=14, ncols=3)

    # ---- Bottom: paired cancers averted per sim vs S&T&T 18% ----
    ax = fig.add_subplot(gs[1])
    all_avr_hi = []
    for cov_idx, cov in enumerate(coverage_levels):
        keys = [f'{strat_key} {cov}' for strat_key, _ in ALL_STRATEGIES]
        stats = [ut.get_paired(paired_df, k, 'cancers') for k in keys]
        med, lo, hi = zip(*stats)
        med, lo, hi = list(med), list(lo), list(hi)
        for k, v, l, h in zip(keys, med, lo, hi):
            print(f'{k}: averted {v:.0f} ({l:.0f}, {h:.0f})')
        all_avr_hi.extend(hi)
        bars = ax.bar(x_base + offsets[cov_idx], med, width=bar_width,
                      color=coverage_colors[cov_idx], label=cov,
                      yerr=ut.yerr(med, lo, hi), capsize=3, ecolor='0.4')
        for bar, h_bar in zip(bars, hi):
            if bar.get_height() > 0:
                ax.text(bar.get_x() + bar.get_width() / 2., max(h_bar, bar.get_height()),
                        f'{sc.sigfig(bar.get_height() / 1e3, 3)}',
                        ha='center', va='bottom', fontsize=text_font)

    ax.set_xticks(x_base)
    ax.set_xticklabels([label for _, label in ALL_STRATEGIES])
    ax.set_title('Cancers averted 2025-2100 (paired diff vs. S&T&T 18%)'); sc.SIticks()
    ax.axhline(0, color='0.4', lw=0.5)
    ax.set_ylim(top=max(all_avr_hi) * 1.15)

    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    sc.savefig(outpath, dpi=100)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--resfolder', default='results/v2.2.6_baseline')
    parser.add_argument('--outpath', default='figures/fig5_comparison.png')
    parser.add_argument('--poster', action='store_true')
    args = parser.parse_args()
    plot_fig5(resfolder=args.resfolder, outpath=args.outpath, poster=args.poster)
    print('Done.')
