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
    _, cum_df = ut.load_scens(resfolder)

    figsize = (14, 12)
    fig = pl.figure(layout="tight", figsize=figsize)
    gs = fig.add_gridspec(2, 1)
    text_font = 14 if not poster else 20

    coverage_levels = ['18%', '35%', '70%']
    coverage_colors = sc.vectocolor(len(coverage_levels)).tolist()

    def cum(scen):
        return ut.get_cum(cum_df, scen, 'cancers')[0]

    baseline_cancers = cum('S&T&T 18%')

    # ---- Top: cumulative cancers ----
    ax = fig.add_subplot(gs[0])
    bar_width = 0.28
    x_base = np.arange(len(ALL_STRATEGIES))
    offsets = [-bar_width, 0, bar_width]

    all_cum = []
    for cov_idx, cov in enumerate(coverage_levels):
        cum_cancers = [cum(f'{strat_key} {cov}') for strat_key, _ in ALL_STRATEGIES]
        for strat_key, v in zip([s for s, _ in ALL_STRATEGIES], cum_cancers):
            print(f'{strat_key} {cov}: {v:.0f} cancers')
        all_cum.extend(cum_cancers)
        bars = ax.bar(x_base + offsets[cov_idx], cum_cancers, width=bar_width,
                      color=coverage_colors[cov_idx], label=cov)
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2., height,
                    f'{sc.sigfig(height / 1e3, 3)}',
                    ha='center', va='bottom', fontsize=text_font)

    ax.set_xticks(x_base)
    ax.set_xticklabels([label for _, label in ALL_STRATEGIES])
    ax.set_title('Cumulative cancers 2025-2100'); sc.SIticks()
    # 15% headroom so the value labels above each bar aren't cropped.
    ax.set_ylim(bottom=0, top=max(all_cum) * 1.15)
    ax.legend(title='Coverage', loc='upper right', frameon=False, fontsize=14, ncols=3)

    # ---- Bottom: cancers averted ----
    ax = fig.add_subplot(gs[1])
    all_avr = []
    for cov_idx, cov in enumerate(coverage_levels):
        averted = [max(baseline_cancers - cum(f'{strat_key} {cov}'), 0)
                   for strat_key, _ in ALL_STRATEGIES]
        for (strat_key, _), v in zip(ALL_STRATEGIES, averted):
            print(f'{strat_key} {cov}: {v:.0f} cancers averted')
        all_avr.extend(averted)
        bars = ax.bar(x_base + offsets[cov_idx], averted, width=bar_width,
                      color=coverage_colors[cov_idx], label=cov)
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax.text(bar.get_x() + bar.get_width() / 2., height,
                        f'{sc.sigfig(height / 1e3, 3)}',
                        ha='center', va='bottom', fontsize=text_font)

    ax.set_xticks(x_base)
    ax.set_xticklabels([label for _, label in ALL_STRATEGIES])
    ax.set_title('Cancers averted 2025-2100 (vs. S&T&T 18%)'); sc.SIticks()
    ax.set_ylim(bottom=0, top=max(all_avr) * 1.15)

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
