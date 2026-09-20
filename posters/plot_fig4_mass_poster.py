"""
Fig 4 (poster variant): one-time mass campaign strategies.

Poster layout drops the HIV+ panel vs. plot_fig4_mass.py.
"""
import argparse
import os

import numpy as np
import pylab as pl
import sciris as sc
from matplotlib.patches import Patch

import utils as ut


def plot_fig4(resfolder='results', outpath='figures/poster/fig4_campaigns.png'):
    ut.set_font(20)
    ts_df, cum_df = ut.load_scens(resfolder)

    fig = pl.figure(layout="tight", figsize=(12, 8))
    gs = fig.add_gridspec(2, 2)

    start_year, end_year = 2016, 2100
    ymax = 25
    vc = sc.vectocolor(3).tolist()
    vc2 = sc.vectocolor(4, cmap='magma').tolist()
    colors = [vc[0], vc2[1], vc2[2]]

    coverage_levels = ['18%', '35%', '70%']
    strategies = {
        'HPV-Faster': 'HPV-Faster',
        'Mass TxV virus': 'Mass TxV 90/0,',
        'Mass TxV lesions': 'Mass TxV 50/90,',
    }
    bar_width = 0.25
    x = np.arange(len(coverage_levels))
    offsets = [-bar_width, 0, bar_width]

    def cum(scen, metric):
        return ut.get_cum(cum_df, scen, metric)[0]

    # ---- Top Left: cumulative cancers ----
    ax = fig.add_subplot(gs[0, 0])
    for idx, (strat_label, strat_key) in enumerate(strategies.items()):
        cum_res = [cum(f'{strat_key} {cov}', 'cancers') for cov in coverage_levels]
        ax.bar(x + offsets[idx], cum_res, width=bar_width, color=colors[idx], label=strat_label)
    ax.set_xticks(x); ax.set_xticklabels(coverage_levels)
    ax.set_title('Cumulative cancers 2025-2100'); sc.SIticks()
    ax.set_xlabel('')
    ax.legend(loc='upper right', frameon=False, fontsize=16)
    ax.set_ylim([0, 100e3])

    # ---- Top Right: combined resource use ----
    ax = fig.add_subplot(gs[0, 1])
    for idx, (strat_label, strat_key) in enumerate(strategies.items()):
        ablations = [cum(f'{strat_key} {cov}', 'ablations') for cov in coverage_levels]
        therapeutics = [cum(f'{strat_key} {cov}', 'txvs') for cov in coverage_levels]
        vaccinations = [cum(f'{strat_key} {cov}', 'vaccinations') for cov in coverage_levels]
        ax.bar(x + offsets[idx], ablations, width=bar_width, color=colors[idx])
        ax.bar(x + offsets[idx], therapeutics, width=bar_width, bottom=ablations,
               color='none', edgecolor=colors[idx], hatch='//', linewidth=1.5)
        bottom_vals = [a + t for a, t in zip(ablations, therapeutics)]
        ax.bar(x + offsets[idx], vaccinations, width=bar_width, bottom=bottom_vals,
               color='none', edgecolor=colors[idx], hatch='\\\\', linewidth=1.5)
    ax.set_xticks(x); ax.set_xticklabels(coverage_levels)
    ax.set_title('Cumulative products 2025-2100'); sc.SIticks()
    ax.set_xlabel('')
    ax.legend(handles=[
        Patch(facecolor='gray', label='Ablations'),
        Patch(facecolor='white', edgecolor='k', hatch='//', label='Therapeutics'),
        Patch(facecolor='white', edgecolor='k', hatch='\\\\', label='Vaccinations'),
    ], loc='upper left', frameon=False, fontsize=16)
    ax.set_ylim([0, 4e6])

    # ---- Bottom Left: time series at 70% ----
    ax = fig.add_subplot(gs[1, 0])
    ax = ut.plot_ts(ax, ts_df, 'S&T&T 18%', 'asr_cancer_incidence', start_year, end_year,
                    color='k', label='Status quo')
    for strat_idx, (strat_label, strat_key) in enumerate(strategies.items()):
        ax = ut.plot_ts(ax, ts_df, f'{strat_key} 70%', 'asr_cancer_incidence',
                        start_year, end_year, color=colors[strat_idx],
                        label=f'{strat_label} 70%')
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('Cervical cancer incidence, 2025-2100')
    ax.legend(loc='upper right', frameon=False, fontsize=16)

    # ---- Bottom Right: products per cancer averted ----
    ax = fig.add_subplot(gs[1, 1])
    baseline_cancers = cum('S&T&T 18%', 'cancers')
    for idx, (strat_label, strat_key) in enumerate(strategies.items()):
        nnt_per_cancer = []
        for cov in coverage_levels:
            scen_key = f'{strat_key} {cov}'
            cancers = cum(scen_key, 'cancers')
            averted = baseline_cancers - cancers
            total = cum(scen_key, 'txvs') + cum(scen_key, 'vaccinations')
            ratio = total / averted if averted > 0 else 0
            nnt_per_cancer.append(ratio)
        ax.bar(x + offsets[idx], nnt_per_cancer, width=bar_width, color=colors[idx], label=strat_label)
    ax.set_xticks(x); ax.set_xticklabels(coverage_levels)
    ax.set_title('Products / cancer averted 2025-2100')
    ax.set_xlabel('')

    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    sc.savefig(outpath, dpi=100)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--resfolder', default='results/v2.2.6_baseline')
    parser.add_argument('--outpath', default='figures/poster/fig4_campaigns.png')
    args = parser.parse_args()
    plot_fig4(resfolder=args.resfolder, outpath=args.outpath)
    print('Done.')
