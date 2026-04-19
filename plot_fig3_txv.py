"""
Fig 3: therapeutic-enhanced screening strategies.

Plots from plot-ready CSVs produced by `run_scenarios.py --run-sim`.
"""
import argparse
import os

import numpy as np
import pylab as pl
import sciris as sc
from matplotlib.patches import Patch

import utils as ut


def plot_fig3(resfolder='results', outpath='figures/fig3_txv.png'):
    ut.set_font(20)
    ts_df, cum_df = ut.load_scens(resfolder)

    fig = pl.figure(layout="tight", figsize=(18, 10))
    gs = fig.add_gridspec(2, 3)

    text_height = [-0.1, 1.2]
    start_year, end_year = 2016, 2100
    ymax = 25
    vc = sc.vectocolor(3).tolist()
    vc2 = sc.vectocolor(4, cmap='magma').tolist()
    colors = [vc[0], vc2[1], vc2[2]]

    coverage_levels = ['18%', '35%', '70%']
    strategies = {
        'Status quo': 'S&T&T',
        'TxV 90/0': 'S&TxV&T&T',
        'TxV 50/90': 'S&TxV',
    }
    bar_width = 0.25
    x = np.arange(len(coverage_levels))
    offsets = [-bar_width, 0, bar_width]

    def cum(scen, metric):
        return ut.get_cum(cum_df, scen, metric)[0]

    # ---- A: Cumulative cancers ----
    ax = fig.add_subplot(gs[0, 0])
    for idx, (strat_label, strat_key) in enumerate(strategies.items()):
        cum_res = []
        for cov in coverage_levels:
            scen_key = f'{strat_key} {cov}'
            v, lb, ub = ut.get_cum(cum_df, scen_key, 'cancers')
            print(f'{scen_key}: {v:.0f} ({lb:.0f}, {ub:.0f}) cancers')
            cum_res.append(v)
        ax.bar(x + offsets[idx], cum_res, width=bar_width, color=colors[idx], label=strat_label)
    ax.set_xticks(x); ax.set_xticklabels(coverage_levels)
    ax.set_title('Cumulative cancers\n2025-2100'); sc.SIticks()
    ax.set_xlabel('')
    ax.legend(loc='upper right', frameon=False, fontsize=16)
    ax.set_ylim([0, 100e3])
    ax.text(*text_height, 'A', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    # ---- B: Cumulative cancers in HIV+ ----
    ax = fig.add_subplot(gs[0, 1])
    for idx, (strat_label, strat_key) in enumerate(strategies.items()):
        hiv_res = [cum(f'{strat_key} {cov}', 'cancers_with_hiv') for cov in coverage_levels]
        ax.bar(x + offsets[idx], hiv_res, width=bar_width, color=colors[idx], label=strat_label)
    ax.set_xticks(x); ax.set_xticklabels(coverage_levels)
    ax.set_title('Cumulative cancers in HIV+ women\n2025-2100'); sc.SIticks()
    ax.set_xlabel('')
    ax.text(*text_height, 'B', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    # ---- C: Combined resource use ----
    ax = fig.add_subplot(gs[0, 2])
    for idx, (strat_label, strat_key) in enumerate(strategies.items()):
        ablations = []
        therapeutics = []
        for cov in coverage_levels:
            scen_key = f'{strat_key} {cov}'
            abl = cum(scen_key, 'ablations')
            txv = cum(scen_key, 'txvs')
            print(f'{scen_key}: {abl:.0f} ablations, {txv:.0f} therapeutics')
            ablations.append(abl); therapeutics.append(txv)
        ax.bar(x + offsets[idx], ablations, width=bar_width, color=colors[idx])
        ax.bar(x + offsets[idx], therapeutics, width=bar_width, bottom=ablations,
               color='none', edgecolor=colors[idx], hatch='//', linewidth=1.5)
    ax.set_xticks(x); ax.set_xticklabels(coverage_levels)
    ax.set_title('Cumulative treatments\n2025-2100'); sc.SIticks()
    ax.set_xlabel('')
    ax.legend(handles=[
        Patch(facecolor='gray', label='Ablations'),
        Patch(facecolor='white', edgecolor='k', hatch='//', label='Therapeutics'),
    ], loc='upper left', frameon=False, fontsize=16)
    ax.text(*text_height, 'C', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    # ---- D: Time series ----
    ax = fig.add_subplot(gs[1, :2])
    ax = ut.plot_ts(ax, ts_df, 'S&T&T 18%', 'asr_cancer_incidence', start_year, end_year,
                    color='k', label='Status quo (S&T&T 18%)')
    for strat_idx, (strat_label, strat_key) in enumerate(strategies.items()):
        scen_key = f'{strat_key} 70%'
        if scen_key == 'S&T&T 18%':
            continue
        ax = ut.plot_ts(ax, ts_df, scen_key, 'asr_cancer_incidence', start_year, end_year,
                        color=colors[strat_idx], label='')
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('ASR cervical cancer incidence, 2025-2100\nComparison of screening strategies')
    legend_labels = ['Status quo', 'SOC algorithm, 70% coverage',
                     '+TxV 90/0, 70% coverage', '+TxV 50/90, 70% coverage']
    thesecolors = ['k'] + colors
    strat_handles = [Patch(facecolor=thesecolors[i], label=legend_labels[i]) for i in range(4)]
    ax.legend(handles=strat_handles, loc='upper right', bbox_to_anchor=(1, 0.8), frameon=False)
    ax.text(-0.05, 1.2, 'D', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    # ---- E: Total treatments per cancer averted ----
    ax = fig.add_subplot(gs[1, 2])
    baseline_cancers = cum('S&T&T 18%', 'cancers')
    for idx, (strat_label, strat_key) in enumerate(strategies.items()):
        treatments_per_cancer = []
        for cov in coverage_levels:
            scen_key = f'{strat_key} {cov}'
            if scen_key == 'S&T&T 18%':
                treatments_per_cancer.append(0)
                continue
            cancers = cum(scen_key, 'cancers')
            averted = baseline_cancers - cancers
            total = cum(scen_key, 'ablations') + cum(scen_key, 'txvs')
            ratio = total / averted if averted > 0 else 0
            treatments_per_cancer.append(ratio)
            print(f'{scen_key}: {ratio:.1f} treatments per cancer averted')

        if strat_key == 'S&T&T':
            x_plot = x[1:] + offsets[idx]
            values_plot = treatments_per_cancer[1:]
        else:
            x_plot = x + offsets[idx]
            values_plot = treatments_per_cancer
        ax.bar(x_plot, values_plot, width=bar_width, color=colors[idx], label=strat_label)

    ax.set_xticks(x); ax.set_xticklabels(coverage_levels)
    ax.set_title('Total treatments per cancer averted\n2025-2100')
    ax.set_xlabel('')
    ax.text(*text_height, 'E', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    sc.savefig(outpath, dpi=100)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--resfolder', default='results/v2.2.6_baseline')
    parser.add_argument('--outpath', default='figures/fig3_txv.png')
    args = parser.parse_args()
    plot_fig3(resfolder=args.resfolder, outpath=args.outpath)
    print('Done.')
