"""
Fig 4: one-time mass campaign strategies.

Plots from plot-ready CSVs produced by `run_scenarios.py --run-sim`.
"""
import argparse
import os

import numpy as np
import pylab as pl
import sciris as sc
from matplotlib.patches import Patch

import utils as ut


def plot_fig4(resfolder='results', outpath='figures/fig4_campaigns.png'):
    ut.set_font(20)
    ts_df, cum_df = ut.load_scens(resfolder)

    fig = pl.figure(layout="tight", figsize=(18, 10))
    gs = fig.add_gridspec(2, 3)

    text_height = [-0.1, 1.2]
    start_year, end_year = 2016, 2100
    ymax = 25
    colors = sc.vectocolor(3).tolist()

    coverage_levels = ['18%', '35%', '70%']
    strategies = {
        'HPV-Faster': 'HPV-Faster',
        'Mass TxV 90/0': 'Mass TxV 90/0,',
        'Mass TxV 50/90': 'Mass TxV 50/90,',
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

    # ---- C: Combined resource use (ablations + therapeutics + vaccinations) ----
    ax = fig.add_subplot(gs[0, 2])
    for idx, (strat_label, strat_key) in enumerate(strategies.items()):
        ablations, therapeutics, vaccinations = [], [], []
        for cov in coverage_levels:
            scen_key = f'{strat_key} {cov}'
            abl = cum(scen_key, 'ablations')
            txv = cum(scen_key, 'txvs')
            vax = cum(scen_key, 'vaccinations')
            print(f'{scen_key}: {abl:.0f} ablations, {txv:.0f} therapeutics, {vax:.0f} vaccinations')
            ablations.append(abl); therapeutics.append(txv); vaccinations.append(vax)
        ax.bar(x + offsets[idx], ablations, width=bar_width, color=colors[idx])
        ax.bar(x + offsets[idx], therapeutics, width=bar_width, bottom=ablations,
               color='none', edgecolor=colors[idx], hatch='//', linewidth=1.5)
        bottom_vals = [a + t for a, t in zip(ablations, therapeutics)]
        ax.bar(x + offsets[idx], vaccinations, width=bar_width, bottom=bottom_vals,
               color='none', edgecolor=colors[idx], hatch='\\\\', linewidth=1.5)
    ax.set_xticks(x); ax.set_xticklabels(coverage_levels)
    ax.set_title('Cumulative products\n2025-2100'); sc.SIticks()
    ax.set_xlabel('')
    ax.legend(handles=[
        Patch(facecolor='gray', label='Ablations'),
        Patch(facecolor='white', edgecolor='k', hatch='//', label='Therapeutics'),
        Patch(facecolor='white', edgecolor='k', hatch='\\\\', label='Vaccinations'),
    ], loc='upper left', frameon=False, fontsize=16)
    ax.set_ylim([0, 4e6])
    ax.text(*text_height, 'C', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    # ---- D: Time series at 70% coverage ----
    ax = fig.add_subplot(gs[1, :2])
    ax = ut.plot_ts(ax, ts_df, 'S&T&T 18%', 'asr_cancer_incidence', start_year, end_year,
                    color='k', label='Status quo (S&T&T 18%)')
    for strat_idx, (strat_label, strat_key) in enumerate(strategies.items()):
        scen_key = f'{strat_key} 70%'
        ax = ut.plot_ts(ax, ts_df, scen_key, 'asr_cancer_incidence', start_year, end_year,
                        color=colors[strat_idx], label=f'{strat_label} 70%')
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('ASR cervical cancer incidence, 2025-2100\nOne-time mass campaigns at 70% coverage')
    ax.legend(loc='upper right', frameon=False)
    ax.text(-0.05, 1.2, 'D', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    # ---- E: Interventions per cancer averted ----
    ax = fig.add_subplot(gs[1, 2])
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
            print(f'{scen_key}: {ratio:.1f} products per cancer averted')
        ax.bar(x + offsets[idx], nnt_per_cancer, width=bar_width, color=colors[idx], label=strat_label)
    ax.set_xticks(x); ax.set_xticklabels(coverage_levels)
    ax.set_title('Interventions per cancer averted\n2025-2100')
    ax.set_xlabel('')
    ax.text(*text_height, 'E', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    sc.savefig(outpath, dpi=100)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--resfolder', default='results/v2.2.6_baseline')
    parser.add_argument('--outpath', default='figures/fig4_campaigns.png')
    args = parser.parse_args()
    plot_fig4(resfolder=args.resfolder, outpath=args.outpath)
    print('Done.')
