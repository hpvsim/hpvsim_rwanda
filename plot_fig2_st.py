"""
Fig 2: screening scenarios with and without VIA triage.

Plots from plot-ready CSVs produced by `run_scenarios.py --run-sim`.
"""
import argparse
import os

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
import sciris as sc

import utils as ut


def plot_fig2(resfolder='results', outpath='figures/fig2_st.png'):
    ut.set_font(20)
    ts_df, cum_df = ut.load_scens(resfolder)

    fig = pl.figure(layout="tight", figsize=(18, 10))
    gs = fig.add_gridspec(2, 3)

    text_height = [-0.1, 1.2]
    start_year, end_year = 2016, 2100
    ymax = 25
    vc = sc.vectocolor(3).tolist()
    colors = vc

    labels = ['18% coverage', '35% coverage', '70% coverage']
    no_triage_keys = [f'S&T {c}' for c in ['18%', '35%', '70%']]
    triage_keys = [f'S&T&T {c}' for c in ['18%', '35%', '70%']]

    bar_width = 0.35
    x = np.arange(len(labels)).astype(float)
    llabels = ['18%', '35%', '70%']
    no_triage_x = x - bar_width / 2
    triage_x = x + bar_width / 2

    # ---- A: Cumulative cancers ----
    ax = fig.add_subplot(gs[0, 0])
    cancers_nt = [ut.get_cum(cum_df, k, 'cancers')[0] for k in no_triage_keys]
    cancers_t = [ut.get_cum(cum_df, k, 'cancers')[0] for k in triage_keys]
    for k in no_triage_keys + triage_keys:
        v, lb, ub = ut.get_cum(cum_df, k, 'cancers')
        print(f'{k}: {v:.0f} ({lb:.0f}, {ub:.0f}) cancers')

    ax.bar(no_triage_x, cancers_nt, width=bar_width, color=colors)
    ax.bar(triage_x, cancers_t, width=bar_width, color=colors,
           edgecolor='k', hatch='//', linewidth=1.5)
    ax.set_xticks(x); ax.set_xticklabels(llabels)
    ax.set_title('Cumulative cancers\n2025-2100'); sc.SIticks()
    ax.set_xlabel('')
    ax.legend(handles=[
        mpatches.Patch(facecolor='gray', label='No triage'),
        mpatches.Patch(facecolor='white', edgecolor='k', hatch='//', label='VIA triage'),
    ], loc='upper right', frameon=False)
    ax.set_ylim([0, 100e3])
    ax.text(*text_height, 'A', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    # ---- B: Cumulative cancers in HIV+ ----
    ax = fig.add_subplot(gs[0, 1])
    hiv_nt = [ut.get_cum(cum_df, k, 'cancers_with_hiv')[0] for k in no_triage_keys]
    hiv_t = [ut.get_cum(cum_df, k, 'cancers_with_hiv')[0] for k in triage_keys]
    ax.bar(no_triage_x, hiv_nt, width=bar_width, color=colors)
    ax.bar(triage_x, hiv_t, width=bar_width, color=colors,
           edgecolor='k', hatch='//', linewidth=1.5)
    ax.set_xticks(x); ax.set_xticklabels(llabels)
    ax.set_title('Cumulative cancers in HIV+ women\n2025-2100'); sc.SIticks()
    ax.set_xlabel('')
    ax.text(*text_height, 'B', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    # ---- C: Cumulative ablations ----
    ax = fig.add_subplot(gs[0, 2])
    abl_nt = [ut.get_cum(cum_df, k, 'ablations')[0] for k in no_triage_keys]
    abl_t = [ut.get_cum(cum_df, k, 'ablations')[0] for k in triage_keys]
    ax.bar(no_triage_x, abl_nt, width=bar_width, color=colors)
    ax.bar(triage_x, abl_t, width=bar_width, color=colors,
           edgecolor='k', hatch='//', linewidth=1.5)
    ax.set_xticks(x); ax.set_xticklabels(llabels)
    ax.set_title('Cumulative ablations\n2025-2100'); sc.SIticks()
    ax.set_xlabel('')
    ax.text(*text_height, 'C', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    # ---- D: Time series ----
    ax = fig.add_subplot(gs[1, :2])
    ax = ut.plot_ts(ax, ts_df, 'S&T&T 18%', 'asr_cancer_incidence', start_year, end_year,
                    color='k', label='Status quo (18% + triage)')
    for cn, slabel in enumerate(labels):
        ax = ut.plot_ts(ax, ts_df, no_triage_keys[cn], 'asr_cancer_incidence',
                        start_year, end_year, color=colors[cn], label=slabel)
    for cn, sname in enumerate(['35% coverage', '70% coverage'], start=1):
        ax = ut.plot_ts(ax, ts_df, triage_keys[cn], 'asr_cancer_incidence',
                        start_year, end_year, color=colors[cn], ls='--', label='')
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('ASR cervical cancer incidence, 2025-2100\nScaled-up screening with/without VIA triage')

    color_handles = [mpatches.Patch(facecolor=colors[i], label=labels[i]) for i in range(3)]
    linestyle_handles = [plt.Line2D([0], [0], color='k', linestyle='-', lw=2),
                         plt.Line2D([0], [0], color='k', linestyle='--', lw=2)]
    legend1 = ax.legend(frameon=False, handles=color_handles, title='',
                        loc='upper right', bbox_to_anchor=(1, 0.7))
    ax.add_artist(legend1)
    ax.legend(linestyle_handles, ['No triage', 'VIA triage'], title='',
              loc='upper right', frameon=False)
    ax.text(-0.05, 1.2, 'D', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    # ---- E: Ablations per cancer averted ----
    ax = fig.add_subplot(gs[1, 2])
    baseline_cancers = ut.get_cum(cum_df, 'S&T&T 18%', 'cancers')[0]

    def ratio(cum_abl, cum_canc):
        averted = baseline_cancers - cum_canc
        return cum_abl / averted if averted > 0 else 0

    nt_ratios = [ratio(abl_nt[i], cancers_nt[i]) for i in range(3)]
    t_ratios = [ratio(abl_t[i], cancers_t[i]) for i in range(3)]
    ax.bar(no_triage_x, nt_ratios, width=bar_width, color=colors)
    ax.bar(triage_x, t_ratios, width=bar_width, color=colors,
           edgecolor='k', hatch='//', linewidth=1.5)
    ax.set_xticks(x); ax.set_xticklabels(llabels)
    ax.set_title('Ablations per cancer averted\n2025-2100')
    ax.set_xlabel('')
    ax.text(*text_height, 'E', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    sc.savefig(outpath, dpi=100)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--resfolder', default='results/v2.2.6_baseline')
    parser.add_argument('--outpath', default='figures/fig2_st.png')
    args = parser.parse_args()
    plot_fig2(resfolder=args.resfolder, outpath=args.outpath)
    print('Done.')
