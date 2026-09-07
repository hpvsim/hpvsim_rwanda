"""
Fig 2 (poster variant): screening scenarios with and without VIA triage.

Poster layout drops the HIV+ panel vs. plot_fig2_st.py.
"""
import argparse
import os

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
import sciris as sc

import utils as ut


def plot_fig2(resfolder='results', outpath='figures/poster/fig2_st.png'):
    ut.set_font(20)
    ts_df, cum_df = ut.load_scens(resfolder)

    fig = pl.figure(layout="tight", figsize=(12, 10))
    gs = fig.add_gridspec(2, 2)

    start_year, end_year = 2016, 2100
    ymax = 25
    colors = sc.vectocolor(3).tolist()

    labels = ['18% coverage', '35% coverage', '70% coverage']
    no_triage_keys = [f'S&T {c}' for c in ['18%', '35%', '70%']]
    triage_keys = [f'S&T&T {c}' for c in ['18%', '35%', '70%']]

    bar_width = 0.35
    x = np.arange(len(labels)).astype(float)
    llabels = ['18%', '35%', '70%']
    no_triage_x = x - bar_width / 2
    triage_x = x + bar_width / 2

    # ---- Top Left: cumulative cancers ----
    ax = fig.add_subplot(gs[0, 0])
    nt = [ut.get_cum(cum_df, k, 'cancers')[0] for k in no_triage_keys]
    t = [ut.get_cum(cum_df, k, 'cancers')[0] for k in triage_keys]
    ax.bar(no_triage_x, nt, width=bar_width, color=colors)
    ax.bar(triage_x, t, width=bar_width, color=colors,
           edgecolor='k', hatch='//', linewidth=1.5)
    ax.set_xticks(x); ax.set_xticklabels(llabels)
    ax.set_title('Cumulative cancers\n2025-2100'); sc.SIticks()
    ax.set_xlabel('')
    ax.legend(handles=[
        mpatches.Patch(facecolor='gray', label='No triage'),
        mpatches.Patch(facecolor='white', edgecolor='k', hatch='//', label='VIA triage'),
    ], loc='upper right', frameon=False)
    ax.set_ylim([0, 100e3])

    # ---- Top Right: cumulative ablations ----
    ax = fig.add_subplot(gs[0, 1])
    abl_nt = [ut.get_cum(cum_df, k, 'ablations')[0] for k in no_triage_keys]
    abl_t = [ut.get_cum(cum_df, k, 'ablations')[0] for k in triage_keys]
    ax.bar(no_triage_x, abl_nt, width=bar_width, color=colors)
    ax.bar(triage_x, abl_t, width=bar_width, color=colors,
           edgecolor='k', hatch='//', linewidth=1.5)
    ax.set_xticks(x); ax.set_xticklabels(llabels)
    ax.set_title('Cumulative ablations\n2025-2100'); sc.SIticks()
    ax.set_xlabel('')

    # ---- Bottom Left: time series ----
    ax = fig.add_subplot(gs[1, 0])
    ax = ut.plot_ts(ax, ts_df, 'S&T&T 18%', 'asr_cancer_incidence', start_year, end_year,
                    color='k', label='Status quo (18% + triage)')
    for cn, slabel in enumerate(labels):
        ax = ut.plot_ts(ax, ts_df, no_triage_keys[cn], 'asr_cancer_incidence',
                        start_year, end_year, color=colors[cn], label=slabel)
    for cn, sname in enumerate(['35% coverage', '70% coverage'], start=1):
        ax = ut.plot_ts(ax, ts_df, triage_keys[cn], 'asr_cancer_incidence',
                        start_year, end_year, color=colors[cn], ls='--', label='')
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('Cervical cancer incidence, 2025-2100\nScaled-up screening with/without VIA triage')

    color_handles = [mpatches.Patch(facecolor=colors[i], label=labels[i]) for i in range(3)]
    ls_handles = [plt.Line2D([0], [0], color='k', linestyle='-', lw=2),
                  plt.Line2D([0], [0], color='k', linestyle='--', lw=2)]
    legend1 = ax.legend(frameon=False, handles=color_handles, title='',
                        loc='upper right', bbox_to_anchor=(1, 0.7))
    ax.add_artist(legend1)
    ax.legend(ls_handles, ['No triage', 'VIA triage'], title='',
              loc='upper right', frameon=False)

    # ---- Bottom Right: ablations per cancer averted ----
    ax = fig.add_subplot(gs[1, 1])
    baseline_cancers = ut.get_cum(cum_df, 'S&T&T 18%', 'cancers')[0]

    def ratio(cum_abl, cum_canc):
        averted = baseline_cancers - cum_canc
        return cum_abl / averted if averted > 0 else 0

    nt_ratios = [ratio(abl_nt[i], nt[i]) for i in range(3)]
    t_ratios = [ratio(abl_t[i], t[i]) for i in range(3)]
    ax.bar(no_triage_x, nt_ratios, width=bar_width, color=colors)
    ax.bar(triage_x, t_ratios, width=bar_width, color=colors,
           edgecolor='k', hatch='//', linewidth=1.5)
    ax.set_xticks(x); ax.set_xticklabels(llabels)
    ax.set_title('Ablations per cancer averted\n2025-2100')
    ax.set_xlabel('')

    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    sc.savefig(outpath, dpi=100)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--resfolder', default='results/v2.2.6_baseline')
    parser.add_argument('--outpath', default='figures/poster/fig2_st.png')
    args = parser.parse_args()
    plot_fig2(resfolder=args.resfolder, outpath=args.outpath)
    print('Done.')
