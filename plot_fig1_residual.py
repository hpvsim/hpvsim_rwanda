"""
Fig 1: residual burden of cervical cancer in Rwanda under ongoing interventions.

Plots from plot-ready CSVs produced by `run_scenarios.py --run-sim`.
"""
import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
import sciris as sc

import utils as ut


def plot_fig1(resfolder='results', outpath='figures/fig1_residual.png', poster=False):
    fs = 20 if not poster else 18
    ut.set_font(fs)

    ts_df, cum_df = ut.load_scens(resfolder)

    figsize = (12, 5) if not poster else (12, 4.5)
    fig = pl.figure(layout="tight", figsize=figsize)
    gs = fig.add_gridspec(1, 3)
    ax = fig.add_subplot(gs[:2])

    start_year, end_year = 2016, 2100

    this_dict = {
        'No interventions': 'No interventions',
        'Status quo': 'Baseline',
    }

    for slabel, scen_key in this_dict.items():
        ls = ':' if slabel == 'No interventions' else '-'
        ax = ut.plot_ts(ax, ts_df, scen_key, 'asr_cancer_incidence',
                        start_year, end_year, color='k', ls=ls, label=slabel,
                        add_bounds=True)
    ax.set_ylim(bottom=0)
    ax.set_title('ASR cervical cancer incidence, 2025-2100')

    linestyle_labels = ['Status quo', 'No interventions'] if not poster else ['Status quo', 'No vax']
    linestyle_handles = [plt.Line2D([0], [0], color='k', linestyle='-', lw=2),
                         plt.Line2D([0], [0], color='k', linestyle=':', lw=2)]
    ax.legend(linestyle_handles, linestyle_labels, title='', loc='upper right', frameon=False)

    # Cumulative cancers 2025-2100
    cum_res = {sname: ut.get_cum(cum_df, scen_key, 'cancers')[0]
               for sname, scen_key in this_dict.items()}
    for sname, scen_key in this_dict.items():
        val, lb, ub = ut.get_cum(cum_df, scen_key, 'cancers')
        print(f'{sname}: {val:.0f} ({lb:.0f}, {ub:.0f}) cancers')

    ax = fig.add_subplot(gs[2])
    bars = list(cum_res.values())
    labels = list(cum_res.keys())
    x = np.arange(len(labels))
    ax.bar(x, bars, color='k')
    ax.set_xticks(x)
    xlabels = ['No\ninterventions', 'Status\nquo'] if not poster else ['No vax', 'Status quo']
    ax.set_xticklabels(xlabels)
    ax.set_title('Cumulative cancers' if not poster else 'Cervical cancers 2025–2100')
    sc.SIticks()

    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    sc.savefig(outpath, dpi=100)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--resfolder', default='results/v2.2.6_baseline')
    parser.add_argument('--outpath', default='figures/fig1_residual.png')
    parser.add_argument('--poster', action='store_true')
    args = parser.parse_args()

    plot_fig1(resfolder=args.resfolder, outpath=args.outpath, poster=args.poster)
    print('Done.')
