"""
Plot 1 for infant vaccination scenarios
"""


import pylab as pl
import sciris as sc
import numpy as np


def plot_single(ax, mres, to_plot, si, ei, color, ls='-', label=None, smooth=True):
    years = mres.year[si:ei]
    best = mres[to_plot][si:ei]
    low = mres[to_plot].low[si:ei]
    high = mres[to_plot].high[si:ei]

    if smooth:
        best = np.convolve(list(best), np.ones(5), "valid")/5
        low = np.convolve(list(low), np.ones(5), "valid")/5
        high = np.convolve(list(high), np.ones(5), "valid")/5
        years = years[4:]

    ax.plot(years, best, color=color, label=label, ls=ls)
    ax.fill_between(years, low, high, alpha=0.1, color=color)
    return ax


def plot_fig1():

    fig, ax = pl.subplots(1, 1, figsize=(10, 6))
    resname = 'asr_cancer_incidence'

    # What to plot
    start_year = 2016
    end_year = 2100

    msim_dict = sc.loadobj('results/st_scens.obj')
    vc = sc.vectocolor(9)
    colors = ['k'] + [vc[0], vc[0], vc[1], vc[1], vc[2], vc[2], vc[3], vc[4], vc[5], vc[6], vc[7], vc[8]]
    mbase = msim_dict['Baseline']
    si = sc.findinds(mbase.year, start_year)[0]
    ei = sc.findinds(mbase.year, end_year)[0]

    cn = 0
    for slabel, mres in msim_dict.items():
        if '50/90' in slabel: ls = '--'
        else: ls = '-'
        ax = plot_single(ax, mres, resname, si, ei, color=colors[cn], ls=ls, label=slabel)
        cn += 1

    ax.set_ylim(bottom=0)  #, top=23)
    ax.legend(loc="lower left")
    if resname == 'asr_cancer_incidence': ax.axhline(y=4, color='k', ls='--')

    fig.tight_layout()
    fig_name = 'figures/fig1_vx_scens_all.png'
    sc.savefig(fig_name, dpi=100)

    return


# %% Run as a script
if __name__ == '__main__':

    # Load scenarios and construct figure
    plot_fig1()




