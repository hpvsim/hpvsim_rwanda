"""
Plot 1 for infant vaccination scenarios
"""


import pylab as pl
import sciris as sc
import numpy as np
import utils as ut


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
    # Add horizontal line at 4

    ax.axhline(4, color='k', ls='--', lw=0.5)
    return ax


def plot_fig1():
    ut.set_font(16)
    fig, axes = pl.subplots(2, 2, figsize=(16, 10))
    axes = axes.ravel()
    resname = 'asr_cancer_incidence'

    # What to plot
    start_year = 2016
    end_year = 2100
    ymax = 25

    msim_dict = sc.loadobj('results/st_scens.obj')
    mbase = msim_dict['Baseline']
    si = sc.findinds(mbase.year, start_year)[0]
    ei = sc.findinds(mbase.year, end_year)[0]

    # Virus-clearing TxV
    ax = axes[0]
    this_dict = {
        'No interventions': msim_dict['No interventions'],
        'Status quo': msim_dict['Baseline'],
        '10% screening': msim_dict['TxV 90/50 with 10% screening'],
        '50% screening': msim_dict['TxV 90/50 with 50% screening'],
        '90% screening': msim_dict['TxV 90/50 with 90% screening'],
    }
    vc = sc.vectocolor(3).tolist()
    colors = ['k'] + vc
    cn = 0
    for slabel, mres in this_dict.items():
        ax = plot_single(ax, mres, resname, si, ei, color=colors[cn], label=slabel)
        cn += 1
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('Virus-clearing TxV')
    ax.legend(loc="upper right", frameon=False)

    # Lesion-regressing TxV
    ax = axes[1]
    this_dict = {
        'No interventions': msim_dict['No interventions'],
        'Status quo': msim_dict['Baseline'],
        '10% screening': msim_dict['TxV 50/90 with 10% screening'],
        '50% screening': msim_dict['TxV 50/90 with 50% screening'],
        '90% screening': msim_dict['TxV 50/90 with 90% screening'],
    }
    vc = sc.vectocolor(3).tolist()
    colors = ['k'] + vc
    cn = 0
    for slabel, mres in this_dict.items():
        ax = plot_single(ax, mres, resname, si, ei, color=colors[cn], label=slabel)
        cn += 1
    ax.set_ylim(bottom=0, top=ymax)
    ax.legend(loc="upper right", frameon=False)
    ax.set_title('Lesion-regressing TxV')

    # Mass vax
    ax = axes[2]
    this_dict = {
        'No interventions': msim_dict['No interventions'],
        'Status quo': msim_dict['Baseline'],
        '10% coverage, 2027': msim_dict['Mass vx 10%'],
        '50% coverage, 2027': msim_dict['Mass vx 50%'],
        '90% coverage, 2027': msim_dict['Mass vx 90%'],
    }
    vc = sc.vectocolor(3).tolist()
    colors = ['k'] + vc
    cn = 0
    for slabel, mres in this_dict.items():
        ax = plot_single(ax, mres, resname, si, ei, color=colors[cn], label=slabel)
        cn += 1
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('Screen, treat & vaccinate 20-25yos')
    ax.legend(loc="upper right", frameon=False)

    # HIV+ vax
    ax = axes[3]
    this_dict = {
        'No interventions': msim_dict['No interventions'],
        'Status quo': msim_dict['Baseline'],
        '10% coverage, 2027': msim_dict['HIV+ vx 10%'],
        '50% coverage, 2027': msim_dict['HIV+ vx 50%'],
        '90% coverage, 2027': msim_dict['HIV+ vx 90%'],
    }
    vc = sc.vectocolor(3).tolist()
    colors = ['k'] + vc
    cn = 0
    for slabel, mres in this_dict.items():
        ax = plot_single(ax, mres, resname, si, ei, color=colors[cn], label=slabel)
        cn += 1
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('Screen, treat & vaccinate WLHIV')

    fig.tight_layout()
    fig_name = 'figures/fig1_vx_scens_all.png'
    sc.savefig(fig_name, dpi=100)

    return


# %% Run as a script
if __name__ == '__main__':

    # Load scenarios and construct figure
    plot_fig1()




