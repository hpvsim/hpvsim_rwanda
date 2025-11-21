"""
Plot 1 for infant vaccination scenarios
"""


import pylab as pl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import sciris as sc
import numpy as np
import utils as ut


def plot_fig1():
    ut.set_font(20)
    fig = pl.figure(layout="tight", figsize=(12, 5))

    gs = fig.add_gridspec(1, 3)
    ax = fig.add_subplot(gs[:2])

    resnames = ['asr_cancer_incidence', 'cancer_incidence_with_hiv']
    colors = ['k', 'r']

    # What to plot
    start_year = 2016
    end_year = 2100
    ymax = 60

    msim_dict = sc.loadobj('results/st_scens.obj')
    mbase = msim_dict['Baseline']
    si = sc.findinds(mbase.year, start_year)[0]
    ei = sc.findinds(mbase.year, end_year)[0]

    this_dict = {
        'No interventions': msim_dict['No interventions'],
        'Status quo': msim_dict['Baseline'],
    }
    for slabel, mres in this_dict.items():
        ls = ':' if slabel == 'No interventions' else '-'
        for resname in resnames:
            color = colors[resnames.index(resname)]
            add_bounds = True if resname == 'asr_cancer_incidence' else False
            ax = ut.plot_single(ax, mres, resname, si, ei, color=color, ls=ls, label=slabel, add_bounds=add_bounds)
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('ASR cervical cancer incidence, 2025-2100')

    # Create handles and labels for the color legend
    # color_handles = [plt.Line2D([0], [0], color='k', lw=2),
    #                  plt.Line2D([0], [0], color='r', lw=2)]
    # color_labels = ['All', 'HIV+']
    circ1 = mpatches.Patch(facecolor='k', label='All')
    circ2 = mpatches.Patch(facecolor='r', label='HIV+')

    # Create handles and labels for the linestyle legend
    linestyle_handles = [plt.Line2D([0], [0], color='k', linestyle='-', lw=2),
                         plt.Line2D([0], [0], color='k', linestyle=':', lw=2)]
    linestyle_labels = ['Status quo', 'No interventions']

    # Create the second legend for linestyles
    legend1 = ax.legend(frameon=False, handles=[circ1, circ2], title='', loc='upper right', bbox_to_anchor=(0.6, 1))
    ax.add_artist(legend1)
    ax.legend(linestyle_handles, linestyle_labels, title='', loc='upper right', frameon=False)

    # Assemble cumulative results
    cum_res = dict()
    resname = 'cancers'
    fi = sc.findinds(mbase.year, 2025)[0]
    for sname, scen in this_dict.items():
        cum_res[sname] = scen[resname].values[fi:].sum()
        lb, ub = scen[resname].low[fi:].sum(), scen[resname].high[fi:].sum()
        print(f'{sname}: {cum_res[sname]} ({lb}, {ub}) cancers')
    # Make a bar plot
    ax = fig.add_subplot(gs[2])
    bars = list(cum_res.values())
    labels = list(cum_res.keys())
    x = np.arange(len(labels))
    ax.bar(x, bars, color='k')
    ax.set_xticks(x)
    ax.set_xticklabels(['No\ninterventions', 'Status\nquo'])
    ax.set_title('Cumulative cancers')
    sc.SIticks()

    fig.tight_layout()
    fig_name = 'figures/fig1_residual.png'
    sc.savefig(fig_name, dpi=100)

    return


# %% Run as a script
if __name__ == '__main__':

    plot_fig1()

    msim_dict = sc.loadobj('results/st_scens.obj')
    mbase = msim_dict['Baseline']
    mno = msim_dict['No interventions']
    start_year = 2016
    end_year = 2100
    si = sc.findinds(mbase.year, start_year)[0]
    ei = sc.findinds(mbase.year, end_year)[0]
    fi = sc.findinds(mbase.year, 2025)[0]

    elim_year = sc.findfirst(msim_dict['Baseline']['asr_cancer_incidence'][si:]<4, die=False)+si

    print(f'Elim year: {mbase.year[elim_year]}')
    print(f'Cancers in 2025: {mbase.cancers[si]} ({mbase.cancers.low[si]}, {mbase.cancers.high[si]})')
    print(f'Cancers in 2100: {mbase.cancers[ei]} ({mbase.cancers.low[ei]}, {mbase.cancers.high[ei]})')
    print(f'Cancers averted: {mno.cancers[fi:].sum()-mbase.cancers[fi:].sum()}')

