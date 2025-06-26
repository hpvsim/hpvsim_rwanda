"""
Plot residul burden of cervical cancer in Rwanda under screening scenarios
"""


import pylab as pl
import sciris as sc
import numpy as np
import utils as ut


def plot_fig3():
    """
    Plot the residual burden of cervical cancer in Rwanda under screening and treatment scenarios.
    """
    ut.set_font(20)
    fig, axes = pl.subplots(3, 1, layout="tight", figsize=(12, 10))
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
    pn = 0
 
    # Virus-clearing TxV
    ax = axes[pn]
    this_dict = {
        'Status quo': msim_dict['Baseline'],
        '10% screening': msim_dict['TxV 90/50 with 10% screening'],
        '50% screening': msim_dict['TxV 90/50 with 50% screening'],
        '90% screening': msim_dict['TxV 90/50 with 90% screening'],
    }
    vc = sc.vectocolor(3).tolist()
    colors = ['k'] + vc
    cn = 0
    for slabel, mres in this_dict.items():
        ax = ut.plot_single(ax, mres, resname, si, ei, color=colors[cn], label=slabel)
        cn += 1
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('Virus-clearing TxV')
    pn += 1

    # Lesion-regressing TxV
    ax = axes[pn]
    this_dict = {
        'Status quo': msim_dict['Baseline'],
        '10% screening': msim_dict['TxV 50/90 with 10% screening'],
        '50% screening': msim_dict['TxV 50/90 with 50% screening'],
        '90% screening': msim_dict['TxV 50/90 with 90% screening'],
    }
    vc = sc.vectocolor(3).tolist()
    colors = ['k'] + vc
    cn = 0
    for slabel, mres in this_dict.items():
        ax = ut.plot_single(ax, mres, resname, si, ei, color=colors[cn], label=slabel)
        cn += 1
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('Lesion-regressing TxV')
    pn += 1

    # Scaled-up S&T
    ax = axes[pn]
    this_dict = {
        'Status quo': msim_dict['Baseline'],
        '10% screening': msim_dict['S&T 10%'],
        '50% screening': msim_dict['S&T 50%'],
        '90% screening': msim_dict['S&T 90%'],
    }
    vc = sc.vectocolor(3).tolist()
    colors = ['k'] + vc
    cn = 0
    for slabel, mres in this_dict.items():
        ls = ':' if cn == 0 else '-'
        ax = ut.plot_single(ax, mres, resname, si, ei, color=colors[cn], label=slabel)
        cn += 1
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('Scaled up screening & POC ablation')
    ax.legend(loc="upper right", frameon=False, bbox_to_anchor=(1, 0.95), fontsize=18)
    pn += 1

    # # Assemble cumulative results
    # cum_res = dict()
    # resname = 'cancers'
    # for sname, scen in msim_dict.items():
    #     cum_res[sname] = scen[resname].values[si:].sum()
    # # Make a bar plot
    # ax = axes[pn]
    # bars = list(cum_res.values())
    # labels = list(cum_res.keys())
    # x = np.arange(len(labels))
    # ax.bar(x, bars, color=sc.vectocolor(len(bars)))
    # ax.set_xticks(x)
    # ax.set_xticklabels(labels, rotation=45, ha='right')
    # ax.set_title('Cumulative cancers')

    fig.tight_layout()
    fig_name = 'figures/fig3_st.png'
    sc.savefig(fig_name, dpi=100)

    return


# %% Run as a script
if __name__ == '__main__':

    # Load scenarios and construct figure
    plot_fig3()




