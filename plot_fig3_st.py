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
        '20% screening': msim_dict['TxV 90/50 with 18% screening'],
        '35% screening': msim_dict['TxV 90/50 with 35% screening'],
        '70% screening': msim_dict['TxV 90/50 with 70% screening'],
    }
    vc = sc.vectocolor(3).tolist()
    colors = ['k'] + vc
    cn = 0
    for slabel, mres in this_dict.items():
        ax = ut.plot_single(ax, mres, resname, si, ei, color=colors[cn], label=slabel)
        cn += 1
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('Virus-clearing therapeutic')
    pn += 1

    # Lesion-regressing TxV
    ax = axes[pn]
    this_dict = {
        'Status quo': msim_dict['Baseline'],
        '20% screening': msim_dict['TxV 50/90 with 18% screening'],
        '35% screening': msim_dict['TxV 50/90 with 35% screening'],
        '70% screening': msim_dict['TxV 50/90 with 70% screening'],
    }
    vc = sc.vectocolor(3).tolist()
    colors = ['k'] + vc
    cn = 0
    for slabel, mres in this_dict.items():
        ax = ut.plot_single(ax, mres, resname, si, ei, color=colors[cn], label=slabel)
        cn += 1
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('Lesion-regressing therapeutic')
    pn += 1

    # Scaled-up S&T
    ax = axes[pn]
    this_dict = {
        'Status quo': msim_dict['Baseline'],
        '20% screening': msim_dict['S&T 18%'],
        '35% screening': msim_dict['S&T 35%'],
        '70% screening': msim_dict['S&T 70%'],
    }
    vc = sc.vectocolor(3).tolist()
    colors = ['k'] + vc
    cn = 0
    for slabel, mres in this_dict.items():
        ls = ':' if cn == 0 else '-'
        ax = ut.plot_single(ax, mres, resname, si, ei, color=colors[cn], label=slabel)
        cn += 1
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('Scaled up screening & ablation')
    ax.legend(loc="upper right", frameon=False, bbox_to_anchor=(1, 0.95), fontsize=18)
    pn += 1

    fig.tight_layout()
    fig_name = 'figures/fig3_st.png'
    sc.savefig(fig_name, dpi=100)

    return


# %% Run as a script
if __name__ == '__main__':

    # Load scenarios and construct figure
    plot_fig3()




