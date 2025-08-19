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
    legend_fontsize = 14
    fig, axes = pl.subplots(2, 2, layout="tight", figsize=(12, 7))
    axes = axes.ravel()
    resname = 'asr_cancer_incidence'
    grid = True

    # What to plot
    start_year = 2016
    end_year = 2100
    ymax = 25

    msim_dict = sc.loadobj('results/st_scens.obj')
    mbase = msim_dict['Baseline']
    si = sc.findinds(mbase.year, start_year)[0]
    ei = sc.findinds(mbase.year, end_year)[0]
    pn = 0
 
    # Virus-clearing TxV within screening
    ax = axes[pn]
    this_dict = {
        'Status quo': msim_dict['Baseline'],
        '18% screening': msim_dict['TxV 90/50 within screening, 18%'],
        '35% screening': msim_dict['TxV 90/50 within screening, 35%'],
        '70% screening': msim_dict['TxV 90/50 within screening, 70%'],
    }
    vc = sc.vectocolor(3).tolist()
    colors = ['k'] + vc
    cn = 0
    for slabel, mres in this_dict.items():
        ax = ut.plot_single(ax, mres, resname, si, ei, color=colors[cn], label=slabel)
        cn += 1
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_ylabel('ASR cervical cancer incidence', fontsize=legend_fontsize)
    ax.set_title('Virus-clearing therapeutic\nwithin screening program')
    if grid: ax.grid(axis='y', linestyle='--', linewidth=0.5)
    pn += 1

    # Lesion-regressing TxV within screening
    ax = axes[pn]
    this_dict = {
        'Status quo': msim_dict['Baseline'],
        '18% screening': msim_dict['TxV 50/90 within screening, 18%'],
        '35% screening': msim_dict['TxV 50/90 within screening, 35%'],
        '70% screening': msim_dict['TxV 50/90 within screening, 70%'],
    }
    vc = sc.vectocolor(3).tolist()
    colors = ['k'] + vc
    cn = 0
    for slabel, mres in this_dict.items():
        ax = ut.plot_single(ax, mres, resname, si, ei, color=colors[cn], label=slabel)
        cn += 1
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('Lesion-regressing therapeutic\nwithin screening program')
    if grid: ax.grid(axis='y', linestyle='--', linewidth=0.5)
    pn += 1

    # Virus-clearing TxV
    ax = axes[pn]
    this_dict = {
        'Status quo': msim_dict['Baseline'],
        '18% screening': msim_dict['Mass TxV 90/50, 18%'],
        '35% screening': msim_dict['Mass TxV 90/50, 35%'],
        '70% screening': msim_dict['Mass TxV 90/50, 70%'],
    }
    vc = sc.vectocolor(3).tolist()
    colors = ['k'] + vc
    cn = 0
    for slabel, mres in this_dict.items():
        ax = ut.plot_single(ax, mres, resname, si, ei, color=colors[cn], label=slabel)
        cn += 1
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_ylabel('ASR cervical cancer incidence', fontsize=legend_fontsize)
    ax.set_title('Virus-clearing therapeutic\nmass delivery')
    if grid: ax.grid(axis='y', linestyle='--', linewidth=0.5)
    pn += 1

    # Lesion-regressing TxV
    ax = axes[pn]
    this_dict = {
        'Status quo': msim_dict['Baseline'],
        '18% coverage': msim_dict['Mass TxV 50/90, 18%'],
        '35% coverage': msim_dict['Mass TxV 50/90, 35%'],
        '70% coverage': msim_dict['Mass TxV 50/90, 70%'],
    }
    vc = sc.vectocolor(3).tolist()
    colors = ['k'] + vc
    cn = 0
    for slabel, mres in this_dict.items():
        ax = ut.plot_single(ax, mres, resname, si, ei, color=colors[cn], label=slabel)
        cn += 1
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('Lesion-regressing therapeutic\nmass delivery')
    if grid: ax.grid(axis='y', linestyle='--', linewidth=0.5)
    pn += 1
    ax.legend(loc="upper right", frameon=False, bbox_to_anchor=(1, 0.95), fontsize=14)

    fig.tight_layout()
    fig_name = 'figures/fig4_st.png'
    sc.savefig(fig_name, dpi=100)

    return


# %% Run as a script
if __name__ == '__main__':

    # Load scenarios and construct figure
    msim_dict = sc.loadobj('results/st_scens.obj')

    # plot_fig3()
    plot_fig3()

