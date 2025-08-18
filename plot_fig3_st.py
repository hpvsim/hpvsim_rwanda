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
    ax.set_title('ASR cervical cancer incidence, 2025-2100\nVirus-clearing therapeutic')
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
    ax.set_title('ASR cervical cancer incidence, 2025-2100\nLesion-regressing therapeutic')
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
    ax.set_title('ASR cervical cancer incidence, 2025-2100\nScaled up screening & ablation')
    ax.legend(loc="upper right", frameon=False, bbox_to_anchor=(1, 0.95), fontsize=18)
    pn += 1

    fig.tight_layout()
    fig_name = 'figures/fig3_st.png'
    sc.savefig(fig_name, dpi=100)

    return


def plot_fig_hiv(msim_dict, start_year=2016, end_year=2100, ymax=50):
    ut.set_font(20)
    fig, axes = pl.subplots(2, 1, layout="tight", figsize=(12, 12))
    mbase = msim_dict['Baseline']
    si = sc.findinds(mbase.year, start_year)[0]
    ei = sc.findinds(mbase.year, end_year)[0]

    # HIV scenarios
    this_dict = {
        'Status quo': msim_dict['Baseline'],
        '20% screening': msim_dict['HIV+ vx 18%'],
        '35% screening': msim_dict['HIV+ vx 35%'],
        '70% screening': msim_dict['HIV+ vx 70%'],
    }

    # Print averted totals and percent averted
    base = msim_dict['Baseline']['cancers_with_hiv'][si:ei].sum()
    optim = msim_dict['HIV+ vx 70%']['cancers_with_hiv'][si:ei].sum()
    print(f'Baseline cancers with HIV: {base:.0f}')
    print(f'Optimized cancers with HIV: {optim:.0f}')
    print(f'Averted cancers with HIV: {base-optim:.0f}')
    print(f'Percent averted: {(base-optim)/base*100:.1f}%')

    vc = sc.vectocolor(3).tolist()
    colors = ['k'] + vc

    cn = 0
    for slabel, mres in this_dict.items():
        ax = axes[0]
        ax = ut.plot_single(ax, mres, 'cancer_incidence_with_hiv', si, ei, color=colors[cn], label=slabel)
        cn += 1
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('Cervical cancer incidence among HIV+ women')
    # ax.legend(loc="upper right", frameon=False, bbox_to_anchor=(1, 0.95), fontsize=18)

    cn = 0
    for slabel, mres in this_dict.items():
        ax = axes[1]
        ax = ut.plot_single(ax, mres, 'cancers_with_hiv', si, ei, color=colors[cn], label=slabel)
        cn += 1
    ax.set_ylim(bottom=0)
    ax.set_title('Cervical cancers among HIV+ women')
    ax.legend(loc="upper right", frameon=False, bbox_to_anchor=(1, 0.95), fontsize=18)

    fig.tight_layout()
    fig_name = 'figures/fig_hiv.png'
    sc.savefig(fig_name, dpi=100)
    return


# %% Run as a script
if __name__ == '__main__':

    # Load scenarios and construct figure
    msim_dict = sc.loadobj('results/st_scens.obj')

    plot_fig3()
    plot_fig_hiv(msim_dict)

    mbase = msim_dict['Baseline']
    strategies = [
        'TxV 50/90 with 18% screening',
        'TxV 90/50 with 18% screening',
        'TxV 50/90 with 35% screening',
        'TxV 90/50 with 35% screening',
        'TxV 50/90 with 70% screening',
        'TxV 90/50 with 70% screening',
        'S&T 18%',
        'S&T 35%',
        'S&T 70%',
        ]

    fi = sc.findinds(mbase.year, 2025)[0]
    for sname in strategies:
        mres = msim_dict[sname]
        elim_year = sc.findfirst(mres['asr_cancer_incidence'][fi:]<4)+fi
        print(f'{sname} elim year: {mbase.year[elim_year]}')


