"""
Plot cervical cancers among WLHIV
"""


import pylab as pl
import sciris as sc
import numpy as np
import utils as ut


def plot_fig3(msim_dict, start_year=2020, end_year=2100, ymax=50):

    # Settings
    ut.set_font(20)
    fig = pl.figure(layout="tight", figsize=(12, 5))

    # Axis and grid
    gs = fig.add_gridspec(1, 3)

    mbase = msim_dict['Baseline']
    si = sc.findinds(mbase.year, start_year)[0]
    ei = sc.findinds(mbase.year, end_year)[0]

    # HIV scenarios
    this_dict = {
        'Status quo': msim_dict['Baseline'],
        '18% coverage': msim_dict['HIV+ vx 18%'],
        '35% coverage': msim_dict['HIV+ vx 35%'],
        '70% coverage': msim_dict['HIV+ vx 70%'],
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

    ax = fig.add_subplot(gs[:2])
    cn = 0
    for slabel, mres in this_dict.items():
        ax = ut.plot_single(ax, mres, 'cancer_incidence_with_hiv', si, ei, color=colors[cn], label=slabel)
        cn += 1
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('Cervical cancer incidence among HIV+ women')
    ax.legend(loc="upper right", frameon=False, bbox_to_anchor=(1, 0.95), fontsize=18)

    # Assemble cumulative results
    cum_res = dict()
    resname = 'cancers_with_hiv'
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
    # Make a bar plot using the colors above for each bar
    ax.bar(x, bars, color=colors)
    ax.set_xticks(x)
    llabels = ['SQ', '18%', '35%', '70%']
    ax.set_xticklabels(llabels)
    ax.set_title('Cumulative cancers\namong HIV+ women\n2025-2100')
    sc.SIticks()


    fig.tight_layout()
    fig_name = 'figures/fig3_hiv.png'
    sc.savefig(fig_name, dpi=100)
    return


# %% Run as a script
if __name__ == '__main__':

    # Load scenarios and construct figure
    msim_dict = sc.loadobj('results/st_scens.obj')

    # plot_fig3()
    plot_fig3(msim_dict)

