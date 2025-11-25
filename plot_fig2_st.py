"""
Plot residul burden of cervical cancer in Rwanda
"""


import pylab as pl
import sciris as sc
import utils as ut
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

def plot_fig2():
    """
    Plot the residual burden of cervical cancer in Rwanda under vaccination scenarios.
    """
    ut.set_font(20)
    fig = pl.figure(layout="tight", figsize=(18, 10))  # Wider for 3 columns
    gs = fig.add_gridspec(2, 3)  # Changed to 2x3 grid
    msim_dict = sc.loadobj('results/st_scens.obj')
    mbase = msim_dict['Baseline']

    text_height = [-0.1, 1.2]

    # What to plot
    start_year = 2016
    end_year = 2100
    ymax = 25
    si = sc.findinds(mbase.year, start_year)[0]
    ei = sc.findinds(mbase.year, end_year)[0]
    vc = sc.vectocolor(3).tolist()
    colors = vc

    this_dict = {
        '18% coverage': msim_dict['S&T 18%'],
        '35% coverage': msim_dict['S&T 35%'],
        '70% coverage': msim_dict['S&T 70%'],
    }
    triage_dict = {
        '18% coverage': msim_dict['S&T&T 18%'],  # This is the actual status quo
        '35% coverage': msim_dict['S&T&T 35%'],
        '70% coverage': msim_dict['S&T&T 70%'],
    }
    triage_list = list(triage_dict.values())

    # Common setup for bar plots
    fi = sc.findinds(mbase.year, 2025)[0]
    bar_width = 0.35
    labels = list(this_dict.keys())
    x = np.arange(len(labels)).astype(float)
    llabels = ['18%', '35%', '70%']

    # Calculate positions for side-by-side bars
    no_triage_x = x - bar_width/2  # Shift no-triage bars left
    triage_x = x + bar_width/2  # Shift triage bars right

    ######################################################
    # Top Left: Cumulative cancers
    ######################################################
    ax = fig.add_subplot(gs[0, 0])

    # Assemble cumulative results
    cum_res = dict()
    resname = 'cancers'
    for sname, scen in this_dict.items():
        cum_res[sname] = scen[resname].values[fi:].sum()
        lb, ub = scen[resname].low[fi:].sum(), scen[resname].high[fi:].sum()
        print(f'{sname}: {cum_res[sname]} ({lb}, {ub}) cancers')

    bars = list(cum_res.values())

    # Assemble VIA triage results
    triage_cum_res = dict()
    for sname, scen in triage_dict.items():
        triage_cum_res[sname] = scen[resname].values[fi:].sum()
        lb, ub = scen[resname].low[fi:].sum(), scen[resname].high[fi:].sum()
        print(f'Triage {sname}: {triage_cum_res[sname]} ({lb}, {ub}) cancers')
    triage_bars = list(triage_cum_res.values())

    # Plot bars
    ax.bar(no_triage_x, bars, width=bar_width, color=colors)
    ax.bar(triage_x, triage_bars, width=bar_width, color=colors, edgecolor='k',
           hatch='//', linewidth=1.5)

    ax.set_xticks(x)
    ax.set_xticklabels(llabels)
    ax.set_title('Cumulative cancers\n2025-2100')
    sc.SIticks()
    ax.set_xlabel('')

    # Create pattern legend for triage vs no triage
    pattern_handles = [
        mpatches.Patch(facecolor='gray', label='No triage'),
        mpatches.Patch(facecolor='white', edgecolor='k', hatch='//', label='VIA triage')
    ]
    ax.legend(handles=pattern_handles, loc='upper right', frameon=False)
    ax.set_ylim([0, 100e3])
    ax.text(*text_height, 'A', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    ######################################################
    # Top Middle: Cumulative cancers in HIV+ women
    ######################################################
    ax = fig.add_subplot(gs[0, 1])

    # Assemble cumulative results for HIV+ cancers
    cum_hiv_res = dict()
    resname = 'cancers_with_hiv'
    for sname, scen in this_dict.items():
        cum_hiv_res[sname] = scen[resname].values[fi:].sum()
        print(f'{sname}: {cum_hiv_res[sname]} cancers in HIV+ women')

    hiv_bars = list(cum_hiv_res.values())

    # Assemble VIA triage results for HIV+
    triage_cum_hiv_res = dict()
    for sname, scen in triage_dict.items():
        triage_cum_hiv_res[sname] = scen[resname].values[fi:].sum()
        print(f'Triage {sname}: {triage_cum_hiv_res[sname]} cancers in HIV+ women')
    triage_hiv_bars = list(triage_cum_hiv_res.values())

    # Plot bars
    ax.bar(no_triage_x, hiv_bars, width=bar_width, color=colors)
    ax.bar(triage_x, triage_hiv_bars, width=bar_width, color=colors, edgecolor='k',
           hatch='//', linewidth=1.5)

    ax.set_xticks(x)
    ax.set_xticklabels(llabels)
    ax.set_title('Cumulative cancers in HIV+ women\n2025-2100')
    sc.SIticks()
    ax.set_xlabel('')
    # Add panel label
    ax.text(*text_height, 'B', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    ######################################################
    # Top Right: Cumulative ablations
    ######################################################
    ax = fig.add_subplot(gs[0, 2])

    # Assemble cumulative ablations
    cum_ablations = dict()
    for sname, scen in this_dict.items():
        cum_ablations[sname] = scen['ablations'][fi:].sum()
        print(f'{sname}: {cum_ablations[sname]} ablations')

    ablation_bars = list(cum_ablations.values())

    # Assemble VIA triage ablations
    triage_cum_ablations = dict()
    for sname, scen in triage_dict.items():
        triage_cum_ablations[sname] = scen['ablations'][fi:].sum()
        print(f'Triage {sname}: {triage_cum_ablations[sname]} ablations')
    triage_ablation_bars = list(triage_cum_ablations.values())

    # Plot bars
    ax.bar(no_triage_x, ablation_bars, width=bar_width, color=colors)
    ax.bar(triage_x, triage_ablation_bars, width=bar_width, color=colors, edgecolor='k',
           hatch='//', linewidth=1.5)

    ax.set_xticks(x)
    ax.set_xticklabels(llabels)
    ax.set_title('Cumulative ablations\n2025-2100')
    sc.SIticks()
    ax.set_xlabel('')
    # Add panel label
    ax.text(*text_height, 'C', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    ######################################################
    # Bottom Left + Middle: Time series (spanning 2 columns)
    ######################################################
    ax = fig.add_subplot(gs[1, :2])  # Bottom row, first two columns

    # Add status quo (18% with triage) first
    ax = ut.plot_single(ax, triage_dict['18% coverage'], 'asr_cancer_incidence', si, ei, color='k', label='Status quo (18% + triage)')

    cn = 0
    for slabel, mres in this_dict.items():
        ax = ut.plot_single(ax, mres, 'asr_cancer_incidence', si, ei, color=colors[cn], label=slabel)
        cn += 1

    # Plot triage scenarios for 35% and 70%
    for cn, sname in enumerate(['35% coverage', '70% coverage'], start=1):
        ax = ut.plot_single(ax, triage_dict[sname], 'asr_cancer_incidence', si, ei, color=colors[cn], ls='--', label='')

    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('ASR cervical cancer incidence, 2025-2100\nScaled-up screening with/without VIA triage')

    # Create handles and labels for the color legend
    circ1 = mpatches.Patch(facecolor=colors[0], label='18% coverage')
    circ2 = mpatches.Patch(facecolor=colors[1], label='35% coverage')
    circ3 = mpatches.Patch(facecolor=colors[2], label='70% coverage')

    # Create handles and labels for the linestyle legend
    linestyle_handles = [plt.Line2D([0], [0], color='k', linestyle='-', lw=2),
                         plt.Line2D([0], [0], color='k', linestyle='--', lw=2)]
    linestyle_labels = ['No triage', 'VIA triage']

    # Create the second legend for linestyles
    legend1 = ax.legend(frameon=False, handles=[circ1, circ2, circ3], title='', loc='upper right', bbox_to_anchor=(1, 0.7))
    ax.add_artist(legend1)
    ax.legend(linestyle_handles, linestyle_labels, title='', loc='upper right', frameon=False)
    # Add panel label
    ax.text(*[-0.05, 1.2], 'D', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    ######################################################
    # Bottom Right: Ablations per cancer averted
    ######################################################
    ax = fig.add_subplot(gs[1, 2])

    # Calculate cancers averted relative to 18% + triage baseline
    baseline_cancers = triage_cum_res['18% coverage']
    cancers_averted = dict()
    for sname in labels:
        cancers_averted[sname] = baseline_cancers - cum_res[sname]
        print(f'{sname}: {cancers_averted[sname]} cancers averted')

    triage_cancers_averted = dict()
    for sname in labels:
        triage_cancers_averted[sname] = baseline_cancers - triage_cum_res[sname]
        print(f'Triage {sname}: {triage_cancers_averted[sname]} cancers averted')

    # Calculate ablations per cancer averted
    ablations_per_cancer = []
    for sname in labels:
        if cancers_averted[sname] > 0:
            ratio = cum_ablations[sname] / cancers_averted[sname]
            ablations_per_cancer.append(ratio)
            print(f'{sname}: {ratio} ablations per cancer averted')
        else:
            ablations_per_cancer.append(0)

    triage_ablations_per_cancer = []
    for sname in labels:
        if triage_cancers_averted[sname] > 0:
            ratio = triage_cum_ablations[sname] / triage_cancers_averted[sname]
            triage_ablations_per_cancer.append(ratio)
            print(f'Triage {sname}: {ratio} ablations per cancer averted')
        else:
            triage_ablations_per_cancer.append(0)

    # Plot bars
    ax.bar(no_triage_x, ablations_per_cancer, width=bar_width, color=colors)
    ax.bar(triage_x, triage_ablations_per_cancer, width=bar_width, color=colors,
           edgecolor='k', hatch='//', linewidth=1.5)

    ax.set_xticks(x)
    ax.set_xticklabels(llabels)
    ax.set_title('Ablations per cancer averted\n2025-2100')
    ax.set_xlabel('')
    # Add panel label
    ax.text(*text_height, 'E', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    fig.tight_layout()
    fig_name = 'figures/fig2_st.png'
    sc.savefig(fig_name, dpi=100)

    return


# %% Run as a script
if __name__ == '__main__':

    # Load scenarios and construct figure
    plot_fig2()

    msim_dict = sc.loadobj('results/st_scens.obj')
