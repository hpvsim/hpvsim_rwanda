"""
Plot residul burden of cervical cancer in Rwanda
"""


import pylab as pl
import sciris as sc
import utils as ut
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

def plot_fig3():
    ut.set_font(20)
    fig = pl.figure(layout="tight", figsize=(18, 10))  # Wider for 3 columns
    gs = fig.add_gridspec(2, 3)  # Changed to 2x3 grid
    msim_dict = sc.loadobj('results/st_scens.obj')  # Updated to load therapeutic scenarios

    text_height = [-0.1, 1.2]

    # What to plot
    start_year = 2016
    end_year = 2100
    ymax = 25
    si = sc.findinds(msim_dict['S&T&T 18%'].year, start_year)[0]
    ei = sc.findinds(msim_dict['S&T&T 18%'].year, end_year)[0]
    vc = sc.vectocolor(3).tolist()
    vc2 = sc.vectocolor(4, cmap='magma').tolist()
    colors = [vc[0], vc2[1], vc2[2]]  # Custom color selection

    # Define the three strategies and coverage levels
    coverage_levels = ['18%', '35%', '70%']
    strategies = {
        'Status quo': 'S&T&T',
        'TxV 90/0': 'S&TxV&T&T',
        'TxV 50/90': 'S&TxV'
    }

    # Common setup for bar plots
    fi = sc.findinds(msim_dict['S&T&T 18%'].year, 2025)[0]
    bar_width = 0.25
    x = np.arange(len(coverage_levels))
    llabels = ['18%', '35%', '70%']

    # Calculate positions for grouped bars (3 bars per coverage level)
    offsets = [-bar_width, 0, bar_width]

    ######################################################
    # Top Left: Cumulative cancers
    ######################################################
    ax = fig.add_subplot(gs[0, 0])

    resname = 'cancers'
    for idx, (strat_label, strat_key) in enumerate(strategies.items()):
        cum_res = []
        for cov in coverage_levels:
            scen_key = f'{strat_key} {cov}'
            val = msim_dict[scen_key][resname].values[fi:].sum()
            lb, ub = msim_dict[scen_key][resname].low[fi:].sum(), msim_dict[scen_key][resname].high[fi:].sum()
            print(f'{scen_key}: {val} ({lb}, {ub}) cancers')
            cum_res.append(val)

        ax.bar(x + offsets[idx], cum_res, width=bar_width, color=colors[idx], label=strat_label)

    ax.set_xticks(x)
    ax.set_xticklabels(llabels)
    ax.set_title('Cumulative cancers\n2025-2100')
    sc.SIticks()
    ax.set_xlabel('')
    # Two columns
    ax.legend(loc='upper right', frameon=False, fontsize=16)
    ax.set_ylim([0, 100e3])

    # Add panel label
    ax.text(*text_height, 'A', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    ######################################################
    # Top Middle: Cumulative cancers in HIV+ women
    ######################################################
    ax = fig.add_subplot(gs[0, 1])

    resname = 'cancers_with_hiv'
    for idx, (strat_label, strat_key) in enumerate(strategies.items()):
        cum_hiv_res = []
        for cov in coverage_levels:
            scen_key = f'{strat_key} {cov}'
            val = msim_dict[scen_key][resname].values[fi:].sum()
            print(f'{scen_key}: {val} cancers in HIV+ women')
            cum_hiv_res.append(val)

        ax.bar(x + offsets[idx], cum_hiv_res, width=bar_width, color=colors[idx], label=strat_label)

    ax.set_xticks(x)
    ax.set_xticklabels(llabels)
    ax.set_title('Cumulative cancers in HIV+ women\n2025-2100')
    sc.SIticks()
    ax.set_xlabel('')

    # Add panel label
    ax.text(*text_height, 'B', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    ######################################################
    # Top Right: Combined resource use (ablations + therapeutics)
    ######################################################
    ax = fig.add_subplot(gs[0, 2])

    for idx, (strat_label, strat_key) in enumerate(strategies.items()):
        ablations = []
        therapeutics = []
        for cov in coverage_levels:
            scen_key = f'{strat_key} {cov}'
            abl = msim_dict[scen_key]['ablations'][fi:].sum()
            txv = msim_dict[scen_key]['txvs'][fi:].sum()
            print(f'{scen_key}: {abl} ablations, {txv} therapeutics')
            ablations.append(abl)
            therapeutics.append(txv)

        # Plot ablations as solid bars
        ax.bar(x + offsets[idx], ablations, width=bar_width, color=colors[idx])
        # Plot therapeutics stacked on top as hatched bars
        ax.bar(x + offsets[idx], therapeutics, width=bar_width, bottom=ablations,
               color='none', edgecolor=colors[idx], hatch='//', linewidth=1.5)

    ax.set_xticks(x)
    ax.set_xticklabels(llabels)
    ax.set_title('Cumulative treatments\n2025-2100')
    sc.SIticks()
    ax.set_xlabel('')

    # Create legend for treatment types
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='gray', label='Ablations'),
        Patch(facecolor='white', edgecolor='k', hatch='//', label='Therapeutics')
    ]
    ax.legend(handles=legend_elements, loc='upper left', frameon=False, fontsize=16)

    # Add panel label
    ax.text(*text_height, 'C', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    ######################################################
    # Bottom Left + Middle: Time series (spanning 2 columns)
    ######################################################
    ax = fig.add_subplot(gs[1, :2])  # Bottom row, first two columns

    # Plot baseline (S&T&T 18%)
    ax = ut.plot_single(ax, msim_dict['S&T&T 18%'], 'asr_cancer_incidence', si, ei,
                        color='k', label='Status quo (S&T&T 18%)')

    # Plot each strategy at different coverage levels
    # line_styles = ['-', '--', ':']
    for cov_idx, cov in enumerate(['70%']):
        for strat_idx, (strat_label, strat_key) in enumerate(strategies.items()):
            scen_key = f'{strat_key} {cov}'
            if scen_key == 'S&T&T 18%':  # Skip baseline, already plotted
                continue
            ax = ut.plot_single(ax, msim_dict[scen_key], 'asr_cancer_incidence', si, ei,
                               color=colors[strat_idx], # ls=line_styles[strat_idx],
                               label=f'{strat_label} {cov}' if cov_idx == 0 and strat_idx == 1 else '')

    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('ASR cervical cancer incidence, 2025-2100\nComparison of screening strategies')
    thesecolors = ['k'] + colors
    labels = ['Status quo', 'SOC algorithm, 70% coverage', '+TxV 90/0, 70% coverage', '+TxV 50/90, 70% coverage']

    # Create legends
    # Coverage legend
    strat_handles = [Patch(facecolor=thesecolors[i], label=labels[i]) for i in range(4)]
    legend1 = ax.legend(handles=strat_handles, loc='upper right', bbox_to_anchor=(1, 0.8), frameon=False)
    ax.add_artist(legend1)

    # # Coverage legend
    # from matplotlib.lines import Line2D
    # strat_handles = [Line2D([0], [0], color='k', linestyle=line_styles[i], lw=2, label=list(strategies.keys())[i])
    #                for i in range(len(coverage_levels))]
    # ax.legend(handles=strat_handles, loc='upper right', frameon=False)

    # Add panel label
    ax.text(-0.05, 1.2, 'D', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    ######################################################
    # Bottom Right: Total treatments per cancer averted
    ######################################################
    ax = fig.add_subplot(gs[1, 2])

    # Get baseline cancers (S&T&T 18%)
    baseline_cancers = msim_dict['S&T&T 18%']['cancers'].values[fi:].sum()

    for idx, (strat_label, strat_key) in enumerate(strategies.items()):
        treatments_per_cancer = []
        for cov in coverage_levels:
            scen_key = f'{strat_key} {cov}'
            if scen_key == 'S&T&T 18%':  # Skip baseline
                treatments_per_cancer.append(0)
                continue

            cancers = msim_dict[scen_key]['cancers'].values[fi:].sum()
            cancers_averted = baseline_cancers - cancers

            ablations = msim_dict[scen_key]['ablations'][fi:].sum()
            therapeutics = msim_dict[scen_key]['txvs'][fi:].sum()
            total_treatments = ablations + therapeutics

            if cancers_averted > 0:
                ratio = total_treatments / cancers_averted
                treatments_per_cancer.append(ratio)
                print(f'{scen_key}: {ratio} treatments per cancer averted')
            else:
                treatments_per_cancer.append(0)

        # Filter out zeros (baseline)
        if strat_key == 'S&T&T':
            # For S&T&T strategy, only plot 35% and 70%
            x_plot = x[1:] + offsets[idx]
            values_plot = treatments_per_cancer[1:]
        else:
            x_plot = x + offsets[idx]
            values_plot = treatments_per_cancer

        ax.bar(x_plot, values_plot, width=bar_width, color=colors[idx], label=strat_label)

    ax.set_xticks(x)
    ax.set_xticklabels(llabels)
    ax.set_title('Total treatments per cancer averted\n2025-2100')
    ax.set_xlabel('')

    # Add panel label
    ax.text(*text_height, 'E', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    fig.tight_layout()
    fig_name = 'figures/fig3_txv.png'
    sc.savefig(fig_name, dpi=100)
    return


# %% Run as a script
if __name__ == '__main__':

    # Load scenarios and construct figure
    plot_fig3()

    msim_dict = sc.loadobj('results/st_scens.obj')
