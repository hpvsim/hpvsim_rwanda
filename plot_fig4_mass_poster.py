"""
Plot residul burden of cervical cancer in Rwanda under campaign scenarios
"""


import pylab as pl
import sciris as sc
import numpy as np
import utils as ut


def plot_fig4():
    """
    Plot the residual burden of cervical cancer in Rwanda under screening and treatment scenarios.
    """
    ut.set_font(20)
    fig = pl.figure(layout="tight", figsize=(12, 8))
    gs = fig.add_gridspec(2, 2)  # Changed to 2x3 grid
    msim_dict = sc.loadobj('results/st_scens.obj')  # Updated to load campaign scenarios

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
        'HPV-Faster': 'HPV-Faster',
        'Mass TxV virus': 'Mass TxV 90/0,',
        'Mass TxV lesions': 'Mass TxV 50/90,'
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
    ax.set_title('Cumulative cancers 2025-2100')
    sc.SIticks()
    ax.set_xlabel('')
    ax.legend(loc='upper right', frameon=False, fontsize=16)
    ax.set_ylim([0, 100e3])


    ######################################################
    # Top Right: Combined resource use (ablations + therapeutics + vaccinations)
    ######################################################
    ax = fig.add_subplot(gs[0, 1])

    for idx, (strat_label, strat_key) in enumerate(strategies.items()):
        ablations = []
        therapeutics = []
        vaccinations = []
        for cov in coverage_levels:
            scen_key = f'{strat_key} {cov}'
            abl = msim_dict[scen_key]['ablations'][fi:].sum()
            txv = msim_dict[scen_key]['txvs'][fi:].sum()
            vax = msim_dict[scen_key]['vaccinations'][fi:].sum()
            print(f'{scen_key}: {abl} ablations, {txv} therapeutics, {vax} vaccinations')
            ablations.append(abl)
            therapeutics.append(txv)
            vaccinations.append(vax)

        # Plot ablations as solid bars
        ax.bar(x + offsets[idx], ablations, width=bar_width, color=colors[idx])
        # Plot therapeutics stacked on top
        ax.bar(x + offsets[idx], therapeutics, width=bar_width, bottom=ablations,
               color='none', edgecolor=colors[idx], hatch='//', linewidth=1.5)
        # Plot vaccinations stacked on top of therapeutics
        bottom_vals = [a + t for a, t in zip(ablations, therapeutics)]
        ax.bar(x + offsets[idx], vaccinations, width=bar_width, bottom=bottom_vals,
               color='none', edgecolor=colors[idx], hatch='\\\\', linewidth=1.5)

    ax.set_xticks(x)
    ax.set_xticklabels(llabels)
    ax.set_title('Cumulative products 2025-2100')
    sc.SIticks()
    ax.set_xlabel('')

    # Create legend for intervention types
    from matplotlib.patches import Patch
    legend_elements = [ 
        Patch(facecolor='gray', label='Ablations'),
        Patch(facecolor='white', edgecolor='k', hatch='//', label='Therapeutics'),
        Patch(facecolor='white', edgecolor='k', hatch='\\\\', label='Vaccinations')
    ]
    ax.legend(handles=legend_elements, loc='upper left', frameon=False, fontsize=16)
    ax.set_ylim([0, 4e6])


    ######################################################
    # Bottom Left + Middle: Time series at 70% coverage (spanning 2 columns)
    ######################################################
    ax = fig.add_subplot(gs[1, 0])  # Bottom row, first two columns

    # Plot baseline (S&T&T 18%)
    ax = ut.plot_single(ax, msim_dict['S&T&T 18%'], 'asr_cancer_incidence', si, ei,
                        color='k', label='Status quo')

    # Plot each strategy at 70% coverage only
    for strat_idx, (strat_label, strat_key) in enumerate(strategies.items()):
        scen_key = f'{strat_key} 70%'
        ax = ut.plot_single(ax, msim_dict[scen_key], 'asr_cancer_incidence', si, ei,
                           color=colors[strat_idx], label=f'{strat_label} 70%')

    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('Cervical cancer incidence, 2025-2100')

    ax.legend(loc='upper right', frameon=False, fontsize=16)


    ######################################################
    # Bottom Right: Number needed to vaccinate/treat per cancer averted
    ######################################################
    ax = fig.add_subplot(gs[1, 1])

    # Get baseline cancers (S&T&T 18%)
    baseline_cancers = msim_dict['S&T&T 18%']['cancers'].values[fi:].sum()

    for idx, (strat_label, strat_key) in enumerate(strategies.items()):
        nnt_per_cancer = []
        for cov in coverage_levels:
            scen_key = f'{strat_key} {cov}'

            cancers = msim_dict[scen_key]['cancers'].values[fi:].sum()
            cancers_averted = baseline_cancers - cancers

            # Calculate total interventions (therapeutics + vaccinations, but not ablations since those aren't mass campaign)
            therapeutics = msim_dict[scen_key]['txvs'][fi:].sum()
            vaccinations = msim_dict[scen_key]['vaccinations'][fi:].sum()
            total_interventions = therapeutics + vaccinations

            if cancers_averted > 0:
                ratio = total_interventions / cancers_averted
                nnt_per_cancer.append(ratio)
                print(f'{scen_key}: {ratio} products per cancer averted')
            else:
                nnt_per_cancer.append(0)

        ax.bar(x + offsets[idx], nnt_per_cancer, width=bar_width, color=colors[idx], label=strat_label)

    ax.set_xticks(x)
    ax.set_xticklabels(llabels)
    ax.set_title('Products / cancer averted 2025-2100')
    ax.set_xlabel('')

    fig.tight_layout()
    fig_name = 'figures/poster/fig4_campaigns.png'
    sc.savefig(fig_name, dpi=100)

    return


# %% Run as a script
if __name__ == '__main__':

    # Load scenarios and construct figure
    msim_dict = sc.loadobj('results/st_scens.obj')

    # plot_fig3()
    plot_fig4()

