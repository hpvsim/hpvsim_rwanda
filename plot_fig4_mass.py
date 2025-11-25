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
    fig = pl.figure(layout="tight", figsize=(18, 10))  # Wider for 3 columns
    gs = fig.add_gridspec(2, 3)  # Changed to 2x3 grid
    msim_dict = sc.loadobj('results/st_scens.obj')  # Updated to load campaign scenarios

    # What to plot
    start_year = 2016
    end_year = 2100
    ymax = 25
    si = sc.findinds(msim_dict['S&T&T 18%'].year, start_year)[0]
    ei = sc.findinds(msim_dict['S&T&T 18%'].year, end_year)[0]
    vc = sc.vectocolor(3).tolist()
    colors = vc

    text_height = [-0.1, 1.2]

    # Define the three strategies and coverage levels
    coverage_levels = ['18%', '35%', '70%']
    strategies = {
        'HPV-Faster': 'HPV-Faster',
        'Mass TxV 90/0': 'Mass TxV 90/0,',
        'Mass TxV 50/90': 'Mass TxV 50/90,'
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
    ax.legend(loc='upper right', frameon=False)

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
    # Top Right: Combined resource use (ablations + therapeutics + vaccinations)
    ######################################################
    ax = fig.add_subplot(gs[0, 2])

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
    ax.set_title('Cumulative interventions\n2025-2100')
    sc.SIticks()
    ax.set_xlabel('')

    # Create legend for intervention types
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='gray', label='Ablations'),
        Patch(facecolor='white', edgecolor='k', hatch='//', label='Therapeutics'),
        Patch(facecolor='white', edgecolor='k', hatch='\\\\', label='Vaccinations')
    ]
    ax.legend(handles=legend_elements, loc='upper left', frameon=False)

    # Add panel label
    ax.text(*text_height, 'C', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    ######################################################
    # Bottom Left + Middle: Time series at 70% coverage (spanning 2 columns)
    ######################################################
    ax = fig.add_subplot(gs[1, :2])  # Bottom row, first two columns

    # Plot baseline (S&T&T 18%)
    ax = ut.plot_single(ax, msim_dict['S&T&T 18%'], 'asr_cancer_incidence', si, ei,
                        color='k', label='Status quo (S&T&T 18%)')

    # Plot each strategy at 70% coverage only
    for strat_idx, (strat_label, strat_key) in enumerate(strategies.items()):
        scen_key = f'{strat_key} 70%'
        ax = ut.plot_single(ax, msim_dict[scen_key], 'asr_cancer_incidence', si, ei,
                           color=colors[strat_idx], label=f'{strat_label} 70%')

    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('ASR cervical cancer incidence, 2025-2100\nOne-time mass campaigns at 70% coverage')

    ax.legend(loc='upper right', frameon=False)

    # Add panel label
    ax.text(-0.05, 1.2, 'D', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    ######################################################
    # Bottom Right: Number needed to vaccinate/treat per cancer averted
    ######################################################
    ax = fig.add_subplot(gs[1, 2])

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
                print(f'{scen_key}: {ratio} interventions per cancer averted')
            else:
                nnt_per_cancer.append(0)

        ax.bar(x + offsets[idx], nnt_per_cancer, width=bar_width, color=colors[idx], label=strat_label)

    ax.set_xticks(x)
    ax.set_xticklabels(llabels)
    ax.set_title('Interventions per cancer averted\n2025-2100')
    ax.set_xlabel('')

    # Add panel label
    ax.text(*text_height, 'E', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    fig.tight_layout()
    fig_name = 'figures/fig4_campaigns.png'
    sc.savefig(fig_name, dpi=100)

    return


# %% Run as a script
if __name__ == '__main__':

    # Load scenarios and construct figure
    msim_dict = sc.loadobj('results/st_scens.obj')

    # plot_fig3()
    plot_fig4()

