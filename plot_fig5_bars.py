"""
Plot residual burden of cervical cancer in Rwanda under vaccination scenarios
"""
import pandas as pd


import pylab as pl
import sciris as sc
import numpy as np
import utils as ut


def plot_fig5(poster=False, end_year=2100):
    """
    Plot the residual burden of cervical cancer in Rwanda under screening and treatment scenarios.
    """
    fs = 20 if not poster else 24
    ut.set_font(fs)
    figsize=(14, 12) if not poster else (14, 12)
    fig = pl.figure(layout="tight", figsize=figsize)
    gs = fig.add_gridspec(2, 1)  # 2 rows, 1 column

    # Load all scenario data from single dictionary
    msim_dict = sc.loadobj('results/st_scens.obj')
    text_font = 14 if not poster else 20

    # Common setup
    fi = sc.findinds(msim_dict['S&T&T 18%'].year, 2025)[0]
    coverage_levels = ['18%', '35%', '70%']

    # Define all strategies with keys and display labels
    all_strategies = [
        ('S&T', 'Screen+\ntreat'),
        ('S&T&T', 'Screen+\ntriage+\ntreat'),
        ('S&TxV&T&T', 'Screen+\nTxV 90/0+\ntriage+treat'),
        ('S&TxV', 'Screen+\nTxV 50/90'),
        ('HPV-Faster', 'HPV-Faster'),
        ('Mass TxV 90/0,', 'Mass TxV\n90/0'),
        ('Mass TxV 50/90,', 'Mass TxV\n50/90'),
    ]

    # Create color palette for coverage levels
    coverage_colors = sc.vectocolor(len(coverage_levels)).tolist()

    # Get baseline cancers
    baseline_cancers = msim_dict['S&T&T 18%']['cancers'].values[fi:].sum()

    ######################################################
    # Top Panel: Cumulative cancers
    ######################################################
    ax = fig.add_subplot(gs[0])

    bar_width = 0.28
    x_base = np.arange(len(all_strategies))
    offsets = [-bar_width, 0, bar_width]

    for cov_idx, cov in enumerate(coverage_levels):
        cum_cancers = []

        for strat_key, strat_label in all_strategies:
            scen_key = f'{strat_key} {cov}'
            val = msim_dict[scen_key]['cancers'].values[fi:].sum()
            cum_cancers.append(val)
            print(f'{scen_key}: {val} cancers')

        bars = ax.bar(x_base + offsets[cov_idx], cum_cancers, width=bar_width,
                      color=coverage_colors[cov_idx], label=cov)

        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                    f'{sc.sigfig(height/1e3, 3)}',
                    ha='center', va='bottom', fontsize=text_font)

    ax.set_xticks(x_base)
    ax.set_xticklabels([label for _, label in all_strategies])
    ax.set_title('Cumulative cancers 2025-2100')
    sc.SIticks()
    ax.set_ylim([0, 100e3])
    ax.legend(title='Coverage', loc='upper right', frameon=False, fontsize=14, ncols=3)

    ######################################################
    # Bottom Panel: Cancers averted
    ######################################################
    ax = fig.add_subplot(gs[1])

    for cov_idx, cov in enumerate(coverage_levels):
        cancers_averted = []

        for strat_key, strat_label in all_strategies:
            scen_key = f'{strat_key} {cov}'
            val = msim_dict[scen_key]['cancers'].values[fi:].sum()
            averted = max(baseline_cancers - val, 0)
            cancers_averted.append(averted)
            print(f'{scen_key}: {averted} cancers averted')

        bars = ax.bar(x_base + offsets[cov_idx], cancers_averted, width=bar_width,
                      color=coverage_colors[cov_idx], label=cov)

        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            if height > 0:  # Only label positive values
                ax.text(bar.get_x() + bar.get_width()/2., height,
                        f'{sc.sigfig(height/1e3, 3)}',
                        ha='center', va='bottom', fontsize=text_font)

    ax.set_xticks(x_base)
    ax.set_xticklabels([label for _, label in all_strategies])
    ax.set_title('Cancers averted 2025-2100 (vs. S&T&T 18%)')
    sc.SIticks()
    ax.set_ylim([0, 38e3])
    # ax.legend(title='Coverage', loc='upper right', frameon=False, fontsize=14)

    fig.tight_layout()
    folder = 'figures/' if not poster else 'figures/poster/'
    fig_name = f'{folder}fig5_comparison.png'
    sc.savefig(fig_name, dpi=100)
    # return perc_averted, new_bars
    return


# %% Run as a script
if __name__ == '__main__':

    # Load scenarios and construct figure
    plot_fig5(poster=False)



