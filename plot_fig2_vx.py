"""
Plot residul burden of cervical cancer in Rwanda under vaccination scenarios
"""


import pylab as pl
import sciris as sc
import numpy as np
import utils as ut
import seaborn as sns


def plot_fig2():
    """
    Plot the residual burden of cervical cancer in Rwanda under vaccination scenarios.
    """
    ut.set_font(20)
    fig = pl.figure(layout="tight", figsize=(12, 5))
    gs = fig.add_gridspec(1, 3)
    msim_dict = sc.loadobj('results/st_scens.obj')
    mbase = msim_dict['Baseline']

    # What to plot
    start_year = 2016
    end_year = 2100
    ymax = 25
    si = sc.findinds(mbase.year, start_year)[0]
    ei = sc.findinds(mbase.year, end_year)[0]
    vc = sc.vectocolor(3).tolist()
    colors = ['k'] + vc

    ax = fig.add_subplot(gs[:2])
    this_dict = {
        'Status quo': msim_dict['Baseline'],
        '20% coverage, 2027': msim_dict['HPV-Faster 18%'],
        '35% coverage, 2027': msim_dict['HPV-Faster 35%'],
        '70% coverage, 2027': msim_dict['HPV-Faster 70%'],
    }

    cn = 0
    for slabel, mres in this_dict.items():
        ax = ut.plot_single(ax, mres, 'asr_cancer_incidence', si, ei, color=colors[cn], label=slabel)
        cn += 1
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('ASR cervical cancer incidence, 2025-2100\nOne-time mass screen, treat & vaccinate 20-50yos')
    ax.legend(loc="upper right", frameon=False)

    # Plot age of causal infection
    ac_colors = sc.gridcolors(3)

    ######################################################
    # Top row: dwelltimes and age of causal
    ######################################################
    ac_df = sc.loadobj('results/age_causal_infection.obj')
    ax = fig.add_subplot(gs[2])
    sns.boxplot(
        x="Health event", y="Age", data=ac_df, ax=ax,
        showfliers=False, palette=ac_colors, hue='Health event', hue_order=['Causal\ninfection', 'HSIL', 'Cancer']
    )
    ax.set_title(f'Age distribution\nof key health events')
    ax.set_xlabel('')
    ax.set_ylim([0, 100])

    fig.tight_layout()
    fig_name = 'figures/fig2_vx.png'
    sc.savefig(fig_name, dpi=100)

    return


# %% Run as a script
if __name__ == '__main__':

    # Load scenarios and construct figure
    plot_fig2()

    msim_dict = sc.loadobj('results/st_scens.obj')
    mbase = msim_dict['Baseline']
    m70 = msim_dict['HPV-Faster 70%']
    fi = sc.findinds(mbase.year, 2025)[0]
    elim_year = sc.findfirst(msim_dict['Baseline']['asr_cancer_incidence'][fi:]<4, die=False)
    if elim_year is None:
        elim_year = len(mbase.year) - 1
    else:
        elim_year = sc.findfirst(msim_dict['Baseline']['asr_cancer_incidence'][fi:]<4)+fi
    print(f'Elim year: {mbase.year[elim_year]}')
    print(f'Cancers averted: {mbase.cancers[fi:].sum()-m70.cancers[fi:].sum()}')

    # Print quantiles of age of causal infection
    ac_df = sc.loadobj('results/age_causal_infection.obj')
    lb = np.quantile(ac_df.loc[ac_df["Health event"] == 'Causal\ninfection']["Age"], 0.9)
    med = np.quantile(ac_df.loc[ac_df["Health event"] == 'Causal\ninfection']["Age"], 0.5)
    ub = np.quantile(ac_df.loc[ac_df["Health event"] == 'Causal\ninfection']["Age"], 0.1)
    print(f'90% of causal infections occur between ages {lb:.1f} and {ub:.1f}, median {med:.1f}')

    # As above but for HSIL
    lb = np.quantile(ac_df.loc[ac_df["Health event"] == 'HSIL']["Age"], 0.9)
    med = np.quantile(ac_df.loc[ac_df["Health event"] == 'HSIL']["Age"], 0.5)
    ub = np.quantile(ac_df.loc[ac_df["Health event"] == 'HSIL']["Age"], 0.1)
    print(f'90% of HSIL occur between ages {lb:.1f} and {ub:.1f}, median {med:.1f}')