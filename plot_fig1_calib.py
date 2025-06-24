"""
Plot the calibration
"""


# Standard imports
import numpy as np  # Add this at the top if not already imported
import sciris as sc
import pylab as pl
import hpvsim as hpv
import pandas as pd
import seaborn as sns

# Imports from this repository
import utils as ut


def plot_calib(calib):
    return


# %% Run as a script
if __name__ == '__main__':

    T = sc.timer()

    calib = sc.loadobj('results/rwanda_calib.obj')
    hiv_df = pd.read_csv('data/rwanda_data.csv').set_index('Unnamed: 0').T
    hiv_df.index = hiv_df.index.astype(int)  # Convert index to integers

    # plot_calib(calib)

    ut.set_font(16)
    fig = pl.figure(layout="tight", figsize=(16, 10))

    gs = fig.add_gridspec(3, 1)
    gs0 = gs[0].subgridspec(1, 4)
    gs1 = gs[1].subgridspec(1, 4)
    gs2 = gs[2].subgridspec(1, 2)

    prev_col = '#5f5cd2'
    canc_col = sc.gridcolors(5)[3]
    ms = 80
    gen_cols = sc.gridcolors(4)
    pn = 0

    # Pull out the results
    analyzer_results = calib.analyzer_results
    sim_results = calib.sim_results
    extra_sim_results = calib.extra_sim_results

    ####################
    # Cancers by age (3 plots)
    ####################
    for rn, res in enumerate(calib.analyzer_results[0].keys()):

        # Data
        if rn == 0: ax = fig.add_subplot(gs0[:2])
        else: ax = fig.add_subplot(gs0[rn+1])

        datadf = calib.target_data[rn]
        if rn == 0:
            age_labels = ['0', '15', '20', '25', '30', '35', '40', '45', '50', '55', '60', '65', '70', '75', '80', '85']
        else:
            age_labels = ['25-35', '35-45', '45-55', '55+']
        x = np.arange(len(age_labels))
        best = datadf.value.values

        # Extract model results
        bins = []
        values = []
        for run_num, run in enumerate(analyzer_results):
            bins += x.tolist()
            yearkey = 2020 if res == 'cancers' else 2017
            values += list(run[res][yearkey])
        modeldf = pd.DataFrame({'bins': bins, 'values': values})

        sns.boxplot(ax=ax, x='bins', y='values', data=modeldf, color=canc_col, showfliers=False)
        # Adjust alpha for the boxes
        for patch in ax.patches:
            patch.set_alpha(0.7)  # Set alpha to 0.7 for the boxes
        ax.scatter(x, best, marker='d', s=ms, color='k')

        ax.set_ylim(bottom=0)
        ax.set_xticks(x, age_labels)
        ax.set_ylabel('')
        ax.set_xlabel('Age')
        ax.set_title('Cancers by age, 2020')
        pn += 1

    ####################
    # Cancer incidence (3 lines on 1 plot)
    ####################
    ax = fig.add_subplot(gs1[:2])
    x = np.arange(1960, 2026)
    si = sc.findfirst(x, 2000)  # Start index for plotting
    x = x[si:]
    rkeys = ['asr_cancer_incidence', 'cancer_incidence_with_hiv',
             'cancer_incidence_no_hiv']
    rlabels = ['Total', 'HIV+', 'HIV-']
    for ai, rkey in enumerate(rkeys):
        bins = []
        values = []
        for run_num, run in enumerate(extra_sim_results):
            bins += x.tolist()
            values += run[rkey][si:].tolist()
        modeldf = pd.DataFrame({'bins': bins, 'values': values})
        sns.lineplot(ax=ax, x='bins', y='values', data=modeldf, errorbar=('pi', 95), label= rlabels[ai])

    datadf = calib.target_data[3]
    ax.scatter(datadf.year.values[0], datadf.value.values[0], marker='d', s=ms, color='k',label='GLOBOCAN')
    ax.legend(loc='upper left', frameon=False, fontsize=12)

    ####################
    # Cancers by genotype
    ####################
    ax = fig.add_subplot(gs1[2])
    # Plot data
    datadf = calib.target_data[-1]
    ydata = datadf.value.values
    x = np.arange(len(ydata))

    bins = []
    values = []
    for run_num, run in enumerate(sim_results):
        bins += x.tolist()
        values += run['cancerous_genotype_dist'].tolist()
    modeldf = pd.DataFrame({'bins': bins, 'values': values})
    sns.boxplot(ax=ax, x='bins', y='values', data=modeldf, hue='bins', palette=gen_cols, showfliers=False)
    ax.scatter(x, ydata, color='k', marker='d', s=ms)

    ax.set_ylim([0, 1])
    ax.set_xticks(np.arange(4), ['16', '18', 'Hi5', 'OHR'])
    ax.legend(loc='upper right', frameon=False, title='', labels=['16', '18', 'Hi5', 'OHR'], fontsize=12)
    ax.set_ylabel('')
    ax.set_xlabel('')
    ax.set_title('Cancers by genotype, 2020')

    ####################
    # CINs by genotype
    ####################
    ax = fig.add_subplot(gs1[3])

    ####################
    # ART coverage over time
    ####################
    ax = fig.add_subplot(gs2[0])
    x = np.arange(1960, 2026)
    si = sc.findfirst(x, 2010)  # Start index for plotting
    x = x[si:]
    rkey = 'art_coverage'
    rlabel = 'ART coverage'
    bins = []
    values = []
    for run_num, run in enumerate(extra_sim_results):
        bins += x.tolist()
        values += run[rkey][si:].tolist()
    modeldf = pd.DataFrame({'bins': bins, 'values': values})
    sns.lineplot(ax=ax, x='bins', y='values', data=modeldf, color=prev_col, errorbar=('pi', 95))
    datadf = hiv_df.loc[(hiv_df.index >= 2000) & (hiv_df.index <= 2025)]
    sns.scatterplot(ax=ax, x=datadf.index, y=datadf[rkey], marker='d', s=ms, color='k', label=rlabel)
    ax.set_ylim(bottom=0)
    ax.set_ylabel('ART coverage (%)')

    ####################
    # HIV prevalence over time
    ####################
    ax = fig.add_subplot(gs2[1])
    rkeys = ['female_hiv_prevalence', 'male_hiv_prevalence']
    rlabels = ['F', 'M']
    colors = sc.gridcolors(len(rkeys))
    for ai, rkey in enumerate(rkeys):
        bins = []
        values = []
        for run_num, run in enumerate(extra_sim_results):
            bins += x.tolist()
            values += run[rkey][si:].tolist()
        modeldf = pd.DataFrame({'bins': bins, 'values': values})
        sns.lineplot(ax=ax, x='bins', y='values', data=modeldf, color=colors[ai], errorbar=('pi', 95), label=rlabels[ai])
        sns.scatterplot(ax=ax, x=datadf.index, y=datadf[rkey], marker='d', s=ms, color=colors[ai])
    ax.legend(loc='upper left', frameon=False)

    fig.tight_layout()
    pl.savefig(f"figures/fig_calib.png", dpi=300)


    T.toc('Done')
 