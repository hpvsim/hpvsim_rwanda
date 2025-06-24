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
    # plot_calib(calib)

    ut.set_font(16)
    fig, axes = pl.subplots(2, 3, figsize=(12, 10))
    axes = axes.flatten()
    prev_col = '#5f5cd2'
    canc_col = '#c1981d'
    ms = 80
    gen_cols = sc.gridcolors(4)

    pn = 0

    # Pull out the results
    analyzer_results = calib.analyzer_results
    sim_results = calib.sim_results
    extra_sim_results = calib.extra_sim_results

    ####################
    # Cancers by age
    ####################
    for rn, res in enumerate(calib.analyzer_results[0].keys()):

        # Data
        ax = axes[pn]
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

        sns.lineplot(ax=ax, x='bins', y='values', data=modeldf, color=canc_col, errorbar=('pi', 95))
        ax.scatter(x, best, marker='d', s=ms, color='k')

        ax.set_ylim(bottom=0)
        ax.set_xticks(x, age_labels)
        ax.set_ylabel('')
        ax.set_xlabel('Age')
        ax.set_title('Cancers by age, 2020')
        pn += 1

    # Cancers by genotype
    rkeys = ['cancerous_genotype_dist']
    rlabels = ['Cancers by genotype']
    for ai, rkey in enumerate(rkeys):
        ax = axes[pn]

        # Plot data
        datadf = calib.target_data[-1]
        ydata = datadf.value.values
        x = np.arange(len(ydata))

        # Extract model results
        bins = []
        values = []
        for run_num, run in enumerate(sim_results):
            bins += x.tolist()
            if sc.isnumber(run[rkey]):
                values += sc.promotetolist(run[rkey])
            else:
                values += run[rkey].tolist()
        modeldf = pd.DataFrame({'bins': bins, 'values': values})

        # Plot model
        sns.boxplot(ax=ax, x='bins', y='values', data=modeldf, palette=gen_cols, showfliers=False)
        ax.scatter(x, ydata, color='k', marker='d', s=ms)

        ax.set_ylim([0, 1])
        ax.set_xticks(np.arange(4), ['16', '18', 'Hi5', 'OHR'])
        ax.set_ylabel('')
        ax.set_xlabel('')
        ax.set_title(rlabels[ai])
        pn += 1

    # New infections by age and sex
    simres = sc.dcp(sim.results)
    years = simres['year']
    year_ind = sc.findinds(years, 1985)[0]
    rsa_df = pd.read_csv('data/RSA_data.csv').set_index('Unnamed: 0').T
    title_dict = dict(
        female_hiv_prevalence='HIV prevalence (%), females 15+',
        hiv_incidence='HIV incidence',
        art_coverage='ART coverage (%)',
    )
    years = years[year_ind:]
    fig, axes = pl.subplots(1, 3, figsize=(10, 4))
    to_plot = ['female_hiv_prevalence', 'art_coverage', ['cancer_incidence_with_hiv', 'cancer_incidence_no_hiv']]
    for iv, ax in enumerate(axes.flatten()):
        val = to_plot[iv]
        if isinstance(val, list):
            for val_to_plot in val:
                label = 'HIV+' if 'with_hiv' in val_to_plot else 'HIV-'
                result = simres[val_to_plot][year_ind:]
                result = np.convolve(list(result), np.ones(5), "valid") / 5
                ax.plot(years[4:], result, label=label)
            ax.legend()
            ax.set_title('Cancer incidence (per 100k)')
        else:

            result = simres[val][year_ind:]
            if iv == 0:
                ax.plot(years, 100 * result, label='HPVsim')
            else:
                ax.plot(years, 100 * result)
            thembisa_val_lb = f'{val}_lb'
            thembisa_val_ub = f'{val}_ub'
            if iv == 0:
                ax.scatter(years, 100 * rsa_df[thembisa_val_lb][:-10], marker='_', label='Thembisa,\n95% uncertainty',
                           color='grey')
                ax.scatter(years, 100 * rsa_df[thembisa_val_ub][:-10], marker='_', color='grey')
            else:
                ax.scatter(years, 100 * rsa_df[thembisa_val_lb][:-10], marker='_', color='grey')
                ax.scatter(years, 100 * rsa_df[thembisa_val_ub][:-10], marker='_', color='grey')
            ax.set_title(title_dict[val])
            if iv == 0:
                ax.legend(title='Source')
        ax.set_ylim(bottom=0)
        sc.SIticks(ax)
    fig.tight_layout()
    fig_name = f'figures/hiv_fit_{location}.png'
    sc.savefig(fig_name, dpi=100)


    fig.tight_layout()
    pl.savefig(f"figures/figS2.png", dpi=300)




    # axes[0].plot(years[year_ind:], extra_sim_results['cancers_with_hiv'][year_ind:], label='HIV+')
    # axes[0].plot(years[year_ind:], extra_sim_results['cancers_no_hiv'][year_ind:], label='HIV-')
    # axes[0].plot(years[year_ind:], extra_sim_results['cancers'][year_ind:], label='Total')
    # axes[0].set_title(f'Cancers over time')
    # axes[0].legend()
    # axes[1].plot(calib.sim.pars['age_bin_edges'][:-1],
    #                  extra_sim_results['cancer_incidence_by_age_with_hiv'][:, -2], label='HIV+')
    # axes[1].plot(calib.sim.pars['age_bin_edges'][:-1],
    #                  extra_sim_results['cancer_incidence_by_age_no_hiv'][:, -2],
    #                  label='HIV-')
    # axes[1].legend()
    #
    # axes[2].plot(years[year_ind:], extra_sim_results['asr_cancer_incidence'][year_ind:])


    T.toc('Done')
 