"""
Plot the calibration
"""


# Standard imports
import numpy as np  # Add this at the top if not already imported
import sciris as sc
import pylab as pl
import hpvsim as hpv

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
    fig, axes = pl.subplot(2, 3, figsize=(12, 10))
    axes = axes.flatten()
    pn = 0

    # Plot the calibration results


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

    fig.show()
    fig.savefig(f'cancers_hiv_calib.png')


    T.toc('Done')
 