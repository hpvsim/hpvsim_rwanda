'''
Utilities
'''

# Imports
import sciris as sc
import hpvsim as hpv
import numpy as np


def set_font(size=None, font='Libertinus Sans'):
    ''' Set a custom font '''
    sc.fonts(add=sc.thisdir(aspath=True) / 'assets' / 'LibertinusSans-Regular.otf')
    sc.options(font=font, fontsize=size)
    return


def shrink_calib(calib, n_results=100):
    plot_indices = calib.df.iloc[0:n_results, 0].values
    calib.analyzer_results = [calib.analyzer_results[i] for i in plot_indices]
    calib.sim_results = [calib.sim_results[i] for i in plot_indices]
    calib.extra_sim_results = [calib.extra_sim_results[i] for i in plot_indices]
    calib.target_data = calib.target_data
    calib.df = calib.df.iloc[0:n_results, ]
    return calib


def plot_single(ax, mres, to_plot, si, ei, color, ls='-', label=None, smooth=True):
    years = mres.year[si:ei]
    best = mres[to_plot][si:ei]
    low = mres[to_plot].low[si:ei]
    high = mres[to_plot].high[si:ei]

    if smooth:
        best = np.convolve(list(best), np.ones(5), "valid")/5
        low = np.convolve(list(low), np.ones(5), "valid")/5
        high = np.convolve(list(high), np.ones(5), "valid")/5
        years = years[4:]

    ax.plot(years, best, color=color, label=label, ls=ls)
    ax.fill_between(years, low, high, alpha=0.1, color=color)

    # Add horizontal line at 4
    ax.axhline(4, color='k', ls='--', lw=0.5)
    return ax


 