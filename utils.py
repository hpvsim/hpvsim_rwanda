'''
Utilities
'''

import numpy as np
import pandas as pd
import sciris as sc


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


# ---------- Scenario CSV loaders (produced by run_scenarios.py) ----------

def load_scens(resfolder='results'):
    ts = pd.read_csv(f'{resfolder}/scens_timeseries.csv')
    cum = pd.read_csv(f'{resfolder}/scens_cumulative.csv')
    return ts, cum


def get_ts(ts_df, scenario, metric):
    sub = ts_df[(ts_df.scenario == scenario) & (ts_df.metric == metric)].sort_values('year')
    return sub.year.values, sub.value.values, sub.low.values, sub.high.values


def get_cum(cum_df, scenario, metric):
    row = cum_df[(cum_df.scenario == scenario) & (cum_df.metric == metric)].iloc[0]
    return float(row.value), float(row.low), float(row.high)


def plot_ts(ax, ts_df, scenario, metric, start_year, end_year,
            color, ls='-', label=None, smooth=True, add_bounds=True):
    years, best, low, high = get_ts(ts_df, scenario, metric)
    mask = (years >= start_year) & (years <= end_year)
    years, best, low, high = years[mask], best[mask], low[mask], high[mask]

    if smooth:
        best = np.convolve(best, np.ones(5), 'valid') / 5
        low = np.convolve(low, np.ones(5), 'valid') / 5
        high = np.convolve(high, np.ones(5), 'valid') / 5
        years = years[4:]

    ax.plot(years, best, color=color, label=label, ls=ls)

    if metric == 'asr_cancer_incidence':
        below = np.where(best < 4)[0]
        if len(below):
            print(f'{label} elim year: {years[below[0]]}')
        else:
            print(f'{label} not eliminated')

    if add_bounds:
        ax.fill_between(years, low, high, alpha=0.1, color=color)
    ax.axhline(4, color='k', ls='--', lw=0.5)
    return ax
