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


# ---------- Scenario CSV loaders (produced by run_scenarios.py) ----------

def load_scens(resfolder='results'):
    """Return (timeseries, cumulative, paired, per_sim) DataFrames.

    Paired holds median/low/high of (baseline_i - scenario_i) per sim.
    Per_sim holds the 2025-2100 sum for every (scenario, sim, metric) triple.
    """
    ts   = pd.read_csv(f'{resfolder}/scens_timeseries.csv')
    cum  = pd.read_csv(f'{resfolder}/scens_cumulative.csv')
    paired  = pd.read_csv(f'{resfolder}/scens_paired.csv')
    per_sim = pd.read_csv(f'{resfolder}/scens_per_sim.csv')
    return ts, cum, paired, per_sim


def get_ts(ts_df, scenario, metric):
    sub = ts_df[(ts_df.scenario == scenario) & (ts_df.metric == metric)].sort_values('year')
    return sub.year.values, sub.value.values, sub.low.values, sub.high.values


def get_cum(cum_df, scenario, metric):
    row = cum_df[(cum_df.scenario == scenario) & (cum_df.metric == metric)].iloc[0]
    return float(row.value), float(row.low), float(row.high)


def get_paired(paired_df, scenario, metric):
    """Median/lo/hi of paired diff (baseline_sum - scenario_sum) per sim."""
    row = paired_df[(paired_df.scenario == scenario) & (paired_df.metric == metric)].iloc[0]
    return float(row.value), float(row.low), float(row.high)


def yerr(med, lo, hi):
    """asymmetric error-bar spec for matplotlib ax.bar(yerr=...)."""
    import numpy as np
    med = np.asarray(med, dtype=float)
    lo  = np.asarray(lo,  dtype=float)
    hi  = np.asarray(hi,  dtype=float)
    return np.vstack([med - lo, hi - med])


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
