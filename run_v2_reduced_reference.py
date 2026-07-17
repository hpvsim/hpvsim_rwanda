"""Patched-v2.3.1 REFERENCE runner (reduced scale) for the v3 migration check.

Runs the SAME reduced config as the v3 run_scenarios.py (n_agents, dt, start,
stop/end, ms_agent_ratio, seeds, scenario subset) on the frozen/patched v2.3.1
engine, and writes the identical CSV schema (scens_timeseries.csv +
scens_cumulative.csv) so compare_baselines.py can put the two side by side.

v2.3.1 already exposes the native results the repo compares on
(asr_cancer_incidence, cancer_incidence_with_hiv/no_hiv, cancers,
cancers_with_hiv/no_hiv, cancer_deaths) and intervention n_products_used.
"""
import os
os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                  NUMEXPR_NUM_THREADS='1', MKL_NUM_THREADS='1')

import numpy as np
import pandas as pd
import sciris as sc
import hpvsim as hpv

import run_sim as rs
from interventions import make_st, make_st_older, make_mv_intvs

# --- Match the v3 reduced config ---
N_AGENTS = 5000
DT = 0.25
START = 1975
END = 2051
MS = 3
N_SEEDS = 2

TS_METRICS = ['asr_cancer_incidence', 'cancer_incidence_with_hiv', 'cancer_incidence_no_hiv']
CUM_METRICS_BOUNDED = ['cancers', 'cancers_with_hiv', 'cancers_no_hiv', 'cancer_deaths']
CUM_METRICS_UNBOUNDED = ['ablations', 'txvs', 'vaccinations', 'screens', 'excisions',
                         'leeps', 'cancer_treatments']
CUM_START_YEAR = 2025

SUBSET = ['No interventions', 'Baseline', 'S&T&T 70%', 'S&T 70%']

_PROGRAMS = {
    'mass_vax': 'vaccinations', 'screening': 'screens', 'ablation': 'ablations',
    'excision': 'leeps', 'radiation': 'cancer_treatments', 'txv': 'txvs',
    'campaign txvx': 'txvs', 'ablation_older': 'ablations',
    'excision_older': 'excisions', 'radiation_older': 'cancer_treatments',
}


END_YEAR = END - 1  # interventions must stay within the reduced horizon


def scenarios():
    d = dict()
    d['No interventions'] = []
    d['Baseline'] = make_st(future_screen_cov=0.18, screen_change_year=2025, end_year=END_YEAR)
    d['S&T&T 70%'] = make_st(screen_change_year=2026, future_screen_cov=0.70, end_year=END_YEAR)
    d['S&T 70%'] = make_st(screen_change_year=2026, future_screen_cov=0.70,
                           tx_assigner_csv='tx_assigner_no_triage', end_year=END_YEAR)
    return {k: d[k] for k in SUBSET}


def build(name, interventions, seed):
    add_vax = name != 'No interventions'
    sim = rs.make_sim(add_st=False, add_vax=add_vax, interventions=interventions,
                      end=END, seed=seed)
    # Override to the reduced config (v2 initializes at run-time).
    sim.pars['n_agents'] = N_AGENTS
    sim.pars['dt'] = DT
    sim.pars['start'] = START
    sim.pars['ms_agent_ratio'] = MS
    sim.pars['verbose'] = 0
    sim.label = name
    return sim


def extract(sim):
    r = sim.results
    years = np.asarray(r['year'][:] if hasattr(r['year'], '__getitem__') else r['year'])
    out = dict(year=np.asarray(years, dtype=float))
    for m in TS_METRICS + CUM_METRICS_BOUNDED:
        out[m] = np.asarray(r[m][:])
    for m in CUM_METRICS_UNBOUNDED:
        out[m] = np.zeros_like(out['year'])
    for iname, metric in _PROGRAMS.items():
        iv = sim.get_intervention(iname, die=False)
        if iv is not None and hasattr(iv, 'n_products_used'):
            out[metric] = out[metric] + np.asarray(iv.n_products_used.values)
    return out


def run():
    scens = scenarios()
    msim_dict = sc.objdict()
    for name, interventions in scens.items():
        per_seed = []
        for seed in range(N_SEEDS):
            sim = build(name, interventions, seed)
            sim.run()
            sim.shrink()
            per_seed.append(extract(sim))
            print(f'  [{name}] seed {seed} done', flush=True)
        years = per_seed[0]['year']
        mres = sc.objdict(year=years)
        for k in TS_METRICS + CUM_METRICS_BOUNDED + CUM_METRICS_UNBOUNDED:
            stack = np.vstack([s[k] for s in per_seed])
            mres[k] = np.nanmean(stack, axis=0)
            mres[f'{k}_low'] = np.nanmin(stack, axis=0)
            mres[f'{k}_high'] = np.nanmax(stack, axis=0)
        msim_dict[name] = mres
    return msim_dict


def save_csvs(msim_dict, resfolder):
    os.makedirs(resfolder, exist_ok=True)
    ts_rows = []
    for scen, mres in msim_dict.items():
        years = np.asarray(mres.year)
        for metric in TS_METRICS:
            for yi, yr in enumerate(years):
                ts_rows.append({'scenario': scen, 'year': float(yr), 'metric': metric,
                                'value': float(mres[metric][yi]),
                                'low': float(mres[f'{metric}_low'][yi]),
                                'high': float(mres[f'{metric}_high'][yi])})
    pd.DataFrame(ts_rows).to_csv(f'{resfolder}/scens_timeseries.csv', index=False)

    cum_rows = []
    for scen, mres in msim_dict.items():
        years = np.asarray(mres.year)
        fi_arr = np.where(years == CUM_START_YEAR)[0]
        fi = int(fi_arr[0]) if len(fi_arr) else 0
        for metric in CUM_METRICS_BOUNDED:
            cum_rows.append({'scenario': scen, 'metric': metric,
                             'value': float(np.nansum(mres[metric][fi:])),
                             'low': float(np.nansum(mres[f'{metric}_low'][fi:])),
                             'high': float(np.nansum(mres[f'{metric}_high'][fi:]))})
        for metric in CUM_METRICS_UNBOUNDED:
            cum_rows.append({'scenario': scen, 'metric': metric,
                             'value': float(np.nansum(mres[metric][fi:])),
                             'low': np.nan, 'high': np.nan})
    pd.DataFrame(cum_rows).to_csv(f'{resfolder}/scens_cumulative.csv', index=False)


if __name__ == '__main__':
    T = sc.timer()
    print('hpvsim', hpv.__version__, hpv.__file__)
    resfolder = 'results/v2.3.1_baseline'
    md = run()
    sc.saveobj(f'{resfolder}/st_scens.obj', md) if os.path.isdir(resfolder) or os.makedirs(resfolder) is None else None
    save_csvs(md, resfolder)
    print(f'Saved to {resfolder}/')
    T.toc('Done')
