"""
Run scenarios (v3 / Starsim).

Two modes:
  python run_scenarios.py --run-sim   # run sims + save plot-ready CSVs (VM)
  python run_scenarios.py             # re-extract CSVs from an existing st_scens.obj

MIGRATION NOTES (v2 -> v3)
--------------------------
* `hpv.MultiSim` is gone; scenarios are just built + run as a list of sims.
* v2 read flat sim results (`sim.results['asr_cancer_incidence']`,
  `['cancer_incidence_with_hiv']`, ...). In v3 those are not native results, so
  each sim carries a `RwandaReport` analyzer (see run_sim.py) that produces the
  v2-faithful annual ASR + HIV-stratified cancer incidence/counts; total cancers
  and cancer deaths are aggregated to annual from the scale-correct
  `sim.results.hpvtotal` flows. Intervention product counts come from the v3
  intervention result objects.
* Multi-seed uncertainty: value = mean across seeds; low/high = min/max.

Reduced scale (see run_sim.py defaults): n_agents, dt, start, stop,
ms_agent_ratio are shared with the patched-v2.3.1 reference runner so the two
engines are compared on an identical config.
"""
import argparse
import os

os.environ.update(
    OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1', MKL_NUM_THREADS='1',
)

import numpy as np
import pandas as pd
import sciris as sc

import run_sim as rs
from interventions import make_st, make_st_older, make_mv_intvs


TS_METRICS = ['asr_cancer_incidence', 'cancer_incidence_with_hiv', 'cancer_incidence_no_hiv']
CUM_METRICS_BOUNDED = ['cancers', 'cancers_with_hiv', 'cancers_no_hiv', 'cancer_deaths']
CUM_METRICS_UNBOUNDED = ['ablations', 'txvs', 'vaccinations', 'screens', 'excisions',
                         'leeps', 'cancer_treatments']
CUM_START_YEAR = 2025

# Settings
debug = 0
n_seeds = [3, 1][debug]

# v3 intervention-name -> (result_key, cum-metric) for product counts.
_PRODUCT_COUNTS = {
    'screening':      ('n_screened', 'screens'),
    'screening_older': ('n_screened', 'screens'),
    'ablation_rx':    ('new_cin_treated', 'ablations'),
    'ablation_older': ('new_cin_treated', 'ablations'),
    'excision_rx':    ('new_cin_treated', 'leeps'),
    'excision_older': ('new_cin_treated', 'excisions'),
    'radiation_rx':   ('new_cancer_treated', 'cancer_treatments'),
    'radiation_older': ('new_cancer_treated', 'cancer_treatments'),
    'txv':            ('new_cin_treated', 'txvs'),
    'campaign txvx':  ('new_cin_treated', 'txvs'),
    'routine_vx':     ('new_doses', 'vaccinations'),
    'mass_vax':       ('new_doses', 'vaccinations'),
}


# Interventions must stay within the (reduced) sim horizon.
END_YEAR = int(np.floor(rs.STOP)) - 1


# %% Scenario definitions (unchanged from v2 except explicit end_year)
def make_baselines():
    scendict = dict()
    scendict['No interventions'] = []
    scendict['Baseline'] = make_st(future_screen_cov=0.18, screen_change_year=2025, end_year=END_YEAR)
    return scendict


def make_campaign_scenarios():
    scendict = dict()
    age_range = [20, 50]
    for cov in [0.18, 0.35, 0.7]:
        scendict[f'Mass TxV 90/0, {int(cov*100)}%'] = make_mv_intvs(txv_pars='precin', campaign_coverage=cov)
        scendict[f'Mass TxV 50/90, {int(cov*100)}%'] = make_mv_intvs(txv_pars='cin', campaign_coverage=cov)
        mass_intvs = make_st_older(screen_cov=cov, age_range=age_range, start_year=2026)
        scendict[f'HPV-Faster {cov*100:.0f}%'] = mass_intvs
    return scendict


def make_st_scenarios():
    scendict = dict()
    start_year = 2026
    for cov_val in [.18, .35, .70]:
        scendict[f'S&T&T {cov_val*100:.0f}%'] = make_st(screen_change_year=start_year, future_screen_cov=cov_val, end_year=END_YEAR)
        scendict[f'S&T {cov_val*100:.0f}%'] = make_st(screen_change_year=start_year, future_screen_cov=cov_val,
                                                      tx_assigner_csv='tx_assigner_no_triage', end_year=END_YEAR)
        scendict[f'S&TxV&T&T {cov_val*100:.0f}%'] = make_st(screen_change_year=start_year, future_screen_cov=cov_val,
                                                           txv_pars='precin', txv=True, end_year=END_YEAR)
        scendict[f'S&TxV {cov_val*100:.0f}%'] = make_st(screen_change_year=start_year, future_screen_cov=cov_val,
                                                       txv_pars='cin', txv=True, end_year=END_YEAR)
    return scendict


# Representative subset for reduced-scale verification (clean single cascades;
# no txv-immunity ambiguity or multi-cascade product overlap).
REDUCED_SUBSET = ['No interventions', 'Baseline', 'S&T&T 70%', 'S&T 70%']


def build_scenarios(subset=None):
    scenarios = sc.mergedicts(make_baselines(), make_st_scenarios(), make_campaign_scenarios())
    if subset is not None:
        scenarios = {k: scenarios[k] for k in subset}
    return scenarios


# %% Run + reduce
def _annual_years(sim):
    rep = next(a for a in sim.analyzers.values() if isinstance(a, rs.RwandaReport))
    return rep.years


def run_one(name, interventions, seed):
    add_vax = name != 'No interventions'
    sim = rs.make_sim(interventions=list(interventions), add_vax=add_vax, add_st=False, seed=seed)
    sim.label = name
    sim.run(verbose=0)
    return sim


def extract_sim(sim):
    """Pull the per-metric annual arrays (v2-faithful) from one finished sim."""
    rep = next(a for a in sim.analyzers.values() if isinstance(a, rs.RwandaReport))
    tab = rep.annual_table()
    years = tab['year']
    out = dict(year=years)
    for m in TS_METRICS:
        out[m] = tab[m]
    # Total cancers / deaths aggregated to annual from scale-correct hpvtotal flows.
    out['cancers'] = rs.annual_from_timevec(sim, 'new_cancers', years)
    out['cancer_deaths'] = rs.annual_from_timevec(sim, 'new_cancer_deaths', years)
    out['cancers_with_hiv'] = tab['cancers_with_hiv']
    out['cancers_no_hiv'] = tab['cancers_no_hiv']
    # Intervention product counts -> annual arrays.
    for m in CUM_METRICS_UNBOUNDED:
        out[m] = np.zeros_like(years)
    tvy = np.floor(np.asarray(sim.results.timevec.years)).astype(int)
    for iname, (rkey, cummetric) in _PRODUCT_COUNTS.items():
        intv = sim.interventions.get(iname) if hasattr(sim.interventions, 'get') else None
        if intv is None and iname in sim.interventions:
            intv = sim.interventions[iname]
        if intv is None or rkey not in getattr(intv, 'results', {}):
            continue
        vals = np.asarray(intv.results[rkey])
        annual = np.array([float(np.sum(vals[tvy == int(y)])) for y in years])
        out[cummetric] = out[cummetric] + annual
    return out


def process(scenarios):
    """Run all seeds for all scenarios; reduce to per-scenario mean/min/max arrays."""
    msim_dict = sc.objdict()
    for name, interventions in scenarios.items():
        per_seed = []
        for seed in range(n_seeds):
            sim = run_one(name, interventions, seed)
            per_seed.append(extract_sim(sim))
            print(f'  [{name}] seed {seed} done', flush=True)
        years = per_seed[0]['year']
        mres = sc.objdict(year=years)
        keys = TS_METRICS + CUM_METRICS_BOUNDED + CUM_METRICS_UNBOUNDED
        for k in keys:
            stack = np.vstack([s[k] for s in per_seed])
            mres[k] = np.nanmean(stack, axis=0)
            mres[f'{k}_low'] = np.nanmin(stack, axis=0)
            mres[f'{k}_high'] = np.nanmax(stack, axis=0)
        msim_dict[name] = mres
    return msim_dict


def save_csvs(msim_dict, resfolder='results'):
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


# %% Run as a script
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-sim', action='store_true')
    parser.add_argument('--resfolder', default='results/v3.0_baseline')
    parser.add_argument('--full', action='store_true', help='Run all 23 scenarios (heavy)')
    parser.add_argument('--seeds', type=int, default=n_seeds)
    args = parser.parse_args()
    n_seeds = args.seeds

    T = sc.timer()
    subset = None if args.full else REDUCED_SUBSET
    scenarios = build_scenarios(subset=subset)

    if args.run_sim:
        print(f'Running {len(scenarios)} scenarios x {n_seeds} seeds '
              f'(n_agents={rs.N_AGENTS}, ms={rs.MS_AGENT_RATIO}, {rs.START}-{rs.STOP})')
        msim_dict = process(scenarios)
        os.makedirs(args.resfolder, exist_ok=True)
        sc.saveobj(f'{args.resfolder}/st_scens.obj', msim_dict)
    else:
        msim_dict = sc.loadobj(f'{args.resfolder}/st_scens.obj')

    save_csvs(msim_dict, resfolder=args.resfolder)
    print(f'Saved scens_timeseries.csv + scens_cumulative.csv to {args.resfolder}/')
    T.toc('Done')
