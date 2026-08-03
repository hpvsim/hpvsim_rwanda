"""
Gap-fill driver: produce the v3 MODEL CSVs for the figS2 calibration diagnostic.

The v3 Rwanda calibration is fixed (tests/regression/rwanda_calib.py, imported by
run_sim.make_sim); there is no saved v3 hpv.Calibration object, so we cannot use
run_calibration.save_figS2_csvs. Instead we run a small multiseed set of the
calibrated natural-history sim (no interventions) and extract the same diagnostic
metrics that save_figS2_csvs writes, matching its exact CSV schemas so
plot_figS2_calib.py renders unchanged.

Outputs -> results/_v3gapS2/:
  figS2_cancers.csv                    (cancers by age, 2020; 16 bins)
  figS2_cancer_incidence_no_hiv.csv    (incidence by age, 2017, HIV-; 4 bins)
  figS2_cancer_incidence_with_hiv.csv  (incidence by age, 2017, HIV+; 4 bins)
  figS2_precin_genotype_dist.csv       (share of precin by genotype, 2020; 4 bins)
  figS2_cancerous_genotype_dist.csv    (share of cancers by genotype, 2020; 4 bins)
  figS2_timeseries.csv                 (asr/inc/HIV metrics, 2000-2025)
  figS2_target_*.csv                   (copied verbatim from v2.2.6_baseline)

Run in the FOREGROUND. Reduced scale: n_agents=15000, ms_agent_ratio=3, n_seeds=3.
"""
import os
os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                  NUMEXPR_NUM_THREADS='1', MKL_NUM_THREADS='1')
import shutil
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss

import run_sim as rs

# ---- Config --------------------------------------------------------------
N_AGENTS = 15000
MS = 3
N_SEEDS = 5
START = 1975
STOP = 2051

# Single-year cancer events at this scale are sub-unit Poisson noise, so the
# by-age SHAPE and the RATE metrics are accumulated over multi-year windows
# (the age distribution / incidence rate are ~stationary over these windows).
CANC_WIN = np.arange(2015, 2026)   # cancers-by-age: report per-year-mean counts
INC_WIN = np.arange(2010, 2020)    # incidence-by-age x HIV: person-year rate
GENO_YEAR = 2020
TS_YEARS = np.arange(2000, 2026)

# cancers-by-age bins (edges match figS2_target_cancers age column: 0,15,20,..,85)
CANC_EDGES = np.array([0, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70,
                       75, 80, 85, 100], dtype=float)
# cancer-incidence-by-age bins (25-35, 35-45, 45-55, 55+)
INC_EDGES = np.array([25, 35, 45, 55, 100], dtype=float)
GENO_KEYS = ['hpv16', 'hpv18', 'hi5', 'ohr']

OUT = 'results/_v3gapS2'
BASELINE = 'results/v2.2.6_baseline'


def _bxp_stats(arr):
    """Boxplot stats across runs (verbatim from run_calibration._bxp_stats)."""
    arr = np.asarray(arr, dtype=float)
    q1, med, q3 = np.percentile(arr, [25, 50, 75])
    iqr = q3 - q1
    lo = float(arr[arr >= q1 - 1.5 * iqr].min())
    hi = float(arr[arr <= q3 + 1.5 * iqr].max())
    return dict(q1=float(q1), med=float(med), q3=float(q3), whislo=lo, whishi=hi)


# ---- Analyzer for by-age cancer metrics ----------------------------------
class FigS2Extra(ss.Analyzer):
    """Accumulate by-age cancer counts (2020) and by-age x HIV incidence (2017)."""

    def init_pre(self, sim):
        super().init_pre(sim)
        from hpvsim.hpv import HPV
        self.hpv_modules = [d for d in sim.diseases.values() if isinstance(d, HPV)]
        tvy = np.floor(np.asarray(sim.results.timevec.years)).astype(int)
        self.canc_by_age = np.zeros(len(CANC_EDGES) - 1)   # summed over CANC_WIN
        nb = len(INC_EDGES) - 1
        self.inc_num_hiv = np.zeros(nb)     # new cancers over INC_WIN (HIV+)
        self.inc_num_nohiv = np.zeros(nb)   # new cancers over INC_WIN (HIV-)
        self.inc_den_hiv = np.zeros(nb)     # person-years over INC_WIN (HIV+)
        self.inc_den_nohiv = np.zeros(nb)   # person-years over INC_WIN (HIV-)
        # final tick index of each INC_WIN year (for the person-year snapshot)
        self._inc_end_ti = set()
        for y in INC_WIN:
            ticks = np.where(tvy == int(y))[0]
            if len(ticks):
                self._inc_end_ti.add(int(ticks[-1]))
        self._canc_win = set(int(y) for y in CANC_WIN)
        self._inc_win = set(int(y) for y in INC_WIN)

    def step(self):
        sim = self.sim
        ti = sim.ti
        yr = int(np.floor(sim.t.now('year')))
        p = sim.people
        ages = p.age.values
        w = p.scale.values
        alive = p.alive.values
        female = p.female.values
        hiv_pos = sim.diseases.hiv.infected.values
        new_c = np.zeros_like(alive)
        for m in self.hpv_modules:
            new_c |= (m.ti_cancerous.values == ti) & m.cancerous.values
        fmask = new_c & alive & female
        if yr in self._canc_win:
            self.canc_by_age += np.histogram(ages[fmask], CANC_EDGES, weights=w[fmask])[0]
        if yr in self._inc_win:
            fh = fmask & hiv_pos
            fn = fmask & ~hiv_pos
            self.inc_num_hiv += np.histogram(ages[fh], INC_EDGES, weights=w[fh])[0]
            self.inc_num_nohiv += np.histogram(ages[fn], INC_EDGES, weights=w[fn])[0]
        if ti in self._inc_end_ti:
            af = alive & female
            afh = af & hiv_pos
            afn = af & ~hiv_pos
            self.inc_den_hiv += np.histogram(ages[afh], INC_EDGES, weights=w[afh])[0]
            self.inc_den_nohiv += np.histogram(ages[afn], INC_EDGES, weights=w[afn])[0]

    def extract(self):
        nb = len(INC_EDGES) - 1
        dh, dn = self.inc_den_hiv, self.inc_den_nohiv
        inc_hiv = np.divide(self.inc_num_hiv, dh, out=np.zeros(nb), where=dh > 0) * 1e5
        inc_nohiv = np.divide(self.inc_num_nohiv, dn, out=np.zeros(nb), where=dn > 0) * 1e5
        # per-year-mean cancer counts (window smoothing; keeps single-year y-scale)
        cancers = self.canc_by_age / len(CANC_WIN)
        return dict(cancers=cancers, inc_hiv=inc_hiv, inc_nohiv=inc_nohiv)


def _annual(vals, tvy, years, how):
    out = np.zeros(len(years))
    for i, y in enumerate(years):
        sel = vals[tvy == int(y)]
        if len(sel):
            out[i] = float(np.sum(sel)) if how == 'sum' else float(np.mean(sel))
    return out


def run_one(seed):
    sim = rs.make_sim(analyzers=[FigS2Extra()], seed=seed, add_vax=False,
                      add_st=False, n_agents=N_AGENTS, start=START, stop=STOP,
                      ms_agent_ratio=MS)
    sim.run(verbose=0)

    extra = next(a for a in sim.analyzers.values() if isinstance(a, FigS2Extra))
    byage = extra.extract()

    report = next(a for a in sim.analyzers.values()
                  if isinstance(a, rs.RwandaReport))
    tab = report.annual_table()

    res = sim.results
    tvy = np.floor(np.asarray(res.timevec.years)).astype(int)

    # genotype distributions at GENO_YEAR (last tick of the year)
    gi = np.where(tvy == GENO_YEAR)[0][-1]
    canc = np.array([float(res[g]['cum_cancers'][gi]) for g in GENO_KEYS])
    precin = np.array([float(res[g]['n_precin'][gi]) for g in GENO_KEYS])
    canc_dist = canc / canc.sum() if canc.sum() > 0 else np.zeros(4)
    precin_dist = precin / precin.sum() if precin.sum() > 0 else np.zeros(4)

    # time series
    ryears = np.asarray(tab['year'], dtype=int)

    def ts_from_report(key):
        out = np.full(len(TS_YEARS), np.nan)
        for i, y in enumerate(TS_YEARS):
            idx = np.where(ryears == int(y))[0]
            if len(idx):
                out[i] = tab[key][int(idx[0])]
        return out

    hiv = res['hiv']
    # HIV module runs on its own (finer) timestep -> use its own timevec.
    htvy = np.floor(np.asarray(hiv['timevec'].years)).astype(int)
    ts = {
        'asr_cancer_incidence': ts_from_report('asr_cancer_incidence'),
        'cancer_incidence_with_hiv': ts_from_report('cancer_incidence_with_hiv'),
        'cancer_incidence_no_hiv': ts_from_report('cancer_incidence_no_hiv'),
        'art_coverage': _annual(np.asarray(hiv['p_on_art']), htvy, TS_YEARS, 'mean'),
        'female_hiv_prevalence': _annual(np.asarray(hiv['prevalence_f']), htvy, TS_YEARS, 'mean'),
        'male_hiv_prevalence': _annual(np.asarray(hiv['prevalence_m']), htvy, TS_YEARS, 'mean'),
        'hiv_infections': _annual(np.asarray(hiv['new_infections']), htvy, TS_YEARS, 'sum'),
        'hiv_deaths': _annual(np.asarray(hiv['new_deaths']), htvy, TS_YEARS, 'sum'),
    }

    return dict(cancers=byage['cancers'], inc_hiv=byage['inc_hiv'],
                inc_nohiv=byage['inc_nohiv'], canc_dist=canc_dist,
                precin_dist=precin_dist, ts=ts)


def main():
    os.makedirs(OUT, exist_ok=True)
    T = sc.timer()
    runs = []
    for s in range(N_SEEDS):
        t0 = sc.timer()
        print(f'--- seed {s} ---', flush=True)
        runs.append(run_one(s))
        t0.toc(f'seed {s} done')

    # ---- by-age boxplot-stat CSVs ----
    def write_bxp(stack, fname):
        stack = np.array(stack)  # (nseeds, nbins)
        rows = [{'bin': bi, **_bxp_stats(stack[:, bi])}
                for bi in range(stack.shape[1])]
        pd.DataFrame(rows).to_csv(f'{OUT}/{fname}', index=False)

    write_bxp([r['cancers'] for r in runs], 'figS2_cancers.csv')
    write_bxp([r['inc_nohiv'] for r in runs], 'figS2_cancer_incidence_no_hiv.csv')
    write_bxp([r['inc_hiv'] for r in runs], 'figS2_cancer_incidence_with_hiv.csv')
    write_bxp([r['precin_dist'] for r in runs], 'figS2_precin_genotype_dist.csv')
    write_bxp([r['canc_dist'] for r in runs], 'figS2_cancerous_genotype_dist.csv')

    # ---- time series (med + pi95 across seeds) ----
    metrics = list(runs[0]['ts'].keys())
    rows = []
    for rkey in metrics:
        stacked = np.array([r['ts'][rkey] for r in runs])  # (nseeds, nyears)
        med = np.nanmedian(stacked, axis=0)
        lo = np.nanpercentile(stacked, 2.5, axis=0)
        hi = np.nanpercentile(stacked, 97.5, axis=0)
        for yi, yr in enumerate(TS_YEARS):
            rows.append({'year': int(yr), 'metric': rkey,
                         'med': float(med[yi]),
                         'pi95_low': float(lo[yi]),
                         'pi95_high': float(hi[yi])})
    pd.DataFrame(rows).to_csv(f'{OUT}/figS2_timeseries.csv', index=False)

    # ---- copy version-independent targets ----
    for f in os.listdir(BASELINE):
        if f.startswith('figS2_target_'):
            shutil.copyfile(f'{BASELINE}/{f}', f'{OUT}/{f}')

    # ---- quick diagnostic printout ----
    print('\n==== SUMMARY (median across seeds) ====')
    canc = np.median([r['cancers'] for r in runs], axis=0)
    print('cancers-by-age win-mean (16 bins):', np.round(canc, 2))
    print('inc_nohiv (25/35/45/55):', np.round(np.median([r['inc_nohiv'] for r in runs], axis=0), 1),
          ' target [3.13,11.67,14.5,12.0]')
    print('inc_hiv   (25/35/45/55):', np.round(np.median([r['inc_hiv'] for r in runs], axis=0), 1),
          ' target [15,76,80,30]')
    print('precin dist (16/18/hi5/ohr):', np.round(np.median([r['precin_dist'] for r in runs], axis=0), 2),
          ' target [0.21,0.11,0.37,0.31]')
    print('cancer dist (16/18/hi5/ohr):', np.round(np.median([r['canc_dist'] for r in runs], axis=0), 2),
          ' target [0.55,0.17,0.25,0.05]')
    asr = np.array([r['ts']['asr_cancer_incidence'] for r in runs])
    winmask = np.isin(TS_YEARS, CANC_WIN)
    asr_win = np.nanmean(asr[:, winmask], axis=1)
    print('ASR window-mean per seed:', np.round(asr_win, 1), ' target ~28.2 (2020 GLOBOCAN)')
    T.toc('ALL DONE')


if __name__ == '__main__':
    main()
