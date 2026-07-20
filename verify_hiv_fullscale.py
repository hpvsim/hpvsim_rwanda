"""Full-scale HIV-HPV verification for the v3 Rwanda migration.

Runs the calibrated Rwanda co-infection sim (incidence-driven HIV, no cancer
interventions) at a fuller scale than the earlier reduced-scale check so the
sparse HIV+ cancer stratum is resolved:

  start=1960 (HPV burn-in), n_agents>=15k, ms_agent_ratio=5, several seeds.

Reports, vs the published Rwanda targets:
  - adult (15-49) HIV prevalence peak (~5%);
  - HIV+ vs HIV- cervical-cancer incidence per 100k, POOLED across seeds
    (sum weighted cancers / sum weighted female-years) over a registry window
    (HIV+ ~33, HIV- ~13).

Everything is scale-weighted by people.scale (grow-multiscale correct); the
sparse HIV+ stratum uses the POOLED estimator, never a per-seed rate average.

Usage:
    .venv/Scripts/python.exe verify_hiv_fullscale.py [n_agents] [n_seeds] [ms] [effects_json]
"""
import os
os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                  NUMEXPR_NUM_THREADS='1', MKL_NUM_THREADS='1')
import sys
import json

import numpy as np
import sciris as sc
import starsim as ss

import run_sim as rs

# --- config (overridable on the CLI) ---------------------------------------
N_AGENTS = int(sys.argv[1]) if len(sys.argv) > 1 else 15_000
N_SEEDS = int(sys.argv[2]) if len(sys.argv) > 2 else 6
MS = int(sys.argv[3]) if len(sys.argv) > 3 else 5
EFFECTS = json.loads(sys.argv[4]) if len(sys.argv) > 4 else None  # override RWANDA_HIV_EFFECTS
START = 1960
STOP = 2020
DT = 0.25

# Registry comparison window for pooled cancer incidence.
WIN_LO, WIN_HI = 2010, 2019


_AGE_EDGES = [25, 35, 45, 55, 200]        # registry by-age bins (calibrate_rwanda)
_AGE_LABELS = ['25', '35', '45', '55']
# by-age registry targets (2017), order matches _AGE_LABELS
TGT_BYAGE = dict(neg=[3.13, 11.67, 14.5, 12.0], pos=[15.0, 76.0, 80.0, 30.0])


class AdultHIVPrev(ss.Analyzer):
    """Record scale-weighted adult (15-49) HIV prevalence per year."""

    def init_pre(self, sim):
        super().init_pre(sim)
        self.rows = []  # (year, w_pos, w_tot)

    def step(self):
        sim = self.sim
        p = sim.people
        yr = int(np.floor(sim.t.now('year')))
        age = p.age.values
        w = p.scale.values
        alive = p.alive.values
        adult = alive & (age >= 15) & (age < 50)
        pos = adult & sim.diseases.hiv.infected.values
        self.rows.append((yr, float(np.sum(w[pos])), float(np.sum(w[adult]))))

    def annual(self):
        out = {}
        agg = {}
        for yr, wp, wt in self.rows:
            a = agg.setdefault(yr, [0.0, 0.0])
            a[0] += wp
            a[1] += wt
        for yr, (wp, wt) in agg.items():
            out[yr] = wp / wt if wt else 0.0
        return out


class ByAgeHIVProbe(ss.Analyzer):
    """Scale-weighted female cancers + female-years by age bin & HIV status.

    Mirrors calibrate_rwanda._WProbe: weighted new cancers and weighted alive
    females per registry age bin, split by current HIV status, plus aggregates
    over all females, females 15+, and females 25+ (denominator sensitivity)."""

    def init_pre(self, sim):
        from hpvsim.hpv import HPV
        self.hpv = [d for d in sim.diseases.values() if isinstance(d, HPV)]
        super().init_pre(sim)
        n = len(sim.t.timevec)
        nb = len(_AGE_LABELS)
        self.yr = np.floor(np.asarray(sim.t.timevec, float)).astype(int)[:n]
        self.canc = {'pos': np.zeros((nb, n)), 'neg': np.zeros((nb, n))}
        self.nf = {'pos': np.zeros((nb, n)), 'neg': np.zeros((nb, n))}
        # aggregate denominators over 3 female age-floors: all / 15+ / 25+
        self.cagg = {'pos': np.zeros(n), 'neg': np.zeros(n)}
        self.nf_all = {'pos': np.zeros(n), 'neg': np.zeros(n)}
        self.nf_15 = {'pos': np.zeros(n), 'neg': np.zeros(n)}
        self.nf_25 = {'pos': np.zeros(n), 'neg': np.zeros(n)}

    def step(self):
        ti = self.sim.ti
        p = self.sim.people
        w = p.scale.values
        al = p.alive.values
        fem = p.female.values & al
        age = p.age.values
        pos = self.sim.diseases.hiv.infected.values
        newc = np.zeros(al.shape, bool)
        for m in self.hpv:
            newc |= (m.cancerous.values & (m.ti_cancerous.values == ti))
        for status, hivmask in (('pos', pos), ('neg', ~pos)):
            fs = fem & hivmask
            self.cagg[status][ti] = (w * (newc & fs)).sum()
            self.nf_all[status][ti] = (w * fs).sum()
            self.nf_15[status][ti] = (w * (fs & (age >= 15))).sum()
            self.nf_25[status][ti] = (w * (fs & (age >= 25))).sum()
            for bi, (lo, hi) in enumerate(zip(_AGE_EDGES[:-1], _AGE_EDGES[1:])):
                ab = (age >= lo) & (age < hi)
                self.canc[status][bi, ti] = (w * (newc & fs & ab)).sum()
                self.nf[status][bi, ti] = (w * (fs & ab)).sum()

    def window(self, lo, hi):
        """Pooled numerators/denominators over calendar-year window [lo,hi]."""
        m = (self.yr >= lo) & (self.yr <= hi)
        out = {}
        for status in ('pos', 'neg'):
            out[f'canc_{status}'] = self.canc[status][:, m].sum(axis=1)
            out[f'nf_{status}'] = self.nf[status][:, m].sum(axis=1)
            out[f'cagg_{status}'] = float(self.cagg[status][m].sum())
            out[f'nf_all_{status}'] = float(self.nf_all[status][m].sum())
            out[f'nf_15_{status}'] = float(self.nf_15[status][m].sum())
            out[f'nf_25_{status}'] = float(self.nf_25[status][m].sum())
        return out


def _run_one(seed, n_agents, ms, effects):
    if effects is not None:
        rs.rc.RWANDA_HIV_EFFECTS = effects  # override in the worker process
    sim = rs.make_sim(add_vax=False, add_st=False, seed=seed,
                      n_agents=n_agents, start=START, stop=STOP, dt=DT,
                      ms_agent_ratio=ms,
                      analyzers=[AdultHIVPrev(), ByAgeHIVProbe()])
    sim.run()
    ahp = next(a for a in sim.analyzers.values() if isinstance(a, AdultHIVPrev))
    probe = next(a for a in sim.analyzers.values() if isinstance(a, ByAgeHIVProbe))
    out = probe.window(WIN_LO, WIN_HI)
    out['adult_hiv_prev'] = ahp.annual()
    return out


def _pool_agg(res_list, status, denom_key):
    """Pooled aggregate incidence per 100k for a HIV status and denominator."""
    num = sum(r[f'cagg_{status}'] for r in res_list)
    den = sum(r[denom_key + f'_{status}'] for r in res_list)
    return (num / den * 1e5 if den else float('nan')), num, den


def _pool_byage(res_list, status):
    """Pooled by-age incidence per 100k (array over _AGE_LABELS)."""
    num = sum(r[f'canc_{status}'] for r in res_list)
    den = sum(r[f'nf_{status}'] for r in res_list)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(den > 0, num / den * 1e5, np.nan), num, den


def report(res_list, tag=''):
    # adult HIV prevalence
    all_years = sorted({y for r in res_list for y in r['adult_hiv_prev']})
    prev_mean = {y: float(np.nanmean([r['adult_hiv_prev'].get(y, np.nan)
                                      for r in res_list])) for y in all_years}
    pk_y = max(prev_mean, key=prev_mean.get)
    print(f'\n=== adult (15-49) HIV prevalence, mean over seeds {tag} ===')
    for y in range(1990, 2021, 5):
        if y in prev_mean:
            print(f'  {y}: {prev_mean[y]*100:.2f}%')
    print(f'  PEAK: {prev_mean[pk_y]*100:.2f}% at {pk_y}   (target ~5%)')

    # by-age HIV-stratified cancer incidence (the trustworthy, age-consistent points)
    print(f'\n=== by-age cervical-cancer incidence /100k, {WIN_LO}-{WIN_HI} '
          f'(pooled, scale-weighted) {tag} ===')
    print(f'{"age":>5} {"HIV- mod":>9} {"HIV- tgt":>9} {"HIV+ mod":>9} {"HIV+ tgt":>9}')
    rn, _, _ = _pool_byage(res_list, 'neg')
    rp, _, _ = _pool_byage(res_list, 'pos')
    for i, lab in enumerate(_AGE_LABELS):
        print(f'{lab:>5} {rn[i]:>9.1f} {TGT_BYAGE["neg"][i]:>9.1f} '
              f'{rp[i]:>9.1f} {TGT_BYAGE["pos"][i]:>9.1f}')

    # aggregate incidence under 3 denominator definitions
    print(f'\n=== aggregate incidence /100k, {WIN_LO}-{WIN_HI} '
          f'(denominator sensitivity) {tag} ===')
    print(f'{"denom":>10} {"HIV- mod":>9} {"HIV+ mod":>9} {"RR":>6}  (tgt HIV- 13.1, HIV+ 33, RR 2.5)')
    for dk, dl in [('nf_all', 'all-fem'), ('nf_15', 'fem 15+'), ('nf_25', 'fem 25+')]:
        hn, _, _ = _pool_agg(res_list, 'neg', dk)
        hp, _, _ = _pool_agg(res_list, 'pos', dk)
        rr = hp / hn if hn else float('nan')
        print(f'{dl:>10} {hn:>9.1f} {hp:>9.1f} {rr:>6.2f}')
    # cancers for resolution context
    cp = sum(r['cagg_pos'] for r in res_list)
    cn = sum(r['cagg_neg'] for r in res_list)
    print(f'  [weighted cancers pooled: HIV+={cp:.0f}, HIV-={cn:.0f}]')


if __name__ == '__main__':
    T = sc.timer()
    import hpvsim as hpv
    print(f'hpvsim {hpv.__version__} {hpv.__file__}')
    print(f'config: n_agents={N_AGENTS} seeds={N_SEEDS} ms={MS} '
          f'{START}-{STOP} dt={DT}')
    print(f'effects override: {EFFECTS}')
    print(f'RWANDA_HIV_EFFECTS in use: {EFFECTS or rs.rc.RWANDA_HIV_EFFECTS}')

    res_list = sc.parallelize(
        _run_one,
        iterkwargs=dict(seed=list(range(N_SEEDS))),
        kwargs=dict(n_agents=N_AGENTS, ms=MS, effects=EFFECTS),
        serial=False,
    )
    report(res_list)
    T.toc('verify done')
