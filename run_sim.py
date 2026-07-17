"""
Define the HPVsim (v3) simulation for Rwanda (HIV-HPV co-infection).

MIGRATION NOTE (v2 -> v3 / Starsim)
-----------------------------------
v2 built the Rwanda sim with a single `pars` dict passing `model_hiv=True`,
`hiv_datafile`, `art_datafile` and an incidence-based HIV module. v3 rebuilds
HIV on Starsim: HIV is a disease module (`hpv.HIV`), the epidemic is imposed by
the incidence-driven importer (`hpv.hiv_incidence_import`, v2-faithful), ART is
the coverage shortcut (`hpv.hiv_art`), and the HIV->HPV effect is applied by the
auto-wired `hpv_hiv_connector` (with the calibrated Rwanda effect strengths).

The full calibrated Rwanda HIV-HPV configuration (genotype pars, network,
init_hpv_dist, HIV effect strengths, severity locus) is the CANONICAL builder in
`tests/regression/rwanda_calib.py` in the hpvsim source repo. We import those
helpers here rather than duplicate them, so this repo always tracks the finalized
calibration. `make_sim` replicates `rwanda_calib.build_rwanda_sim` but adds the
screening / TxV / vaccination interventions and the reporting analyzer.
"""
import os
import sys

import numpy as np
import sciris as sc
import starsim as ss
import hpvsim as hpv

# --- Locate the canonical v3 Rwanda calibration builders -------------------
# These live in the hpvsim source repo under tests/regression and are the
# single source of truth for the finalized Rwanda HIV-HPV calibration. They
# must NOT be modified or copied; we import them. Override the location with
# HPVSIM_RWANDA_REF if the source repo is elsewhere.
_REF = os.environ.get(
    'HPVSIM_RWANDA_REF',
    r'C:\Users\ryanhu\PycharmProjects\hpvsim_claudecontrol\tests\regression',
)
if _REF not in sys.path:
    sys.path.insert(0, _REF)
import rwanda_calib as rc  # noqa: E402
from hpvsim.cross_genotype import CrossImmunity  # noqa: E402
from hpvsim.hiv import hpv_hiv_connector  # noqa: E402

from interventions import make_st, make_vx  # noqa: E402


# %% Reduced-scale defaults (shared by run_scenarios.py and the v2 reference)
# Modest ms_agent_ratio so HIV+ cancer (sparse) is still resolved, small pop +
# few seeds so each scenario runs in a few minutes.
N_AGENTS = 5000
DT = 0.25
START = 1975
STOP = 2051          # +1 over the last fully-covered annual year (2050)
MS_AGENT_RATIO = 3


# %% Reporting analyzer -----------------------------------------------------
# WHO standard-population weights (same age edges as the v2 standard_pop).
_ASR_EDGES = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65,
                       70, 75, 80, 85, 100], dtype=float)
_ASR_WEIGHTS = np.array([.12, .10, .09, .09, .08, .08, .06, .06, .06, .06,
                         .05, .04, .04, .03, .02, .01, .005, .005, 0])[:-1]


class RwandaReport(ss.Analyzer):
    """Annual, scale-weighted, v2-faithful cancer reporting for Rwanda.

    Per calendar year produces (matching the v2 result definitions this repo
    compares against):
      - asr_cancer_incidence: WHO-standardized cervical-cancer incidence.
        numerator = full-year new cancers by age (accumulated over all dt
        sub-steps); denominator = year-end alive non-cancerous females by age;
        x1e5; dotted with WHO weights. (v2's built-in AgeResults captures only
        the final dt sub-step, undercounting ~4x at dt<1 -- see M10.)
      - cancer_incidence_with_hiv / _no_hiv: crude per-100k incidence among
        alive HIV+/HIV- women (v2: cancers_[no_]hiv / n_females_[no_]hiv_alive).
      - cancers_with_hiv / cancers_no_hiv: annual scale-weighted cancer counts.
    All histograms/counts use people.scale weights (grow-multiscale correct).
    """

    def __init__(self, **kw):
        super().__init__(**kw)
        self.edges = _ASR_EDGES
        self.weights = _ASR_WEIGHTS
        self.nbins = len(self.edges) - 1

    def init_pre(self, sim):
        super().init_pre(sim)
        from hpvsim.hpv import HPV
        self.hpv_modules = [d for d in sim.diseases.values() if isinstance(d, HPV)]
        tvy = np.floor(np.asarray(sim.results.timevec.years)).astype(int)
        self.years = sorted(int(y) for y in np.unique(tvy)
                            if int(y) >= int(np.ceil(sim.t.start.years))
                            and int(y) < int(np.floor(sim.t.stop.years)))
        z = {y: 0.0 for y in self.years}
        self.num_age = {y: np.zeros(self.nbins) for y in self.years}   # ASR numerator
        self.denom_age = {y: None for y in self.years}                 # ASR denominator
        self.num_hiv = dict(z); self.num_nohiv = dict(z)               # cancer counts
        self.den_hiv = {y: None for y in self.years}
        self.den_nohiv = {y: None for y in self.years}
        # final tick index of each requested year
        self._year_end_ti = {}
        for y in self.years:
            ticks = np.where((tvy >= y) & (tvy < y + 1))[0]
            if len(ticks):
                self._year_end_ti[int(ticks[-1])] = y

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
        canc = np.zeros_like(alive)
        for m in self.hpv_modules:
            new_c |= (m.ti_cancerous.values == ti) & m.cancerous.values
            canc |= m.cancerous.values
        if yr in self.num_age:
            fmask = new_c & alive & female
            self.num_age[yr] += np.histogram(ages[fmask], self.edges, weights=w[fmask])[0]
            self.num_hiv[yr] += float(np.sum(w[fmask & hiv_pos]))
            self.num_nohiv[yr] += float(np.sum(w[fmask & ~hiv_pos]))
        if ti in self._year_end_ti:
            y = self._year_end_ti[ti]
            at_risk = alive & female & ~canc
            self.denom_age[y] = np.histogram(ages[at_risk], self.edges, weights=w[at_risk])[0]
            # v2 HIV denom = all alive HIV+/HIV- females (not restricted to non-cancerous)
            af = alive & female
            self.den_hiv[y] = float(np.sum(w[af & hiv_pos]))
            self.den_nohiv[y] = float(np.sum(w[af & ~hiv_pos]))

    def annual_table(self):
        """Return dict of year-indexed arrays for the reported metrics."""
        years = np.array(self.years, dtype=float)
        asr = np.full(len(years), np.nan)
        inc_hiv = np.full(len(years), np.nan)
        inc_nohiv = np.full(len(years), np.nan)
        canc_hiv = np.zeros(len(years))
        canc_nohiv = np.zeros(len(years))
        for i, y in enumerate(self.years):
            d = self.denom_age[y]
            if d is not None:
                asi = np.divide(self.num_age[y], d, out=np.zeros(self.nbins), where=d > 0) * 1e5
                asr[i] = float(np.dot(asi, self.weights))
            if self.den_hiv[y]:
                inc_hiv[i] = self.num_hiv[y] / self.den_hiv[y] * 1e5
            if self.den_nohiv[y]:
                inc_nohiv[i] = self.num_nohiv[y] / self.den_nohiv[y] * 1e5
            canc_hiv[i] = self.num_hiv[y]
            canc_nohiv[i] = self.num_nohiv[y]
        return dict(year=years, asr_cancer_incidence=asr,
                    cancer_incidence_with_hiv=inc_hiv,
                    cancer_incidence_no_hiv=inc_nohiv,
                    cancers_with_hiv=canc_hiv, cancers_no_hiv=canc_nohiv)


def annual_from_timevec(sim, result_key, years):
    """Aggregate a per-timestep hpvtotal flow result to annual sums."""
    tvy = np.floor(np.asarray(sim.results.timevec.years)).astype(int)
    vals = np.asarray(sim.results.hpvtotal[result_key])
    out = np.zeros(len(years))
    for i, y in enumerate(years):
        out[i] = float(np.sum(vals[tvy == int(y)]))
    return out


# %% Simulation creation ----------------------------------------------------
def make_sim(interventions=None, analyzers=None, seed=0, add_vax=True,
             add_st=False, n_agents=N_AGENTS, start=START, stop=STOP, dt=DT,
             ms_agent_ratio=MS_AGENT_RATIO, end=None):
    """Build the calibrated Rwanda HIV-HPV v3 sim, with interventions/analyzers.

    Mirrors `rwanda_calib.build_rwanda_sim` (incidence-driven HIV, v2-faithful)
    but injects the screening/TxV/vaccination interventions and the RwandaReport
    analyzer. `end` is accepted as a v2-compatibility alias for `stop`.
    """
    if end is not None:
        stop = end
    # Interventions must stay within [start, stop); last usable year is stop-1.
    last_year = int(np.floor(stop)) - 1
    interventions = sc.autolist(interventions)
    if add_vax:
        interventions += make_vx(end_year=last_year)
    if add_st:
        interventions += make_st(end_year=last_year)

    analyzers = sc.autolist(analyzers)
    analyzers += RwandaReport()

    connectors = [
        CrossImmunity(rel_sev_loc=rc.REL_SEV_LOC),
        hpv_hiv_connector(effects=rc.RWANDA_HIV_EFFECTS),
    ]
    hiv = hpv.HIV.from_location('rwanda', beta_m2f=0.0, init_prev_data=0.0)
    hiv_intvs = [
        hpv.hiv_incidence_import.from_location('rwanda'),
        hpv.hiv_art.from_location('rwanda'),
    ]

    sim = hpv.Sim(
        location='rwanda',
        rand_seed=seed,
        n_agents=n_agents,
        start=start,
        stop=stop,
        dt=dt,
        ms_agent_ratio=ms_agent_ratio,
        genotypes=rc.GENOTYPES,
        genotype_pars=rc.rwanda_genotype_pars(),
        init_hpv_dist=rc.RWANDA_INIT_HPV_DIST,
        networks=[rc.make_rwanda_network()],
        connectors=connectors,
        diseases=[hiv],
        interventions=list(interventions) + hiv_intvs,
        analyzers=list(analyzers),
        verbose=0,
    )
    return sim


def run_sim(interventions=None, analyzers=None, seed=0, verbose=0,
            do_save=False, end=None, add_vax=True, add_st=True):
    sim = make_sim(interventions=interventions, analyzers=analyzers, seed=seed,
                   add_vax=add_vax, add_st=add_st, end=end)
    sim.label = f'Sim--{seed}'
    sim.run(verbose=verbose)
    if do_save:
        sim.save('results/rwanda.sim')
    return sim


# %% Run as a script
if __name__ == '__main__':
    T = sc.timer()
    print('hpvsim', hpv.__version__, hpv.__file__)
    sim = run_sim(add_vax=True, add_st=True, seed=0)
    rep = next(a for a in sim.analyzers.values() if isinstance(a, RwandaReport))
    tab = rep.annual_table()
    for y in [2020, 2030, 2040]:
        idx = np.where(tab['year'] == y)[0]
        if len(idx):
            i = int(idx[0])
            print(f'{y}: ASR={tab["asr_cancer_incidence"][i]:.1f} '
                  f'inc_hiv={tab["cancer_incidence_with_hiv"][i]:.1f} '
                  f'inc_nohiv={tab["cancer_incidence_no_hiv"][i]:.1f}')
    T.toc('Done')
