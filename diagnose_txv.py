"""Categorize the 50K cancers in S&TxV 70% by mode-of-failure.

Buckets:
  1. Never screened
  2. Screened but always HPV-DNA negative (false negatives at every screen)
  3. Positive at screen but never given TxV (pre-2030 or refused/LTFU)
  4. TxV given but failed to clear (subcategorized by state at last admin: precin vs cin)
  5. Effective TxV, then post-treatment infection progressed to cancer

Notes:
  - We only count cancers with onset ti in [2025-01-01, 2100-01-01].
  - "Effective TxV" = the tx.administer efficacy draw succeeded for at least
    one genotype on that dose.
  - Post-TxV reinfection is judged by: the cancer's causative infection's
    ti_infected > any successful TxV ti.
"""
import argparse
import collections
import numpy as np
import sciris as sc
import starsim as ss
import hpvsim as hpv
from hpvsim.utils import iter_hpv_modules

import run_scenarios as rsc
import run_sim as rs
import interventions as intv


# ---- Patched hpv.txvx.administer -------------------------------------------

_tracker = {
    'txv_precin_attempts':  collections.Counter(),
    'txv_cin_attempts':     collections.Counter(),
    'txv_precin_successes': collections.Counter(),
    'txv_cin_successes':    collections.Counter(),
    'ti_last_effective':    {},  # uid -> ti of latest successful clearance
    'ti_first_txv':         {},  # uid -> first ti of any TxV admin
}

_original_administer = hpv.txvx.administer


def _patched_administer(self, uids, return_format='dict'):
    # Snapshot pre-treatment state
    uids_np = np.asarray(uids)
    uid_set = set(int(u) for u in uids_np)
    was_precin = set()
    was_cin = set()
    for mod in iter_hpv_modules(self.sim):
        pre_uids = set(int(u) for u in np.asarray(mod.precin.uids))
        cin_uids = set(int(u) for u in np.asarray(mod.cin.uids))
        was_precin |= (pre_uids & uid_set)
        was_cin |= (cin_uids & uid_set)

    result = _original_administer(self, uids, return_format='dict')
    successful = set(int(u) for u in np.asarray(result['successful']))

    ti = self.sim.ti
    for u in uid_set:
        if u not in _tracker['ti_first_txv']:
            _tracker['ti_first_txv'][u] = ti
        # cin priority — if a woman was CIN on any genotype, count as CIN encounter
        if u in was_cin:
            _tracker['txv_cin_attempts'][u] += 1
            if u in successful:
                _tracker['txv_cin_successes'][u] += 1
        elif u in was_precin:
            _tracker['txv_precin_attempts'][u] += 1
            if u in successful:
                _tracker['txv_precin_successes'][u] += 1
        if u in successful:
            _tracker['ti_last_effective'][u] = ti

    return result if return_format == 'dict' else result['successful']


hpv.txvx.administer = _patched_administer


# ---- Analyzer to track screening flags -------------------------------------

class ScreenPositiveTracker(ss.Analyzer):
    """Track positive screens AND cancer transitions per uid."""
    def __init__(self):
        super().__init__(name='screen_pos_tracker')
        self.n_positive_screens = collections.Counter()  # uid -> count
        self.ti_first_positive = {}                       # uid -> ti
        self.cancer_uid_to_ti = {}                        # uid -> ti of first cancer
        self.cancer_uid_to_genotype = {}                  # uid -> genotype str
        self.cancer_uid_to_age = {}                       # uid -> age at cancer
        self.cancer_uid_to_hiv = {}                       # uid -> HIV+ at cancer
        self.cancer_uid_to_age_causal = {}                # uid -> age at causal infection
        self.cancer_uid_to_age_cin = {}                   # uid -> age at CIN onset
        self.cancer_uid_to_scale = {}                     # uid -> people.scale at transition

    def step(self):
        sim = self.sim
        # 1) Screen-positive tracking
        scr = sim.interventions.get('screening')
        if scr is not None:
            pos = scr.outcomes.get('positive', ss.uids())
            if len(pos):
                ti = sim.ti
                for u in np.asarray(pos):
                    u = int(u)
                    self.n_positive_screens[u] += 1
                    if u not in self.ti_first_positive:
                        self.ti_first_positive[u] = ti

        # 2) Cancer transitions — diff against prior-step cancerous set per module.
        if not hasattr(self, '_prev_cancerous'):
            self._prev_cancerous = {m.genotype: set() for m in iter_hpv_modules(sim)}
        hivm = sim.diseases.get('hiv') if hasattr(sim.diseases, 'get') else None
        dt_year = sim.t.dt_year
        for m in iter_hpv_modules(sim):
            curr = set(int(u) for u in np.asarray(m.cancerous.uids))
            new = curr - self._prev_cancerous[m.genotype]
            for u in new:
                if u not in self.cancer_uid_to_ti:
                    age_now = float(sim.people.age.raw[u])
                    ti_inf = float(m.ti_infected.raw[u])
                    ti_cin = float(m.ti_cin.raw[u])
                    self.cancer_uid_to_ti[u] = sim.ti
                    self.cancer_uid_to_genotype[u] = m.genotype
                    self.cancer_uid_to_age[u] = age_now
                    self.cancer_uid_to_hiv[u] = bool(hivm.infected.raw[u]) if hivm is not None else False
                    self.cancer_uid_to_age_causal[u] = (
                        age_now - (sim.ti - ti_inf) * dt_year if np.isfinite(ti_inf) else np.nan)
                    self.cancer_uid_to_age_cin[u] = (
                        age_now - (sim.ti - ti_cin) * dt_year if np.isfinite(ti_cin) else np.nan)
                    self.cancer_uid_to_scale[u] = float(sim.people.scale.raw[u])
            self._prev_cancerous[m.genotype] = curr


# ---- Run one sim and categorize --------------------------------------------

def run_one(top_par, seed_idx, end=2100):
    """Run one S&TxV 70% sim with the given (calibrated pars + seed)."""
    intvs = intv.make_st(future_screen_cov=0.7, txv=True, txv_pars='cin', end_year=end)
    sim = rs.make_sim(add_st=False, interventions=intvs,
                       analyzers=[ScreenPositiveTracker()],
                       stop=end, calib_pars=top_par)
    sim.run(verbose=0)
    # sim copies analyzers on init; grab the populated copy
    tracker = sim.analyzers['screen_pos_tracker']
    return sim, tracker


_ns_profile = {'ages': [], 'hiv': [], 'ti_cancer': [], 'alive_at_end': [],
               'age_in_2020': [], 'age_at_intv_start_or_death': []}


def categorize(sim, tracker):
    """Return {bucket: scaled count} for all cancers in [2025, 2100].
    Each cancer transition is weighted by people.scale × pop_scale so the
    total matches sim.results.all_hpv.new_cancers summed over the window.
    """
    _ns_profile['ages'].clear()
    _ns_profile['hiv'].clear()
    _ns_profile['ti_cancer'].clear()
    _ns_profile['alive_at_end'].clear()
    buckets = collections.Counter()
    pop_scale = float(sim.pars.pop_scale)
    hpv_modules = list(iter_hpv_modules(sim))

    scr = sim.interventions['screening']
    txv_intv = sim.interventions.get('txv')

    # ti bounds for cancer counting
    yearvec = sim.t.yearvec
    ti_2025 = int(np.searchsorted(yearvec, 2025.0))
    ti_2100 = int(np.searchsorted(yearvec, 2100.0))

    # Iterate over EVERY cancer transition recorded during the sim
    all_cancers = tracker.cancer_uid_to_ti  # uid -> ti_cancer
    print(f'Cancer transitions recorded: {len(all_cancers)}')

    for uid, ti_cancer in sorted(all_cancers.items()):
        # Filter to 2025-2100 window
        if ti_cancer < ti_2025 or ti_cancer >= ti_2100:
            continue
        # Genotype the cancer was on
        gt = tracker.cancer_uid_to_genotype[uid]
        cancerous_gts = [m for m in hpv_modules if m.genotype == gt]

        # Per-agent weight: matches sim.results.new_cancers scaling.
        w = tracker.cancer_uid_to_scale.get(uid, 1.0) * pop_scale

        # Bucket 1: never screened — profile these
        n_screens = float(scr.screens.raw[uid])
        if n_screens == 0:
            buckets['1_never_screened'] += w
            _ns_profile['ages'].append(tracker.cancer_uid_to_age.get(uid, np.nan))
            _ns_profile['hiv'].append(tracker.cancer_uid_to_hiv.get(uid, False))
            _ns_profile['ti_cancer'].append(ti_cancer)
            _ns_profile['alive_at_end'].append(bool(sim.people.alive.raw[uid]))
            continue

        # Bucket 2: screened but always negative
        n_pos = tracker.n_positive_screens.get(uid, 0)
        if n_pos == 0:
            buckets['2_screened_always_negative'] += w
            continue

        # Bucket 3: positive but never given TxV
        if txv_intv is None or not bool(txv_intv.tx_vaccinated.raw[uid]):
            buckets['3_positive_but_no_txv'] += w
            continue

        # Bucket 4: TxV given but failed. Subcategorize by state at last admin.
        n_precin_att = _tracker['txv_precin_attempts'].get(uid, 0)
        n_cin_att = _tracker['txv_cin_attempts'].get(uid, 0)
        n_precin_ok = _tracker['txv_precin_successes'].get(uid, 0)
        n_cin_ok = _tracker['txv_cin_successes'].get(uid, 0)
        n_total_success = n_precin_ok + n_cin_ok
        if n_total_success == 0:
            if n_cin_att > 0:
                buckets['4a_txv_failed_at_cin'] += w
            elif n_precin_att > 0:
                buckets['4b_txv_failed_at_precin'] += w
            else:
                buckets['4c_txv_admin_no_lesion'] += w
            continue

        # Bucket 5: effective TxV, then reinfection led to cancer
        ti_last_effective = _tracker['ti_last_effective'].get(uid, -1)
        ti_infections = [int(m.ti_infected.raw[uid]) for m in cancerous_gts
                         if not np.isnan(m.ti_infected.raw[uid])]
        latest_infection_ti = max(ti_infections) if ti_infections else -1
        if latest_infection_ti > ti_last_effective:
            buckets['5_effective_txv_then_reinfection'] += w
        else:
            buckets['6_effective_txv_but_current_infection_still_progressed'] += w

    return buckets


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--rep', type=int, default=0, help='calibration ensemble index (0-9)')
    ap.add_argument('--end', type=int, default=2100)
    args = ap.parse_args()

    top_pars = rsc._top_pars(10)
    top_par = top_pars[args.rep]
    print(f'Using calibration rep {args.rep} (pars: {list(top_par.items())[:3]}...)')

    T = sc.timer()
    sim, tracker = run_one(top_par, args.rep, end=args.end)
    T.toc('sim run')

    buckets = categorize(sim, tracker)
    total = sum(buckets.values())

    # Sanity check vs sim's own counter
    r = sim.results.all_hpv.new_cancers
    tv = np.asarray(r.timevec.years if hasattr(r.timevec, 'years') else r.timevec)
    mask = (tv >= 2025) & (tv < args.end)
    sim_total = float(np.asarray(r.values)[mask].sum())

    print()
    print(f'=== Cancer decomposition (S&TxV 70%, rep {args.rep}, 2025-{args.end}) ===')
    print(f'Tracker total (scaled): {total:,.0f}')
    print(f'sim.results.new_cancers window sum: {sim_total:,.0f}')
    for k in sorted(buckets):
        pct = 100 * buckets[k] / total if total else 0
        print(f'  {k:60s} {buckets[k]:>10,.0f}  ({pct:5.1f}%)')

    # Age-at-stage stratified by screening status
    print()
    print('=== Age at causal infection / CIN / cancer, by screening status ===')
    yearvec = sim.t.yearvec
    ti_2020 = int(np.searchsorted(yearvec, 2020.0))
    ti_2025 = int(np.searchsorted(yearvec, 2025.0))
    ti_2100 = int(np.searchsorted(yearvec, 2100.0))
    scr = sim.interventions['screening']

    for label, screened_flag in (('never-screened', 0), ('ever-screened', 1)):
        uids = []
        ages_causal = []; ages_cin = []; ages_cancer = []
        for uid, ti_c in tracker.cancer_uid_to_ti.items():
            if ti_c < ti_2025 or ti_c >= ti_2100:
                continue
            n_scr = float(scr.screens.raw[uid])
            is_scr = (n_scr > 0)
            if is_scr != bool(screened_flag):
                continue
            uids.append(uid)
            ages_causal.append(tracker.cancer_uid_to_age_causal.get(uid, np.nan))
            ages_cin.append(tracker.cancer_uid_to_age_cin.get(uid, np.nan))
            ages_cancer.append(tracker.cancer_uid_to_age.get(uid, np.nan))
        if not uids:
            continue
        ac = np.asarray(ages_causal)
        acin = np.asarray(ages_cin)
        acan = np.asarray(ages_cancer)
        print(f'\n  {label} (n={len(uids)}):')
        for name, arr in [('age at causal HPV', ac), ('age at CIN', acin), ('age at cancer', acan)]:
            valid = np.isfinite(arr)
            if valid.sum():
                v = arr[valid]
                print(f'    {name:22s} median={np.median(v):5.1f}  p10={np.percentile(v,10):5.1f}  p90={np.percentile(v,90):5.1f}')
        # Was she in 30-50 during the intervention era (i.e., age 30-50 at any ti >= ti_2020)?
        # Approximation: her age at cancer minus (current_ti - ti_cancer). But screening runs from 2020+.
        # Age at 2020 = age_at_cancer - (year_at_cancer - 2020)
        cancer_years_arr = np.array([yearvec[int(tracker.cancer_uid_to_ti[u])] for u in uids])
        age_at_2020 = acan - (cancer_years_arr - 2020)
        # Was she in [30, 50) at any point during 2020 - min(2100, death)?
        # Her age at 2020 is age_at_2020. If age_at_2020 < 50 AND age_at_2020 + years_alive_post_2020 > 30, she overlapped.
        # Simpler: overlap = (age_at_2020 < 50) AND (age_at_cancer > 30)  → she was <50 in 2020 and >30 at cancer, so had 30-50 overlap
        overlap = (age_at_2020 < 50) & (acan > 30)
        # Refine: age_at_2020 must be < 50 AND she must have reached 30 by 2020 or later (i.e., age_at_2020 >= 30, OR she reached 30 during intervention)
        # Actually: she was 30-50 during 2020+ iff there exists year Y in [2020, min(2100, year_at_cancer)] with 30 <= age_at_year(Y) < 50
        # age_at_year(Y) = age_at_cancer - (year_at_cancer - Y)
        # We need Y such that 30 <= acan - (year_at_cancer - Y) < 50 AND Y in [2020, year_at_cancer]
        # Solving: Y >= year_at_cancer - (acan - 30) AND Y < year_at_cancer - (acan - 50)
        # And Y in [2020, year_at_cancer]
        lo_y = np.maximum(cancer_years_arr - (acan - 30), 2020)
        hi_y = np.minimum(cancer_years_arr - (acan - 50), cancer_years_arr)
        window_overlap = (hi_y > lo_y) & (hi_y > 2020) & (lo_y < 2100)
        pct = 100 * window_overlap.sum() / len(uids)
        print(f'    ever in 30-50 window during 2020-2100: {window_overlap.sum()} ({pct:.1f}%)')
