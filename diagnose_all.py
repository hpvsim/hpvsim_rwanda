"""Multi-scenario cancer decomposition + birth-cohort analysis.

Runs each of the six flagship scenarios once per calibration rep (default 3),
records every cancer transition, and categorizes cancers by:
  * mode-of-failure bucket (never-screened, screened-negative, positive-no-ablate,
    ablate/txv failed, ablate/txv success then reinfection)
  * 10-year birth cohort
  * screening reach at any age

Writes CSVs to results/diagnostic/ so downstream analysis is fast.
"""
import argparse
import collections
import os
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import hpvsim as hpv
from hpvsim.utils import iter_hpv_modules

import run_scenarios as rsc
import run_sim as rs
import interventions as intv


# ---- Global patch trackers (per-run reset via .clear()) ----------------

TXV_TRK = {
    'ti_last_effective': {},  # uid -> ti of latest successful clearance
    'ti_first_txv':      {},  # uid -> first ti of any TxV admin
    'attempts':          collections.Counter(),
    'successes':         collections.Counter(),
    'was_precin':        collections.Counter(),  # attempts at precin state
    'was_cin':           collections.Counter(),  # attempts at cin state
}

ABL_TRK = {
    'attempts':          collections.Counter(),
    'successes':         collections.Counter(),
    'ti_last_effective': {},  # uid -> ti of most recent successful ablation
}

_original_txv_admin = hpv.txvx.administer
_original_ablation_admin = None  # set at import time below

def _reset_trackers():
    for d in (TXV_TRK, ABL_TRK):
        for k, v in d.items():
            v.clear() if hasattr(v, 'clear') else d.__setitem__(k, {})
    TXV_TRK['ti_last_effective'] = {}
    TXV_TRK['ti_first_txv'] = {}
    ABL_TRK['ti_last_effective'] = {}


def _patched_txv_admin(self, uids, return_format='dict'):
    uids_np = np.asarray(uids)
    uid_set = set(int(u) for u in uids_np)
    was_precin = set(); was_cin = set()
    for mod in iter_hpv_modules(self.sim):
        pre_uids = set(int(u) for u in np.asarray(mod.precin.uids))
        cin_uids = set(int(u) for u in np.asarray(mod.cin.uids))
        was_precin |= (pre_uids & uid_set)
        was_cin |= (cin_uids & uid_set)
    result = _original_txv_admin(self, uids, return_format='dict')
    successful = set(int(u) for u in np.asarray(result['successful']))
    ti = self.sim.ti
    for u in uid_set:
        TXV_TRK['ti_first_txv'].setdefault(u, ti)
        TXV_TRK['attempts'][u] += 1
        if u in was_cin:
            TXV_TRK['was_cin'][u] += 1
        elif u in was_precin:
            TXV_TRK['was_precin'][u] += 1
        if u in successful:
            TXV_TRK['successes'][u] += 1
            TXV_TRK['ti_last_effective'][u] = ti
    return result if return_format == 'dict' else result['successful']

hpv.txvx.administer = _patched_txv_admin


# Ablation tracker: intercept hpv.tx.administer for ablation products so we
# can score ablation attempts and successes per uid.
_original_tx_admin = hpv.tx.administer

def _patched_tx_admin(self, uids, return_format='dict'):
    """Track ablation attempts / successes. Only for products named 'ablation'.
    Skip if the concrete class is hpv.txvx (its .administer calls super, which
    would double-count if we tracked it here as well)."""
    uids_np = np.asarray(uids)
    result = _original_tx_admin(self, uids, return_format='dict')
    if isinstance(self, hpv.txvx):
        return result if return_format == 'dict' else result['successful']
    # Only score ablations here (the mass_older path uses a differently-named
    # module but same product base). We inspect the module_name.
    mn = getattr(self, 'module_name', getattr(self, 'name', ''))
    is_ablation = ('ablation' in mn.lower())
    if is_ablation:
        successful = set(int(u) for u in np.asarray(result['successful']))
        ti = self.sim.ti
        for u in uids_np:
            u = int(u)
            ABL_TRK['attempts'][u] += 1
            if u in successful:
                ABL_TRK['successes'][u] += 1
                ABL_TRK['ti_last_effective'][u] = ti
    return result if return_format == 'dict' else result['successful']

hpv.tx.administer = _patched_tx_admin


# ---- Analyzer: track cancers + screen positives ------------------------

class ScenTracker(ss.Analyzer):
    def __init__(self):
        super().__init__(name='scen_tracker')
        self.n_positive_screens = collections.Counter()
        self.n_screens = collections.Counter()
        self.ti_first_positive = {}
        self.ti_first_screen = {}
        self.cancer = {}  # uid -> dict of per-cancer attributes

    def step(self):
        sim = self.sim
        # Screen tracking
        for scr_name in ('screening', 'screening_older'):
            scr = sim.interventions.get(scr_name)
            if scr is None:
                continue
            # positives
            pos = scr.outcomes.get('positive', ss.uids()) if scr.outcomes else ss.uids()
            if len(pos):
                ti = sim.ti
                for u in np.asarray(pos):
                    u = int(u)
                    self.n_positive_screens[u] += 1
                    self.ti_first_positive.setdefault(u, ti)
            # negatives + inadequates count as screens too via ti_screened; use screens flag.
            # We use the intervention's own per-uid screens counter for total,
            # so we only need to record 'first screen ti' here.
            # ti_screened is BoolArr for latest ti; per-uid count in scr.screens.
            if hasattr(scr, 'ti_screened'):
                # Anyone whose ti_screened equals current ti got screened now.
                ti = sim.ti
                just = (scr.ti_screened == ti).uids
                for u in np.asarray(just):
                    u = int(u)
                    self.n_screens[u] += 1
                    self.ti_first_screen.setdefault(u, ti)

        # Cancer transition tracking
        if not hasattr(self, '_prev'):
            self._prev = {m.genotype: set() for m in iter_hpv_modules(sim)}
        hivm = sim.diseases.get('hiv') if hasattr(sim.diseases, 'get') else None
        dt_year = sim.t.dt_year
        for m in iter_hpv_modules(sim):
            curr = set(int(u) for u in np.asarray(m.cancerous.uids))
            new = curr - self._prev[m.genotype]
            for u in new:
                if u in self.cancer:
                    continue  # first cancer only per uid
                age_now = float(sim.people.age.raw[u])
                ti_inf = float(m.ti_infected.raw[u])
                ti_cin = float(m.ti_cin.raw[u])
                self.cancer[u] = {
                    'ti': sim.ti,
                    'genotype': m.genotype,
                    'age': age_now,
                    'hiv': bool(hivm.infected.raw[u]) if hivm is not None else False,
                    'age_causal': age_now - (sim.ti - ti_inf) * dt_year if np.isfinite(ti_inf) else np.nan,
                    'age_cin':    age_now - (sim.ti - ti_cin) * dt_year if np.isfinite(ti_cin) else np.nan,
                    'scale':      float(sim.people.scale.raw[u]),
                    'fine':       bool(sim.people.fine.raw[u]) if hasattr(sim.people, 'fine') else False,
                }
            self._prev[m.genotype] = curr


# ---- Run one scenario --------------------------------------------------

def build_scenario_intvs(name, end=2100):
    """Return (interventions list, name) for a named scenario."""
    if name == 'baseline':
        return intv.make_st(future_screen_cov=0.18, end_year=end)
    if name == 'sTT18':
        return intv.make_st(future_screen_cov=0.18, end_year=end)
    if name == 'sTT70':
        return intv.make_st(future_screen_cov=0.70, end_year=end)
    if name == 'sT70':
        return intv.make_st(future_screen_cov=0.70,
                            tx_assigner_csv='tx_assigner_no_triage', end_year=end)
    if name == 'sTxV70':
        return intv.make_st(future_screen_cov=0.70, txv=True, txv_pars='cin', end_year=end)
    if name == 'sTxV18':
        return intv.make_st(future_screen_cov=0.18, txv=True, txv_pars='cin', end_year=end)
    if name == 'hpvfaster70':
        return intv.make_st_older(screen_cov=0.70, age_range=[20, 50], end_year=end)
    if name == 'hpvfaster70_ablate_only':
        # HPV-Faster minus the mass adult vax.
        intvs = intv.make_st_older(screen_cov=0.70, age_range=[20, 50], end_year=end)
        return [i for i in intvs if getattr(i, 'name', '') != 'mass_vax']
    if name == 'hpvfaster70_vax_only':
        # Baseline S&T&T + a standalone mass adult prophylactic vax at 2027 for 20-50yo.
        base = intv.make_st(future_screen_cov=0.18, end_year=end)
        mass_vx = hpv.campaign_vx(
            name='mass_vax',
            product='nonavalent',
            age_range=[20, 50],
            prob=0.70,
            years=[2027],
        )
        return base + [mass_vx]
    if name == 'masstxv70':
        return intv.make_mv_intvs(campaign_coverage=0.70, txv_pars='cin', end_year=end)
    raise ValueError(f'unknown scenario {name}')


def build_normalized_scenarios(intv_start_year=2030, end_year=2100):
    """Return the dict of normalized scenario -> intv list. Thin wrapper
    around rsc.make_normalized_scenarios so this file owns nothing that
    interventions.py doesn't."""
    return rsc.make_normalized_scenarios(intv_start_year=intv_start_year,
                                         end_year=end_year)


def run_one(name, top_par, seed_idx, end=2100, intvs=None):
    """intvs, if supplied, overrides the scenario-name dispatch. Used by
    the normalized-set entrypoint to hand the pre-built interventions in
    directly."""
    _reset_trackers()
    if intvs is None:
        intvs = build_scenario_intvs(name, end=end)
    sim = rs.make_sim(add_st=False, interventions=intvs,
                      analyzers=[ScenTracker()],
                      stop=end, calib_pars=top_par)
    sim.run(verbose=0)
    tracker = sim.analyzers['scen_tracker']
    return sim, tracker


# ---- Categorization ----------------------------------------------------

def categorize(sim, tracker, txv_used=False, hpv_faster=False, end=2100,
               accounting_start=2025):
    """Return list of per-cancer rows for accounting_start..end window."""
    yearvec = sim.t.yearvec
    ti_start = int(np.searchsorted(yearvec, float(accounting_start)))
    ti_end = int(np.searchsorted(yearvec, float(end)))
    pop_scale = float(sim.pars.pop_scale)
    hpv_modules = list(iter_hpv_modules(sim))

    scr = sim.interventions.get('screening')
    scr_older = sim.interventions.get('screening_older')
    txv_intv = sim.interventions.get('txv') or sim.interventions.get('campaign_txvx')
    mass_vx = sim.interventions.get('mass_vax')

    rows = []
    for uid, info in tracker.cancer.items():
        ti_c = info['ti']
        if ti_c < ti_start or ti_c >= ti_end:
            continue
        # weight per sim.results.new_cancers scaling
        w = info['scale'] * pop_scale
        n_screens = int(scr.screens.raw[uid]) if scr is not None else 0
        n_older = 0
        if scr_older is not None:
            # scr_older is a campaign_screening — has screens FloatArr
            if hasattr(scr_older, 'screens'):
                n_older = int(scr_older.screens.raw[uid])
        n_pos = tracker.n_positive_screens.get(uid, 0)
        n_txv_att = TXV_TRK['attempts'].get(uid, 0)
        n_txv_ok = TXV_TRK['successes'].get(uid, 0)
        n_txv_precin = TXV_TRK['was_precin'].get(uid, 0)
        n_txv_cin = TXV_TRK['was_cin'].get(uid, 0)
        n_abl_att = ABL_TRK['attempts'].get(uid, 0)
        n_abl_ok = ABL_TRK['successes'].get(uid, 0)
        ever_scr = (n_screens + n_older) > 0
        ever_pos = n_pos > 0
        # For mass vax reach flag
        got_mass_vx = False
        if mass_vx is not None and hasattr(mass_vx, 'vaccinated'):
            got_mass_vx = bool(mass_vx.vaccinated.raw[uid])
        # For routine vax reach
        rvx = sim.interventions.get('routine_vx')
        got_routine_vx = False
        if rvx is not None and hasattr(rvx, 'vaccinated'):
            got_routine_vx = bool(rvx.vaccinated.raw[uid])

        # Bucket assignment
        if not ever_scr:
            bucket = '1_never_screened'
        elif not ever_pos:
            bucket = '2_screened_always_negative'
        elif n_txv_att == 0 and n_abl_att == 0:
            bucket = '3_positive_no_treatment'
        else:
            # Some treatment attempt occurred
            if n_txv_ok == 0 and n_abl_ok == 0:
                bucket = '4_treatment_failed'
            else:
                # Effective at some point; did reinfection follow?
                gt = info['genotype']
                cancerous_gts = [m for m in hpv_modules if m.genotype == gt]
                ti_inf_this = [int(m.ti_infected.raw[uid]) for m in cancerous_gts
                                if not np.isnan(m.ti_infected.raw[uid])]
                latest_infection = max(ti_inf_this) if ti_inf_this else -1
                latest_success = max(
                    TXV_TRK['ti_last_effective'].get(uid, -1),
                    ABL_TRK['ti_last_effective'].get(uid, -1),
                )
                if latest_infection > latest_success:
                    bucket = '5_treated_then_reinfected'
                else:
                    bucket = '6_treated_but_progressed_anyway'

        # Year of cancer, birth year
        year_c = float(yearvec[ti_c])
        birth_year = year_c - info['age']
        birth_cohort_10 = int(np.floor(birth_year / 10.0) * 10)
        age_in_2027 = info['age'] - (year_c - 2027)
        age_in_2030 = info['age'] - (year_c - 2030)

        rows.append({
            'uid': uid,
            'ti': ti_c,
            'year_cancer': year_c,
            'age_cancer': info['age'],
            'age_causal': info['age_causal'],
            'age_cin':    info['age_cin'],
            'genotype':   info['genotype'],
            'hiv':        info['hiv'],
            'weight':     w,
            'fine':       info['fine'],
            'birth_year': birth_year,
            'birth_cohort_10': birth_cohort_10,
            'age_in_2027': age_in_2027,
            'age_in_2030': age_in_2030,
            'n_screens':      n_screens,
            'n_screens_older': n_older,
            'n_positive':     n_pos,
            'n_txv_attempts': n_txv_att,
            'n_txv_success':  n_txv_ok,
            'n_txv_precin':   n_txv_precin,
            'n_txv_cin':      n_txv_cin,
            'n_abl_attempts': n_abl_att,
            'n_abl_success':  n_abl_ok,
            'got_mass_vx':    got_mass_vx,
            'got_routine_vx': got_routine_vx,
            'ever_screened':  ever_scr,
            'ever_positive':  ever_pos,
            'bucket':         bucket,
        })
    return rows


def sim_totals(sim, end=2100, accounting_start=2025):
    r = sim.results.all_hpv.new_cancers
    tv = np.asarray(r.timevec.years if hasattr(r.timevec, 'years') else r.timevec)
    mask = (tv >= float(accounting_start)) & (tv < float(end))
    return float(np.asarray(r.values)[mask].sum())


# ---- Main --------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--scenarios', nargs='+', default=None,
                    help='Scenario names. If omitted with --normalized, runs all '
                         'normalized scenarios; otherwise defaults to the flagship '
                         'short-name list.')
    ap.add_argument('--reps', type=int, default=3)
    ap.add_argument('--end', type=int, default=2100)
    ap.add_argument('--outdir', default='results/diagnostic')
    ap.add_argument('--parallel', action='store_true')
    ap.add_argument('--normalized', action='store_true',
                    help='Use the normalized scenario set (all interventions '
                         'start in --intv-start-year); accounting window '
                         'defaults to that year.')
    ap.add_argument('--intv-start-year', type=int, default=2030)
    ap.add_argument('--accounting-start', type=int, default=None,
                    help='Cumulative accounting start year. Defaults to 2025 '
                         '(default set) or --intv-start-year (normalized).')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    top_pars = rsc._top_pars(args.reps)

    # Resolve scenarios + accounting window
    if args.normalized:
        norm = build_normalized_scenarios(intv_start_year=args.intv_start_year,
                                          end_year=args.end)
        scen_names = args.scenarios if args.scenarios else list(norm.keys())
        prebuilt = {n: norm[n] for n in scen_names}
        accounting_start = args.accounting_start or args.intv_start_year
    else:
        scen_names = args.scenarios or [
            'sTT18', 'sTT70', 'sT70', 'sTxV70', 'sTxV18', 'hpvfaster70',
        ]
        prebuilt = None
        accounting_start = args.accounting_start or 2025

    print(f'Running {len(scen_names)} scenarios x {args.reps} reps  '
          f'(accounting {accounting_start}..{args.end})')

    all_rows = []
    summary_rows = []
    for scen in scen_names:
        for rep, tp in enumerate(top_pars):
            T = sc.timer()
            print(f'\n===== {scen} rep {rep} =====')
            intvs = prebuilt[scen] if prebuilt is not None else None
            sim, tracker = run_one(scen, dict(tp), rep, end=args.end, intvs=intvs)
            rows = categorize(sim, tracker,
                              txv_used=('TxV' in scen or scen in ('sTxV70', 'sTxV18', 'masstxv70')),
                              hpv_faster=('Faster' in scen or scen.startswith('hpvfaster')),
                              end=args.end,
                              accounting_start=accounting_start)
            for r in rows:
                r['scenario'] = scen
                r['rep'] = rep
            all_rows.extend(rows)
            total_tracked = sum(r['weight'] for r in rows)
            total_sim = sim_totals(sim, end=args.end, accounting_start=accounting_start)
            print(f'  Tracked: {total_tracked:,.0f}  vs sim: {total_sim:,.0f}  '
                  f'(delta {total_tracked - total_sim:+,.0f})')
            summary_rows.append({
                'scenario': scen, 'rep': rep,
                'tracked': total_tracked, 'sim_total': total_sim,
                'n_cancer_uids': len(rows),
                'accounting_start': accounting_start,
            })
            T.toc(f'{scen} rep {rep}')

    df = pd.DataFrame(all_rows)
    df.to_csv(f'{args.outdir}/cancer_rows.csv', index=False)
    pd.DataFrame(summary_rows).to_csv(f'{args.outdir}/summary.csv', index=False)
    print(f'\nWrote {len(df)} cancer rows to {args.outdir}/cancer_rows.csv')


if __name__ == '__main__':
    main()
