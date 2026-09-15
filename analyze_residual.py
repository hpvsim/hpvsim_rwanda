"""Residual cancer analysis: who is still getting cancer, and can we reach them?

Reads results/diagnostic_normalized/cancer_rows.csv (produced by
diagnose_all.py --normalized) and produces the tables that fill in
docs/residual_analysis.md.

Everything is medians across reps unless stated otherwise. Weights
account for people.scale x pop_scale so cancer counts match sim
totals (verified against summary.csv `sim_total`).
"""
import argparse
import os

import numpy as np
import pandas as pd


ROUTINE_VAX_START = 2011          # bivalent, age 11-12
INTV_START = 2030                 # normalized start year
SCREEN_AGE = (30, 50)             # routine screening window
MASS_ADULT_AGE = (20, 50)         # HPV-Faster / Mass TxV window


def load(path):
    df = pd.read_csv(path)
    df['birth_year'] = df['birth_year'].astype(float)
    df['birth_cohort_10'] = df['birth_cohort_10'].astype(int)
    return df


def scenario_totals(df):
    """Table 1 of the residual doc: cumulative cancers per scenario."""
    per_rep = df.groupby(['scenario', 'rep'])['weight'].sum().unstack('rep')
    out = pd.DataFrame({
        'median': per_rep.median(axis=1),
        'p10':    per_rep.quantile(0.10, axis=1),
        'p90':    per_rep.quantile(0.90, axis=1),
    }).round(0)

    if 'No interventions' in out.index:
        baseline = out.loc['No interventions', 'median']
        out['averted_vs_none'] = (baseline - out['median']).round(0)
    return out


def by_birth_cohort(df, scenario):
    """§2a: cancers by 10-year birth cohort for one scenario."""
    sub = df[df.scenario == scenario].copy()
    per_rep = (sub.groupby(['rep', 'birth_cohort_10'])['weight']
                  .sum().unstack('birth_cohort_10', fill_value=0))
    med = per_rep.median(axis=0).round(0)
    lo  = per_rep.quantile(0.10, axis=0).round(0)
    hi  = per_rep.quantile(0.90, axis=0).round(0)
    total = med.sum()
    out = pd.DataFrame({'median': med, 'p10': lo, 'p90': hi,
                        'share_pct': (med / total * 100).round(1)})
    out.index.name = 'birth_cohort_10'
    return out


def by_age_causal(df, scenario, bins=(0, 20, 30, 40, 50, 60, 200)):
    """§2b: cancers by age at causal infection."""
    sub = df[df.scenario == scenario].copy()
    sub['age_causal_bin'] = pd.cut(sub['age_causal'], bins=bins, right=False)
    per_rep = (sub.groupby(['rep', 'age_causal_bin'], observed=True)['weight']
                  .sum().unstack('age_causal_bin', fill_value=0))
    med = per_rep.median(axis=0).round(0)
    total = med.sum()
    return pd.DataFrame({'median': med,
                         'share_pct': (med / total * 100).round(1)})


def age_causal_novax_vs_baseline(df, bins=(0, 20, 30, 40, 50, 60, 200)):
    """Side-by-side comparison of age at causal HPV infection under
    No interventions (no vax, no screening) vs Baseline (vax + status-quo
    screening). Shows that routine vaccination shifts the age
    distribution of causal HPV infection later - i.e., protects the
    younger cohorts."""
    scens = ['No interventions', 'Baseline']
    scens = [s for s in scens if s in df.scenario.unique()]
    if len(scens) < 2:
        return None
    out = {}
    for scen in scens:
        sub = df[df.scenario == scen].copy()
        sub['age_causal_bin'] = pd.cut(sub['age_causal'], bins=bins, right=False)
        per_rep = (sub.groupby(['rep', 'age_causal_bin'], observed=True)['weight']
                      .sum().unstack('age_causal_bin', fill_value=0))
        med = per_rep.median(axis=0)
        total = med.sum()
        out[f'{scen}: n'] = med.round(0).astype(int)
        out[f'{scen}: %'] = (med / total * 100).round(1)
    return pd.DataFrame(out)


def by_hiv(df, scenario):
    """§2c: share of residual that is HIV+."""
    sub = df[df.scenario == scenario].copy()
    per_rep = (sub.groupby(['rep', 'hiv'])['weight'].sum()
                  .unstack('hiv', fill_value=0))
    per_rep.columns = ['hiv_negative' if not c else 'hiv_positive'
                       for c in per_rep.columns]
    med = per_rep.median(axis=0).round(0)
    total = med.sum()
    return pd.DataFrame({'median': med,
                         'share_pct': (med / total * 100).round(1)})


def by_bucket(df, scenario):
    """§3: bucket decomposition of the residual."""
    sub = df[df.scenario == scenario].copy()
    per_rep = (sub.groupby(['rep', 'bucket'])['weight'].sum()
                  .unstack('bucket', fill_value=0))
    med = per_rep.median(axis=0).round(0)
    total = med.sum()
    return pd.DataFrame({'median': med,
                         'share_pct': (med / total * 100).round(1)})


def reachability(df, scenario, intv_start=INTV_START):
    """§4: partition residual into 'reached in principle' vs 'unreachable'.

    A cancer is 'reached in principle' if the woman was in *any*
    intervention's addressable age window at any point from intv_start
    to her cancer diagnosis:
      - routine vax:     age 11-12 at any t >= 2011 (already done by
                         intv_start=2030 for anyone born pre-2018)
      - routine screen:  age 30-50 at any t >= intv_start
      - mass adult:      age 20-50 at intv_start (one-off campaign)
      - TxV eligibility: same as routine screen (must be screen-positive)

    We approximate reachability from the birth year alone:
      - reachable_by_screen = birth_year <= (intv_start - SCREEN_AGE[0])
                              AND birth_year >= (2100 - SCREEN_AGE[1])
                              (born early enough to hit the window
                              before 2100 and late enough not to have
                              aged out entirely by intv_start)
    """
    sub = df[df.scenario == scenario].copy()

    by = sub['birth_year']
    # In screening window for at least one year between intv_start & 2100
    age_lo, age_hi = SCREEN_AGE
    reached_screen = ((by >= 2100 - age_hi) & (by <= intv_start - age_lo)) | \
                     ((by >= intv_start - age_hi) & (by <= intv_start - age_lo)) | \
                     ((by <= intv_start - age_lo) & (by >= 2100 - age_hi))
    # Simpler: was 30-50 at any point between intv_start and cancer_year
    # A woman aged 30-50 in [intv_start, cancer_year] means
    #   cancer_year - 50 <= birth_year <= intv_start - 30
    #   (i.e., she reached 30 by cancer_year and hadn't passed 50 by intv_start)
    cyear = sub['year_cancer']
    reached_screen = (cyear - age_hi <= by) & (by <= intv_start - age_lo)

    # Mass adult campaign: age 20-50 in intv_start (one-off in intv_start)
    a_lo, a_hi = MASS_ADULT_AGE
    reached_mass = (intv_start - a_hi <= by) & (by <= intv_start - a_lo)

    # Routine prophylactic vax: age 11-12 at any t >= 2011
    # equivalent to birth_year >= 1999 (turned 12 in 2011)
    reached_vax = by >= (ROUTINE_VAX_START - 12)

    reached_any = reached_screen | reached_mass | reached_vax
    sub['reached'] = reached_any

    per_rep = (sub.groupby(['rep', 'reached'])['weight'].sum()
                  .unstack('reached', fill_value=0))
    per_rep.columns = ['unreachable' if not c else 'reached_in_principle'
                       for c in per_rep.columns]
    med = per_rep.median(axis=0).round(0)
    total = med.sum()
    return pd.DataFrame({'median': med,
                         'share_pct': (med / total * 100).round(1)})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in', dest='inpath',
                    default='results/diagnostic_normalized/cancer_rows.csv')
    ap.add_argument('--outdir', default='results/diagnostic_normalized')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = load(args.inpath)

    # §1: totals per scenario
    totals = scenario_totals(df)
    totals.to_csv(f'{args.outdir}/residual_by_scenario.csv')
    print('\n=== §1: Cumulative cancers per scenario (2030-2100) ===')
    print(totals.to_string())

    # §2b bonus: age-at-infection under No interventions vs Baseline
    cmp_age = age_causal_novax_vs_baseline(df)
    if cmp_age is not None:
        cmp_age.to_csv(f'{args.outdir}/residual_age_causal_novax_vs_baseline.csv')
        print('\n=== §2b bonus: age at causal HPV infection - vax vs no vax ===')
        print(cmp_age.to_string())

    # Pick the best-performing intervention scenario as the "residual"
    intv_scens = [s for s in totals.index if s not in ('No interventions',)]
    best = totals.loc[intv_scens, 'median'].idxmin()
    print(f'\nBest intervention scenario: {best}')
    print(f'Residual: {int(totals.loc[best, "median"]):,} '
          f'[{int(totals.loc[best, "p10"]):,}, {int(totals.loc[best, "p90"]):,}]')

    # §2-4 for the best scenario, Baseline, and any TxV scenario present
    # (S&TxV's 8% bucket-3 is a specific policy talking point: TxV is
    # delivered directly at the screening visit, so LTFU protection is
    # what TxV effectively buys vs. HPV-Faster's triage-then-ablate chain.)
    focus = [best]
    for s in ['Baseline', 'S&TxV 70%', 'No interventions']:
        if s in totals.index and s not in focus:
            focus.append(s)

    for scen in focus:
        tag = scen.replace(' ', '_').replace('&', 'and').replace('/', '_').replace(',', '')

        by_birth = by_birth_cohort(df, scen)
        by_birth.to_csv(f'{args.outdir}/residual_birth_cohort__{tag}.csv')

        by_age = by_age_causal(df, scen)
        by_age.to_csv(f'{args.outdir}/residual_age_causal__{tag}.csv')

        by_h = by_hiv(df, scen)
        by_h.to_csv(f'{args.outdir}/residual_hiv__{tag}.csv')

        by_bu = by_bucket(df, scen)
        by_bu.to_csv(f'{args.outdir}/residual_bucket__{tag}.csv')

        reach = reachability(df, scen)
        reach.to_csv(f'{args.outdir}/residual_reachability__{tag}.csv')

        print(f'\n=== [{scen}] Birth cohort ===')
        print(by_birth.to_string())
        print(f'\n=== [{scen}] Age at causal infection ===')
        print(by_age.to_string())
        print(f'\n=== [{scen}] HIV status ===')
        print(by_h.to_string())
        print(f'\n=== [{scen}] Bucket ===')
        print(by_bu.to_string())
        print(f'\n=== [{scen}] Reachability ===')
        print(reach.to_string())


if __name__ == '__main__':
    main()
