"""Analyze the diagnostic sweep output. Reads results/diagnostic/cancer_rows.csv."""
import numpy as np
import pandas as pd
import os

DIAG = 'results/diagnostic/cancer_rows.csv'
OUT = 'results/diagnostic'
os.makedirs(OUT, exist_ok=True)


def load():
    df = pd.read_csv(DIAG)
    return df


def summarize(df):
    """Print per-scenario bucket totals (median across reps of the scaled sum)."""
    # per rep totals per bucket
    per_rep = df.groupby(['scenario', 'rep', 'bucket'])['weight'].sum().unstack('bucket', fill_value=0)
    per_rep_total = per_rep.sum(axis=1)
    per_rep['total'] = per_rep_total

    print('\n=== Per-scenario median scaled bucket totals (2025-2100) ===')
    med = per_rep.groupby('scenario').median()
    print(med.round(0).to_string())
    return per_rep, med


def per_rep_totals(df):
    """Return per (scenario, rep) totals."""
    return df.groupby(['scenario', 'rep'])['weight'].sum().unstack('rep')


def birth_cohort_table(df):
    """Median scaled cancers by scenario × 10yr birth cohort."""
    per_rep = df.groupby(['scenario', 'rep', 'birth_cohort_10'])['weight'].sum().unstack('birth_cohort_10', fill_value=0)
    med = per_rep.groupby('scenario').median()
    return med


def age_in_2027_table(df):
    """Median scaled cancers by scenario × age band in 2027."""
    df = df.copy()
    bins = [-100, 0, 10, 20, 30, 50, 70, 200]
    labels = ['unborn_at_end', 'unborn', '<10', '10-20', '20-30', '30-50', '50-70', '70+']
    df['age_band_2027'] = pd.cut(df['age_in_2027'], bins=bins, labels=labels)
    per_rep = df.groupby(['scenario', 'rep', 'age_band_2027'], observed=True)['weight'].sum().unstack('age_band_2027', fill_value=0)
    med = per_rep.groupby('scenario').median()
    return med


def mass_vax_reach(df):
    """For HPV-Faster: how many cancers occurred in women who got mass vax vs not?"""
    hf = df[df['scenario'].str.startswith('hpvfaster')].copy()
    per_rep = hf.groupby(['scenario', 'rep', 'got_mass_vx'])['weight'].sum().unstack('got_mass_vx', fill_value=0)
    return per_rep.groupby('scenario').median()


def screen_status_for_hpvfaster(df):
    """Compare 'never_screened' bucket size across scenarios — key for Q4."""
    ns = df[df['bucket'] == '1_never_screened']
    per_rep = ns.groupby(['scenario', 'rep'])['weight'].sum().unstack('rep')
    med = per_rep.median(axis=1)
    print('\nNever-screened cancer count (median across reps):')
    print(med.round(0).sort_values().to_string())
    return med


if __name__ == '__main__':
    df = load()
    print(f'\nLoaded {len(df):,} cancer rows across {df["scenario"].nunique()} scenarios × {df["rep"].nunique()} reps')

    per_rep, med = summarize(df)
    med.to_csv(f'{OUT}/bucket_medians.csv')

    print('\n=== Per-(scenario,rep) tracked totals ===')
    print(per_rep_totals(df).round(0).to_string())

    print('\n=== Birth-cohort medians ===')
    bc = birth_cohort_table(df)
    print(bc.round(0).to_string())
    bc.to_csv(f'{OUT}/birth_cohort_medians.csv')

    print('\n=== Age-in-2027 medians ===')
    ab = age_in_2027_table(df)
    print(ab.round(0).to_string())
    ab.to_csv(f'{OUT}/age_in_2027_medians.csv')

    print('\n=== HPV-Faster mass vax reach ===')
    mv = mass_vax_reach(df)
    print(mv.round(0).to_string())

    screen_status_for_hpvfaster(df)

    # Cancer decomposition ratios
    print('\n=== Fraction of cancers in each bucket per scenario ===')
    frac = med.copy()
    total_col = frac.sum(axis=1)
    for c in frac.columns:
        frac[c] = (100 * frac[c] / total_col).round(1)
    print(frac.to_string())
    frac.to_csv(f'{OUT}/bucket_fraction.csv')

    print('\nDone.')
