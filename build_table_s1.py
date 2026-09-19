"""Build Table S1: prior distributions + posterior summary + description/source
for the 17 calibrated parameters (R1.1).

Reads raw_results/rwanda_calib.obj. Priors come from calib.calib_pars;
posterior comes from the N-best trials of calib.df. Description and prior
source come from PARAM_METADATA below.

Writes results/table_s1.csv.
"""
import argparse
import os

import numpy as np
import pandas as pd
import sciris as sc


# Per-parameter metadata for R1.1: what each calibrated parameter represents
# and where the prior range came from. Posterior for every row is from this
# study (calibrated against Rwanda cancer + HIV cancer + HPV genotype targets).
PARAM_METADATA = {
    'beta': (
        'HPV transmission probability per act',
        'HPVsim defaults; range widened to accommodate Rwanda-specific fit',
    ),
    'hpv16.cin_fn.k': (
        'HPV16 dose-response shape (dysplasia progression rate)',
        'Prior from HPVsim natural-history calibration (Stuart et al. 2024)',
    ),
    'hpv18.cin_fn.k': (
        'HPV18 dose-response shape (dysplasia progression rate)',
        'Prior from HPVsim natural-history calibration (Stuart et al. 2024)',
    ),
    'hi5.cin_fn.k': (
        'Pooled Hi5 (31/33/45/52/58) dose-response shape',
        'Prior from HPVsim natural-history calibration (Stuart et al. 2024)',
    ),
    'ohr.cin_fn.k': (
        'Pooled OHR (35/39/51/56/59) dose-response shape',
        'Prior from HPVsim natural-history calibration (Stuart et al. 2024)',
    ),
    'age_risk.risk': (
        'Age-dependent multiplier on HPV acquisition risk',
        'Introduced in this study to fit the observed age distribution of causal infection',
    ),
    'imm_init.low': (
        'Lower bound of initial humoral immunity after natural clearance',
        'Introduced in this study; prior from HPVsim natural-history defaults',
    ),
    'cell_imm_init.low': (
        'Lower bound of initial cell-mediated immunity after lesion regression',
        'Introduced in this study; prior from HPVsim natural-history defaults',
    ),
    'network.m_cross_layer': (
        'Male probability of a concurrent partnership across marital/casual layers',
        'Prior from HPVsim network calibration (Stuart et al. 2024)',
    ),
    'network.f_cross_layer': (
        'Female probability of a concurrent partnership across marital/casual layers',
        'Prior from HPVsim network calibration (Stuart et al. 2024)',
    ),
    'network.m_partners_casual': (
        'Male mean casual partners per year (Poisson rate)',
        'Prior from HPVsim network calibration (Stuart et al. 2024)',
    ),
    'network.f_partners_casual': (
        'Female mean casual partners per year (Poisson rate)',
        'Prior from HPVsim network calibration (Stuart et al. 2024)',
    ),
    'hiv.rel_sus_lo': (
        'HIV+ relative HPV susceptibility at low CD4',
        'Prior informed by Liu et al. 2018 meta-analysis (ref 25)',
    ),
    'hiv.rel_sus_hi': (
        'HIV+ relative HPV susceptibility at high CD4',
        'Prior informed by Liu et al. 2018 meta-analysis (ref 25)',
    ),
    'hiv.rel_sev_lo': (
        'HIV+ relative multiplier on HPV disease progression at low CD4',
        'Prior informed by Liu et al. 2018 meta-analysis (ref 25)',
    ),
    'hiv.rel_sev_hi': (
        'HIV+ relative multiplier on HPV disease progression at high CD4',
        'Prior informed by Liu et al. 2018 meta-analysis (ref 25)',
    ),
    'hiv.p_effective_art': (
        'Probability that ART fully suppresses HIV-driven HPV effects',
        'Prior from UNAIDS ART effectiveness estimates',
    ),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in', dest='inpath', default='raw_results/rwanda_calib.obj')
    ap.add_argument('--n_top', type=int, default=50)
    ap.add_argument('--out', default='results/table_s1.csv')
    args = ap.parse_args()

    calib = sc.loadobj(args.inpath)
    priors = calib.calib_pars
    df = calib.df

    print(f'Loaded {len(df)} trials; using top {args.n_top} by mismatch.')
    top = df.nsmallest(args.n_top, 'mismatch')

    rows = []
    for pname, prior in priors.items():
        vals = top[pname].values.astype(float)
        desc, source = PARAM_METADATA.get(pname, ('', ''))
        rows.append({
            'parameter':       pname,
            'description':     desc,
            'prior_source':    source,
            'prior_low':       prior['low'],
            'prior_high':      prior['high'],
            'prior_guess':     prior['guess'],
            'posterior_mean':  float(np.mean(vals)),
            'posterior_median': float(np.median(vals)),
            'posterior_ci_lo': float(np.quantile(vals, 0.025)),
            'posterior_ci_hi': float(np.quantile(vals, 0.975)),
            'posterior_min':   float(np.min(vals)),
            'posterior_max':   float(np.max(vals)),
        })

    out = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out) or '.', exist_ok=True)
    out.to_csv(args.out, index=False, float_format='%.4f')
    print(f'Wrote {args.out}')
    print()
    print(out.to_string(index=False))


if __name__ == '__main__':
    main()
