"""Build Table S1: prior distributions + posterior summary (top 50 trials).

Reviewer 1 comment 1: 'a full table of prior distributions and posterior
mean/95% CIs across the 50 best-fitting parameter sets'.

Reads raw_results/rwanda_calib.obj. Priors come from calib.calib_pars;
posterior comes from the N-best trials of calib.df.

Writes results/table_s1.csv (columns: parameter, prior_low, prior_high,
posterior_mean, posterior_ci_lo, posterior_ci_hi, posterior_median).
"""
import argparse
import os

import numpy as np
import pandas as pd
import sciris as sc


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
        rows.append({
            'parameter':       pname,
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
