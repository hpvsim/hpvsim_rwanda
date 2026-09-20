"""
Fig S3 (revision): natural history under baseline Rwanda programs.

Two-panel figure showing (A) age at causal HPV infection and (B) dwell
times from causal infection to cancer, both under the baseline scenario
(routine prophylactic vaccination + status-quo 18% screen-and-treat).

Reads analyzer arrays saved by run_age_causal.py at
results/age_causal_rwanda.obj.
"""
import argparse
import os

import numpy as np
import matplotlib.pyplot as pl
import sciris as sc
from scipy.stats import gaussian_kde

import utils as ut


SCEN = 'Baseline'
COLOR = '#27ae60'
DARK = '#0e6251'


def _weighted_percentile(a, w, q):
    order = np.argsort(a)
    a, w = a[order], w[order]
    c = np.cumsum(w) - 0.5 * w
    c /= w.sum()
    return np.interp(q / 100.0, c, a)


DWELL_ORDER = ['precin', 'cin', 'total']
DWELL_LABELS = {
    'precin': 'HPV to CIN2+',
    'cin':    'CIN2+ to cancer',
    'total':  'HPV to cancer',
}


def plot_figS2(inpath='results/age_causal_rwanda.obj',
               outpath='figures/figS3_natural_history.png'):
    ut.set_font(16)
    data = sc.loadobj(inpath)[SCEN]

    fig, (ax_age, ax_dwell) = pl.subplots(1, 2, figsize=(13, 4.5))

    # ---------- Panel A: age at causal infection ----------
    bins = np.arange(10, 81)
    bin_centers = bins[:-1]
    ages = data['age_causal']
    w = data['weights']
    m = (ages >= bins[0]) & (ages <= bins[-1])
    ages, w = ages[m], w[m]

    counts, _ = np.histogram(ages, bins=bins, weights=w)
    pcts = counts / counts.sum() * 100
    ax_age.bar(bin_centers, pcts, color=COLOR, edgecolor='none',
               width=1.0, alpha=0.5)

    try:
        kde = gaussian_kde(ages, bw_method=0.15, weights=w)
    except TypeError:
        kde = gaussian_kde(ages, bw_method=0.15)
    x_smooth = np.linspace(bins[0], bins[-1], 500)
    kde_vals = kde(x_smooth)
    kde_scaled = kde_vals / kde_vals.sum() * pcts.sum() * (500 / (bins[-1] - bins[0]))
    ax_age.plot(x_smooth, kde_scaled, color=DARK, linewidth=2)

    q25 = _weighted_percentile(ages, w, 25)
    med = _weighted_percentile(ages, w, 50)
    q75 = _weighted_percentile(ages, w, 75)
    print(f'{SCEN} age at causal infection: 25th={q25:.1f}, median={med:.1f}, 75th={q75:.1f}')

    ymax = max(pcts) if len(pcts) else 1.0
    for val, label, yoff in [(q25, '25th', 1.18), (med, 'Median', 1.08), (q75, '75th', 1.18)]:
        color = '#e74c3c' if label == 'Median' else '#555555'
        lw = 2 if label == 'Median' else 1.5
        ax_age.axvline(val, color=color, linestyle='--', linewidth=lw, zorder=5)
        ax_age.text(val, ymax * yoff, label, ha='center', fontsize=9,
                    color=color, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                              edgecolor='none', alpha=1.0), zorder=10)

    ax_age.set_xticks(np.arange(10, 81, 10))
    ax_age.set_xlim(9.5, 80.5)
    ax_age.set_ylim(0, ymax * 1.30)
    ax_age.set_xlabel('Age')
    ax_age.set_ylabel('% of causal infections')
    ax_age.yaxis.set_major_formatter(pl.FuncFormatter(lambda x, _: f'{x:.0f}%'))
    ax_age.set_title('A. Age at causal HPV infection', loc='left')
    ax_age.spines['top'].set_visible(False)
    ax_age.spines['right'].set_visible(False)

    # ---------- Panel B: dwell times ----------
    dw = data['weights']
    stats = []
    for k in DWELL_ORDER:
        x = data['dwelltime'][k]
        if len(x) == 0:
            continue
        wk = dw if len(dw) == len(x) else np.ones_like(x)
        stats.append(dict(
            med=_weighted_percentile(x, wk, 50),
            q1=_weighted_percentile(x, wk, 25),
            q3=_weighted_percentile(x, wk, 75),
            whislo=_weighted_percentile(x, wk, 5),
            whishi=_weighted_percentile(x, wk, 95),
            fliers=[],
            label=DWELL_LABELS[k],
        ))

    bp = ax_dwell.bxp(stats, showfliers=False, patch_artist=True)
    for patch in bp['boxes']:
        patch.set_facecolor(COLOR)
        patch.set_alpha(0.6)
        patch.set_edgecolor('#333')
    for m in bp['medians']:
        m.set_color('#e74c3c')
        m.set_linewidth(2)

    ymax_data = max(s['whishi'] for s in stats)
    for i, s in enumerate(stats, start=1):
        ax_dwell.text(i, s['med'], f' {s["med"]:.1f}y', va='center', ha='left',
                      fontsize=9, color='#333')
        print(f'{SCEN} {DWELL_ORDER[i-1]}: median {s["med"]:.1f} years')

    ax_dwell.set_ylabel('Dwell time (years)')
    ax_dwell.set_ylim(0, max(ymax_data * 1.15, 5))
    ax_dwell.set_title('B. Dwell times to cervical cancer', loc='left')
    ax_dwell.spines['top'].set_visible(False)
    ax_dwell.spines['right'].set_visible(False)
    ax_dwell.grid(axis='y', alpha=0.3, linestyle='--')

    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    sc.savefig(outpath, dpi=200)
    print(f'Saved {outpath}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--inpath', default='results/age_causal_rwanda.obj')
    parser.add_argument('--outpath', default='figures/figS3_natural_history.png')
    args = parser.parse_args()
    plot_figS2(inpath=args.inpath, outpath=args.outpath)
    print('Done.')
