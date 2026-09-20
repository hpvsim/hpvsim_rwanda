"""
Fig S5 (revision): workforce capacity sensitivity (R2.5).

Cancers averted 2030-2100 vs S&T&T 18% baseline as a function of workforce
cap on treat_num.max_capacity, for each S&T-family strategy at 70% coverage.

Reads paired-diff CSVs from:
 - results/sens_workforce/scens_paired.csv  (workforce-capped runs)
 - results/scens_paired.csv                  (no-cap main sweep, for reference)
"""
import argparse
import os

import numpy as np
import pylab as pl
import sciris as sc

import utils as ut
from run_sensitivity_workforce import CAPS


ARMS = [
    ('S&T&T 70%', 'S&T&T (with VIA triage)'),
    ('S&T 70%',   'S&T (no triage)'),
    ('S&TxV&T&T 70%', 'S&TxV+T&T (virus-clearing)'),
    ('S&TxV 70%',     'S&TxV (lesion-regressing)'),
]


def plot_figS5(cap_resfolder='results/sens_workforce',
               nocap_resfolder='results',
               outpath='figures/figS5_workforce.png'):
    ut.set_font(16)
    _, _, cap_paired, _ = ut.load_scens(cap_resfolder)
    _, _, nocap_paired, _ = ut.load_scens(nocap_resfolder)

    fig, ax = pl.subplots(figsize=(10, 7), layout='tight')
    colors = sc.gridcolors(len(ARMS))

    x_caps = np.array(CAPS)

    for (arm_key, arm_label), color in zip(ARMS, colors):
        med, lo, hi = [], [], []
        for cap in CAPS:
            scen = f'{arm_key} (cap{cap})'
            m, l, h = ut.get_paired(cap_paired, scen, 'cancers')
            med.append(m); lo.append(l); hi.append(h)
            print(f'{scen}: averted {m:.0f} ({l:.0f}, {h:.0f})')
        med = np.asarray(med); lo = np.asarray(lo); hi = np.asarray(hi)
        ax.plot(x_caps, med, marker='o', color=color, label=arm_label, lw=2)
        ax.fill_between(x_caps, lo, hi, color=color, alpha=0.15)

        # Reference no-cap value (horizontal dashed line at right edge)
        nocap_m, _, _ = ut.get_paired(nocap_paired, arm_key, 'cancers')
        ax.axhline(nocap_m, color=color, ls='--', lw=1, alpha=0.5)
        ax.annotate(f'no cap: {nocap_m/1e3:.0f}K',
                    xy=(x_caps[-1], nocap_m), xytext=(6, 0),
                    textcoords='offset points', color=color, fontsize=11,
                    va='center', ha='left')

    ax.axhline(0, color='0.4', lw=0.5)
    ax.set_xticks(CAPS)
    ax.set_xticklabels([f'{c}\n(~{c*4000/1000:.0f}K/yr)' for c in CAPS])
    ax.set_xlabel('Workforce cap on treatments (agents per timestep)\n'
                  '(approximate annual population-scale equivalent shown below)')
    ax.set_ylabel('Cancers averted 2030-2100\n(vs S&T&T 18%)')
    ax.set_title('Workforce capacity sensitivity, 70% coverage scenarios')
    ax.legend(frameon=False, loc='upper left', fontsize=12)
    sc.SIticks(ax=ax)

    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    sc.savefig(outpath, dpi=150)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cap-resfolder', default='results/sens_workforce')
    parser.add_argument('--nocap-resfolder', default='results')
    parser.add_argument('--outpath', default='figures/figS5_workforce.png')
    args = parser.parse_args()
    plot_figS5(cap_resfolder=args.cap_resfolder,
               nocap_resfolder=args.nocap_resfolder,
               outpath=args.outpath)
    print('Done.')
