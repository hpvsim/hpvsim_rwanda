"""
Fig S3 (revision): TxV introduction-year sensitivity.

Cancers averted 2030-2100 vs S&T&T 18% baseline, plotted as a function
of TxV introduction year for each of the four TxV-carrying arms.

Reads paired-diff CSV produced by `run_sensitivity.py --run-sim`
(default resfolder: results/sens_txv_year).
"""
import argparse
import os

import numpy as np
import pylab as pl
import sciris as sc

import utils as ut
from run_sensitivity import TXV_INTRO_YEARS, TXV_COV


ARMS = [
    ('S&TxV&T&T', 'S&TxV+T&T (virus-clearing 90/0)'),
    ('S&TxV',     'S&TxV (lesion-regressing 50/90)'),
    ('Mass TxV 90/0',  'Mass TxV (virus-clearing 90/0)'),
    ('Mass TxV 50/90', 'Mass TxV (lesion-regressing 50/90)'),
]


def plot_figS3(resfolder='results/sens_txv_year',
               outpath='figures/figS3_txv_intro.png'):
    ut.set_font(16)
    _, _, paired_df, _ = ut.load_scens(resfolder)

    fig, ax = pl.subplots(figsize=(10, 7), layout='tight')
    colors = sc.gridcolors(len(ARMS))
    cov_lbl = f'{TXV_COV*100:.0f}%'

    for (arm_key, arm_label), color in zip(ARMS, colors):
        med, lo, hi = [], [], []
        for yr in TXV_INTRO_YEARS:
            scen = f'{arm_key} {cov_lbl} (TxV {yr})'
            m, l, h = ut.get_paired(paired_df, scen, 'cancers')
            med.append(m); lo.append(l); hi.append(h)
            print(f'{scen}: averted {m:.0f} ({l:.0f}, {h:.0f})')
        med = np.asarray(med); lo = np.asarray(lo); hi = np.asarray(hi)
        ax.plot(TXV_INTRO_YEARS, med, marker='o', color=color, label=arm_label, lw=2)
        ax.fill_between(TXV_INTRO_YEARS, lo, hi, color=color, alpha=0.15)

    ax.axhline(0, color='0.4', lw=0.5)
    ax.set_xticks(TXV_INTRO_YEARS)
    ax.set_xlabel('TxV introduction year')
    ax.set_ylabel('Cancers averted 2030-2100\n(vs S&T&T 18%)')
    ax.set_title(f'TxV introduction-year sensitivity ({cov_lbl} coverage)')
    ax.legend(frameon=False, loc='upper right', fontsize=13)
    sc.SIticks(ax=ax)

    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    sc.savefig(outpath, dpi=150)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--resfolder', default='results/sens_txv_year')
    parser.add_argument('--outpath', default='figures/figS3_txv_intro.png')
    args = parser.parse_args()
    plot_figS3(resfolder=args.resfolder, outpath=args.outpath)
    print('Done.')
