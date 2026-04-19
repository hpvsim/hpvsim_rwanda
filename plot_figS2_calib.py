"""
Fig S2: calibration diagnostic.

Plots from plot-ready CSVs produced by `run_calibration.py` (extract step).
"""
import argparse
import os

import numpy as np
import pandas as pd
import pylab as pl
import sciris as sc

import utils as ut


def plot_calib(resfolder='results', outpath='figures/fig_calib.png',
               hiv_datafile='data/rwanda_data.csv'):
    ut.set_font(16)
    fig = pl.figure(layout="tight", figsize=(16, 10))

    gs = fig.add_gridspec(3, 1)
    gs0 = gs[0].subgridspec(1, 4)
    gs1 = gs[1].subgridspec(1, 4)
    gs2 = gs[2].subgridspec(1, 4)

    prev_col = '#5f5cd2'
    canc_col = sc.gridcolors(5)[3]
    ms = 80
    gen_cols = sc.gridcolors(4)

    # ---- Row 0: cancers by age ----
    age_metrics = [
        ('cancers', ['0', '15', '20', '25', '30', '35', '40', '45', '50', '55',
                     '60', '65', '70', '75', '80', '85'],
         'Cancers by age, 2020', gs0[:2]),
        ('cancer_incidence_no_hiv', ['25-35', '35-45', '45-55', '55+'],
         'Cancer incidence by age,\n2017, HIV- women', gs0[2]),
        ('cancer_incidence_with_hiv', ['25-35', '35-45', '45-55', '55+'],
         'Cancer incidence by age,\n2017, HIV+ women', gs0[3]),
    ]
    for rkey, age_labels, title, spec in age_metrics:
        ax = fig.add_subplot(spec)
        model_df = pd.read_csv(f'{resfolder}/figS2_{rkey}.csv').sort_values('bin')
        target_df = pd.read_csv(f'{resfolder}/figS2_target_{rkey}.csv')

        stats = [dict(med=r.med, q1=r.q1, q3=r.q3, whislo=r.whislo, whishi=r.whishi,
                      fliers=[], label=str(int(r.bin)))
                 for _, r in model_df.iterrows()]
        bp = ax.bxp(stats, positions=model_df['bin'].values,
                    patch_artist=True, showfliers=False, manage_ticks=False)
        for patch in bp['boxes']:
            patch.set_facecolor(canc_col); patch.set_alpha(0.7)
        x = np.arange(len(age_labels))
        ax.scatter(x, target_df.value.values, marker='d', s=ms, color='k')
        ax.set_xticks(x); ax.set_xticklabels(age_labels)
        ax.set_ylim(bottom=0); ax.set_ylabel(''); ax.set_xlabel('')
        ax.set_title(title)

    # ---- Row 1, left: time series (asr + by HIV status) ----
    ts_df = pd.read_csv(f'{resfolder}/figS2_timeseries.csv')
    ax = fig.add_subplot(gs1[:2])
    rkeys = ['asr_cancer_incidence', 'cancer_incidence_with_hiv', 'cancer_incidence_no_hiv']
    rlabels = ['Total', 'HIV+', 'HIV-']
    for rkey, rlabel in zip(rkeys, rlabels):
        sub = ts_df[ts_df.metric == rkey].sort_values('year')
        ax.plot(sub.year, sub.med, label=rlabel)
        ax.fill_between(sub.year, sub.pi95_low, sub.pi95_high, alpha=0.2)

    target = pd.read_csv(f'{resfolder}/figS2_target_asr_cancer_incidence.csv')
    ax.scatter(target.year.values[0], target.value.values[0], marker='d', s=ms,
               color='k', label='GLOBOCAN')
    ax.legend(loc='upper left', frameon=False, fontsize=12)
    ax.set_title('Cancer incidence, 2020')
    ax.set_ylabel(''); ax.set_xlabel('')

    # ---- Row 1, right: genotype distributions ----
    geno_metrics = [
        ('precin_genotype_dist', 'Share of LSILs\nby genotype, 2020', gs1[2], gen_cols[0]),
        ('cancerous_genotype_dist', 'Share of cancers\nby genotype, 2020', gs1[3], gen_cols[1]),
    ]
    for rkey, title, spec, color in geno_metrics:
        ax = fig.add_subplot(spec)
        model_df = pd.read_csv(f'{resfolder}/figS2_{rkey}.csv').sort_values('bin')
        target_df = pd.read_csv(f'{resfolder}/figS2_target_{rkey}.csv')

        stats = [dict(med=r.med, q1=r.q1, q3=r.q3, whislo=r.whislo, whishi=r.whishi,
                      fliers=[], label=str(int(r.bin)))
                 for _, r in model_df.iterrows()]
        bp = ax.bxp(stats, positions=model_df['bin'].values,
                    patch_artist=True, showfliers=False, manage_ticks=False)
        for patch in bp['boxes']:
            patch.set_facecolor(color); patch.set_alpha(0.5)

        ax.scatter(np.arange(len(target_df)), target_df.value.values, color='k', marker='d', s=ms)
        ax.set_xticks(np.arange(4)); ax.set_xticklabels(['16', '18', 'Hi5', 'OHR'])
        ax.set_ylabel(''); ax.set_xlabel('')
        ax.set_title(title)

    # ---- Row 2: ART, HIV prevalence, HIV infections, HIV deaths ----
    hiv_df = pd.read_csv(hiv_datafile).set_index('Unnamed: 0').T
    hiv_df.index = hiv_df.index.astype(int)
    hiv_window = hiv_df.loc[(hiv_df.index >= 2000) & (hiv_df.index <= 2025)]

    def plot_ts_panel(ax, rkey, color, title, scale=1.0, legend=False, label=None):
        sub = ts_df[ts_df.metric == rkey].sort_values('year')
        ax.plot(sub.year, sub.med * scale, color=color, label=label)
        ax.fill_between(sub.year, sub.pi95_low * scale, sub.pi95_high * scale,
                        alpha=0.2, color=color)
        if rkey in hiv_window.columns:
            ax.scatter(hiv_window.index, hiv_window[rkey] * scale, marker='d',
                       s=ms, color='k')
        ax.set_xlabel(''); ax.set_ylabel(''); ax.set_title(title)
        if not legend and ax.get_legend() is not None:
            ax.get_legend().remove()

    ax = fig.add_subplot(gs2[0])
    plot_ts_panel(ax, 'art_coverage', prev_col, 'ART coverage (%)', scale=100)
    ax.set_ylim(bottom=0)

    ax = fig.add_subplot(gs2[1])
    sex_colors = sc.gridcolors(2)
    plot_ts_panel(ax, 'female_hiv_prevalence', sex_colors[0], 'HIV prevalence (%)', scale=100, label='Female')
    plot_ts_panel(ax, 'male_hiv_prevalence', sex_colors[1], 'HIV prevalence (%)', scale=100, label='Male')
    ax.legend(loc='upper right', frameon=False)

    ax = fig.add_subplot(gs2[2])
    plot_ts_panel(ax, 'hiv_infections', prev_col, 'HIV infections')
    sc.SIticks()

    ax = fig.add_subplot(gs2[3])
    plot_ts_panel(ax, 'hiv_deaths', prev_col, 'HIV deaths')
    sc.SIticks()

    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    pl.savefig(outpath, dpi=300)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--resfolder', default='results/v2.2.6_baseline')
    parser.add_argument('--outpath', default='figures/fig_calib.png')
    args = parser.parse_args()

    T = sc.timer()
    plot_calib(resfolder=args.resfolder, outpath=args.outpath)
    T.toc('Done')
