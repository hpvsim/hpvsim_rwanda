"""
Fig 4: one-time mass campaign strategies.

Plots from plot-ready CSVs produced by `run_scenarios.py --run-sim`.
"""
import argparse
import os

import numpy as np
import pylab as pl
import sciris as sc
from matplotlib.patches import Patch

import utils as ut


def plot_fig4(resfolder='results', outpath='figures/fig4_campaigns.png'):
    ut.set_font(20)
    ts_df, cum_df, paired_df, per_sim_df = ut.load_scens(resfolder)

    fig = pl.figure(layout="tight", figsize=(18, 10))
    gs = fig.add_gridspec(2, 3)

    text_height = [-0.1, 1.2]
    start_year, end_year = 2016, 2100
    colors = sc.vectocolor(3).tolist()

    coverage_levels = ['18%', '35%', '70%']
    strategies = {
        'HPV-Faster': 'HPV-Faster',
        'Mass TxV 90/0': 'Mass TxV 90/0,',
        'Mass TxV 50/90': 'Mass TxV 50/90,',
    }
    bar_width = 0.25
    x = np.arange(len(coverage_levels))
    offsets = [-bar_width, 0, bar_width]

    def cum(scen, metric):
        return ut.get_cum(cum_df, scen, metric)[0]

    def stack(scen_keys, metric):
        med, lo, hi = zip(*(ut.get_cum(cum_df, k, metric) for k in scen_keys))
        return list(med), list(lo), list(hi)

    def total_stats(scen_keys, metrics):
        """Per-sim sum across `metrics`, then median/10/90 across sims."""
        med, lo, hi = [], [], []
        for k in scen_keys:
            sub = per_sim_df[(per_sim_df.scenario == k) & (per_sim_df.metric.isin(metrics))]
            tot = sub.groupby('sim')['value'].sum().values
            med.append(float(np.median(tot))); lo.append(float(np.quantile(tot, 0.1))); hi.append(float(np.quantile(tot, 0.9)))
        return med, lo, hi

    # ---- A: Cumulative cancers ----
    ax = fig.add_subplot(gs[0, 0])
    for idx, (strat_label, strat_key) in enumerate(strategies.items()):
        keys = [f'{strat_key} {cov}' for cov in coverage_levels]
        med, lo, hi = stack(keys, 'cancers')
        for k, v, l, h in zip(keys, med, lo, hi):
            print(f'{k}: {v:.0f} ({l:.0f}, {h:.0f}) cancers')
        ax.bar(x + offsets[idx], med, width=bar_width, color=colors[idx], label=strat_label,
               yerr=ut.yerr(med, lo, hi), capsize=3, ecolor='0.4')
    ax.set_xticks(x); ax.set_xticklabels(coverage_levels)
    ax.set_title('Cumulative cancers\n2030-2100'); sc.SIticks()
    ax.set_xlabel('')
    ax.legend(loc='upper right', frameon=False, fontsize=16)
    ax.set_ylim(bottom=0)
    ax.text(*text_height, 'A', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    # ---- B: Cumulative cancers in HIV+ ----
    ax = fig.add_subplot(gs[0, 1])
    for idx, (strat_label, strat_key) in enumerate(strategies.items()):
        keys = [f'{strat_key} {cov}' for cov in coverage_levels]
        med, lo, hi = stack(keys, 'cancers_with_hiv')
        ax.bar(x + offsets[idx], med, width=bar_width, color=colors[idx], label=strat_label,
               yerr=ut.yerr(med, lo, hi), capsize=3, ecolor='0.4')
    ax.set_xticks(x); ax.set_xticklabels(coverage_levels)
    ax.set_title('Cumulative cancers in HIV+ women\n2030-2100'); sc.SIticks()
    ax.set_xlabel('')
    ax.set_ylim(bottom=0)
    ax.text(*text_height, 'B', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    # ---- C: Combined resource use (stacked; error bars on total) ----
    ax = fig.add_subplot(gs[0, 2])
    for idx, (strat_label, strat_key) in enumerate(strategies.items()):
        keys = [f'{strat_key} {cov}' for cov in coverage_levels]
        m_abl, _, _ = stack(keys, 'ablations')
        m_tx,  _, _ = stack(keys, 'txvs')
        m_vx,  _, _ = stack(keys, 'vaccinations')
        for k, a, t, v in zip(keys, m_abl, m_tx, m_vx):
            print(f'{k}: {a:.0f} ablations, {t:.0f} therapeutics, {v:.0f} vaccinations')
        med, lo, hi = total_stats(keys, ('ablations', 'txvs', 'vaccinations'))
        ax.bar(x + offsets[idx], m_abl, width=bar_width, color=colors[idx])
        ax.bar(x + offsets[idx], m_tx, width=bar_width, bottom=m_abl,
               color='none', edgecolor=colors[idx], hatch='//', linewidth=1.5)
        bottom_vals = [a + t for a, t in zip(m_abl, m_tx)]
        ax.bar(x + offsets[idx], m_vx, width=bar_width, bottom=bottom_vals,
               color='none', edgecolor=colors[idx], hatch='\\\\', linewidth=1.5,
               yerr=ut.yerr(med, lo, hi), capsize=3, ecolor='0.4')
    ax.set_xticks(x); ax.set_xticklabels(coverage_levels)
    ax.set_title('Cumulative products\n2030-2100'); sc.SIticks()
    ax.set_xlabel('')
    ax.legend(handles=[
        Patch(facecolor='gray', label='Ablations'),
        Patch(facecolor='white', edgecolor='k', hatch='//', label='Therapeutics'),
        Patch(facecolor='white', edgecolor='k', hatch='\\\\', label='Vaccinations'),
    ], loc='upper left', frameon=False, fontsize=16)
    ax.set_ylim(bottom=0)
    ax.text(*text_height, 'C', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    # ---- D: Time series at 70% coverage ----
    ax = fig.add_subplot(gs[1, :2])
    ax = ut.plot_ts(ax, ts_df, 'S&T&T 18%', 'asr_cancer_incidence', start_year, end_year,
                    color='k', label='Status quo (S&T&T 18%)')
    for strat_idx, (strat_label, strat_key) in enumerate(strategies.items()):
        scen_key = f'{strat_key} 70%'
        ax = ut.plot_ts(ax, ts_df, scen_key, 'asr_cancer_incidence', start_year, end_year,
                        color=colors[strat_idx], label=f'{strat_label} 70%')
    ax.set_ylim(bottom=0)
    ax.set_title('ASR cervical cancer incidence, 2030-2100\nOne-time mass campaigns at 70% coverage')
    ax.legend(loc='upper right', frameon=False)
    ax.text(-0.05, 1.2, 'D', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    # ---- E: Interventions per cancer averted (per-sim ratio, median +/- 10/90) ----
    ax = fig.add_subplot(gs[1, 2])
    def ratio_stats(scen_key, num_metrics=('txvs', 'vaccinations')):
        pv = per_sim_df[per_sim_df.scenario == scen_key].set_index(['sim', 'metric'])['value']
        base = per_sim_df[per_sim_df.scenario == 'S&T&T 18%'].set_index(['sim', 'metric'])['value']
        sims = pv.index.get_level_values('sim').unique()
        ratios = []
        for s in sims:
            averted = base.loc[(s, 'cancers')] - pv.loc[(s, 'cancers')]
            if averted > 0:
                total = sum(pv.loc[(s, m)] for m in num_metrics if (s, m) in pv.index)
                ratios.append(total / averted)
        if not ratios:
            return 0.0, 0.0, 0.0
        arr = np.asarray(ratios)
        return float(np.median(arr)), float(np.quantile(arr, 0.1)), float(np.quantile(arr, 0.9))

    for idx, (strat_label, strat_key) in enumerate(strategies.items()):
        keys = [f'{strat_key} {cov}' for cov in coverage_levels]
        med, lo, hi = zip(*(ratio_stats(k) for k in keys))
        for k, m in zip(keys, med):
            print(f'{k}: {m:.1f} products per cancer averted')
        ax.bar(x + offsets[idx], med, width=bar_width, color=colors[idx], label=strat_label,
               yerr=ut.yerr(med, lo, hi), capsize=3, ecolor='0.4')
    ax.set_xticks(x); ax.set_xticklabels(coverage_levels)
    ax.set_title('Interventions per cancer averted\n2030-2100')
    ax.set_xlabel('')
    ax.set_ylim(bottom=0)
    ax.text(*text_height, 'E', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    sc.savefig(outpath, dpi=100)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--resfolder', default='results/v2.2.6_baseline')
    parser.add_argument('--outpath', default='figures/fig4_campaigns.png')
    args = parser.parse_args()
    plot_fig4(resfolder=args.resfolder, outpath=args.outpath)
    print('Done.')
