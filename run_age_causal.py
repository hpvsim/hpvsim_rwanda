"""
Age at causal HPV infection & dwell-time distributions (Rwanda).

Runs two scenarios to 2050 with the ``age_causal_infection`` analyzer
(start=2020) attached:
    1. No interventions
    2. Baseline (routine vax + status-quo screen-and-treat, cov=0.18)

Produces two figures under figures/:
    - age_causal_bar_rwanda.png    (age-at-causal-infection, per scenario)
    - dwelltimes_rwanda.png        (precin/cin/total dwell time, per scenario)

Analyzer arrays are pickled to raw_results/age_causal_rwanda.obj so plots
can be regenerated without rerunning the sims.
"""

import os

os.environ.update(
    OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1', MKL_NUM_THREADS='1',
)

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import sciris as sc
import hpvsim as hpv
from scipy.stats import gaussian_kde

import run_sim as rs
from interventions import make_st


SCENARIOS = ['No interventions', 'Baseline']
END = 2050
START_ANALYZER = 2020


# ---------- run ----------

def _build_sim(scen_label, end=END):
    analyzers = [hpv.age_causal_infection(start=START_ANALYZER)]
    if scen_label == 'No interventions':
        return rs.make_sim(add_vax=False, add_st=False, interventions=[],
                           analyzers=analyzers, stop=end, use_calib=True)
    if scen_label == 'Baseline':
        return rs.make_sim(add_vax=True, add_st=False,
                           interventions=make_st(future_screen_cov=0.18, end_year=end),
                           analyzers=analyzers, stop=end, use_calib=True)
    raise ValueError(scen_label)


def _extract(sim):
    a = sim.analyzers['age_causal_infection']
    return dict(
        age_causal=np.asarray(a.age_causal),
        age_cin=np.asarray(a.age_cin),
        age_cancer=np.asarray(a.age_cancer),
        weights=np.asarray(a.weights),
        dwelltime={k: np.asarray(v) for k, v in a.dwelltime.items()},
    )


def run(end=END, verbose=1/12, out='raw_results/age_causal_rwanda.obj'):
    data = {}
    for scen in SCENARIOS:
        print(f'\n=== Running: {scen} (to {end}) ===')
        sim = _build_sim(scen, end=end)
        sim.pars.verbose = verbose
        sim.run()
        data[scen] = _extract(sim)
        print(f'  {scen}: n_cancers={len(data[scen]["age_cancer"])}  '
              f'weighted={data[scen]["weights"].sum():.1f}')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    sc.saveobj(out, data)
    print(f'\nSaved analyzer arrays to {out}')
    return data


# ---------- plot ----------

SCEN_COLORS = {'No interventions': '#3498db', 'Baseline': '#27ae60'}


def _weighted_percentile(a, w, q):
    """Weighted percentile; q in [0, 100]."""
    order = np.argsort(a)
    a, w = a[order], w[order]
    c = np.cumsum(w) - 0.5 * w
    c /= w.sum()
    return np.interp(q / 100.0, c, a)


def plot_age_causal_bar(data, outpath='figures/age_causal_bar_rwanda.png'):
    """Kenya-style causal-age bar+KDE, one panel per scenario."""
    fig, axes = pl.subplots(1, len(SCENARIOS), figsize=(9, 3.4), sharey=False)
    if len(SCENARIOS) == 1:
        axes = [axes]

    bins = np.arange(10, 81)
    bin_centers = bins[:-1]

    for ax, scen in zip(axes, SCENARIOS):
        ages = data[scen]['age_causal']
        w = data[scen]['weights']
        if len(ages) == 0:
            ax.set_title(f'{scen}\n(no causal infections)')
            continue
        # keep only reasonable range (age_causal analyzer can spill outside
        # the plot window; mirror the Kenya <50 filter loosely by clipping
        # to the bin range)
        m = (ages >= bins[0]) & (ages <= bins[-1])
        ages, w = ages[m], w[m]

        counts, _ = np.histogram(ages, bins=bins, weights=w)
        pcts = counts / counts.sum() * 100

        ax.bar(bin_centers, pcts, color=SCEN_COLORS[scen], edgecolor='none',
               width=1.0, alpha=0.5)

        # KDE for smooth overlay (unweighted; weights are ~constant here)
        try:
            kde = gaussian_kde(ages, bw_method=0.15, weights=w)
        except TypeError:  # older scipy
            kde = gaussian_kde(ages, bw_method=0.15)
        x_smooth = np.linspace(bins[0], bins[-1], 500)
        kde_vals = kde(x_smooth)
        kde_scaled = kde_vals / kde_vals.sum() * pcts.sum() * (500 / (bins[-1] - bins[0]))
        dark = '#1a5276' if scen == 'No interventions' else '#0e6251'
        ax.plot(x_smooth, kde_scaled, color=dark, linewidth=2)

        q25 = _weighted_percentile(ages, w, 25)
        med = _weighted_percentile(ages, w, 50)
        q75 = _weighted_percentile(ages, w, 75)
        print(f'{scen}: 25th={q25:.1f}, Median={med:.1f}, 75th={q75:.1f}')

        ymax = max(pcts) if len(pcts) else 1.0
        for val, label, yoff in [(q25, '25th', 1.18), (med, 'Median', 1.08), (q75, '75th', 1.18)]:
            color = '#e74c3c' if label == 'Median' else '#555555'
            lw = 2 if label == 'Median' else 1.5
            ax.axvline(val, color=color, linestyle='--', linewidth=lw, zorder=5)
            ax.text(val, ymax * yoff, label, ha='center', fontsize=9,
                    color=color, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                              edgecolor='none', alpha=1.0), zorder=10)

        ax.set_xticks(np.arange(10, 81, 10))
        ax.set_xlim(9.5, 80.5)
        ax.set_ylim(0, ymax * 1.30)
        ax.set_xlabel('Age', fontsize=12)
        ax.set_ylabel('% of causal infections', fontsize=12)
        ax.yaxis.set_major_formatter(pl.FuncFormatter(lambda x, _: f'{x:.0f}%'))
        ax.tick_params(axis='y', labelsize=10)
        ax.set_title(scen, fontsize=13)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.suptitle('Age at causal HPV infection (Rwanda, cancers 2020-2050)',
                 fontsize=13, y=1.02)
    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    sc.savefig(outpath, dpi=200)
    print(f'Figure saved to: {outpath}')
    return fig


DWELL_ORDER = ['precin', 'cin', 'total']
DWELL_LABELS = {
    'precin': 'Causal infection\n→ CIN2+',
    'cin':    'CIN2+ → cancer',
    'total':  'Causal infection\n→ cancer',
}


def plot_dwelltimes(data, outpath='figures/dwelltimes_rwanda.png'):
    """Weighted boxplot of dwell times (precin, cin, total) per scenario."""
    fig, axes = pl.subplots(1, len(SCENARIOS), figsize=(9, 3.6), sharey=True)
    if len(SCENARIOS) == 1:
        axes = [axes]

    def _wq(x, w, q):
        return _weighted_percentile(x, w, q)

    for ax, scen in zip(axes, SCENARIOS):
        d = data[scen]
        w = d['weights']
        stats = []
        for k in DWELL_ORDER:
            x = d['dwelltime'][k]
            if len(x) == 0:
                stats.append(None); continue
            stats.append(dict(
                med=_wq(x, w, 50),
                q1=_wq(x, w, 25),
                q3=_wq(x, w, 75),
                whislo=_wq(x, w, 5),
                whishi=_wq(x, w, 95),
                fliers=[],
                label=DWELL_LABELS[k],
            ))
        valid = [s for s in stats if s is not None]
        if not valid:
            ax.set_title(f'{scen}\n(no cancers)'); continue

        bp = ax.bxp(valid, showfliers=False, patch_artist=True)
        for patch in bp['boxes']:
            patch.set_facecolor(SCEN_COLORS[scen])
            patch.set_alpha(0.6)
            patch.set_edgecolor('#333')
        for med in bp['medians']:
            med.set_color('#e74c3c')
            med.set_linewidth(2)

        # Annotate median years
        ymax_data = max(s['whishi'] for s in valid)
        for i, s in enumerate(valid, start=1):
            ax.text(i, s['med'], f' {s["med"]:.1f}y', va='center', ha='left',
                    fontsize=9, color='#333')

        ax.set_title(scen, fontsize=13)
        ax.set_ylabel('Dwell time (years)', fontsize=12)
        ax.set_ylim(0, max(ymax_data * 1.15, 5))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', alpha=0.3, linestyle='--')

    fig.suptitle('Dwell times to cervical cancer (Rwanda, cancers 2020-2050)',
                 fontsize=13, y=1.02)
    fig.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    sc.savefig(outpath, dpi=200)
    print(f'Figure saved to: {outpath}')

    # Also print a compact summary table for the log
    rows = []
    for scen in SCENARIOS:
        d = data[scen]; w = d['weights']
        for k in DWELL_ORDER:
            x = d['dwelltime'][k]
            if len(x) == 0: continue
            rows.append(dict(scenario=scen, transition=k,
                             median=_wq(x, w, 50),
                             q25=_wq(x, w, 25),
                             q75=_wq(x, w, 75),
                             n=int(len(x))))
    print('\nDwell-time summary (years):')
    print(pd.DataFrame(rows).to_string(index=False, float_format='%.2f'))
    return fig


# ---------- entry ----------

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-sim', action='store_true',
                        help='Run the two sims (heavy). Otherwise reload cached arrays.')
    parser.add_argument('--end', type=int, default=END)
    parser.add_argument('--cache', default='raw_results/age_causal_rwanda.obj')
    args = parser.parse_args()

    T = sc.timer()
    if args.run_sim or not os.path.exists(args.cache):
        data = run(end=args.end, out=args.cache)
    else:
        print(f'Loading cached arrays from {args.cache}')
        data = sc.loadobj(args.cache)

    plot_age_causal_bar(data)
    plot_dwelltimes(data)
    T.toc('Done')
