# HPVsim model for Rwanda

Code for modelling HPV transmission, cervical cancer burden, and prevention strategies (screening, therapeutic vaccination, and HPV-Faster) in Rwanda, and for reproducing the figures of the accompanying manuscript.

**Results in this repository were produced with HPVsim v2.2.6.** The plot-ready baseline lives in [`results/v2.2.6_baseline/`](results/v2.2.6_baseline/); all plot scripts default to reading from this folder.

## Installation

```bash
pip install hpvsim==2.2.6 seaborn optuna
```

Python 3.9+.

## Workflow: heavy sims on VM, plots locally

Every heavy script has two modes: `--run-sim` runs the simulation and saves lightweight CSVs under `results/`; running without the flag loads the CSVs and produces the figure. The intended flow is:

1. **VM:** run `python run_scenarios.py --run-sim`, commit and push the resulting `results/*.csv` files.
2. **Local:** pull, run `python plot_fig1_residual.py` etc. to render figures from the CSVs.

Transient top-level `results/*.obj` / `.csv` / `.sim` / `.msim` / `.zip` are gitignored. Plot-ready CSVs frozen under `results/v<version>_baseline/` **are** committed — these are the durable artifacts.

### Heavy (VM) scripts

| Script | What it produces |
|---|---|
| `run_scenarios.py --run-sim` | `scens_timeseries.csv`, `scens_cumulative.csv` (23 scenarios × annual metrics + cumulative 2025–2100 sums) |
| `run_calibration.py --run-sim` | `rwanda_calib.obj` (heavy) + `figS2_*.csv` plot-ready summaries |
| `run_sim.py` | Single sim for debugging (no persistent outputs) |

### Plot scripts (local)

| Script | Manuscript figure |
|---|---|
| `plot_fig1_residual.py` | Fig 1 — residual cancer burden under ongoing interventions |
| `plot_fig2_st.py` | Fig 2 — screening ± VIA triage at 18/35/70% coverage |
| `plot_fig3_txv.py` | Fig 3 — therapeutic-vaccine-enhanced screening |
| `plot_fig4_mass.py` | Fig 4 — one-time mass campaigns (HPV-Faster + Mass TxV) |
| `plot_fig5_bars.py` | Fig 5 — all seven strategies side-by-side |
| `plot_figS2_calib.py` | Fig S2 — calibration diagnostic |

Poster variants (`plot_fig2_st_poster.py`, `plot_fig3_txv_poster.py`, `plot_fig4_mass_poster.py`) emit to `figures/poster/`.

### Other scripts

- `run_sim.py` — sim builder shared by all runners
- `interventions.py` — screening, TxV, and vaccination intervention factories
- `compare_baselines.py` — cross-version comparison (see below)
- `utils.py` — CSV loaders and `plot_ts` helper

## Reproducing the manuscript figures

```bash
# Clone and install
git clone git@github.com:hpvsim/hpvsim_rwanda.git
cd hpvsim_rwanda
pip install hpvsim==2.2.6 seaborn

# Render all figures from the committed v2.2.6 baseline
python plot_fig1_residual.py
python plot_fig2_st.py
python plot_fig3_txv.py
python plot_fig4_mass.py
python plot_fig5_bars.py
python plot_figS2_calib.py
```

## Cross-version comparison (v2.2.6 → v2.3 → v3.0)

When a new HPVsim version ships, regenerate the baseline in a clean env and compare side-by-side:

```bash
# 1. On a VM, in a clean env pinned to the new version
conda create -n hpvsim230 python=3.11 -y && conda activate hpvsim230
pip install hpvsim==2.3.0 seaborn optuna

# 2. Re-run the heavy scripts
python run_calibration.py --run-sim          # calibration (slowest)
python run_scenarios.py --run-sim            # all 23 scenarios

# 3. Freeze the fresh CSVs into a versioned baseline dir
mkdir -p results/v2.3.0_baseline
cp results/scens_*.csv results/figS2_*.csv results/v2.3.0_baseline/

# 4. Commit + push, then locally compare
python compare_baselines.py --baselines v2.2.6_baseline v2.3.0_baseline
```

`compare_baselines.py` extends trivially to v3.0 by appending the new baseline name to `--baselines`. The comparison prints a cumulative-cancers + elimination-year table and emits `figures/compare_baselines.png` (ASR trajectories overlaid per scenario).

## Inputs

- `data/` — HIV incidence, mortality, ART coverage, and cancer incidence targets (IARC, GLOBOCAN)
- `hpvdna.csv`, `tx_assigner*.csv`, `txvx_pars_*.csv` — screening + TxV intervention specifications

## Citation

If you use this code, please cite the accompanying Rwanda manuscript (in review).

## Further information

See [hpvsim.org](https://hpvsim.org) and [docs.hpvsim.org](https://docs.hpvsim.org).
