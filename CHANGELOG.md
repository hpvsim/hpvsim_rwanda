# Changelog

## 2026-09 — v3.2 revision

- Added sensitivity analyses for TxV introduction year (`run_sensitivity_txv.py` + `plot_figS4_txv_intro.py`) and workforce cap (`run_sensitivity_workforce.py` + `plot_figS5_workforce.py`).
- Added threshold cost-effectiveness reporting driven by `build_table_s3_costing.py` (Table S3) alongside the existing Table S1/S2 builders.
- Wired the DALY analyzer into `run_sim.py` / `run_scenarios.py` so DALY outputs propagate through to the paired-scenario CSVs.
- Fixed the excision counter double-counting in `interventions.py`.
- Renumbered references and added Spencer, Pan, and Canfell citations throughout MS and SM.
- Repo tidy for the revision: dropped one-off diagnostic scripts (analyze_*, diagnose_*, compare_baselines) and their `results/diagnostic*` CSVs; retired the `v2.2.6_baseline` / `v2.3.0_baseline` snapshots and pointed plot script `--resfolder` defaults at `results/`; retained the 1500-trial calibration backup in `raw_results/` as a fallback for the 10k re-run.

## 2026-04-19 — HPVsim v2.2.6 lift

- Split the workflow so heavy simulations run on a VM via `--run-sim` and produce plot-ready CSVs; local plot scripts load the CSVs (no pickles).
- Added `run_scenarios.py --run-sim` to emit `scens_timeseries.csv` (year × scenario × metric with 95% CIs) and `scens_cumulative.csv` (2025–2100 sums).
- Added `run_calibration.py` CSV extraction (`figS2_*.csv`: boxplot stats, med/pi95 time series, targets).
- Refactored all paper plots (`plot_fig1_residual`, `plot_fig2_st`, `plot_fig3_txv`, `plot_fig4_mass`, `plot_fig5_bars`, `plot_figS2_calib`) to read from CSVs and accept `--resfolder` / `--outpath`.
- Refactored matching poster variants (`plot_fig2_st_poster`, `plot_fig3_txv_poster`, `plot_fig4_mass_poster`).
- Froze plot-ready baseline under `results/v2.2.6_baseline/`; plot scripts default to this folder.
- Dropped orphan `plot_fig3_hiv.py` (references scenarios no longer generated) and `plot_fig2_vx.py` (not in the manuscript).
- Removed `age_causal_infection` analyzer hook from `run_sim.py` (was consumed only by the deleted `plot_fig2_vx.py`).
- Added `.gitignore` excluding transient `.obj/.pkl/.sim/.msim/.csv` outputs in `results/` (plot-ready CSVs under `results/v<version>_baseline/` are tracked).

## Earlier

- Mass-campaign + therapeutic-vaccine scenario suite (branch `cleanup`, merged into `main`): added `plot_fig2_st.py`, `plot_fig3_txv.py`, `plot_fig4_mass.py`, and poster variants; reorganised `interventions.py`.
- Initial calibration and ongoing S&T&T screening + therapeutic vaccine analyses for the Rwanda manuscript.
