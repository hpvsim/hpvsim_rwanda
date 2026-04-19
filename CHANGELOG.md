# Changelog

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
