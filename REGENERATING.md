# Regenerating the v3 figures (hpvsim_rwanda)

Which script produces which figure for the v2→v3 review. Run from this repo dir with the v3
venv (`.venv`, hpvsim 3.0.0); confirm `hpvsim.__version__ == '3.0.0'` first. This is the HIV repo;
sims keep `ms_agent_ratio=3` and `n_agents≈15000` so the sparse HIV+ cancer stratum resolves.
Run foreground (no fragile long background jobs).

| Figure | Driver → output | Then plot |
|---|---|---|
| fig1_residual, fig2_st, fig3_txv, fig4_campaigns, fig5_comparison | `python run_scenarios.py --run-sim --full --resfolder results/_v3gap` — runs all 23 scenarios (3 seeds) and writes `scens_{timeseries,cumulative}.csv` | `python plot_fig1_residual.py --resfolder results/_v3gap --outpath ...` (and `plot_fig2_st.py` … `plot_fig5_bars.py`) |
| figS2 (calibration diagnostic: cancer incidence by age × HIV, HPV-type shares, ART/HIV prevalence) | `python _gap_figS2.py` — runs the calibrated natural-history sim (5 seeds) and extracts the 6 model CSVs into `results/_v3gapS2/`, then copies the version-independent `figS2_target_*.csv` there | `python plot_figS2_calib.py --resfolder results/_v3gapS2 --outpath ...` |

Notes:
- The full 23-scenario set runs at the default `STOP=2051` now (the campaign builders propagate
  `end_year` — earlier it needed `STOP≥2101`).
- Two fixes were required for the full set to run: the engine TxV-delivery fix (hpvsim branch
  `fix/txvx-administer-dispatch`) and the vx-module naming fix (here). Both are committed.
- Absolute cancer counts are far below the v2.2.6 national-scaled baseline (`total_pop=n_agents`);
  the aggregate ASR runs ~3× low vs GLOBOCAN even though the by-age HIV− registry rates match — a
  known by-age-vs-aggregate calibration tension, not a rendering bug. See the review repo
  `hpvsim_v23_migration_review` for details.
