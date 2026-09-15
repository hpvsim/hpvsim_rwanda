# Reviewer response — concrete plan

Structured against `docs/Reviewer Comments.md` (R1 1-4, R2 1-11).
Grouping the 15 comments by what they actually require (see the
revision plan step 4 table).

**Dependencies:** all new runs use the normalized scenario set
(everything in 2030, accounting 2030-2100) built on the v3.2 pipeline.
Sensitivity knobs plumbed into `make_st`:

- `txv_start_year=2030` — for R1.3 / R2.1
- `treat_capacity=None` — for R2.5
- `future_treat_cov=0.75` — already parametric, for R2.6
- `primary=None` (uses `make_hpv_test` by default) — swap for VIA product for R2.6
- `txv_pars` CSV path — swap for a scaled-efficacy version for R2.4

---

## Table of comments → actions

| # | Comment (short) | Type | Action |
|---|---|---|---|
| R1.1 | Table S1 of prior + posterior | New table | **DONE** — `build_table_s1.py`; `results/table_s1.csv` (14 params × top 50 trials) |
| R1.2 | Sensitivity around calibration uncertainty | Sensitivity | Already captured in the paired-uncertainty bars (top 50 pars + seeds). Text-only response pointing to the intervals |
| R1.3 | Vary TxV introduction year 2030 → 2050 | Sensitivity | New run: S&TxV 70% and Mass TxV 50/90 70% at intro years {2030, 2035, 2040, 2045, 2050} |
| R1.4 | Discussion framing | Text | Text-only |
| R2.1 | 2030 TxV is implausible; who is this for? | Text + sensitivity | Same run as R1.3; text-only otherwise |
| R2.2 | Cost / affordability | New calc | See R2.8 |
| R2.3 | Framing / motivation | Text | Text-only |
| R2.4 | TxV assumptions need references + sensitivity | Text + sensitivity | New run: scale TxV efficacy CSV by ×{0.8, 1.0, 1.2}; cite the two vaccine profile sources |
| R2.5 | Workforce implications | Sensitivity + external | New run: `treat_capacity` sweep. External: Rwanda ablation volume + provider count — needs local data source, may only appear as a caveat |
| R2.6 | VIA sensitivity + 25% LTFU | Sensitivity | New run: replicate S&T&T 70% with (a) VIA product (b) `future_treat_cov=0.75` reduced to 0.5 (25% LTFU beyond current) |
| R2.7 | Result tables of numbers requiring intervention | New tables | Extract from existing scen runs: cumulative screens, treatments, TxV doses per scenario |
| R2.8 | Threshold cost per DALY averted | New calc | DALYs from cancer deaths × life-years lost (unweighted); solve for cost that meets Rwanda's CET (~$225 per DALY at 1× GDP per capita) |
| R2.9 | Framing | Text | Text-only |
| R2.10 | Framing | Text | Text-only |
| R2.11 | Framing | Text | Text-only |

---

## New runs required (5 sweeps)

Each runs the same 10 top-calibration reps for tight uncertainty.

### Sweep A — TxV intro year (R1.3 + R2.1)

Scenarios: `S&TxV 70%` and `Mass TxV 50/90, 70%` at 5 intro years.
`(2 scenarios × 5 years) × 10 reps = 100 sims`.

Command sketch:
```bash
python run_scenarios.py --run-sim --scenario-set sensitivity_txv_year --end 2100
```
Requires a new scenario factory `make_txv_year_scenarios(intv_start_year=2030, txv_years=[2030,2035,2040,2045,2050])` in run_scenarios.py.

Analysis: for each intro year, cumulative cancers 2030-2100 and elimination year. Reviewer question is what delay we can tolerate before losing the effect.

### Sweep B — Workforce cap (R2.5)

Scenarios: `S&T&T 70%` at max_capacity ∈ {∞, 5, 3, 1} agents/timestep.
At Rwanda pop_scale × 4 timesteps/year, these map to roughly {∞, 26k, 15.6k, 5.2k} treatments/year.
`4 caps × 10 reps = 40 sims`.

Analysis: sensitivity of cancers averted to treatment capacity. If real-world Rwanda capacity is ~5k/yr, that maps to `max_capacity=1` — the elimination pathway needs stakeholder discussion of workforce expansion.

### Sweep C — VIA + LTFU (R2.6)

Requires a new `via.csv` product with realistic VIA per-stage sensitivity
(precin ~15%, cin ~55%, cancer ~80%, specificity ~95%).

Scenarios: `S&T&T 70%` (baseline) vs `S&T&T 70%, VIA` vs `S&T&T 70%, 25% LTFU`.
`3 scenarios × 10 reps = 30 sims`.

Analysis: how sensitive is the 70% headline to primary-test choice and to treatment LTFU.

### Sweep D — TxV efficacy (R2.4)

Scenarios: `S&TxV 70%` and `Mass TxV 50/90, 70%` at efficacy scale ∈ {0.8, 1.0, 1.2}.
`2 scenarios × 3 scales × 10 reps = 60 sims`.

Requires `txvx_pars_cin_scaled_0.8.csv` etc, or a `txv_efficacy_mult` knob in `make_st`.

Analysis: reviewer wants to see the TxV numbers aren't a knife-edge on optimistic efficacy assumptions.

### Sweep E — Elimination push (§5 of residual analysis)

Scenarios that stack the strongest levers: `S&TxV 90%` + workforce=∞ + `treat_cov=0.9`.
Comparator: current `S&TxV 70%` baseline.
`~2 × 10 reps = 20 sims`.

Analysis: what fraction of the residual cancers close when we push everything to plausible upper bounds?

Total new sims across A-E: ~250. At 90s/sim on 32 cores that's ~12 minutes VM time.

---

## Text-only responses (7 comments)

For R1.2, R1.4, R2.1 (framing part), R2.3, R2.9, R2.10, R2.11:
draft in a single markdown file `docs/reviewer_response_text.md`
once the new runs are done and we know what numbers to cite.

---

## R2.5 external data (workforce)

Ablation providers and current volume in Rwanda — not a model run.
Sources to try:
- Rwanda MoH cervical cancer strategic plan
- WHO country profile
- Recent Rwanda-based publications (Fadhil, Umuhoza)

If not obtainable, treat as a caveat in the discussion: "we do not
model workforce constraints, but the treatment volumes implied by the
S&T 70% scenario (X per year) exceed current national capacity of Y,
so scale-up is contingent on training and infrastructure investment
outside the scope of this analysis."

---

## R2.8 cost threshold — mechanics

For each scenario, per rep:
1. Cumulative new_cancers, new_cancer_deaths 2030-2100.
2. DALYs = (cancer_deaths × YLL) + (cancer_incidence × disability_weight × avg_dur_disease).
   Use Rwanda life expectancy at each age of death, WHO disability weights.
3. Threshold cost per person treated = (baseline_DALYs - scenario_DALYs) × CET / n_persons_treated_incremental
   where CET is Rwanda's cost-effectiveness threshold (~$225 per DALY at 1× GDP p.c.; more conservative $75 at 0.3× GDP).

Not a full CEA — just the "must be cheaper than $X per dose" number the
reviewer asks for. Deliverable: `results/cost_thresholds.csv` (scenario ×
threshold at 0.3×/1×/3× GDP CET).

---

## Suggested execution order

1. **DONE**: Table S1 (R1.1) — 5 min
2. **Now**: Sweep A (TxV year, R1.3+R2.1) — 100 sims
3. Sweep B (workforce, R2.5) — 40 sims
4. Sweep C (VIA + LTFU, R2.6) — 30 sims + `via.csv` construction
5. Sweep D (TxV efficacy, R2.4) — 60 sims + efficacy multiplier
6. Sweep E (elimination push) — 20 sims
7. R2.7 tables — script extract from existing scenario CSVs
8. R2.8 cost thresholds — script + DALY assumptions doc
9. Draft `docs/reviewer_response_text.md` — 15 numbered responses
10. Update the manuscript
