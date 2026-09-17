# Handoff — hpvsim_rwanda revision cycle

Written 2026-09-17 for the next agent picking this up. Read this
whole file before touching anything.

---

## TL;DR

- The paper is at **step 4** of `docs/revision_plan.md` — reviewer
  response. Steps 1-3 (engineering uplift, v2→v3 migration,
  calibration) are done.
- Working branch is **`age-causal`** off `migrate-v3.2`. It is
  **ahead of origin, not pushed** — the user has been holding pushes
  while iterating on modelling decisions.
- The latest change is a **structural refactor of `make_st`** (see
  §"What just changed"). A **scenario sweep is currently running**
  in the background — do not start a competing one until it
  completes; check status first.
- The next task once the sweep + figures land is the **reviewer
  response plan** (`docs/reviewer_response_plan.md`, 15 numbered
  comments across R1 and R2, most requiring new sensitivity sweeps).

---

## Repo layout you need to know

Working directory: `/home/robyn/hpvsim_rwanda`.

Sibling deps (editable installs):
- `/home/robyn/hpvsim` — hpvsim core, on branch `restore-hpv-latency`.
- `/home/robyn/starsim` — starsim core.
- `/home/robyn/stisim` — stisim (for the HIV module).

Env: `/home/robyn/miniconda/envs/rwanda/bin/python` (aliased in this
project's default `python`). VS Code pre-activates
`/home/robyn/hpvsim/.conda` via PATH — **never `pip install` without
an absolute path**, or you will silently install into hpvsim's env
and break the sibling repo. See auto-memory for details
(`feedback_conda_env_prefix.md`).

Key files:
- `interventions.py` — intervention constructors. Recently
  refactored; §"What just changed" explains.
- `run_scenarios.py` — 23-scenario harness that runs everything via
  `ss.MultiSim`. CLI: `python run_scenarios.py --run-sim --end 2100`.
- `run_calibration.py` — Optuna calibration (~1500 trials, top-50
  extracted). Only run on a VM; do not run locally.
- `run_sim.py` — `make_sim`. This is where per-agent pars live
  (debut ages, layer_probs, HIV pars).
- `plot_fig{1,2,3,4,5}_*.py` + `plot_figS2_calib.py` — figure
  scripts, all reading `results/scens_*.csv` + `results/figS2_*.csv`.
- `build_table_s1.py` — Reviewer 1 comment 1 (prior + posterior
  table).
- `diagnose_all.py` — per-cancer decomposition harness with an
  `ScenTracker` analyzer for the residual analysis. Slow (~90s/sim).
- `analyze_residual.py` — aggregates the diagnose_all output into
  the tables in `docs/residual_analysis.md`.

Docs to read in order:
1. `docs/Reviewer Comments.md` — the 15 numbered comments.
2. `docs/revision_plan.md` — the four-step plan (we're at step 4).
3. `docs/reviewer_response_plan.md` — the concrete plan for step 4.
4. `docs/residual_analysis.md` — Phase 3 write-up (stale — the
   numbers in it are from the pre-refactor sweep; the structure and
   framing are still useful).
5. `docs/txv_investigation.md` + `docs/txv_investigation_review.md`
   — earlier debugging record for the S&TxV ordering question; still
   useful reference for the mechanics.

---

## What just changed

The user asked me to fix two things flagged during figure review:

### 1. Double-LTFU in `make_st` (fixed)

The S&T variant was chaining `routine_triage prob=0.9` × `treat_num
prob=0.75` = 68% clearance instead of the paper's single-visit 90%
attendance. Cause: `treat_num` applies its own acceptance filter on
top of the triage attendance step, which is correct for the
two-visit S&T&T model but a double count for the same-visit S&T.

Fix: when `tx_assigner_csv == 'tx_assigner_no_triage'`, `make_st`
force-overrides `treat_prob = 1.0` in the intv-era `treat_num`
interventions. The single-event 90% LTFU is supplied by the
`routine_triage prob=0.9` step, matching the paper. Chain per
lesion for S&T: `0.9 × 1.0 × 1.0 × 0.936 = 84%`. Same fix
propagated to `make_st_older` where `campaign_triage prob` was `1`
(no LTFU by assumption) — now `0.9` to match the paper.

### 2. SQ 2020-2027 + intv 2028+ era split (fixed)

Per user spec: "ALL scenarios [should] be the same over 2020-2027.
They should all run with the status-quo screening program, i.e.
18% screening coverage, S&T&T with VIA triage and 75% treatment
probability. From 2028 we do the screening coverage switch."

`make_st` now builds two eras with unique intervention names:

- **SQ era** (`start_year..coverage_change_year-1`, 2020-2027):
  hardcoded S&T&T at 18% coverage, VIA triage, 75% treatment.
  Intervention names have `_sq` suffix (`tx_assigner_sq_intv`,
  `ablation_sq_intv`, etc.); products carry unique `module_name`s
  (`ablation_sq_prod`, `excision_sq_prod`, `radiation_sq_prod`) so
  they don't collide with the intv-era's default-named products.
- **Intv era** (`coverage_change_year..end_year`, 2028+): the
  variant-specific path. For S&T&T that's identical to SQ era but
  at scaled coverage. For S&T it's no-triage + treat_prob=1.0. For
  S&TxV{,&T&T} it adds TxV from `txv_start_year=2030`.

Screening is one intervention with a time-varying `prob` covering
both eras (avoids the eligibility-across-two-screenings problem).

Removed: the deprecated `screen_change_year` param (replaced by
`coverage_change_year=2028`), `prev_screen_cov` (hardcoded 18% for
SQ era), `primary` param (was for VIA test injection; put back if
a sensitivity needs it), `make_normalized_scenarios` in
`run_scenarios.py` and its wiring in `diagnose_all.py`, and the
`q1_test_change_year.py` scratch script.

### 3. HPV-Faster campaign LTFU (fixed)

Per user: "For HPV-faster, change the campaign_triage prob to 0.9."
Done. `campaign_triage(name='tx_assigner_older', prob=0.9)`.
`treat_cov=1` default is kept so ablation adds no further LTFU.
Chain matches S&T's `0.9 × 1.0 × 1.0 × 0.936 = 84%`. Also shifted
`make_st_older(start_year=2028)` (was 2027) to align with the
S&T-family coverage-change year.

---

## Current state

### Git

```
Branch: age-causal
Ahead of origin by 6 commits (not pushed - the user has been
withholding push while iterating).

Recent commits:
  b53e00d  Refit Rwanda calibration with 5y cancer detection lag
  31f5dbe  Investigate age at causal HPV infection; switch debut to normal
  41daadc  Residual analysis rewrite: address the four review points
  8f26403  diagnose_all: drop routine vax for 'No interventions' scenario
  b2788e2  Residual analysis: who gets the 26K residual cancers under HPV-Faster 70%
  474f147  Phase 4 primitives: VIA product + TxV efficacy multiplier

Uncommitted (the SQ+intv+LTFU refactor described above):
  M interventions.py     - SQ+intv era split, LTFU fix, make_st_older
  M run_scenarios.py     - dropped make_normalized_scenarios + CLI flag
  M diagnose_all.py      - dropped build_normalized_scenarios + --normalized
  M plot_fig{1..5}_*.py  - title strings 2025->2030
  M utils.py             - docstring 2025->2030
  M results/table_s1.csv - regenerated on new calibration
  D q1_test_change_year.py

DO NOT push without the user's explicit go-ahead. They have
repeatedly declined pushes while modelling decisions are open.
```

### Sim runs

**A scenario sweep is currently running** (background,
`b8b3rxrwc`). It rebuilds `results/scens_{timeseries,cumulative,
paired,per_sim}.csv` + `results/st_scens.obj` with the new SQ+intv
era `make_st`. Expected duration ~12 minutes on this 160-core box
via `ss.MultiSim`. Check status with the notification event or:

```bash
tail -3 results/scen_run.log | grep -v clipped | grep -v betamap
ps -p $(pgrep -f 'run_scenarios.py --run-sim') -o pid,pcpu,etime 2>/dev/null
```

**Once it completes**, regenerate figures 1-5:

```bash
cd /home/robyn/hpvsim_rwanda
for f in fig1_residual fig2_st fig3_txv fig4_mass fig5_bars; do
  python plot_${f}.py --resfolder results 2>&1 | tail -6
done
```

The user will want to eyeball the new fig 5 bar ordering to see if
the S&T ↔ S&T&T ↔ S&TxV comparison now looks sensible under the
paper-consistent LTFU treatment (S&T should now avert ~45K instead
of the previous 37K).

### Calibration

Fresh calibration at `raw_results/rwanda_calib.obj` (Sep 17 11:38).
Best mismatch 1.08. 17 parameters (added `age_risk.risk`,
`imm_init.low`, `cell_imm_init.low` on top of the 14 already in
the base calibration). Median age at causal HPV infection now
~27.5 (target ~28) — this was the main driver of the recent
recalibration.

`results/table_s1.csv` is already regenerated against the new
calib.

---

## What's next

### Immediate (post-sweep)

1. Regenerate figures 1-5 (commands above).
2. Present the new numbers to the user for review. Key comparison:
   the S&T ↔ S&T&T bar heights in fig 5. Previously S&T averted 37K
   vs S&T&T 13.6K; under the LTFU fix S&T should be higher (~45K),
   and the S&TxV bars should also shift because the SQ era for
   2020-2027 is now identical across scenarios.
3. If the user approves the numbers, they'll say so. Only then
   consider pushing.

### After that: the reviewer response (step 4)

Read `docs/reviewer_response_plan.md` for the concrete grouping of
the 15 numbered reviewer comments into 5 sensitivity sweeps + 3
new tables + a cost-threshold calc + 7 text-only responses. Table
S1 (R1.1) is already delivered.

Sensitivity infrastructure is already plumbed into `make_st`:

- `txv_start_year=2030` — for R1.3 / R2.1 (vary 2030→2050).
- `treat_capacity=None` — for R2.5 (workforce cap).
- `future_treat_cov=0.75` — for R2.6 (LTFU sensitivity).
- `txv_efficacy_mult=1.0` — for R2.4 (±20% TxV efficacy).
- `via.csv` exists — for R2.6 (primary VIA vs HPV DNA). Pass the
  product via a new `primary=` kwarg (need to re-add — I removed it
  in the SQ+intv refactor; put back if the R2.6 sweep uses it).

The five sweeps as sized in the plan doc total ~250 sims (~15
minutes on this box). Do them via `run_scenarios.py` with new
scenario factories, not by hand-editing sims.

### Text-only responses

7 comments (R1.2, R1.4, R2.1 framing, R2.3, R2.9, R2.10, R2.11)
need a written response only, no runs. Draft in
`docs/reviewer_response_text.md`. Wait until the sensitivity sweeps
land so you can cite specific numbers.

---

## Standing user preferences (from prior sessions)

- **Domain fluency**: assume the user is an HPVsim core dev and
  first author on the paper. Skip textbook explanations.
- **Terse responses**: no trailing summaries, no headers on simple
  answers. State results and decisions directly.
- **Do not push without explicit approval.** The user has withheld
  pushes repeatedly. Commit locally, wait.
- **Do not run calibration locally.** Only on a VM. See
  `feedback_do_shrink_calibration.md` — `do_shrink=False` in bulk
  runs will OOM the box.
- **Simpler is usually right on calibration.** See
  `feedback_calibration_knob_creep.md` — widening bounds or adding
  knobs tends to worsen fit rather than help.
- **Auto-memory system** at
  `/home/robyn/.claude/projects/-home-robyn-hpvsim/memory/` — read
  MEMORY.md at session start. Update it when learning new
  cross-session facts. `feedback_*.md` files are the load-bearing
  ones.

---

## Known open concerns

- The user is worried the age-at-causal-infection was previously
  too old (fed the "unreachable older cohorts" framing in the
  residual analysis harder than the biology warranted). The
  recalibration on this branch should fix that; new fig 5 will
  confirm. If age-at-infection still looks off (median much >30)
  after the new sweep, this is unfinished.
- `results/diagnostic_normalized/cancer_rows.csv` (90MB, gitignored)
  and the residual analysis in `docs/residual_analysis.md` are
  built on the PRE-refactor scenario dynamics. After the current
  sweep finishes, the residual analysis needs a rerun via
  `diagnose_all.py` + `analyze_residual.py` if we want its numbers
  updated. The user hasn't asked for this yet — they may just want
  the paper's figures to look right first.
- Two secondary artifacts in the S&TxV chain I flagged but did not
  fix: (a) `linked_txvx prob=0.9` and `routine_triage prob=0.9`
  are independent coins (0.81 joint attendance) but the paper
  implies same-visit (0.9 joint). Second-order. (b) TxV starts in
  2030 while ablation runs from 2028 — 2-year window with only
  ablation before TxV kicks in. Consistent with the paper spec but
  causes some CIN backlog in the S&TxV family. Not blocking.

---

## Gotchas / footguns

- **Env**: as noted, VS Code pre-activates hpvsim's `.conda`.
  Absolute `python`/`pip` paths only.
- **Background sims**: `run_in_background=true` starts a new shell
  in `/home/robyn/hpvsim` (the primary CWD), not the rwanda repo.
  Always `cd /home/robyn/hpvsim_rwanda && ...` inside the bash
  command.
- **90MB CSVs**: `results/diagnostic*/cancer_rows.csv` is
  gitignored. The small aggregates (`summary.csv`, `residual_*.csv`)
  are tracked.
- **Beta clip warnings** flood the sim logs. Grep `-v clipped -v
  betamap` when tailing.
- **Analyzer copies**: `sim.analyzers['name']` returns a *copy* of
  the analyzer instance you passed in, not the original. Reference
  `sim.analyzers['name']` after `sim.run()` to inspect state, not
  the instance you constructed.
- **Fig 5 baseline**: fig 5's paired-diff baseline is `S&T&T 18%`
  (not `No interventions`). All "averted" numbers are vs
  `S&T&T 18%`.
