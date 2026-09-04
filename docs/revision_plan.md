# Revision plan: hpvsim_rwanda

Four steps, each on its own branch. Steps 1–3 are engineering; step 4 is the science response.

Current state: results produced with **HPVsim v2.2.6**, frozen in `results/v2.2.6_baseline/`. Target: **v3.2.0**.

---

## Step 1 — Engineering uplift

**Branch:** `eng-uplift` → PR → you merge.

Pure hygiene on the v2.2.6 code. No behaviour changes, so the frozen baseline must still reproduce byte-identically at the end.

1. Run `/idm-eng-plugin:eng-quality-checker` to get a scored report, then apply fixes via the fixer skill.
2. Known gaps visible without running anything:
   - No `requirements.txt` / `pyproject.toml` — the hpvsim version is pinned only in prose in the README. Add a real pin (`hpvsim==2.2.6`).
   - No tests, no CI.
   - `utils.py` and `interventions.py` carry the shared logic but have no module docstrings.
   - Six near-duplicate `plot_fig*_poster.py` variants duplicate their non-poster siblings.
3. Verification gate: re-run every plot script against `results/v2.2.6_baseline/` and confirm the figures are unchanged.

**Deliverable:** PR that changes only code quality, with the "figures unchanged" check stated in the PR body.

---

## Step 2 — Migration to hpvsim v3.2

**Branch:** `migrate-v3.2`. You push, pull onto the VMs, recalibrate.

**Success criterion (yours):** all figures tell more or less the same story. Parameter values may change; overall results must not.

This is not a version bump — v2.2.6 → v3.2.0 crosses a ground-up Starsim rebuild plus four sets of documented regressions. Expect to recalibrate from scratch, not to port parameters.

### 2a. Blocker to resolve first: HIV results no longer exist

The Rwanda model calibrates to, and Figure S2 plots, HIV results that v3.2 removed. Checked against `hpvsim/hiv.py`:

| Used by Rwanda | Status in v3.2 |
|---|---|
| `female_hiv_prevalence`, `male_hiv_prevalence` | **Gone.** v3.2 exposes `prevalence` only — no sex disaggregation |
| `hiv_infections`, `hiv_deaths` | **Gone.** Now `new_infections` / `new_deaths` on the HIV module |
| `art_coverage` (result) | **Gone.** Now `p_on_art`; `art_coverage` survives only as an input |
| `n_females_with_hiv_alive`, `n_males_with_hiv_alive` | **Gone** with the ~104 deleted sex-by-age strata |
| `cancers_by_age_with_hiv` / `_no_hiv` | Rebuild via `hpv.by_age` |
| `cancers_with_hiv`, `cancer_incidence_with_hiv` | Survive, but **redefined** — see 2b |

v3.2 deliberately ships only the 24 HIV results it recomputes with per-agent weighting. The sex-stratified ones were dropped because they shared the `count_nonzero` defect.

**Decision needed:** where does sex-stratified HIV prevalence go?

- **Recommended — hpvsim (layer 2).** Add scale-weighted `prevalence_f` / `prevalence_m` to `hpv.HIV` alongside the existing recomputed results. Any HIV–HPV model wants these, the weighting logic already exists in `_rescale_stisim_results`, and it is a small addition to code you are already editing on `restore-hpv-latency`.
- Alternative — a Rwanda-local analyzer (layer 3). Faster, but every other localization re-solves it.

Sex-stratified prevalence is generically useful, so it belongs upstream. Flagging rather than assuming, since it means a second PR against hpvsim.

### 2b. Regressions that will legitimately move the numbers

Do not treat these as migration bugs. Rwanda runs `dt=0.25`, `ms_agent_ratio=100`, `n_agents=10e3`, so all of them bite:

| Change | Version | Effect on Rwanda |
|---|---|---|
| `ablation`/`excision` now clear precin | 3.1.0 | **Screen-and-treat averts more cancers.** Directly inflates the Fig 2/3/4 effect sizes — the paper's headline |
| `cancer_incidence_with_hiv` / `_no_hiv` now annual, female denominator | 3.2.0 | At `dt=0.25` these read several times too low in v2.2.6. Fig 1 and Fig S2 panels change scale; calibration targets must be refitted |
| HIV stocks scale-weighted | 3.2.0 | v2.2.6 over-reported ~6x at `ms_agent_ratio=100`. The HIV calibration targets were fitted against wrong values |
| `transm2f` default 3.69 → 2.0 | 3.1.0 | Recalibrate |
| `layer_probs` / `f_cross_layer` / `m_cross_layer` → annual probabilities | 2.3.0 | `run_calibration.py` overrides both cross-layer pars. Convert with `1 - (1 - p)**dt` or transmission collapses |
| Vaccine immunity now sterilizing | 2.3.0 | Vaccine efficacy differs |
| `pop_scale > 1` by default when `location` is set | 3.1.0 | Absolute case counts shift unless `total_pop=n_agents`. Decide which convention the paper reports |
| Interventions default to `sex='f'` | 3.1.0 | Check `interventions.py` |
| Starsim RNG, no shared stream with v2 | 3.0.0 | Never expect bit-identical results; compare on overlapping intervals |

### 2c. API port

- `hpv.Sim`: `end=` → `stop=`; no positional pars dict; `'hr'` → `hi5`/`ohr`.
- `hpv.Calibration`: `datafiles=` / `genotype_pars=` / `hiv_pars=` are removed. Signature is now `(sim, calib_pars, *, data=None, ...)` with nested `calib_pars`. `run_calibration.py` needs rewriting, not patching.
- `hpv.AgeResults` → `hpv.by_age`.
- Results are per-module (`sim.results.all_hpv.*`) rather than one flat dict.
- HIV setup collapses to `hpv.Sim(model_hiv=True)`; `pip install hpvsim[hiv]` now that stisim is optional.
- `hpv.MultiSim` / `hpv.save` / `hpv.load` / `sim.short_summary` are gone.
- Read `hpvsim/docs/migration.qmd` before starting.

### 2d. Making the success criterion testable

`compare_baselines.py` already exists for exactly this and takes `--baselines v2.2.6_baseline v3.2.0_baseline`. Use it as the gate rather than eyeballing figures.

Proposed pass conditions, per scenario:

1. **Ordering preserved** — the ranking of the 7 strategies by cumulative cancers averted is unchanged. This is the paper's actual claim.
2. **Elimination years** within ~5 years of v2.2.6, and the qualitative gaps hold: lesion-regressing TxV at 70% remains the earliest, status quo remains ~2080.
3. **Cumulative cancers** overlap on 10–90% intervals; where they do not, the shift is attributable to a row in the 2b table.
4. **Calibration quality** at least as good as v2.2.6 on the shared targets.

Condition 1 is the one that matters. Expect 3 to fail in a specific direction — the precin-clearing fix means screen-and-treat scenarios should avert *more*, so the paper's conclusions get stronger, not weaker. Worth confirming rather than assuming.

### 2e. Sequence

1. Port `run_sim.py` first; get a single sim running to 2025.
2. Resolve 2a, upstream if we take the hpvsim route.
3. Rewrite `run_calibration.py` for the new API. Convert the cross-layer pars.
4. Recalibrate on the VM. Do not run this locally.
5. Port `run_scenarios.py` + plot scripts; freeze `results/v3.2.0_baseline/`.
6. Run `compare_baselines.py`; write the comparison into the PR body.

---

## Step 3 — Verify calibration, merge migration

**Branch:** same `migrate-v3.2` → PR → merge.

1. Check the recalibrated fit against Figure S2 targets: cancer by age, cancer by HIV status, genotype distributions in LSIL and cancer, HIV prevalence/infections/deaths/ART.
2. Run the step 2d gate; record which of the four conditions passed.
3. Fill in **Table S1** — prior distributions and posterior mean/95% intervals across the 50 best-fitting sets. It is currently an empty skeleton in the SM, and it is Reviewer 1's comment 1. The recalibration produces exactly this, so generate it here rather than as a separate task later.
4. Update README: version pin, baseline folder, install line.

---

## Step 4 — Reviewer response and revision plan

**Branch:** `reviewer-response`.

Produce a point-by-point response against the numbered comments in `docs/Reviewer Comments.md` (R1 1–4, R2 1–11), then a revision plan.

Grouping the 15 comments by what they actually require:

| Theme | Comments | Work |
|---|---|---|
| **Sensitivity analysis** | R1 2, R1 3, R2 6 | The dominant ask. R1 3 (vary TxV introduction 2030→2050) is cheapest against existing scenario machinery. R2 6 (VIA sensitivity, 25% LTFU) reuses the same harness |
| **Parameter table** | R1 1 | Table S1, delivered in step 3 |
| **Cost / affordability** | R2 2, R2 8 | R2 8 asks only for a threshold cost per DALY averted, not a full CEA. Tractable |
| **Result tables** | R2 7 | Appendix tables of numbers requiring intervention |
| **Workforce implications** | R2 5 | Current ablation volume in Rwanda + provider count. Needs a local data source, not a model run |
| **Framing / text only** | R1 4, R2 1, R2 3, R2 9, R2 10, R2 11 | No new runs |
| **TxV assumptions need references** | R2 4 | Citations for the two vaccine profiles |

Notes:

- **Sequencing.** Everything in the sensitivity row is new scenario runs, so it should be built on the v3.2 code after step 3 — not against v2.2.6. That is the main reason to do the migration before the revisions.
- **R1 3 and R2 1 are the same underlying objection** (a 2030 TxV is implausible; who is this for?). One sensitivity analysis plus a framing change answers both.
- **R2 5** may not be answerable from the model at all; treat as a caveat in the discussion if the Rwanda ablation-volume data is not available.
- Decide up front which comments get a *response only* versus *new analysis*. Current read: 6 text-only, 4 new runs, 3 new tables, 2 needing external data.

---

## Open questions

1. Sex-stratified HIV prevalence — hpvsim (recommended) or Rwanda-local? Gates step 2.
2. `total_pop` convention: keep v2.2.6 absolute counts (`total_pop=n_agents`) or move to real-population scale? Affects every headline number in the abstract.
3. Turn on HPV latency (`hpv_control_prob`, new in v3.2)? It changes cancer burden substantially and has never been fitted to data. Recommend leaving off for this paper, and noting it as a limitation.
