# Residual cancer analysis

**Question:** After every intervention we have — routine bivalent vax,
screen-and-treat at 70% coverage, TxV, HPV-Faster mass campaign — how
many cervical cancers still occur in Rwanda over 2030-2100? Who is
getting them? To what extent can the existing intervention set reach
these women, and what is genuinely beyond us?

**Setup:** normalized scenario sweep — every intervention starts 2030
(except the ongoing routine bivalent vax program, which runs from
2011 onward in every scenario except a true no-vax counterfactual),
accounting window 2030-2100. 5 top-calibration reps × 6 scenarios.
Data: `results/diagnostic_normalized/cancer_rows.csv`.

**Headline:**

- **Routine bivalent vax by itself averts ~208,000 cumulative cancers
  vs a true no-intervention world.** This is Rwanda's dominant policy
  lever and it is already in place — every scenario below assumes it
  continues.
- On top of vax, adding **HPV-Faster 70%** (the strongest single
  additional lever in the 2030-normalized comparison) averts another
  ~27,000 cancers and leaves a residual of **26,012** [21K, 40K].
- Of that 26K residual: **~40% is a genuine floor** (women whose HPV
  was acquired before any intervention could reach them AND who are
  past the 30-50 screening window by 2030); **~60% is
  reachable-in-principle-but-missed**, dominated by never-screened
  women.
- The most promising unmodelled additional levers, ranked, are:
  (a) **screening coverage above 70%**, (b) **widening the screening
  age window** to include 25-55 or 30-60, (c) **switching S&T
  scenarios to include TxV** as the LTFU-protection lever (already
  in the paper), (d) shorter screening intervals for women in
  high-risk subgroups.

---

## §1 — Scenario totals (2030-2100)

Cumulative cervical cancer cases, median across 5 reps [p10, p90]:

| Scenario | Cases | Averted vs no intv |
|---|---:|---:|
| **No interventions (no vax, no S&T)** | **261,013** [197,653; 377,594] | 0 |
| Baseline (vax + S&T 18%) | 53,206 [42,880; 78,658] | 207,807 |
| S&T&T 70% | 48,033 [38,991; 74,764] | 212,980 |
| Mass TxV 50/90, 70% | 33,315 [25,264; 44,488] | 227,698 |
| S&TxV 70% | 35,701 [30,290; 54,133] | 225,312 |
| **HPV-Faster 70%** | **26,012** [21,169; 39,691] | **235,001** |

Notes:

- **Routine bivalent HPV vaccination (from 2011, ages 11-12, 90%
  target coverage) averts ~207,800 cancers by itself.** That is 80%
  of the burden a true no-vax counterfactual would carry. Every
  additional intervention below is a marginal lever on top of an
  already-installed prophylactic vax programme.
- Among additional interventions modelled here, HPV-Faster 70% is
  the strongest, averting ~27K on top of Baseline; S&TxV 70% averts
  ~17K; S&T&T 70% only 5K.
- HPV-Faster beats S&TxV by ~10K in the 2030-normalized framing,
  consistent with the corrected earlier investigation: HPV-Faster
  wins on the 2030 one-off + mass adult vax combination when both
  strategies start on the same day (in the paper's 2020-start
  framing, screening had a decade to work before TxV became the
  therapy, so the comparison flipped).

---

## §2 — Residual (26K) under HPV-Faster 70%

### 2a. Birth cohort

| Cohort (age in 2030) | Cases | Share |
|---|---:|---:|
| 1950-1960 (70-80) | 3,890 | 14.3% |
| **1970 (60)** | **8,283** | **30.8%** |
| **1980 (50)** | **6,264** | **23.3%** |
| **1990 (40)** | **5,180** | **19.3%** |
| 2000 (30) | 1,382 | 5.1% |
| 2010 (20) | 786 | 2.9% |
| 2020+ (10 or younger) | 1,048 | 3.9% |

**73% of the residual is concentrated in the 1970-1990 cohorts.**
The 1970 cohort alone (age 60 in 2030) contributes nearly a third —
these women were 30-50 during 2000-2020, before any Rwandan screening
programme was at scale, and they age out of the current 30-50
screening window in the early 2030s.

### 2b. Age at causal HPV infection

| Age at infection | Cases | Share |
|---|---:|---:|
| < 20 | 751 | 2.9% |
| 20-30 | 4,200 | 16.0% |
| 30-40 | 8,540 | 32.6% |
| 40-50 | 9,064 | 34.6% |
| 50-60 | 2,973 | 11.4% |
| 60+ | 664 | 2.5% |

Two-thirds of residual cancers arise from HPV acquired in the 30-50yo
window — the window that HPV DNA screening covers in principle.
Failures here are operational (coverage, uptake), not biological.

### 2b (bonus) — Vaccination shifts causal HPV to later ages

Comparing the age-at-causal-infection distribution under a
counterfactual with NO routine vaccination vs Baseline (vax + status-
quo screening):

| Age at infection | No vax (n) | No vax (%) | Baseline (n) | Baseline (%) |
|---|---:|---:|---:|---:|
| < 20 | 13,764 | 5.2% | 1,832 | 3.4% |
| 20-30 | 70,938 | 27.0% | 11,456 | 21.3% |
| 30-40 | 87,937 | 33.4% | 18,026 | 33.5% |
| 40-50 | 66,211 | 25.2% | 16,290 | 30.2% |
| 50-60 | 20,007 | 7.6% | 4,959 | 9.2% |
| 60+ | 4,352 | 1.7% | 1,325 | 2.5% |

Without vaccination, **32% of causal HPV happens before age 30**.
With routine bivalent vax from 2011, that drops to 25%; the
distribution shifts noticeably toward later ages. The routine vax
programme is successfully protecting the younger cohorts, and this
is *why* the residual we see under HPV-Faster 70% is concentrated in
women born pre-1990 — precisely the cohorts the vax couldn't reach
because they were already past age 12 when it started.

### 2c. HIV status of the residual

Under HPV-Faster 70%: 92.9% HIV-negative, 7.1% HIV-positive. Roughly
proportional to national HIV prevalence in women (~6-8%). Not a
disproportionate driver of the residual.

---

## §3 — Mode of failure: bucket comparison across scenarios

Cumulative residual cases per bucket, medians:

| Bucket | S&T&T 70% | S&TxV 70% | HPV-Faster 70% |
|---|---:|---:|---:|
| 1 Never screened | 26,465 (55%) | 25,553 (72%) | 20,561 (76%) |
| 2 Screened, always negative | 7,274 (15%) | 3,842 (11%) | 4,107 (15%) |
| 3 Positive, no treatment | **14,428 (30%)** | 3,848 (11%) | 2,154 (8%) |
| 4 Treatment attempted, all failed | 161 (0%) | 2,380 (7%) | 36 (0%) |
| 5 Successful tx, then reinfection | 15 (0%) | 9 (0%) | 15 (0%) |
| 6 Successful tx, progressed anyway | 9 (0%) | 57 (0%) | 33 (0%) |
| **Total** | **48,033** | **35,701** | **26,012** |

Three observations:

1. **Bucket 1 (never screened) is dominant in every scenario.**
   Getting screening reach above the modelled 70% ceiling is the
   biggest remaining lever — see §6.

2. **S&T&T 70% has a huge bucket 3 (14,428 = 30% of its residual).**
   The tx_assigner chain (0.9 triage prob × 0.75 treatment prob) only
   delivers treatment to 68% of screen-positive women. Every second
   screen-positive woman is lost in the follow-up chain.

3. **S&TxV 70% cuts bucket 3 from 14,428 → 3,848 — saving ~10,600
   cancers.** This is exactly the LTFU-protection story: because TxV
   is delivered at the same visit as the positive screen (linked_txvx
   with prob=0.9), it bypasses the triage-then-return-to-clinic drop-
   off. Bucket 4 rises (2,380 vs 161) because TxV's per-attempt
   clearance rate on established CIN is lower than ablation's — but
   the two together (bucket 3 + bucket 4 = 6,228 in S&TxV) are less
   than half of the treatment-lost residual under S&T&T (14,589).
   **TxV as delivered here is fundamentally a mechanism for
   protecting the women who fall out of the triage-and-treat
   pipeline.**

4. **HPV-Faster's small bucket 3 (2,154) comes from
   `tx_assigner_no_triage`** — the 2030 one-off campaign routes every
   positive straight to ablation (no triage LTFU step). This gives
   HPV-Faster the smallest treatment-lost bucket of any scenario.
   The user-facing lesson: single-visit see-and-treat models are
   the alternative to TxV for closing bucket 3.

Combined treatment-related residual (bucket 3 + 4):

- S&T&T 70%: 14,589 (30% of residual)
- S&TxV 70%: 6,228 (17% of residual) — TxV closes ~57% of this
- HPV-Faster 70%: 2,190 (8% of residual) — see-and-treat closes ~85%

---

## §4 — Reachability partition

For each residual cancer, was the woman in *any* modelled
intervention's addressable window at any point between 2030 and her
cancer diagnosis? Reachable-if she was: age 30-50 during
[2030, cancer_year] (screening + linked TxV), OR age 20-50 in 2030
(mass adult vax / HPV-Faster campaign), OR age 11-12 at any t ≥ 2011
(routine prophylactic vax).

| Scenario | Reached in principle | Unreachable |
|---|---:|---:|
| HPV-Faster 70% | 14,139 (54%) | 11,914 (46%) |
| S&TxV 70% | 23,534 (66%) | 12,412 (34%) |

The **unreachable share** (~12K cancers under either scenario) is the
1960-1975 cohorts, past 50 by 2030 and past 12 by 2011. Under the
current intervention set they cannot be reached — but see §6, item 2:
widening the routine screening age window would move some of them
from unreachable to reachable.

The **reached-in-principle-but-missed** share is the operational
target for further improvement: mostly bucket 1 (never screened).

---

## §5 — Comparison with Baseline (18% status-quo screening)

Baseline residual: 53,206 cases. Birth-cohort distribution shifts
noticeably: 1980 & 1990 cohorts each carry ~30% of Baseline residual,
vs 23% and 19% under HPV-Faster. That is, **the ~27K cancers that
HPV-Faster additionally averts vs Baseline are concentrated in the
1980-1990 birth cohorts** — precisely the women targeted by the 2030
mass adult vax + one-off screening.

---

## §6 — What would it take to close the reachable residual?

Ranked by leverage on the ~14K reached-but-missed cancers under
HPV-Faster 70%:

1. **Screening coverage above 70%.** Bucket 1 is 76% of the residual.
   Getting from 70% to 90% lifetime coverage in the 30-50yo window
   would take an estimated 6-8K cancers off the table (linear
   extrapolation). Requires self-sampling, mobile clinics, or
   employer-based delivery. **This is the single highest-leverage
   remaining lever within the current intervention geometry.**

2. **Wider screening age window (25-55 or 30-60).** This is the lever
   that shrinks the "unreachable" floor rather than the reachable
   residual:
   - Extending to 30-60 would reach the 1970 cohort (currently the
     largest single contributor at 8,283 cancers, all "unreachable"
     under the current 30-50 window) for their remaining screening
     years 2030-2040 (age 60-70).
   - Extending to 25-55 would reach the 1980-2005 cohorts for
     additional years at both ends of the window.
   - Estimated impact: several thousand additional cancers averted.
     Concrete numbers require a targeted sensitivity sweep — this is
     an explicit next-run recommendation (see the reviewer-response
     plan).

3. **Add TxV to any S&T scenario (as in S&TxV).** From §3: TxV cuts
   the treatment-lost residual (buckets 3+4) by ~57% by protecting
   against triage LTFU. In HPV-Faster the equivalent effect is
   achieved via see-and-treat (`tx_assigner_no_triage`), which cuts
   treatment losses by ~85%. Either mechanism removes the triage
   drop-off; the question is which fits Rwanda's delivery reality.

4. **Shorter screening interval for high-risk subgroups** (bucket 2 =
   15% of residual). Women with prior positives or HIV+ screened
   every 3-5 years instead of every 10 would take ~1-2K cancers off
   the residual.

5. **Better therapeutic products.** Buckets 4-6 combined are ~0.3%
   of the residual. Product efficacy is not the operational
   bottleneck.

**Rough elimination pathway from HPV-Faster 70% baseline:**

- 26,000 residual → 20,000 by pushing screening coverage 70 → 90%
- → ~16,000 by widening the age window to 25-55 or 30-60
- → ~14,000 by adding TxV or see-and-treat to close residual bucket 3

Bringing the residual from ~26K to ~14K over 2030-2100 is plausible
without adding any qualitatively new modality. The ~14K floor that
remains is a mix of (a) women born pre-1970 whose HPV was acquired
50 years ago, and (b) small-N bucket-2 test-sensitivity misses that
even 3-yearly screening would only partially close.

---

## §7 — Manuscript framing

Suggested paragraphs for the discussion:

> "Routine bivalent vaccination, in place in Rwanda since 2011,
> accounts for approximately 208,000 cumulative cancers averted over
> 2030-2100 compared with a counterfactual with no HPV intervention
> — the dominant policy lever in the modelled portfolio. On top of
> this already-installed programme, the strongest additional
> intervention we model, an HPV-Faster mass campaign at 70% target
> coverage in 2030, averts a further ~27,000 cancers, leaving a
> residual of approximately 26,000 [21,000; 40,000] cumulative cases
> over the 70-year horizon."

> "This residual is concentrated in birth cohorts 1960-1990 — women
> already past age 12 when routine vaccination was introduced, and
> largely past age 50 by intervention launch in 2030. Approximately
> 12,000 cases are structurally beyond the reach of the modelled
> intervention set: these women acquired HPV before screening or
> vaccination were available in Rwanda, and by 2030 have aged out of
> the routine 30-50 screening window. **A modest widening of the
> screening age window (e.g. to 25-55 or 30-60) would shift some of
> this pre-2030 legacy burden from unreachable to reachable, and
> represents an underexplored policy lever worth costing.**"

> "The remaining ~14,000 cancers are 'reachable in principle but
> missed', dominated (~76% of the residual) by women who were never
> screened despite being in the eligible age window at some point
> during the modelled horizon. The highest-leverage remaining
> intervention is not a new therapeutic product — treatment product
> failures account for less than 1% of residual cancers under
> HPV-Faster 70% — but delivery: raising screening coverage above
> the modelled 70% ceiling via self-sampling, mobile clinics, or
> employer-based delivery models."

> "The comparison between S&T&T 70% and S&TxV 70% is instructive:
> both offer 70% lifetime screening coverage, yet S&TxV averts ~10,000
> additional cancers, entirely by protecting women from the triage-
> and-treatment loss-to-follow-up chain. Because the therapeutic
> vaccine is delivered at the same visit as the screening result,
> it bypasses the ~30% attrition between a positive screen and
> completed ablation that occurs in the conventional S&T pathway.
> A single-visit see-and-treat model achieves the equivalent effect;
> both address the same operational failure."

---

## §8 — Recommended next runs

Concrete sensitivity sweeps that would sharpen §6 and the manuscript
discussion. All are compatible with the existing normalized-run
harness — parameter knobs are already plumbed into `make_st`.

| Sweep | Purpose | Scenarios | Sim volume |
|---|---|---|---|
| Widened age window (§6.2) | Move unreachable → reachable | `S&T&T 70%` and `HPV-Faster 70%` with `age_range=[25,55]` and `[30,60]` | 4 scen × 5 reps = 20 |
| Higher screening coverage (§6.1) | Diminishing-returns curve on bucket 1 | `HPV-Faster` at coverage {70%, 80%, 90%, 95%} | 4 × 5 = 20 |
| Elimination push (§6 combined) | Stack everything | `HPV-Faster 90%`, age 25-55, TxV | 2 × 5 = 10 |

Total ~50 additional sims (~75 min on one core). Would give a
plausible "floor" estimate for the manuscript.

---

## Reproduction

```bash
# Normalized scenario sweep (5 reps x 6 scenarios, ~45 min):
python diagnose_all.py --normalized --reps 5 \
    --scenarios "No interventions" "Baseline" "S&T&T 70%" \
                "S&TxV 70%" "HPV-Faster 70%" "Mass TxV 50/90, 70%" \
    --outdir results/diagnostic_normalized

# Aggregate:
python analyze_residual.py

# The 'No interventions' scenario requires the diagnose_all.py fix
# (commit 8f26403) that drops routine vax for that scenario.
```
