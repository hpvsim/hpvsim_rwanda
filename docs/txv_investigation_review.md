# Adversarial review of txv_investigation.md

Reviewer: Claude Code (Opus 4.7), 2026-09-14.
Method: recomputed every rep-0 bucket count via `ugrep -c` against
`cancer_rows.csv`; cross-checked mechanism claims against `interventions.py`,
`tx_assigner*.csv`, `txvx_pars_cin.csv`, `hpvsim/data/products_tx.csv`, and
`run_scenarios.py`.

## Verdict summary

- **Q1 (5K averted, screen_change_year=2027 as dominant cause):** Arithmetic
  and mechanism verified. The sensitivity at change_year=2020 only recovers
  the 3-rep median to 7.4K, still short of v2's 9K — the ramp date is a
  partial explanation, not the sole cause.
- **Q2 (VIA-triage LTFU explains the 7K delta):** Mechanism correct, but
  the specific arithmetic "bucket-3 difference fully accounts for the Q2 gap"
  is **overstated**. Bucket 3 explains 76 % of rep-0 delta (4,265 of 5,629
  uids); buckets 1 and 2 contribute the remaining 24 %.
- **Q3 (`txv_pars='cin'` shuts off ablation post-2030):** Verified in
  `interventions.py:79`. Efficacy arithmetic has a minor error (omits the
  `triage_prob=0.9` on `assign_treatment`), but the qualitative rank
  survives.
- **Q4 (HPV-Faster decomposition):** Two problems. **(a) A stale number:**
  the report cites the vax_only median as "~12K (2-rep interim)"; the sweep
  has since finished — rep 2 is present in `summary.csv` and the correct
  3-rep median averted is **10.8K, not 12.3K**. **(b) The birth-cohort
  narrative is wrong:** HPV-Faster's advantage over S&T at rep 0 comes
  primarily from the **1980 (~3.0K uids) and 1990 (~1.9K uids) birth
  cohorts** — women aged 37-47 in 2027, i.e. *inside* the sustained
  screening age range. The 20-29-in-2027 (2000 cohort) contributes only
  ~350 uids of the differential. The report's "sustained screening
  never touches" story is not the dominant driver.
- **Overall:** 3-rep sample is thin; rep 0 is a high-cancer outlier
  (~2× reps 1 and 2 in sim totals) that inflates every bucket claim
  quoted. The report's mechanism story is largely correct; the
  quantitative decomposition should be treated as directional, not point.

---

## Per-question critique

### Q1: S&T&T 70% averts only 5K (v2 was 9K); driven by 2027 ramp

**Verified/false:** Mechanism verified; magnitude partially explained.

**Independent recomputation.**

`sTT18` rep 0 sim total = 81,991; `sTT70` rep 0 sim total = 76,666. Δ = 5,325.
Matches `q1_change_year.csv` `avert_2027` for rep 0 exactly.

Rep-0 bucket recount (my `ugrep -c` numbers all match the report table exactly):

| bucket | sTT18 | sTT70 | Δ |
|---|---|---|---|
| 1_never_screened | 22,595 | 11,883 | −10,712 |
| 3_positive_no_treatment | 3,016 | 9,032 | **+6,016** |
| 4_treatment_failed | 22 | 106 | +84 |

Ramp-date sensitivity (`q1_change_year.csv`) — I recomputed the medians:

| change_year | rep 0 | rep 1 | rep 2 | median |
|---|---|---|---|---|
| 2027 | 5,326 | 3,637 | 2,803 | 3,637 |
| 2025 | 8,718 | 4,024 | 7,816 | 7,816 |
| 2020 | 12,984 | 5,358 | 7,375 | 7,375 |

**Missing considerations.**

1. **The ramp date is not sufficient to reach v2's 9K.** Even at change_year
   2020, the 3-rep median is 7.4K averted — still 1.6K short of v2's 9K.
   The report's TL;DR framing "change_year=2027 accounts for ~2K vs older
   baseline (2025)" is arithmetic-correct (3.6K → 7.8K = +4.2K, or 3.6K →
   7.4K at 2020 = +3.8K on the median), but does not close the whole gap.
2. The report acknowledges the coverage-formula quirk (89 %/91 %
   ever-screened at "70 %") but does not integrate it into the causal
   story. If the intent were to model 70 % lifetime coverage per woman,
   the per-year rate should be higher (or vice versa: the "70 %" arm is
   already over-delivering, so the incremental benefit ceiling is
   ~9 %-of-population-not-yet-reached, which is a much stricter bound
   than the report acknowledges).
3. Rep 0 is a **high-outlier**. sTT18 sim totals: 81,991 / 42,469 / 55,714
   across reps 0/1/2. Rep 0 is 47 % higher than rep 1. Bucket claims
   presented for rep 0 are not typical of the ensemble.
4. **Sample-size issue for the sensitivity:** the median across only 3
   reps is fragile. The 2025 median (7,816) is higher than the 2020
   median (7,375) — implying moving the ramp *earlier* averts *fewer*
   cancers on rep 2 (7,816 → 7,375), which is unphysical and points to
   noise dominating signal in individual reps.

**Remaining uncertainty.** What was `screen_change_year` in v2.2.6 is
still unknown. Until that is verified, "2027 ramp date accounts for 2K of
the 4K gap" is a conjecture. The report flags this appropriately as
residual uncertainty — it should not be presented in the TL;DR as
resolved.

### Q2: Dropping VIA triage saves 7K → attributed to bucket 3

**Verified/false:** Direction correct. Attribution "bucket-3 fully
accounts" is overstated.

**Independent recomputation.**

Rep-0 bucket-by-bucket sTT70 − sT70 delta (in uids):

| bucket | sTT70 | sT70 | Δ |
|---|---|---|---|
| 1_never_screened | 11,883 | 11,179 | +704 |
| 2_screened_always_negative | 4,563 | 3,850 | +713 |
| 3_positive_no_treatment | 9,032 | 4,767 | +4,265 |
| 4_treatment_failed | 106 | 146 | −40 |
| 5_treated_then_reinfected | 8 | 15 | −7 |
| 6_treated_but_progressed_anyway | 7 | 13 | −6 |
| **Total uids** | 25,599 | 19,970 | **+5,629** |

Bucket 3 accounts for **4,265 / 5,629 = 76 %** of the rep-0 uid delta;
buckets 1 and 2 add another 25 %. The report claims "bucket difference
fully accounts for the Q2 gap" and cites 4,265 × 2.98 ≈ 12.7K weighted
cancers vs a sim-level gap of 16,835 (sTT70 76,666 − sT70 59,831).
That ratio is also 76 %, not 100 %.

The remaining +704 in bucket 1 (never_screened) is unexpected — the two
scenarios use the same screening intervention with the same seed within
a rep, so never-screened counts should be identical modulo secondary
knock-on effects (e.g. mortality differences shifting the at-risk pool).
Not investigated by the report.

**Multi-rep check.** Bucket 3 delta (sTT70 − sT70):
- rep 0: 9,032 − 4,767 = 4,265
- rep 1: 4,279 − 2,208 = 2,071
- rep 2: 5,850 − 2,944 = 2,906
- Median: 2,906 uids, i.e. **8,655 weighted cancers, not 12,700**.

Report's headline "4,265 × 2.98 ≈ 12.7K" is a rep-0-only calculation.
The median across 3 reps gives ~8.7K, closer to the actual
median-of-medians cancer delta (52,530 − 59,583 = 7,053 in the 10-rep
sweep). So the effect size lines up better once you average, but the
report's specific number is inflated by taking rep 0 alone.

**Mechanism.** The `tx_assigner.csv` vs `tx_assigner_no_triage.csv`
rows are exactly as the report describes (precin 0.3 vs 1.0; cin 0.6
vs 1.0). Verified.

**Missing considerations.**

- The report ignores the `triage_prob=0.9` on `assign_treatment`
  (`interventions.py:81`). Even in S&T (no-triage), 10 % of HPV+ women
  fail the triage-arrival step and land in bucket 3. That's likely why
  sT70 still has 4,767 in bucket 3 rather than ~0.
- The report says "S&T doesn't have triage-stage LTFU" — this is wrong.
  S&T has *identical* triage LTFU (via `triage_prob=0.9`); what it
  removes is the *routing* step (VIA judging most precin women "not
  eligible"). The distinction matters because the "hidden" 10 %
  triage-arrival LTFU is invariant across scenarios.
- Extra ablations per averted cancer = 164 verified (Δablations
  1,159,113 / Δcancers 7,053).

**Remaining uncertainty.** v2's 12.5K S&T advantage vs v3's 7.1K
advantage is not resolved. The report speculates about `treat_num`
prob=1.0 in v2 vs 0.75 in v3 without confirming from v2.2.6 source.

### Q3: Adding TxV to S&T adds nothing at 70 %

**Verified/false:** Mechanism claim (`triage_end_year = min(2030,
end_year) if txv_pars=='cin'`) verified in `interventions.py:79`. TxV
efficacy (precin 0.5, cin 0.936) verified in `txvx_pars_cin.csv`.

**Independent recomputation.**

Rep-0 bucket data as in the report — I re-checked bucket 4:

| scenario | 4_treatment_failed | 3_positive_no_treatment |
|---|---|---|
| sT70 | 146 | 4,767 |
| sTxV70 | 1,012 | 3,938 |
| Δ (sTxV − sT) | **+866** | **−829** |

Confirmed: the bucket-4 excess in sTxV (+866) nearly cancels the
bucket-3 saving (−829). Net cancer delta at rep 0 is small.

**Multi-rep check.** Bucket-4 excess of sTxV over sT:
- rep 0: 1,012 − 146 = 866
- rep 1: 458 − 66 = 392
- rep 2: 707 − 107 = 600

TxV consistently produces more treatment-failure cancers across all
reps. Median = 600, not 866. Consistent direction but ~30 % lower
magnitude than rep 0.

**Missing considerations.**

- The efficacy arithmetic in the report (bullet 3: "S&T ablation …
  0.75 × 0.93 = 0.70") **omits the `triage_prob=0.9` factor** on
  `assign_treatment`. Correct arithmetic:
  - S&T (post-2020, precin): 0.9 × 1.0 × 0.75 × 0.936 = **0.632**
  - S&TxV (post-2030, precin): 0.9 × 0.5 = **0.45**
  - S&T (cin): 0.9 × 1.0 × 0.75 × 0.936 = **0.632**
  - S&TxV (cin): 0.9 × 0.936 = **0.842**

  Qualitative rank (S&T-precin > TxV-precin; S&T-cin < TxV-cin) is
  preserved, so the conclusion "precin loss dominates" still stands.
  But the specific ratios in the report (0.70 vs 0.45) understate S&T's
  strength.
- The report **does** flag the confounder in bullet 6: "S&T has a
  10-year head start of 'everyone gets ablated' while S&TxV is still
  LTFU'ing on VIA triage". This is a real disadvantage for S&TxV that
  is separate from the post-2030 efficacy question, and the report
  correctly identifies it. Good.
- The report does not quantify what fraction of the 2K sTxV70 vs sT70
  delta comes from the 2020-2030 head start vs the post-2030 efficacy
  swap. A cheap counterfactual (TxV on from 2020) would decompose
  this, and the report recommends running it.
- The report does not check whether the TxV-cin efficacy of 0.936
  applies *per dose* or *per course*. This matters — if it's per dose
  and multiple doses are administered, effective per-course efficacy
  could be higher. Unverified.

**Remaining uncertainty.** Report itself says medium-high confidence.
The mechanism story is solid; the exact numeric decomposition should be
treated as illustrative until the "TxV from 2020" and "TxV + parallel
ablation" sensitivity runs land.

### Q4: HPV-Faster averts 33K = ablate_only 17.6K + vax_only 12K, minus overlap

**Verified/false:** Mechanism verified (three structural asymmetries:
2027 one-shot, 20-50 age range, adult catch-up prophylactic vax).
Quantitative decomposition has one **stale figure** and the
**birth-cohort narrative is wrong**.

**Independent recomputation.**

Vax_only medians. The report says "rep2 pending" and quotes an interim
median of 12,265 (mean of reps 0-1 differences). But `summary.csv`
shows rep 2 for `hpvfaster70_vax_only` is done: sim total = 44,929.
Correct 3-rep averted values vs baseline:
- rep 0: 81,991 − 66,366 = **15,625**
- rep 1: 42,469 − 33,563 = **8,906**
- rep 2: 55,714 − 44,929 = **10,786**
- **Median = 10,786** (not the report's "interim ~12,265")

Full HPV-Faster averted (3-rep median):
- rep 0: 81,991 − 42,632 = 39,359
- rep 1: 42,469 − 22,494 = 19,975
- rep 2: 55,714 − 28,967 = 26,747
- Median = 26,747. Matches report.

Ablate_only averted (3-rep median):
- rep 0: 81,991 − 54,076 = 27,915
- rep 1: 42,469 − 28,928 = 13,541
- rep 2: 55,714 − 38,084 = 17,630
- Median = 17,630. Matches report.

**Sum-of-parts vs whole per-rep:**

| rep | ablate_only | vax_only | sum | full | overlap |
|---|---|---|---|---|---|
| 0 | 27,915 | 15,625 | 43,540 | 39,359 | 4,181 (10.6 %) |
| 1 | 13,541 | 8,906 | 22,447 | 19,975 | 2,472 (12.4 %) |
| 2 | 17,630 | 10,786 | 28,416 | 26,747 | 1,669 (6.2 %) |

Overlap is 6-12 % — the report's "10 % overlap" claim is fine.

Updated decomposition (report's "two-thirds mass ablation, one-third
adult catch-up vax"): 17,630 (ablate) / 26,747 (full) = **66 % ablation
contribution**, 10,786 (vax) / 26,747 (full) = **40 % vax
contribution**. The "one-third vax" claim should really be "one-third
to two-fifths vax". Not a big shift, but the specific 15,625/12K
number quoted for vax_only is a stale interim.

**Median-of-medians vs paired median mismatch (still unexplained by
report).** The 10-rep paired sweep gives 33,177 averted; the 3-rep
diagnostic gives 26,747. If the 3 reps in the diagnostic were the
top-3-mismatch calibration trials (i.e., a subset of the 10 the main
sweep used), then reps 3-9 in the main run must average ~35-40K
averted to lift the 10-rep median to 33K. Rep 0 is 39K, so this is
plausible but unconfirmed. The report handwaves this ("rep 4 in the
10-rep run averted 40K") without evidence.

**Birth-cohort claim is wrong.** The report says HPV-Faster's advantage
is that it "reaches 20-29yo cohorts sustained screening never touches"
and that "the ~1997-2007 birth cohort is 20-30 in 2027 and thus in
HPV-Faster's target but not (yet) in S&T's target. This cohort is at
peak HPV incidence."

Recount of rep-0 cancers by birth cohort (unweighted uids):

| birth cohort | age@2027 | sTT18 | sT70 | hpvfaster70 | HF−S&T |
|---|---|---|---|---|---|
| 1960 | 67 | 1,761 | 1,747 | 1,743 | −4 |
| 1970 | 57 | 4,672 | 4,146 | 3,524 | −622 |
| 1980 | 47 | 7,879 | 6,170 | 3,161 | **−3,009** |
| 1990 | 37 | 7,848 | 4,861 | 2,989 | **−1,872** |
| 2000 | 27 | 2,147 | 1,325 | 972 | −353 |
| 2010 | 17 | 1,015 | 501 | 299 | −202 |

Sum of HF-vs-S&T advantage: ~6,062 uids (× 2.98 = ~18,065 weighted
cancers). That matches the sim-level delta (59,831 − 42,632 = 17,199)
within ~5 % rounding.

**Where the delta comes from:**
- 1980 cohort (aged 47 in 2027): **−3,009 uids = 50 % of HF's advantage**.
  These women *are* in the sustained-screening 30-50 range for the full
  20-year window.
- 1990 cohort (aged 37 in 2027): **−1,872 uids = 31 % of the advantage**.
  Also fully inside the sustained-screening window.
- 2000 cohort (aged 27 in 2027): only −353 uids = 6 % of the advantage.
  This is the "20-29 in 2027" cohort the report claims is the driver.
- 2010 cohort (aged 17 in 2027): −202 uids = 3 % of the advantage.
  These get HPV-Faster at 2027 but would eventually enter routine
  screening in 2033+.

**So the dominant driver of HPV-Faster's edge over S&T is *not*
reaching the un-reached 20-29 birth cohort.** It is one-shot, no-LTFU
treatment of the 30-49-year-olds who *are* in sustained screening's
range but who S&T only sees over 20 years of annual lottery (with
25 % ablation LTFU each visit).

The report's mechanism claim (a) "2027 one-off captures 20-29yo whom
sustained 30-50 screening never touches" contributes only ~6-9 % of the
actual HF-vs-S&T advantage in rep 0. The dominant mechanism is the
100 % treatment probability applied one-shot to the entire eligible
population, which recovers years of accumulated ablation LTFU in the
40s cohort.

**3.4M adult vax number verified.** Baseline vaccinations = 14,339,768;
HPV-Faster 70 % = 17,725,071. Δ = 3,385,303 doses. Correct.

**Excision count sanity check flagged, not resolved.** HPV-Faster
excisions = 26.7M, ablations = 3.1M. Report speculates it's a retry
loop but does not resolve. If the excision counter is inflated by
per-timestep re-eligibility of failed ablation cases (via
`excision_eligible: triage_out | abl_fail`), then the excision total
cannot be interpreted as discrete women treated. This is a real code
concern that should be traced. Report correctly flags but does not
resolve.

**Missing considerations.**

- The report does not verify whether `hpvfaster70_vax_only` correctly
  models "no ablation" (it may still have baseline S&T ablations
  running underneath, since HPV-Faster is layered on top of `make_st`
  per line 246-249). Without reading `diagnose_all.py`'s
  `vax_only`/`ablate_only` variant definitions, one cannot rule out
  that the components include baseline S&T noise, inflating both
  components.
- No control for MS's specific quantitative claim: MS says
  "vaccinating and screening women aged 20-50 could avert 22,700
  cases from a single-year intervention" — v3 HPV-Faster averts
  33,000. The v3 number is ~50 % higher than the MS abstract. The
  report notes this in the "v2 vs v3 delta explanation" but does not
  reconcile whether the MS number is v2.2.6 or an earlier revision.

**Remaining uncertainty.** The 3-rep decomposition is fragile (rep 0
gives ~30 % higher averted than rep 2 for every scenario). The
component "vax_only" rep 2 that the report treats as pending is
actually done, and the corrected 3-rep median (10.8K) is 12 % lower
than the report's interim (12.3K).

---

## Methodological concerns

1. **3 reps is too few.** Diagnostic-sweep cancer counts vary by ~2×
   across reps 0 vs 1 (rep 0 is systematically the highest). All bucket
   claims in the report tables are rep-0 values, presented without
   equivalent rep 1/2 recounts. The multi-rep summary I ran shows the
   median is often 30-40 % below rep 0 for the same quantity. **The
   report should either report all three per-rep values or use the
   median explicitly.** Right now the numbers read as authoritative but
   are single-rep snapshots.

2. **Median-of-medians ≠ paired median ≠ per-rep median.** Report
   correctly notes this early (para 55) but then blurs it in Q1
   (comparing "5K" from paired to "5,326" from rep 0 to "5.3K" from
   uid-weighted rep 0). Values within 5 % of each other doesn't mean
   they measure the same thing.

3. **The report has an internal inconsistency.** TL;DR says "the 5K
   averted matches Fig 5 exactly". But Fig 5's paired median is 5,436;
   rep 0's sim-total delta is 5,326; the report's uid-weighted rep-0
   estimate is 5,290. These three numbers happen to be close, but the
   "matches exactly" claim is a coincidence of noise, not a proof that
   the mechanism is fully explained.

4. **No sensitivity on triage_prob.** The 10 % LTFU at the
   `assign_treatment` step (via `triage_prob=0.9`) affects all scenarios
   equally in principle, but its interaction with the post-2030 TxV
   switchover in S&TxV is unquantified.

5. **Coverage formula quirk (89-91 % actual vs 70 % intended) is
   acknowledged but not used to reframe results.** If "70 %" over-
   delivers to 91 % ever-screened, then the S&T&T 70 % arm's *headroom*
   for additional averted cancers is tiny (~9 % of women), which is a
   stronger explanation for Q1 than the 2027 ramp. The report gives it
   third billing.

6. **The `hpvfaster70_ablate_only` and `hpvfaster70_vax_only` variant
   definitions are not audited.** Both are cited as evidence of
   decomposition, but without reading `diagnose_all.py`'s treatment of
   these variants, I cannot confirm they cleanly separate the two
   effects. A run that "removes ablation" but leaves in the mass vax
   might still have baseline S&T running underneath.

---

## Claims you cannot verify from the available evidence

- What `screen_change_year` was set to in v2.2.6 (report acknowledges
  this).
- Whether v2 `treat_num` used `prob=1.0` or something else.
- Whether the 26.7M HPV-Faster excisions represent 26.7M discrete
  women-events or repeated events on the same women. The
  `treat_num.step()` code visible does not obviously double-count, but
  the `excision_eligible: triage_out | abl_fail` logic could keep a
  woman eligible across timesteps until she clears — this requires
  reading the underlying starsim `treat_num` implementation, which the
  environment does not allow me to access.
- Whether the MS's "22,700 averted from HPV-Faster" and "16K averted
  from TxV" numbers were generated from v2.2.6 code exactly, or from
  an earlier v2 revision. If different, some of the v2/v3 deltas in
  the report's cross-comparison table are apples-to-oranges.
- Whether the 3 diagnostic reps are actually the same top-3 trials as
  the first 3 of the 10-rep main sweep. If they are, the median-of-3
  should sit somewhere in the middle of the 10-rep distribution — but
  since we see 3-rep median averted = 27K vs 10-rep median = 33K, the
  10-rep run's other 7 reps must include some very-high-avert draws.
  This is not implausible but is not shown.

---

## What the researcher should trust vs re-check

**Trust:**

- Cumulative-metrics numbers in the tables sourced from
  `scens_cumulative.csv`. I re-checked several and they all match.
- Interventions.py logic descriptions: `triage_end_year = min(2030,
  end_year) if txv_pars=='cin'`; per-year screening-probability
  formula; `future_treat_cov=0.75`; HPV-Faster older-cohort `prob=1`.
  All verified against source.
- TxV efficacy shape (precin 0.5, cin 0.936, latent 0). Verified.
- Ablation efficacy (0.936 for both precin and cin). Verified against
  `hpvsim/data/products_tx.csv`.
- 3.4M extra vaccinations in HPV-Faster. Verified.
- Extra-ablations-per-averted-cancer for Q2 (164). Verified.
- Coverage-formula arithmetic ("70 %" → 91 % ever-screened over
  30-50). Verified.
- q1_change_year sensitivity numbers. Verified from CSV.

**Re-check before using in publication:**

- **Vax_only median for Q4:** the correct 3-rep median averted is
  **10,786, not 12,265**. Every downstream percentage ("one-third",
  "40 %") should be recomputed from the corrected number.
- **Q2 attribution:** bucket-3 explains 76 % of the Q2 gap in rep 0
  and about the same in the median, not "fully".
- **Birth-cohort mechanism narrative in Q4:** the "20-29 in 2027" story
  is not supported by the actual birth-cohort distribution. The
  dominant contributor is the 1980 and 1990 cohorts (37-47 in 2027) —
  i.e. women who *are* in sustained screening's range but benefit from
  HPV-Faster's one-shot 100 %-treatment. The report should rewrite
  this section around "one-shot no-LTFU treatment of the 30-49-yo
  population" rather than "reaching cohorts sustained screening
  misses".
- **Q3 efficacy arithmetic:** include `triage_prob=0.9` factor to get
  correct per-course probabilities. Direction unchanged but the
  numbers "0.70 vs 0.45" understate S&T's effective coverage.
- **Statement that S&T has no triage-stage LTFU** (Q2, Q3): incorrect
  — S&T has identical 10 % `triage_prob` LTFU as S&T&T, it only
  removes the VIA routing. The distinction matters if anyone asks
  what "removing VIA" changes vs "removing LTFU".

**Investigate further before publication:**

- Excision retry loop (Q4 caveat): 26.7M excisions on 3.1M ablations
  suggests repeated events. Not resolved. Read starsim `treat_num` and
  the eligibility logic on failed-ablation carryover.
- Confirm v2.2.6 `screen_change_year` and `future_treat_cov` values
  before claiming those are the sources of v2/v3 divergence.
- Run 10-rep diagnostic sweep before quoting any percentage-attribution
  figures.
- Run `analyze_diag.py` to generate `birth_cohort_medians.csv` and
  `bucket_medians.csv` — currently these files are referenced in the
  report but do not exist.
- Confirm that the `hpvfaster70_ablate_only` and
  `hpvfaster70_vax_only` variants in `diagnose_all.py` truly isolate
  the components (a run that only removes the mass-vax add-on, versus
  a run that only removes the older-cohort ablation).

---

## One-line takeaway

The mechanisms identified in the report are structurally right, but
the specific numbers to quote back at reviewers (especially the Q4
decomposition and birth-cohort narrative) are wrong in ways that
matter. Rerun with the finished vax_only rep 2 data, produce the
missing `analyze_diag.py` outputs, correct the birth-cohort story, and
soften "fully accounts" language to "primarily accounts (~76 %)" and
the report is publishable as an internal diagnostic. As a manuscript
appendix or reviewer response, it needs one more pass.
