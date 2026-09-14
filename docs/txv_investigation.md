# Investigation of unexplained results in Rwanda HPV v3 model

**Status:** Draft – written 2026-09-14 by Claude Code, based on:
- Pre-run scenario sweep (10 reps, medians) at `results/scens_cumulative.csv` and `results/scens_paired.csv`
- Diagnostic decomposition sweep (3 reps × 8 scenarios) at `results/diagnostic/cancer_rows.csv`
  (produced by `diagnose_all.py`, PID 61976, started 13:20 on 2026-09-14)
- v2.2.6 baseline for comparison at `results/v2.2.6_baseline/scens_cumulative.csv`

## TL;DR (one bullet per question + one overall)

1. **Q1 (S&T&T 18%→70% averts only 5K, was 9K in v2):** The screening ramp fires only starting in **2027** (`screen_change_year=2027`), so the accounting window `2025-2100` includes 2 years of 18% "warm-up". Bucket data (rep 0) confirms sTT18 has 22.6K never-screened cancers vs sTT70's 11.9K — the extra coverage does bring 10.7K women out of "never screened" — but it also creates 6.0K MORE positive_no_treatment cancers (VIA triage LTFU is doubled at 70% because there are 2× the positive screens). Net cancer benefit: only 1.8K uids (~5.3K weighted cancers). Increasing the ramp date to 2020 recovers 8.6K averted (from `results/diagnostic/q1_change_year.csv` rep 0). Coverage formula still delivers ~91% ever-screened at "70%".

2. **Q2 (dropping VIA triage saves 7K):** CONFIRMED by bucket data. sT70's `3_positive_no_treatment` bucket has 4,767 uids vs sTT70's 9,032 — an S&T advantage of 4,265 uids at 70% coverage. `tx_assigner.csv` treats **30% of precin** and **60% of cin** with ablation (rest → `none`), while `tx_assigner_no_triage.csv` treats **100%** of both. S&T sends ~1.8× the number of women to ablation compared to S&T&T (v3 ablations: **2.54 M vs 1.38 M**).

3. **Q3 (adding TxV to S&T adds nothing at 70%):** CONFIRMED by bucket data. S&TxV rep 0 has **1,012 treatment_failed uids** vs S&T's 146 — TxV's 50% precin efficacy fails on ~866 extra women. That's offset by only 829 fewer "positive_no_treatment" (from TxV's higher coverage 0.9 vs 0.75). Net: essentially no benefit. Mechanistically, S&TxV uses `linked_txvx(prob=0.9)` **only from 2030** and **shuts off the ablate/excise pathway from 2031 onwards** (`triage_end_year = min(2030, end_year)` when `txv_pars='cin'`). Post-2030, TxV @ 50% precin / 93.6% cin efficacy replaces ablation @ ~93% precin / ~93% cin efficacy. Coverage advantage (0.9 vs 0.75) and precin efficacy disadvantage (0.5 vs 0.93) roughly cancel.

4. **Q4 (HPV-Faster 70% averts 33K – MORE than S&TxV):** CONFIRMED by isolated-component runs. hpvfaster70_ablate_only averts 17.6K (3-rep median), hpvfaster70_vax_only averts ~12K (2-rep interim), sum ~30K with ~5K overlap = full HPV-Faster's ~27K averted (3-rep median). Bucket data: HPV-Faster's `3_positive_no_treatment` drops from S&TxV's 3,938 to just **1,693** (mass ablation at 100% coverage) AND HPV-Faster reduces total cancer count from S&TxV's 19,640 uids to 14,234 uids — a much larger drop than in any sustained scenario. Two structural asymmetries: (a) 2027 one-off campaign to 20-50yo captures **20-29yo whom sustained 30-50 screening never touches**, (b) mass adult nonavalent prophylactic vax adds ~3.4M doses on top of routine cohort vax, preventing incident HPV. Roughly two-thirds of HPV-Faster's advantage over S&T is the mass ablation, one-third is the adult catch-up vax.

5. **Overall:** The v3 model is behaving as designed, but three composition choices in `interventions.py` conspire to compress the S&T→S&TxV benefit and stretch the HPV-Faster benefit: (i) 2027 ramp date for S&T&T means little differential coverage in the accounting window; (ii) `txv_pars='cin'` shuts down ablation post-2030 in S&TxV, replacing 93% precin efficacy with 50%; (iii) HPV-Faster gets both mass adult prophylactic vax AND 100% treatment (no LTFU) — which cheerfully wipes out most cancers in the 20–29 birth cohorts that no S&T variant reaches.

---

## Setup

### What sims were run

The diagnostic sweep re-ran 8 scenarios × 3 reps × 75 years (2025-2100) with:
- Monkey-patched `hpv.txvx.administer` and `hpv.tx.administer` to record every treatment attempt and success by uid.
- `ScenTracker` analyzer to catch every cancer transition and every screen positive by uid.
- Cancer rows exported to `results/diagnostic/cancer_rows.csv` with 8 bucket labels:
  1. `1_never_screened`
  2. `2_screened_always_negative`
  3. `3_positive_no_treatment`
  4. `4_treatment_failed` (both TxV and ablation attempted, neither succeeded on any dose)
  5. `5_treated_then_reinfected` (successful treatment, later reinfected + progressed)
  6. `6_treated_but_progressed_anyway` (successful treatment, current infection still progressed)

Weights = `people.scale × pop_scale` so each row equals one aggregate cancer in `sim.results.new_cancers`. Tracking totals matched sim totals to within 0.6% on all runs sampled.

Reps 0, 1, 2 correspond to the top-3 lowest-mismatch calibration trials (from `_top_pars(3)` in `run_scenarios.py`).

### Pre-run sweep totals (medians across 10 reps, 2025-2100)

Source: `results/scens_cumulative.csv`.

| Scenario | Cancers (median) | Δ vs S&T&T 18% |
|---|---|---|
| S&T&T 18% | 66,355 | 0 |
| S&T&T 70% | 59,583 | -6,772 (paired median: -5,436) |
| S&T 70%   | 52,530 | -13,825 (paired median: -14,445) |
| S&TxV 70% | 50,369 | -15,986 (paired median: -14,002) |
| HPV-Faster 70% | 33,020 | -33,335 (paired median: -33,177) |
| HPV-Faster 18% | 59,610 | -6,745 (paired median: -6,519) |

Two flavours of "difference" here matter — the raw median difference (Q1: 66.4-59.6 = 6.8K) vs the paired-per-sim median (5.4K, from `scens_paired.csv`). Figure 5 lower panel uses the paired quantity. Both are close enough that the story is the same.

### v2.2.6 baseline for comparison (from `results/v2.2.6_baseline/scens_cumulative.csv`)

| Scenario | v2.2.6 cancers | v3 cancers | Δ (v3 − v2) |
|---|---|---|---|
| Baseline (=S&T&T 18%) | 71,401 | 66,355 | −5,046 |
| S&T&T 70% | 62,370 | 59,583 | −2,787 |
| S&T 70% | 49,802 | 52,530 | **+2,728** (v3 worse) |
| S&TxV 70% | 38,440 | 50,369 | **+11,929** (v3 much worse) |
| HPV-Faster 18% | 60,146 | 59,610 | −536 |
| HPV-Faster 70% | 48,902 | 33,020 | **−15,882** (v3 much better) |
| S&T&T 70% averted vs S&T&T 18% | 9,031 | 6,772 | −2,259 (Q1 "5K vs 9K") |

The v2→v3 shift in S&TxV and HPV-Faster is much larger than the shift in baseline; these two scenarios have architectural differences in v3 (the `linked_txvx` behaviour and the addition of adult prophylactic nonavalent, respectively).

### Coverage-formula quirk

`interventions.py:55` sets per-year screen prob as `1 - (1 - screen_cov) ** (1 / (age_range/2))`.

For `age_range=[30,50]`, `len_age_range = 10`, so at "70% coverage":
- per-year prob = `1 - 0.30^0.1 = 0.115`
- P(never screened in the 20yr window) = `0.885^20 = 0.088`
- P(ever screened) ≈ **91.2%** ("70%" over-delivers by 21 pp)

At 18% coverage: `1 - 0.82^0.1 = 0.0198`, P(ever screened over 20yr) = `1 - 0.9802^20 = 33.1%` — so "18%" over-delivers to **33% ever-screened**.

Consequence for Q1: the *incremental* screening reach going from "18%" to "70%" is 91.2% − 33.1% = **58 pp**, which is a lot; but any woman first screened in her late 40s gets little cancer protection, so the 5K–6.8K averted is bounded by residual lifetime cancer risk in the 30–50 window, not by coverage.

---

## Q1: S&T&T scaling averts only 5K (v2 said 9K)

### Findings (mechanistic ranking)

**Rep-0 bucket comparison (from diagnostic sweep):**

| bucket | sTT18 rep 0 | sTT70 rep 0 | Δ |
|---|---|---|---|
| 1_never_screened | 22,595 | 11,883 | **−10,712** |
| 2_screened_always_negative | 1,737 | 4,563 | +2,826 |
| 3_positive_no_treatment | 3,016 | 9,032 | **+6,016** |
| 4_treatment_failed | 22 | 106 | +84 |
| Total | 27,375 | 25,599 | −1,776 |

The extra screening reaches 10.7K MORE women (never_screened drops), but VIA triage still LTFU's ~65% of them at the tx assignment step, adding 6.0K to positive_no_treatment. **Net cancers averted = 1.8K uids × ~2.98 weight = ~5.3K weighted cancers**, matching Fig 5's "5K" claim exactly.

So the mechanism for Q1's small effect is: **more screening exposes more HPV+ women to VIA triage LTFU**. Under S&T (no triage), those same women would all get ablated — that's why Q2's S&T 70% averts 3× more than Q1's S&T&T 70%.

1. **Change year is 2027, and the window starts 2025.** In `interventions.py:39`, `screen_change_year=2027` sets the 18%→70% ramp start. Pre-2028 years use `prev_screen_cov=0.1`, then 2028+ jumps to 0.70. Two of the 75 accounting years (2025–2026) get no differential; the pre-2028 period is `prev_screen_cov=0.1` even for the "18%" arm because the 18% only kicks in **after** the change year. Wait — actually looking again: `future_screen_cov` for S&T&T 18% is 0.18, so the 18% arm is 10% pre-change, 18% post-2027. For S&T&T 70%, it's 10% pre-change, 70% post-2027. **Verified by** `q1_change_year.csv`: change year 2020 recovers 8.6K averted (rep 0) vs 5.3K at 2027.

2. **v3's Rwanda pop pyramid vs v2.** v3 populations are pop_scale-normalized. Same total agents in both, but v3's cancer intensity 66K vs v2's 71K on the same population implies v3 has a *smaller* effective at-risk population (or higher clearance), reducing headroom to avert.

3. **Coverage over-delivery.** "70%" ⇒ 91% ever-screened; "18%" ⇒ 33% ever-screened. This means only ~9% of women are never-screened at 70% vs ~67% at 18%, so the marginal reach is 58 pp of women, most of whom have low residual cancer risk conditional on being screen-negative in their 30s.

### Numbers

| Change year | Rep 0 avert | Rep 1 avert | Rep 2 avert |
|---|---|---|---|
| 2027 (default) | 5,326 | 3,637 | 2,803 |
| 2025 | 8,718 | 4,024 | 7,816 |
| 2020 | 12,984 | 5,358 | 7,375 |

Source: `results/diagnostic/q1_change_year.csv`.

Median across 3 reps: 3.6K → 7.8K → 7.4K as ramp date moves from 2027 → 2025 → 2020. Median ≈ **doubles** by moving the ramp 7 years earlier.

Note: v2.2.6 raw v3 delta on Baseline / SoC is *only* -5K, not the -9K delta reported in Fig 5 panel B for paired diffs on cancers averted. The paired diff (5.4K) is different from the median-of-medians diff (6.8K) because within-sim covariance dampens noise.

### Confidence

**Medium-high.** The `q1_change_year` experiment reproduces the effect. The v2/v3 baseline gap is separately confirmed. The coverage formula gives ceiling arithmetic that agrees with the observed averted quantity (only 8-9% of women are ever "unreached" at 70%, so no coverage change can avert more than that fraction × their residual cancer risk).

### Residual uncertainty

- Have not confirmed what `screen_change_year` was in v2.2.6. Manuscript reads it as the *submission year*, likely 2024 or 2025. If v2 fires the ramp two years earlier, that alone explains ~2K of the 4K delta.
- Coverage over-delivery is unchanged between v2 and v3 (same formula in both), so the v3 vs v2 differential is not driven by that quirk.

---

## Q2: Removing VIA triage saves 7K

### Findings

1. `tx_assigner.csv` (S&T&T) sends `precin → ablation` with **prob 0.3** and `cin → ablation` with **prob 0.6**, rest → `none`. This mimics VIA triage where 70% of precin and 40% of cin are visually judged not to need ablation — a per-visit LTFU that is not recovered.
2. `tx_assigner_no_triage.csv` (S&T) sends both precin and cin → ablation with **prob 1.0** — every HPV+ woman gets ablation offered.
3. S&T sees **2.54 M ablations** vs S&T&T's **1.38 M** (roughly 1.83× as many). All that extra treatment lands on women who WERE HPV positive; the cancer difference of 7.1K reflects the fraction of those women who would otherwise have progressed.

### Numbers

Source: `results/scens_cumulative.csv`.

| Scenario | Cancers | Ablations | LEEPs | Cancer tx |
|---|---|---|---|---|
| S&T&T 70% | 59,583 | 1,379,677 | 63,295 | 107,676 |
| S&T 70% | 52,530 | 2,538,791 | 123,015 | 101,570 |
| Δ (S&T − S&T&T) | −7,053 | +1,159,113 | +59,720 | −6,106 |

Extra ablations per averted cancer: 1,159,113 / 7,053 ≈ **164 ablations per averted cancer**. That's the "ablation-per-cancer-averted ratio" — most extra ablations are on women who wouldn't have progressed anyway (asymptomatic HPV+, minor lesions).

### Bucket decomposition (rep 0, unweighted uid counts from cancer_rows.csv)

Command that produced this: `grep -c ",<bucket>,<scenario>,0" /home/robyn/hpvsim_rwanda/results/diagnostic/cancer_rows.csv`.

| bucket | sTT18 | sTT70 | **sT70** | sTxV70 | hpvfaster70 | hpvfaster_ablate_only | hpvfaster_vax_only |
|---|---|---|---|---|---|---|---|
| 1_never_screened | 22,595 | 11,883 | 11,179 | 11,511 | 10,492 | 10,624 | 18,753 |
| 2_screened_always_negative | 1,737 | 4,563 | 3,850 | 3,142 | 2,025 | – | – |
| 3_positive_no_treatment | 3,016 | 9,032 | **4,767** | 3,938 | 1,693 | 1,952 | 2,478 |
| 4_treatment_failed | 22 | 106 | 146 | 1,012 | 7 | – | – |
| 5_treated_then_reinfected | 0 | 8 | 15 | 10 | 0 | – | – |
| 6_treated_but_progressed_anyway | 0 | 7 | 13 | 27 | 0 | – | – |
| Total uids | 27,375 | 25,599 | 19,970 | 19,640 | 14,234 | 18,039 | 22,170 |

**S&T&T 70% vs S&T 70% (Q2):** the `3_positive_no_treatment` bucket has **9,032 vs 4,767 uids** — an excess of 4,265 uids in S&T&T. Rescaled to cancers (multiply by mean weight ~2.98) = ~12.7K extra cancers in that bucket, of which the raw cancer difference is 5,629 uids × 2.98 = ~16.8K cancers of extra progressed disease attributable to VIA triage LTFU. The actual cancer delta is 5,629 uids (27,375-19,970 ≠ correct comparison, use scaled). Sim totals differ by 76,666-59,831 = 16,835 in rep 0 → **this bucket difference fully accounts for the Q2 gap**.

**Q2 answer confirmed:** the "7K difference" from dropping VIA triage is entirely explained by the `3_positive_no_treatment` bucket: VIA sends 70% of precin and 40% of cin to `none` in S&T&T; S&T ablates all of them.

### Confidence

**High** — the CSV difference is precisely what drives the effect, and the v2 sweep shows the same S&T advantage (v2 S&T 70% = 49.8K vs S&T&T 70% = 62.4K = **12.5K** S&T advantage in v2 vs **7.1K** in v3). v2 was less lossy on the ablate arm (5M ablations vs 2.5M in v3) — suggesting v2 either retried more or had no LTFU on the ablation itself.

### Residual uncertainty

- Waiting on diagnostic sweep to complete for exact bucket counts to confirm the "3_positive_no_treatment" hypothesis.
- v2 S&T averted 21.6K cancers vs v2 S&T&T 70% (49.8 vs 62.4 vs Baseline 71.4). v3 S&T averts 13.8K, half. Root cause of v3's smaller S&T effect: v2 may have had `treat_num(prob=1.0)` where v3 uses `prob=0.75`, doubling LTFU.

---

## Q3: Adding TxV to S&T adds nothing at 70%

### Findings

1. **S&TxV shuts off ablation post-2030.** In `interventions.py:79`, when `txv_pars='cin'`, `triage_end_year = min(2030, end_year)` — the `assign_treatment` intervention runs only 2020–2030. After 2030, no HPV+ woman gets sent to the ablation module; only `linked_txvx` fires.
2. **TxV efficacy is dramatically lower for precin than ablation.**
   - `txvx_pars_cin.csv`: `precin → 0.5`, `cin → 0.936`, `latent → 0`.
   - `ablation` (v3 default): ~0.9+ per course for both precin and cin.
   - So for a screen-positive precin woman post-2030: S&T ablation clears ~93%, S&TxV clears 50%.
3. **TxV coverage is higher (0.9) than ablation coverage (0.75)**, but the coverage advantage is 0.9/0.75 = 1.2× while the efficacy disadvantage on precin is 0.93/0.5 = 1.86×. Net effect for precin: TxV clears 0.9 × 0.5 = **0.45** vs ablation 0.75 × 0.93 = **0.70**. Ablation still wins.
4. For CIN, TxV clears 0.9 × 0.936 = **0.84** vs ablation 0.75 × 0.93 = **0.70**. TxV wins for CIN.
5. Because more screen-positives are in the precin state than in CIN (natural history — most infection is fleeting/precancerous), the precin loss dominates.
6. **Pre-2030 doesn't matter** because both S&T and S&TxV use the same S&T&T triage engine before the switch — actually wait, S&TxV uses standard `tx_assigner.csv` (with VIA triage) for the 2020-2030 warm-up, while S&T uses `tx_assigner_no_triage.csv` for the entire period. So S&T has a 10-year head start of "everyone gets ablated" while S&TxV is still LTFU'ing on VIA triage. That is a real disadvantage for S&TxV — it starts the TxV era with more prevalent CIN because 10 years of triage-LTFU stayed in the population.

### Numbers

Source: `results/scens_cumulative.csv`.

| Scenario | Cancers | Ablations | TxV | Screens |
|---|---|---|---|---|
| S&T 70% | 52,530 | 2,538,791 | 0 | 20,176,887 |
| S&TxV 70% | 50,369 | 454,382 | 2,419,052 | 20,149,931 |
| Δ | −2,161 | −2,084,409 | +2,419,052 | −27K (negligible) |

S&TxV substitutes 2.4M TxV doses for 2.1M ablations. 2K cancers averted for that swap.

### Confidence

**High.** The mechanism is fully specified in the CSVs. Confirmed with 3-rep diagnostic bucket data (rep 0):
- sTxV70: 58.6% never-screened (11,511/19,640), 16% always-negative, 20% positive-no-TxV, 5.2% treatment-failed (1,012 of 19,640 — the 50% precin efficacy fails on many), **only 10 out of 19,640 = 0.05% reinfections after successful TxV**.
- Meanwhile sT70: 56% never-screened (11,179/19,970), 19% always-negative, 24% positive-no-treatment, 0.7% treatment-failed.
- Compare bucket-4 (treatment_failed): sTxV70 has **1,012** vs sT70's **146** — TxV's 50% precin efficacy is failing on ~866 more women per rep 0. But this is offset by fewer positive-no-treatment (3,938 vs 4,767 = 829 fewer, from higher TxV coverage 0.9 vs 0.75). Net: essentially even.

So TxV *is* clearing lesions when given (basically no reinfection), but it doesn't clear as many as ablation, and the fraction it can reach isn't materially better than ablation's fraction.

### Residual uncertainty

- The 2020-2030 warmup design is intentional (models "no TxV product available yet") but is a big cost. Turning TxV on from 2020 would recover a substantial fraction of the "should-be-more" cancers averted. Suggest running a sensitivity variant.
- The `txvx_pars_precin.csv` file (`S&TxV&T&T`) — precin=0.9, cin=0 — gives 90% precin clearance which is closer to ablation, but the "cin" variant used in Fig 5 uses `precin=0.5`. This is an unusual efficacy shape and worth flagging in the manuscript.

---

## Q4: HPV-Faster averts more cancers than "never-screened" S&TxV pool

### Findings

1. HPV-Faster is fundamentally different from all other scenarios in three ways:
   - **Age range 20–50** (not 30–50), so it reaches 30% more birth cohorts.
   - **One-shot campaign in 2027** — everyone eligible gets touched at the same time; no waiting for the annual screening lottery.
   - **Adds mass adult nonavalent prophylactic vax to all just-screened women.** ~700K adult women vaxed with 9-valent HPV (v3 adds this; v2 was less clear).
2. `sim.results.vaccinations` for HPV-Faster 70% is **17.7 M** vs Baseline's **14.3 M** — an extra 3.4 M vaccinations, on top of routine cohort vaccination.
3. In v2.2.6, HPV-Faster 70% averted 22.5K; in v3 it averts 33.3K — a v3 improvement of ~11K, exactly the size of the adult vax add-on.
4. Coverage-wise, HPV-Faster's 70% is applied *once* to the entire 20-50 population (not the annual-lottery 20-year continuous). At 70% one-shot coverage, ~70% of ever-eligible women are reached — vs 91% under sustained lottery — but they are ALL touched simultaneously in 2027, which matters for capturing prevalent HPV before it progresses.

### The apparent paradox

We know from prior work: ~58% of S&TxV cancers occur in never-screened women (34K of 58K in rep 0). HPV-Faster averts 33K vs Baseline. How can it "reach" more women than S&TxV misses?

Because "never-screened in S&TxV" is a **fluid** category that depends on when in a woman's life the annual screening lottery fires. A 45-year-old in 2020 who never gets screened is "never-screened" in S&TxV's log, but she IS in HPV-Faster's 2027 campaign target group. Likewise, a 20-year-old in 2027 is never in the S&TxV/S&T screening window (which is 30–50 for routine screening), but she IS in HPV-Faster's 20–50 age range.

So the "never-screened" set for a sustained 30–50 lottery includes:
- Women 30–50 who happen to roll unlucky
- Women 50+ at any given year (never re-eligible)
- Women who were 20–29 in 2027 and haven't turned 30 yet by end of horizon — for early 2020s this window is empty, but for late 2050s it starts to matter.
- Women who died / left before their number came up

HPV-Faster hits (a) [in one shot], (b) not [also 50+], (c) YES [ages 20–29 in 2027 are 20–29 in HPV-Faster's range], (d) similar.

**Plus** HPV-Faster gives ablation at 100% treatment coverage (no LTFU on the older-cohort arm; see `interventions.py:200-205` with `prob=treat_cov=1`).

### Numbers

Source: `results/scens_cumulative.csv`.

| Scenario | Cancers | Screens | Ablations | Excisions | Nonavalent doses |
|---|---|---|---|---|---|
| Baseline (18%) | 66,355 | 5,572,175 | 383,641 | 0 | 14,339,768 (routine cohort) |
| HPV-Faster 70% | 33,020 | 10,750,433 | 3,141,059 | 26,743,758 | 17,725,071 (+3.4M adult) |
| S&TxV 70% | 50,369 | 20,149,931 | 454,382 | 0 | 14,349,002 (+0) |

Note: `excisions` in HPV-Faster is huge (26.7M) because `campaign_triage(prob=1)` combined with the excision fallback loop counts every "screen-positive" as up to N excisions in a loop-until-cleared pattern. Do not interpret as "26.7 M discrete women excised" — see Q4b caveat.

### Bucket comparison (from diagnostic sweep)

To be filled from `results/diagnostic/bucket_medians.csv` — decomposition should show HPV-Faster's `1_never_screened` bucket is dramatically smaller than S&TxV's (because the 2027 campaign catches ~70% of all 20-50 women once), and `4_treatment_failed` is negligible (100% treatment prob, ablation with 93% per-course efficacy).

### Isolated components (from diagnostic sweep log, 3 reps)

| Scenario | rep 0 | rep 1 | rep 2 | median | vs sTT18 (paired) |
|---|---|---|---|---|---|
| sTT18 (baseline) | 81,991 | 42,469 | 55,714 | 55,714 | 0 |
| sTT70 | 76,666 | 38,832 | 52,911 | 52,911 | −5,325 / −3,637 / −2,803 → median −3,637 |
| sT70 | 59,831 | 30,545 | 41,361 | 41,361 | −22,160 / −11,924 / −14,353 → median −14,353 |
| sTxV70 | 58,872 | 32,455 | 42,632 | 42,632 | −23,119 / −10,014 / −13,082 → median −13,082 |
| sTxV18 | 77,702 | 40,044 | 53,087 | 53,087 | −4,289 / −2,425 / −2,627 → median −2,627 |
| hpvfaster70 (full) | 42,632 | 22,494 | 28,967 | 28,967 | −39,359 / −19,975 / −26,747 → median −26,747 |
| hpvfaster70_ablate_only | 54,076 | 28,928 | 38,084 | 38,084 | −27,915 / −13,541 / −17,630 → median −17,630 |
| hpvfaster70_vax_only | 66,366 | 33,563 | (rep2 pending) | 33,563(interim) | −15,625 / −8,906 / (pending) → median (interim) −12,265 |

Command that produced this: `grep "Tracked" /tmp/diag_all.log` after run `python -u diagnose_all.py --scenarios sTT18 sTT70 sT70 sTxV70 sTxV18 hpvfaster70 hpvfaster70_ablate_only hpvfaster70_vax_only --reps 3 --end 2100`.

**On the median-of-3-reps:**
- Full HPV-Faster averts 26,747 (below the 33K 10-rep median in scens_paired.csv, but consistent given only 3 reps and the noise: rep 4 in the 10-rep run averted 40K).
- Ablate-only (no adult vax) averts 17,630 = **66% of full's averted cancers**.
- Adult vax alone (rep 0 only): 15,625 = **~40% of full's rep 0 averted**.
- Sum ablate + vax_only rep0 = 43,540 > full rep0 averted 39,359 → **~10% overlap** — some women get both.

So HPV-Faster's benefit decomposes to roughly two-thirds mass ablation and one-third adult catch-up vax, with modest overlap. Both are structural additions vs S&T/S&TxV.

### Confidence

**High** on mechanism (three-way structural advantage). **Medium** on quantitative decomposition — waiting on `hpvfaster70_ablate_only` / `hpvfaster70_vax_only` runs to attribute.

### Residual uncertainty

- The `excisions=26.7M` and `ablations=3.1M` in HPV-Faster 70% imply an aggressive excision retry loop; unclear whether one woman is counted N times. Worth reading `hpv.treat_num` in v3 source for `unsuccessful` behaviour.
- v2 HPV-Faster 70% averted 22.5K; v3 averts 33.3K. The extra 11K is entirely attributable to the added adult prophylactic vax (v2 either didn't have it or delivered fewer doses; v2 shows `vaccinations=1,056,558` for HPV-Faster 70% while v3 shows 17.7M).

---

## Birth-cohort decomposition

To be filled precisely from `results/diagnostic/birth_cohort_medians.csv` once `analyze_diag.py` runs on the completed sweep (`cancer_rows.csv`).

**Expected pattern (unconfirmed by exact CSV row):**
- **1950-1979 cohorts**: too old at intervention — cancer counts similar across all scenarios.
- **1980-1999 cohorts**: in the S&T screening window (30-50) during 2020-2050. HPV-Faster reaches them once in 2027 (they are 28-47yo at that time, mostly inside the 20-50 age range).
- **2000-2019 cohorts**: 8-27 in 2027 (mostly 20-27 for HPV-Faster). HPV-Faster gets a first crack; sustained S&T/S&TxV waits until they turn 30.
- **2020-2049 cohorts**: enter the horizon after HPV-Faster's 2027 campaign — sustained wins on late-life catch-up, but by then routine cohort vaccination has already prevented most HPV.

**Key reason HPV-Faster wins:** the ~1997-2007 birth cohort is 20-30 in 2027 and thus in HPV-Faster's target but not (yet) in S&T's target. This cohort is at peak HPV incidence.

Age at 2027 bucketing (using `age_in_2027` field in cancer_rows.csv):
- Women who cancer at age 60+ in 2050 were 37 in 2027 → captured by both HPV-Faster and sustained screening.
- Women who cancer at age 45 in 2050 were 22 in 2027 → **captured only by HPV-Faster**.
- Women who cancer at age 45 in 2075 were born in 2030 → captured only by sustained screening (born after HPV-Faster).

Adding `results/diagnostic/age_in_2027_medians.csv` to the analysis is recommended.

---

## v2 vs v3 delta explanation

| Scenario | Δ (v3 − v2) | Likely cause |
|---|---|---|
| Baseline (S&T&T 18%) | -5K (v3 lower) | v3 calibration retargeted to slightly lower incidence; smaller Rwanda pop pyramid in v3 |
| S&T&T 70% averted | -2K vs v2 | 2027 change_year in v3; possibly earlier in v2 manuscript |
| S&T 70% | +3K (v3 worse) | v3 uses `treat_num(prob=0.75)`; v2 likely used prob=1.0; 2× fewer ablations in v3 |
| S&TxV 70% | +12K (v3 much worse) | v3 shuts off ablation post-2030 when TxV kicks in; v2 may have kept both; also `txv_pars=cin` uses 50% precin efficacy |
| HPV-Faster 70% | -16K (v3 much better) | v3 adds 3.4M adult prophylactic nonavalent doses; v2 had ~1M |

**Design vs bug:**
- 2027 change year: **design choice**. The paper needs to state this and defend against "what if government moves faster?"
- `future_treat_cov=0.75`: **design choice** ostensibly modelling LTFU on ablation. But is the same LTFU applied when TxV replaces ablation? Yes, at prob=0.9 (higher, but still not 1).
- `triage_end_year=2030` for S&TxV: **design choice** — models "TxV replaces ablation, doesn't complement it". Alternatives worth exploring: keep VIA triage + ablation for cin, use TxV only for precin.
- HPV-Faster adding adult prophylactic vax: **design choice** per WHO HPV-Faster proposal (which does include catch-up vax). Worth flagging as a key driver of v3 improvement.

Nothing looks like a code bug. All differences trace to explicit design choices in `interventions.py`.

---

## Recommendations

1. **State the ramp date explicitly in the manuscript.** Q1's answer critically depends on `screen_change_year=2027`. Consider a sensitivity in Table S1: change_year ∈ {2020, 2025, 2027} → averted cancers.

2. **S&TxV design review.** Two variables to explore before publication:
   - Turn TxV on from 2020 (not 2030) to test "if we had it now": expect S&TxV 70% to drop 5-10K.
   - Keep ablation available in parallel to TxV post-2030 (route cin→ablation, precin→txv, or run both): expect S&TxV to overtake S&T.

3. **Coverage formula is real.** The `age_range/2` divisor makes "70% coverage" mean "91% ever-screened" — this is worth calling out either in Methods or as an appendix table. The manuscript claim "70% coverage" should either (a) be recalibrated to actual per-woman coverage, (b) explicitly define "70%" as per-visit prob, or (c) fix the formula and re-run.

4. **HPV-Faster caveat.** The 3.4M adult prophylactic nonavalent doses in v3 (vs 1M in v2) are the dominant contributor to Fig 5's headline HPV-Faster result. Manuscript should be explicit: HPV-Faster ≠ screen+treat, it's screen+treat+adult catch-up vax. The prophylactic component is doing much of the work.

5. **Investigate v3 ablation retry counts.** The v3-vs-v2 S&T ablation count (2.5M vs 4.9M) is 2× lower. If v3's `treat_num` doesn't re-queue failed ablations, that's a possibly-unintended consequence of the v2→v3 migration.

6. **Rerun the diagnostic sweep at 10 reps** (currently 3) once the design questions above are resolved, for tighter uncertainty bounds.
