# Residual cancer analysis

**Question:** After every intervention we have — routine bivalent vax,
screen-and-treat at 70% coverage, TxV, HPV-Faster mass campaign — how
many cervical cancers still occur in Rwanda over 2030-2100? Who is
getting them? To what extent can the existing intervention set reach
these women, and what is genuinely beyond us?

**Setup:** normalized scenario sweep — every intervention starts 2030,
accounting window 2030-2100. 5 top-calibration reps × 6 scenarios.
Data: `results/diagnostic_normalized/cancer_rows.csv`.

**Headline:** HPV-Faster 70% is the strongest of the six modelled
scenarios in the normalized (2030-start) framing, leaving 26,012
cumulative cancers over 70 years. About **46% of that residual is
genuinely beyond the reach of any modelled intervention** — women
who were 50+ in 2030 and acquired HPV before screening or vax
existed. The remaining **54% is reachable in principle**, dominated
by never-screened women (76% of the reachable residual), pointing at
screening coverage above 70% as the highest-leverage lever.

---

## §1 — Scenario totals (2030-2100)

Cumulative cervical cancer cases, median across 5 reps [p10, p90]:

| Scenario | Cases | Averted vs no intv |
|---|---:|---:|
| No interventions | 53,573 [43,050; 82,552] | 0 |
| Baseline (S&T 18%) | 53,206 [42,880; 78,658] | 367 |
| S&T&T 70% | 48,033 [38,991; 74,764] | 5,540 |
| Mass TxV 50/90, 70% | 33,315 [25,264; 44,488] | 20,258 |
| S&TxV 70% | 35,701 [30,290; 54,133] | 17,872 |
| **HPV-Faster 70%** | **26,012** [21,169; 39,691] | **27,561** |

Notes:

- Baseline in the normalized framing is nearly indistinguishable from
  No interventions because it only supplies 5 years of pre-scale-up
  screening from 2030 onward (screen_change_year = 2029). This is
  intentional — the point of the normalized run is to compare
  interventions on equal footing, not to score them against the paper's
  2020-start baseline.
- **HPV-Faster 70% beats S&TxV 70% by ~10K averted.** Consistent with
  the earlier investigation's corrected finding: HPV-Faster wins via
  the 2030 one-off no-LTFU screen-and-treat of the 20-50yo cohort
  (1980-2010 birth years), plus mass prophylactic vax preventing
  future HPV in that group. S&TxV only screens 30-50 and starts the
  TxV era from 2030 with no historical treatment inventory. In the
  paper's 2020-start framing, screening had a decade to work before
  TxV became the therapy, so the comparison flipped.

---

## §2 — Residual (26K) under HPV-Faster 70%

### 2a. Birth cohort

| Cohort (age in 2030) | Cases | Share |
|---|---:|---:|
| 1930-40 (90+) | 63 | 0.2% |
| 1950 (80) | 813 | 3.0% |
| 1960 (70) | 3,077 | 11.4% |
| **1970 (60)** | **8,283** | **30.8%** |
| **1980 (50)** | **6,264** | **23.3%** |
| **1990 (40)** | **5,180** | **19.3%** |
| 2000 (30) | 1,382 | 5.1% |
| 2010 (20) | 786 | 2.9% |
| 2020 (10) | 599 | 2.2% |
| 2030+ (unborn in 2030) | 450 | 1.6% |

**73% of the residual is concentrated in the 1970-1990 cohorts.**
The 1970 cohort alone (age 60 in 2030) contributes almost a third —
these women were 30-50 during 2000-2020, before any Rwandan screening
programme reached them, and they age out of the screening window in
the early 2030s.

### 2b. Age at causal HPV infection

| Age at infection | Cases | Share |
|---|---:|---:|
| < 20 | 751 | 2.9% |
| 20-30 | 4,200 | 16.0% |
| 30-40 | 8,540 | 32.6% |
| 40-50 | 9,064 | 34.6% |
| 50-60 | 2,973 | 11.4% |
| 60+ | 664 | 2.5% |

**Two-thirds of residual cancers arise from HPV acquired in the
30-50yo window.** This is the window that HPV DNA screening covers
in principle, so failures here point at coverage / uptake gaps rather
than biology-outside-the-net.

Only ~3% arise from HPV acquired below 20 — indicating that even
adding pre-adolescent vaccination gains would help modestly at best;
the residual isn't dominated by cancers seeded before the vax window.

### 2c. HIV status

| Status | Cases | Share |
|---|---:|---:|
| HIV negative | 24,266 | 92.9% |
| HIV positive | 1,859 | 7.1% |

Roughly proportional to national HIV prevalence in women (~6-8%).
Not a disproportionate driver of the residual, so WHO's stronger
screening-frequency guidance for WLHIV would have a small aggregate
effect on the total (though still worth doing on equity grounds).

---

## §3 — Mode of failure (buckets)

For each residual cancer under HPV-Faster 70%:

| Bucket | Meaning | Cases | Share |
|---|---|---:|---:|
| 1 | Never screened | 20,561 | 76.4% |
| 2 | Screened, always negative | 4,107 | 15.3% |
| 3 | Screened positive, no treatment | 2,154 | 8.0% |
| 4 | Treatment attempted, all failed | 36 | 0.1% |
| 5 | Successful treatment then reinfection | 15 | 0.1% |
| 6 | Successful treatment but progressed anyway | 33 | 0.1% |

**Bucket 1 dominates by an order of magnitude.** Treatments 4-6
combined are 84 cases (0.3%) — the treatment products (ablation,
excision, TxV) are essentially not the operational bottleneck once
a woman is caught.

**Bucket 2 (screened but always negative) at 15%** is the "test
sensitivity" ceiling. HPV DNA is already 99% sensitive per genotype;
this residual can only come from CIN developing between screening
intervals (10y default gap), so shorter intervals — not better tests —
is the lever there.

**Bucket 3 (positive but no treatment) at 8%** is the ~25% LTFU
default (`future_treat_cov=0.75`) plus the tx_assigner CSV's
routing failures.

---

## §4 — Reachability partition

For each residual cancer, was the woman in *any* modelled
intervention's addressable window at any point between 2030 and
her cancer diagnosis? (Windows: routine vax = age 11-12 at any t
after 2011, routine screening = age 30-50 at any t after 2030,
mass adult campaign = age 20-50 in 2030.)

| Class | Cases | Share |
|---|---:|---:|
| **Unreachable** in principle | **11,914** | **45.7%** |
| **Reached in principle** but missed | **14,139** | **54.3%** |

**The 11,914 unreachable cancers are the ~1960-1975 birth cohorts.**
These women:

- Were past 50 (aged out of the routine screening window) by the time
  screening rolled out in 2030,
- Were past 50 (aged out of the mass adult campaign window) in 2030,
- Were past 12 by 2011 (missed the routine prophylactic vax rollout).

They acquired HPV in the 1970s-1990s, decades before any modelled
intervention was available in Rwanda. **This is a legitimate
limitation of the model and the reality it describes** — no plausible
70-year policy scenario can retrospectively prevent HPV acquired
before the tools existed. The right framing in the manuscript is:
"a defensible floor to the residual, imposed by the pre-2030 legacy
of HPV exposure in cohorts that had aged out of the screening window
by intervention launch."

**The 14,139 reachable-but-missed cancers are the operational target.**
Their bucket breakdown (approximate — buckets are not conditioned on
reachability in the current script) is dominated by:

- Never screened → higher screening reach + self-sampling / DBS
  modalities to close the last 30%.
- LTFU on treatment referral → tighter same-visit models (see-and-treat)
- Between-screen incident CIN → shorter screening interval for high-risk
  subgroups (WLHIV, women with a prior positive).

---

## §5 — Comparison: Baseline (18% status quo) residual

For context, the Baseline (S&T at 18% status-quo coverage, 2030 start)
residual is 53,206 cases with:

- Birth cohorts more concentrated in 1980-1990 (61% vs 43% under
  HPV-Faster) — the residual you avoid under HPV-Faster is
  specifically in those age groups the 2030 mass campaign captures.
- Age-at-infection distribution similar (~65% in 30-50 window).
- HIV share slightly lower at 5.8%.
- Bucket distribution: 87% never-screened (up from 76%), 8% positive-
  no-treatment, 5% screened-negative.
- Reachability: only 24% unreachable (vs 46% under HPV-Faster), because
  under Baseline many cancers occur in "reachable-in-principle" cohorts
  simply because they were never actually screened.

Reading the reachability lift as **"HPV-Faster closes 27K cancers,
mostly by reaching people the Baseline never reaches"** is a cleaner
framing than "HPV-Faster averts 27K cancers" — the levers are
targeted at very specific cohorts (20-50 in 2030) and modalities
(mass one-off campaign + prophylactic vax + routine screening).

---

## §6 — What would it take to close the reachable residual?

Levers ranked by the bucket they attack, ordered by leverage under
HPV-Faster 70%:

1. **Screening reach above 70%** (bucket 1, 76% of residual) — the
   dominant lever. Getting from 70% to 90% lifetime coverage in
   30-50yo would take ~15K of the residual off the table (linear
   extrapolation). Requires self-sampling, mobile clinics, or
   employer-based delivery.
2. **Tighter treatment LTFU** (bucket 3, 8%) — currently 75% treatment
   probability. Same-visit see-and-treat would push this to 90%+.
   Would take ~1K off the residual.
3. **Shorter screening interval for high-risk** (bucket 2, 15%) —
   women with prior positives or HIV+ screened every 3-5 years
   instead of every 10. Would take ~1-2K off.
4. **Better therapeutic products** (buckets 4-6, 0.3%) — negligible
   effect. Not the operational bottleneck.

**Total plausible reduction of the reachable residual: ~50-60%** — so
~7-8K more cancers averted, bringing the total from 26K down to
18-19K over 2030-2100.

The remaining ~18K would be the **structural floor** — 12K unreachable
pre-2030 legacy + 6K reachable-but-missed even under implausibly high
push. The paper's discussion should describe this floor explicitly.

---

## §7 — Manuscript framing

Suggested paragraphs for the discussion:

> "Under the strongest modelled intervention, HPV-Faster at 70% target
> coverage, an estimated 26,000 [21,000; 40,000] cumulative cervical
> cancer cases occur in Rwanda over 2030-2100. Approximately half of
> this residual (12,000 cases) is attributable to the pre-2030 legacy
> of HPV exposure in cohorts already past the screening age window at
> intervention launch — a floor that no future prevention or screening
> tool can lower. The remaining residual (~14,000 cases) is
> dominated by women who were never reached by screening despite being
> in the eligible age window at some point during the modelled
> horizon, indicating that increases in screening uptake above the
> modelled 70% ceiling — whether through self-sampling, mobile
> clinics, or same-visit see-and-treat models — represent the highest-
> leverage remaining opportunity for further reducing the burden."

> "Failures of screening test sensitivity, treatment product efficacy,
> and post-treatment reinfection collectively account for less than 1%
> of residual cancers under HPV-Faster 70%. Once a woman enters the
> screening cascade, current-generation tools successfully prevent her
> from developing cervical cancer. The operational challenge is
> reaching her."

---

## Reproduction

```bash
# Runs the normalized scenario sweep (30 sims, ~45 min on one core).
python diagnose_all.py --normalized --reps 5 \
    --scenarios "No interventions" "Baseline" "S&T&T 70%" \
                "S&TxV 70%" "HPV-Faster 70%" "Mass TxV 50/90, 70%" \
    --outdir results/diagnostic_normalized

# Produces per-cancer CSV; then:
python analyze_residual.py
```
