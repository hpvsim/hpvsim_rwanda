# Residual cancer analysis

**Question (paraphrased from the user):** After every intervention we have —
routine bivalent vax, screen-and-treat at 70% coverage, TxV, HPV-Faster
mass campaign — roughly 50-70K cumulative cervical cancers still occur
in Rwanda over 2030-2100. Who is getting them? To what extent can the
existing intervention set reach these women, and what is genuinely
beyond us?

The goal is to either (a) name what additional interventions would push
Rwanda close to elimination or (b) account credibly for why we cannot.

**Data:** `results/diagnostic_normalized/cancer_rows.csv` — every cancer
diagnosed 2030-2100 across the normalized scenario sweep (5 reps ×
{No interventions, Baseline, S&T&T 70%, S&TxV 70%, HPV-Faster 70%,
Mass TxV 50/90 70%}), tagged with:

- birth_year, age at cancer, age at causal infection, age at first CIN
- HIV status, genotype
- screening reach (n_screens, n_positive), TxV attempts & successes,
  ablation attempts & successes
- mass adult vax reach, routine vax reach
- bucket (1_never_screened / 2_screened_always_negative /
  3_positive_no_treatment / 4_treatment_failed /
  5_treated_then_reinfected / 6_treated_but_progressed_anyway)

Accounting window is 2030-2100. All interventions in the sweep begin in
2030 so scenarios are directly comparable.

---

## Structure

### 1. The residual for each scenario

For each scenario, cumulative cancer count over 2030-2100 (median +
[10, 90] across the 5 reps):

| Scenario | Cumulative cancers | Averted vs No intv |
|---|---|---|
| No interventions | _fill_ | 0 |
| Baseline (S&T 18%) | _fill_ | _fill_ |
| S&T&T 70% | _fill_ | _fill_ |
| S&TxV 70% | _fill_ | _fill_ |
| Mass TxV 50/90, 70% | _fill_ | _fill_ |
| HPV-Faster 70% | _fill_ | _fill_ |

The "residual" for each intervention is what's left, and that's the
denominator for the reachability analysis below.

### 2. Who is in the residual?

For the best-performing scenario (candidates: HPV-Faster 70% or the
"stack them all" combo — decide from §1), decompose the residual by:

#### 2a. Birth cohort

10-year birth-cohort bins. For each: cancer count, share of residual,
age at cancer (median), and whether the woman was reachable by any
2030-onward intervention at any point in her life.

Expected pattern: older cohorts (born pre-1980) are past-30 in 2030,
so screening is available but they carry the full history of untreated
HPV. Younger cohorts (born 2000+) get routine prophylactic vax at 11-12,
and are late to enter the screening age window — but their cancers
should be very few.

#### 2b. Age at causal event

For each cancer:
- Age at HPV acquisition (age_causal)
- Age at first CIN (age_cin)
- Age at cancer diagnosis (age_cancer)

Cross-tab against birth cohort. Identifies the biological windows the
intervention set fails to catch:

- Acquired HPV before 25 → prophylactic vax could reach if born
  ≥ 2000 (routine); mass adult vax could reach 20-50yo in 2030.
- Acquired HPV 30-45 → routine screening window covers this in
  principle; failure modes are the buckets in §3.
- Acquired HPV after 50 → out of the screening window entirely.

#### 2c. HIV status

Share of residual that is HIV-positive; compare against HIV prevalence
in the modelled population. If disproportionate, that's a specific
policy lever (higher screening frequency for WLHIV, per WHO 2021).

### 3. Bucket decomposition of the residual

For the best-performing scenario, split the residual into the 6
buckets. This is the operational bottleneck decomposition:

| Bucket | Meaning | Count | Share | Reachable how? |
|---|---|---|---|---|
| 1 | Never screened | Woman never touched by S&T | _fill_ | _fill_ | ↑ screening reach / new modality (self-sampling) |
| 2 | Screened, always negative | Test missed her infection | _fill_ | _fill_ | ↑ test sensitivity (nothing to do — HPV DNA is already 99%) |
| 3 | Positive, no treatment | Screened+, but no treatment attempt | _fill_ | _fill_ | ↑ LTFU recovery, ↓ tx_assigner LTFU |
| 4 | Treatment failed | Treatment attempted, all failed | _fill_ | _fill_ | Better treatment product |
| 5 | Treated then reinfected | Successful treatment, later reinfection → cancer | _fill_ | _fill_ | Prophylactic vaccine post-treatment |
| 6 | Treated but progressed anyway | Cancer was already inevitable at treatment time | _fill_ | _fill_ | Earlier detection |

### 4. The reachable pool

For each residual cancer, ask: at any point in this woman's life
between 2030 and her cancer diagnosis, could *any* of the modelled
interventions have plausibly reached her? Combine:

- Was she in the routine-vax age window (11-12) at any point 2030+?
- Was she in the screening age window (30-50) at any point 2030+?
- Was she in the mass-adult-vax / HPV-Faster window (20-50) in 2030?
- Was she in the TxV eligibility window (screened positive 2030+)?

"Unreachable in principle" = **none** of the above. That is:
- Born pre-1980 (already 50+ in 2030) → aged out of screening
- OR born post-2020 (only 10 in 2030, past the mass campaign) but
  still gets routine vax at 11-12 → shouldn't be in the residual
  much; if she is, prophylactic vax failed for her genotype.

Expected finding: most of the "unreachable" residual is older women
(born 1960-1980) whose HPV was acquired 1985-2010, well before any
screening was available. This is a legacy of Rwanda not having had
screening at population scale until the 2020s.

### 5. What would it take to eliminate the reachable residual?

For the residual cancers we *could* reach, quantify what push is
needed:

- If bucket 1 dominates: what screening coverage (%) closes bucket 1?
- If bucket 3 dominates: what treatment coverage (currently 75%) closes
  bucket 3?
- If bucket 5 dominates: catch-up prophylactic vax for treated women.
- If bucket 6 dominates: shorten the screening interval (currently 10y).

Concretely: rerun 1-2 "eliminative" scenarios with the levers turned
up as far as they can plausibly go — e.g., 90% screening coverage,
90% treatment coverage, no LTFU. See what the residual looks like then.

### 6. The unreachable residual → the paper's caveat

The remainder — cancers in women aged 50+ in 2030 who acquired HPV
before any intervention was available — is the paper's legitimate
limitation. Frame in the discussion as:

> "By 2100, roughly N cases remain in each modelled scenario. These are
> concentrated in birth cohorts pre-1980, whose lifetime HPV exposure
> predates the availability of screening or vaccination in Rwanda.
> Elimination of this legacy burden is beyond the reach of any
> currently modelled tool set."

If N is small (< 20K over 70 years for a 13M-population country), this
is a defensible caveat. If N is large (>50K), we need a stronger story.

---

## Analysis script

To be added: `analyze_residual.py`, reads
`results/diagnostic_normalized/cancer_rows.csv` and emits:

- `residual_by_scenario.csv` (§1)
- `residual_by_birth_cohort.csv` (§2a) — per (scenario, cohort)
- `residual_by_age_bucket.csv` (§2b) — per (scenario, age_causal_band)
- `residual_by_hiv.csv` (§2c)
- `residual_by_bucket.csv` (§3)
- `residual_reachability.csv` (§4)

All medians across reps; the paired-diff quantiles are not meaningful
here because we're conditioning on cancer having occurred, not
comparing paired arms.
