# Response to reviewers, npj Vaccines resubmission

*Manuscript: "Strategies to accelerate cervical cancer elimination: a modeling study in Rwanda"*

## To the editor

Thank you for allowing us to submit this revised manuscript. We are grateful for the two thoughtful reviewer reports; between them, they identified several substantive areas where the manuscript could be strengthened. In response, we have refit our calibration, corrected an inconsistency between the model and its written description, added five new sensitivity analyses, produced a threshold cost-effectiveness calculation, added supplementary tables of resource requirements, and updated our discussion and recommendations. A full point-by-point response to each numbered comment follows.

Text quoted verbatim from the reviewer reports is shown in italics; our responses and pointers to the revised manuscript follow each.

## Reviewer 1

### R1.1: Table of parameter values and sources

> *For modeling studies, it is helpful to the reader when the authors provide a table of parameter values and their sources. In the methods section, the authors provide the source of parameter values (Globocan 2020), and a list the parameters for the model. They describe the process of selecting the 50 best parameter sets, but they do not list the parameter values used in their simulations. It would help the reader if the authors would provide a table listing the parameters and the average and range of parameter values in their 50 parameter sets, along with the specific source.*

**Response:** We thank the reviewer for this suggestion, which we have adopted. We have added Table S1 to the supplementary materials, which lists all 17 calibrated parameters together with (i) a brief description of what each parameter represents, (ii) the source of the prior range (HPVsim defaults, prior HPVsim natural-history and network calibrations reported in Stuart et al. 2024, or literature-informed for the HIV-related parameters), (iii) the prior range (low, high, guess), and (iv) posterior summaries across the 50 best-fitting parameter sets (mean, median, 95% credible interval, min, max).

Non-calibrated model inputs (Rwanda demography from UN World Population Prospects, HIV incidence and ART coverage from UNAIDS, the Globocan 2020 cancer targets, and the HPV natural-history structure inherited from HPVsim) are cited inline in the Methods as before. We have added a pointer to Table S1 in the calibration paragraph of the Methods.

---

### R1.2: Sensitivity analyses

> *Modeling studies normally use sensitivity analyses to determine the sensitivity of the model results to variation in parameter values. It would help if the investigators would provide sensitivity analyses for time-to-elimination and for CC cases prevented for the key parameters in their agent-based model. This will give the reader more confidence in the study results.*

**Response:** We agree that sensitivity analysis is a core part of any modeling study. Our approach treats parameter uncertainty as a first-class output rather than a secondary check: every scenario is run once per each of the 50 best-fitting parameter sets from calibration (Table S1), and every result reported in the main text, including cumulative cancers, cancers averted, and elimination years, carries a 10-90% interval computed across those 50 runs. The width of each interval directly quantifies how each result would change as the underlying model parameters vary across the range compatible with the observed data.

Beyond this parameter uncertainty, we explore a set of targeted structural sensitivities around specific points of uncertainty that the reviewers highlighted. The most consequential of these are the assumed characteristics of the therapeutic vaccine, which we address both by comparing two very different vaccine profiles side-by-side (virus-clearing and lesion-regressing, main text §Results) and by varying the year the product becomes available across {2030, 2035, 2040, 2045, 2050} for every therapeutic-vaccine-carrying scenario (Figure S3; see also response to R1.3). We also test the sensitivity of the screen-treat impact to the assumption that women accept same-visit ablation without a triage confirmation of disease (main text §Results; see also response to R2.3).

---

### R1.3: Therapeutic vaccine availability date

> *The authors indicate that therapeutic HPV vaccines do not exist currently, but their model assumes that in 2030, they will exist and be deployable. Having effective therapeutic vaccines developed, evaluated, licensed, and deployed four years from now seems overly ambitious. It will help the reader if some sensitivity analysis could be conducted by varying the time to widespread availability and adoption of therapeutic vaccines, perhaps varying the therapeutic vaccine start year from 2030 to 2050 or so.*

**Response:** We agree that 2030 is an ambitious assumption and that a sensitivity analysis over later intro years is informative. We have added such an analysis, sweeping the therapeutic vaccine introduction year across {2030, 2035, 2040, 2045, 2050} for all four TxV-carrying strategies at 70% coverage. Figure S3 in the supplementary materials shows the resulting cumulative cancers averted 2030-2100 relative to status quo.

The finding is that all four strategies remain net-beneficial even under the most delayed intro year (2050), averting between 17,000 (Mass TxV 50/90) and 26,000 (S&TxV) cases relative to status quo. However, the sensitivity of impact to intro year varies substantially between strategies. The one-off Mass TxV 50/90 campaign loses the most benefit with delay (58% reduction in cases averted between the 2030 and 2050 intro years), because older cohorts age past the treatable window before the campaign takes place. In contrast, the linked screening-plus-TxV strategies retain most of their benefit (25 to 30% reduction across the same delay window), because screening continues to catch cases in the years before the TxV product becomes available.

This reinforces the reviewer's underlying point: policymakers should not count on a specific 2030 availability date, but the qualitative case for TxV as a complementary tool holds up across the plausible delay range.

---

### R1.4: Recommendation is incomplete

> *The recommendation of the authors ("investing in improved screening and treatment strategies to complement existing prevention efforts and achieve cervical cancer elimination decades sooner") is incomplete, since optimal strategies are dependent on availability of therapeutic vaccines. I think that the authors should include investing in therapeutic vaccine development into their recommendation. Such an investment seems necessary based on the data the authors present.*

**Response:** We agree, and have broadened the recommendation. The abstract's concluding statement now reads: "These findings support continued investment in improved screening and treatment strategies to complement existing prevention efforts and achieve cervical cancer elimination decades sooner, and in the continued development of HPV therapeutic vaccines as a complementary future option." The final Discussion paragraph has been updated in the same spirit, noting that while the largest near-term gains in Rwanda come from scaling existing screening and mass campaign approaches, therapeutic vaccine development remains a valuable long-term investment for settings where screening infrastructure is difficult to scale.

---

## Reviewer 2

### R2.1: Target audience

> *It was not clear to me who the target audience is for this paper. As neither of the 2 therapeutic vaccines are apparently available (or indeed close to being available) and we do not therefore know their potential efficacy (including over the long-term), delivery requirements, or cost then the analyses of these interventions does not seem aimed at policy makers. So is the intended audience research funders?*

**Response:** We appreciate this challenge and have addressed it in two ways.

First, on audience: the paper is intended primarily for two groups. (1) Country planners in settings that, like Rwanda, have already achieved high prophylactic vaccination coverage and are now considering how to allocate resources against the residual burden. For these planners, the paper quantifies the timing and magnitude of impact from several distinct pathways, including options that do not rely on any new product (scaled screen-and-treat, HPV-Faster). (2) Research funders and product developers who need modeling evidence to prioritize investment in therapeutic vaccine development and clinical trials. For these audiences, the paper quantifies the potential impact of two very different TxV profiles across a range of delivery contexts, showing where a TxV would add the most value and how sensitive that value is to the year the product becomes available. We have added a sentence to the introduction clarifying this intended audience.

Second, on the specific concern about 2030 availability: we have added Figure S3 showing the sensitivity of the TxV-carrying strategies to intro years 2030 through 2050 (see also response to R1.3). All four TxV-carrying strategies remain net-beneficial even under a 2050 intro, showing that the paper's qualitative conclusions do not depend on TxV availability by 2030.

---

### R2.2: Novelty and cost

> *The analyses demonstrating the advantages of improving implementation of screening and treatment strategies that are known to be effective provide some potentially useful estimates of the burden of disease that might be prevented if the triage/treat aspect of cervical cancer care was radically improved. The possible effect of a mass screening campaign that is linked to radically improved triage/treat aspects of care are also potentially interesting. However, I am not sure that these results are very surprising in principle. If access to effective treatment is low then improving access to treatment should improve outcomes. If these analyses are aimed at policy makers then I suspect the critical issue for them is what might it cost to adopt either of these strategies and what would be the relative benefits since policy makers are faced with difficult resource allocation decisions between alternative effective interventions for cancer care and in other areas.*

**Response:** We agree that a policy-facing analysis of this space benefits from a cost lens, and we have added one in the form of a threshold cost-effectiveness analysis (new Table S3 in the supplementary materials, described in Methods and Results). We solve for the maximum therapeutic vaccine per-dose price at which each therapeutic-vaccine-carrying strategy would be cost-effective, at three willingness-to-pay tiers (an opportunity-cost anchor of $130/DALY averted, and Rwanda GDP-per-capita anchors of $450 and $900/DALY), against two comparators (status quo screen-triage-treat at 18% coverage, and the strongest non-therapeutic-vaccine alternative of scaled screen-and-treat at 70% coverage). Unit costs for the non-therapeutic ingredients (HPV DNA screens, ablation, LEEP, radiation, prophylactic doses) are drawn from published literature ranges. All streams are discounted at 3%/yr from 2030. The headline finding, described in Results, is that a lesion-regressing therapeutic vaccine could be priced at around $43/dose (range $38-$55 across cost inputs) and still be cost-effective vs status quo at $450/DALY, whether delivered as a mass campaign or embedded in existing screening; a therapeutic vaccine added on top of triage plus ablative treatment is dominated at every willingness-to-pay tier we examined.

We also note the reviewer's related point about novelty. Our aim was not to produce a surprising finding in the mechanistic sense (we agree that improving treatment access improves outcomes), but rather to provide Rwanda-specific quantitative estimates that translate the general principle into a decision-relevant comparison: how much residual burden remains under status quo, how large the marginal gains are from each of seven distinct strategies, how similar in scale several of these strategies actually turn out to be, and what per-dose price a therapeutic vaccine would need to hit to be a competitive addition to Rwanda's toolkit. To our knowledge, this Rwanda-specific comparison has not been reported in prior modeling work.

We do, however, want to be transparent about the limits of the threshold analysis. It uses published literature ranges for unit costs rather than Rwanda-specific micro-costing, does not account for the fixed capital and workforce investments required to scale screening or deliver a mass campaign, and captures only the health system cost perspective (not the broader societal costs of illness). A full cost-effectiveness study will be needed once the therapeutic vaccine candidates are further along in development and Rwanda-specific delivery cost data can be collected. We have added a caveat to this effect in the Discussion (limitations paragraph).

---

### R2.3: Screen-treat vs screen-triage-treat

> *Please clarify what the actual procedure is for 'screen – treat' in comparison to 'screen-triage-treat'. In the 'screen – treat' process direct visualisation of the cervix is still required so the only change is not using acetic acid with the consequence that all HPV+ve women are treated. Will the women be asked if they accept treatment even if the provider is unsure they have 'disease' – it would seem women would need to be carefully consented? Linked to this should any risks or harms of ablation be included in the models as the screen-treat approach will massively increase the number of women getting a procedure which the majority may not need.*

**Response:** These are important points and we address each in turn.

On procedural clarity: we have expanded the Methods description of the screen-treat scenario to make explicit that direct visualisation of the cervix is still required. The only change relative to screen-triage-treat is that VIA is not used to filter treatment eligibility, and all HPV DNA-positive women are offered same-visit ablation. The 90% same-visit attendance and 25% LTFU-at-treatment assumptions are the same in both scenarios, matching the source data (Muhimpundu et al. 2021, ref 28).

On consent and acceptance: the reviewer raises an important concern that women may be less willing to accept ablation without a triage confirmation of visible disease. In response, we have added a sensitivity analysis in which 50% of eligible women decline treatment under screen-treat (versus 25% in the base case). Even under this pessimistic acceptance assumption, screen-treat at 70% coverage would avert 23,000 cases relative to status quo, substantially more than the equivalent screen-triage-treat scenario (16,000 averted). This result is now presented in the main text (Results §Screening scenarios).

On ablation harms: we agree that ablation carries risk of complications and that the screen-treat approach substantially increases the volume of procedures. We have not modeled ablation-specific harms directly. However, Figure 2C reports the total number of ablations required over 2030-2100 under each scenario (2.1M to 5.5M depending on coverage and triage strategy). We have added text to the Discussion (§Limitations) acknowledging this limitation and noting that scenarios without triage require 55 to 90% more ablations depending on coverage level, with implications for both resources and women's wellbeing. A full accounting of ablation harms would require country-specific data on complication rates that is beyond the scope of this analysis.

---

### R2.4: Basis for therapeutic vaccine assumptions

> *What is the basis for these assumptions? There are no references in this section? [followed by quoted paragraph on TxV efficacy assumptions]*

**Response:** The reviewer is correct that this section previously lacked references, and we have substantially revised the Methods paragraph on therapeutic vaccine assumptions to trace each choice to source material. Our assumptions draw on three primary sources: (i) the 2024 WHO Preferred Product Characteristics (PPCs) for therapeutic HPV vaccines and the underlying consultation report, which explicitly distinguish two product archetypes that correspond directly to our virus-clearing and lesion-regressing scenarios and provide indicative efficacy ranges [refs 29 and 30]; (ii) the trial-design consensus document by Dull et al., which documents endpoint definitions, spontaneous HPV clearance rates, and the 90% and 50% target efficacies now standard in trial power calculations [ref 15]; and (iii) prior modeling studies that examined near-identical product profiles in Uganda [ref 12], across multiple settings [refs 31 and 17], and in China [ref 16]. The revised Methods paragraph traces each numeric value to one or more of these sources.

We also want to be transparent about three specific parameter choices where the reviewer's concern is well-founded and where we have adjusted the framing in the manuscript:

- **90% efficacy against high-grade lesions** (lesion-regressing archetype) exceeds efficacies observed to date in phase II/III trials, which have been modest [ref 14]. The revised Methods now describes this explicitly as a target-product-profile assumption in the sense of the WHO PPCs, consistent with the "0-90" and "50-50" product profiles used in the Cohen and Canfell modeling reports [refs 31 and 17].
- **3-month delay to action** is not directly observed. The revised Methods describes this as a deliberately conservative simplification representing the time for a CD8-mediated cellular response to develop and act. Trial endpoints are assessed at longer timepoints (36 weeks for lesion regression, at least 24 months for virologic clearance), but these windows incorporate confirmatory testing rather than measuring time-to-effect [ref 15]. Our choice of 3 months avoids an unrealistic instantaneous-action assumption while still allowing the vaccine to influence cases arising in cohorts vaccinated near the end of the projection horizon.
- **50% efficacy against infection in the lesion-regressing product** differs from the Daffodil report's UC2 assumption of 0% [ref 17]. We retain a partial infection-clearing effect on the basis that an E6/E7-directed cellular response would be expected to act on productive infection as well as on established lesions; the revised Methods now says so explicitly.

We have not run a formal sensitivity analysis over the 3-month delay parameter given the length of the current revision, but note that shorter delays would improve TxV impact modestly and longer delays would reduce it, without changing the overall ordering of the strategies.

---

### R2.5: Ablation workforce implications

> *It is reported that 2.5-4M additional ablations would be required under the screen-treat scenario. It would be useful to know how many ablations are currently conducted and based on even a simple calculation how many additional providers would be needed as the assumption is that there is an expansion in the workforce to accommodate this.*

**Response:** We agree that workforce capacity is a critical practical constraint on any scaled cervical cancer program. We have added a workforce sensitivity analysis (Figure S5) that caps the annual number of ablative and excisional procedures at three levels approximating 1.5x, 2x, and 5x the current national throughput at 18% coverage (about 34,700 ablations per year in 2028, as estimated by the model). The analysis is applied to the four S&T-family scenarios at 70% coverage.

The findings show a large divergence between strategies that rely on ablation and strategies that leverage a therapeutic vaccine:

For ablation-dependent strategies (S&T&T 70% and S&T 70%), even the most permissive cap does not fully absorb the demand generated by 70% screening coverage. Under all three caps, the incremental benefit from scaling coverage from 18% to 70% is essentially eliminated (medians near zero cases averted for both scenarios, with wide 10-90% intervals crossing zero), because women identified as needing treatment queue up and some progress to cancer before receiving it.

For strategies that use a therapeutic vaccine, the workforce cap has a much smaller effect. S&TxV 70% averts around 37,000 cases regardless of the ablation cap, because the therapeutic vaccine (delivered by injection rather than a surgical procedure) is not subject to the same workforce constraint. S&TxV+T&T 70% loses about half its no-cap benefit (from +31,000 to +15,000) because it still uses ablation as a secondary treatment path.

The practical implication is that scaling screening beyond current national workforce capacity delivers its intended benefit only in strategies that reduce the burden on ablation providers, either by substituting a therapeutic vaccine (S&TxV) or by making substantial workforce investment. In the absence of such investment, the incremental benefit of scaled screening is largely captured by the treatment queue rather than by patients.

We were unable to obtain Rwanda-specific data on current ablation throughput or provider counts, and instead take the model's estimate of current 18%-coverage throughput as a baseline reference.

---

### R2.6: Sensitivity on VIA sensitivity and loss to follow up

> *It is stated that "Two key gaps in the delivery of these programs are (1) the poor sensitivity of VIA, and (2) the 25% loss to follow up between the initial determination of treatment eligibility (either on the basis of a positive HPV test or VIA triage) and treatment." Does this not argue for a sensitivity analysis to explore the effects of addressing either or both of these determinants?*

**Response:** We agree that these two gaps warrant sensitivity analysis, and both are addressed in the revised paper.

On VIA sensitivity: the comparison between screen-triage-treat and screen-treat is itself a sensitivity analysis on the contribution of VIA. Removing VIA triage at 70% coverage increases the cumulative cancers averted from 16,000 (with VIA) to 37,000 (without VIA), a greater-than-twofold improvement attributable to VIA's poor per-lesion sensitivity. This finding is discussed in Results §Screening scenarios.

On loss to follow up: we agree that acceptance of same-visit ablation without a triage confirmation of visible disease may be substantially lower than the 75% we assume in the base case, and this is a genuine implementation risk that could materially reduce the benefit of removing triage. To quantify this, we added a sensitivity scenario at S&T 70% with 50% LTFU. In this scenario, screen-treat averts 23,000 cases relative to status quo, compared with 37,000 in the base case. The marginal gain from bypassing VIA triage (relative to the equivalent screen-triage-treat scenario at 16,000 averted) is thus reduced from around 21,000 to around 7,000 additional averted cancers. This sensitivity is now presented in Results §Screening scenarios.

---

### R2.7: Tables of results

> *It would be good to have some tables of results as appendices to accompany the figures to illustrate the numbers requiring interventions.*

**Response:** We agree, and have added Table S2 to the supplementary materials showing cumulative intervention counts required over 2030-2100 for each of the 23 scenarios. Table S2 reports, per scenario, cumulative cancers, screens, and precancerous treatments for all scenarios. The table complements the figure-only summaries of impact by making the operational scale of each strategy explicit, and is also cited alongside our discussion of workforce implications (see response to R2.5).

---

### R2.8: Cost-effectiveness threshold

> *It is appreciated that this is not a cost-effectiveness study but it should be possible to estimate what the cost per life-saved or cost per-DALY averted would need to be to suggest that therapeutic vaccines would be a good buy when compared with other cost-effective interventions given that many LMIC including Rwanda are making difficult choices on what to invest in.*

**Response:** We agree, and have added this analysis (see also the response to R2.2). New Table S3 in the supplementary materials reports, for every therapeutic-vaccine-carrying strategy, the maximum per-dose price at which the strategy would be cost-effective against two comparators (status-quo screen-triage-treat at 18% coverage, and the strongest non-therapeutic-vaccine alternative of scaled screen-and-treat at 70% coverage), at three willingness-to-pay tiers ($130, $450, and $900 per DALY averted; the latter two correspond to 0.5 and 1 × Rwanda's GDP per capita). Cumulative DALYs 2030-2100 are computed via the HPVsim incidence-based DALY analyzer using GBD 2017 disability weights; costs use published unit-cost ranges from the sub-Saharan African cervical cancer costing literature; both streams are discounted at 3%/yr from 2030. Full methodology is described in the Methods section of the revised manuscript.

The main findings: at $450/DALY (0.5 × GDP per capita), a lesion-regressing therapeutic vaccine would be cost-effective at around $43/dose whether delivered as a mass campaign or embedded in existing screening (Table S3). Against the more demanding comparator of a fully scaled-up screen-and-treat program at 70% coverage, the threshold rises substantially (to around $72-$90/dose for the same two strategies), reflecting the operational savings from replacing a large fraction of the ablative-treatment workload with a therapeutic dose. In contrast, layering a therapeutic vaccine on top of both triage and treatment (the screen-triage-treat-plus-therapeutic-vaccine family) is dominated at every willingness-to-pay tier we examined, because the additional procedures required outweigh the marginal health gains.

We are careful to flag in Methods, Results, and Discussion that this is a threshold analysis rather than a full cost-effectiveness study, and that it should be revisited once therapeutic vaccine candidates are further along in development and Rwanda-specific delivery cost data are available.

---

### R2.9: Same-visit VIA and treatment

> *Is it true that VIA and treatment cannot be done in the same visit – especially as HPV testing is the primary screening tool?*

**Response:** We appreciate the reviewer highlighting this. To clarify: our model does assume that VIA triage and treatment are performed in the same visit. The 25% loss to follow up we apply represents women who do not accept or receive treatment at that same visit (for reasons including provider availability, procedural anxiety, or clinical judgment), rather than a between-visit dropout. This matches the source data from the Rwanda program (Muhimpundu et al. 2021, ref 28), where the shortfall was attributed to low availability of trained providers rather than to a multi-visit protocol. We have edited the Discussion (paragraph on screening challenges) to remove the phrasing that suggested VIA and treatment are performed at separate visits, and to reflect the same-visit but incomplete-uptake framing that we actually model.

---

### R2.10: Multi-purpose screening infrastructure

> *There are challenges establishing screening infrastructure of course – but such infrastructure is often multi-purpose and might be developed as part of wider women's health initiatives or primary care services including breast cancer screening and other services etc.*

**Response:** We agree that screening infrastructure is often multi-purpose and that investment in cervical cancer screening can be shared with, and support, other women's health services (breast cancer screening, cervical HPV DNA testing, contraception counselling, general reproductive-health outreach). We have added a sentence to the Discussion acknowledging this and noting that it strengthens the practical and economic case for investing in the screening pathways we model, since the fixed costs of establishing and operating the infrastructure are shared across multiple health domains.

---

### R2.11: Missing denominator

> *Introduction P3, an estimated age-standardized rate (ASR) of cervical cancer incidence of 28 - should this have a denominator?*

**Response:** Yes, thank you for catching this. The Introduction now reads "an estimated age-standardized rate (ASR) of cervical cancer incidence of 28 per 100,000 women (18)."
