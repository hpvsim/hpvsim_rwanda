# Supplementary materials for Strategies to accelerate cervical cancer elimination: a modeling study in Rwanda

## Figure S1: Rwanda's screening and treatment algorithm

![Figure S1: Rwanda's screening and treatment algorithm](../figures/v1_submission/figS1_screening_algorithm.png)

*Caption: Rwanda's screening and treatment algorithm.*

## Figure S2: model calibration

![Figure S2: model calibration](../figures/figS2_calib.png)

*Caption: model calibration to data. Black markers indicate data, with sources listed below. The first row shows the fit to cancer by age for the whole population (data from GLOBOCAN) and disaggregated by HIV status. The second row shows the fit to additional GLOBOCAN data, including the age-standardized rate of cancer incidence (left panel), and the distribution of LSILs and cancers by genotype. The bottom row shows the model fit to HIV data, which is sourced from UNAIDS. For each indicator, we plot the medians and 10-90% ranges of the 50 best-fitting parameter sets.*

## Figure S3: natural history under baseline programs

![Figure S3: natural history under baseline programs](../figures/figS3_natural_history.png)

*Caption: natural history of HPV infection and cervical cancer in Rwandan women under the baseline scenario (routine prophylactic vaccination plus status-quo 18% screen-and-treat), across cancers occurring 2020-2050. (A) Age at causal HPV infection, showing the distribution and 25th/50th/75th percentiles. (B) Dwell times (in years) from causal infection to CIN2+, from CIN2+ to cancer, and total time from causal infection to cancer.*

## Figure S4: TxV introduction-year sensitivity

![Figure S4: TxV introduction-year sensitivity](../figures/figS4_txv_intro.png)

*Caption: cumulative cancers averted 2030-2100 relative to the S&T&T 18% baseline as a function of the year the therapeutic vaccine becomes available, across the four TxV-carrying strategies at 70% coverage. Shaded regions show 10-90% intervals across the 50 best-fitting parameter sets.*

## Figure S5: workforce capacity sensitivity

![Figure S5: workforce capacity sensitivity](../figures/figS5_workforce.png)

*Caption: cumulative cancers averted 2030-2100 relative to the S&T&T 18% baseline as a function of workforce capacity for ablation and excision procedures, expressed as a multiplier of the current (2028) national throughput at 18% coverage (~34,700 ablations per year). Lines show the four S&T-family strategies at 70% coverage; dashed lines show the no-cap value for each. Shaded regions show 10-90% intervals across the 50 best-fitting parameter sets.*

## Table S1: calibrated parameters — priors and posteriors

| parameter | description | prior_source | prior_low | prior_high | posterior_median | posterior_ci_lo | posterior_ci_hi |
|---|---|---|---|---|---|---|---|
| beta | HPV transmission probability per act | HPVsim defaults; range widened to accommodate Rwanda-specific fit | 0.020 | 0.700 | 0.460 | 0.280 | 0.680 |
| hpv16.cin_fn.k | HPV16 dose-response shape (dysplasia progression rate) | Prior from HPVsim natural-history calibration (Stuart et al. 2024) | 0.250 | 0.350 | 0.330 | 0.260 | 0.350 |
| hpv18.cin_fn.k | HPV18 dose-response shape (dysplasia progression rate) | Prior from HPVsim natural-history calibration (Stuart et al. 2024) | 0.200 | 0.300 | 0.240 | 0.220 | 0.278 |
| hi5.cin_fn.k | Pooled Hi5 (31/33/45/52/58) dose-response shape | Prior from HPVsim natural-history calibration (Stuart et al. 2024) | 0.150 | 0.250 | 0.200 | 0.150 | 0.240 |
| ohr.cin_fn.k | Pooled OHR (35/39/51/56/59) dose-response shape | Prior from HPVsim natural-history calibration (Stuart et al. 2024) | 0.150 | 0.250 | 0.170 | 0.150 | 0.250 |
| age_risk.risk | Age-dependent multiplier on HPV acquisition risk | Introduced in this study to fit the observed age distribution of causal infection | 1.500 | 3.500 | 3.000 | 1.500 | 3.500 |
| imm_init.low | Lower bound of initial humoral immunity after natural clearance | Introduced in this study; prior from HPVsim natural-history defaults | 0.300 | 0.850 | 0.600 | 0.500 | 0.850 |
| cell_imm_init.low | Lower bound of initial cell-mediated immunity after lesion regression | Introduced in this study; prior from HPVsim natural-history defaults | 0.200 | 0.700 | 0.250 | 0.200 | 0.689 |
| network.m_cross_layer | Male probability of a concurrent partnership across marital/casual layers | Prior from HPVsim network calibration (Stuart et al. 2024) | 0.350 | 0.950 | 0.850 | 0.711 | 0.900 |
| network.f_cross_layer | Female probability of a concurrent partnership across marital/casual layers | Prior from HPVsim network calibration (Stuart et al. 2024) | 0.200 | 0.950 | 0.800 | 0.750 | 0.939 |
| network.m_partners_casual | Male mean casual partners per year (Poisson rate) | Prior from HPVsim network calibration (Stuart et al. 2024) | 0.100 | 0.600 | 0.250 | 0.100 | 0.527 |
| network.f_partners_casual | Female mean casual partners per year (Poisson rate) | Prior from HPVsim network calibration (Stuart et al. 2024) | 0.100 | 0.600 | 0.350 | 0.111 | 0.550 |
| hiv.rel_sus_lo | HIV+ relative HPV susceptibility at low CD4 | Prior informed by Liu et al. 2018 meta-analysis (ref 25) | 2.000 | 5.000 | 3.000 | 2.000 | 5.000 |
| hiv.rel_sus_hi | HIV+ relative HPV susceptibility at high CD4 | Prior informed by Liu et al. 2018 meta-analysis (ref 25) | 2.000 | 4.000 | 2.250 | 2.000 | 3.500 |
| hiv.rel_sev_lo | HIV+ relative multiplier on HPV disease progression at low CD4 | Prior informed by Liu et al. 2018 meta-analysis (ref 25) | 1.500 | 5.000 | 2.500 | 1.500 | 4.694 |
| hiv.rel_sev_hi | HIV+ relative multiplier on HPV disease progression at high CD4 | Prior informed by Liu et al. 2018 meta-analysis (ref 25) | 1.500 | 5.000 | 1.750 | 1.500 | 2.194 |
| hiv.p_effective_art | Probability that ART fully suppresses HIV-driven HPV effects | Prior from UNAIDS ART effectiveness estimates | 0.700 | 0.950 | 0.840 | 0.760 | 0.950 |

*Caption: Prior distributions and posterior estimates for the 17 parameters of the model that were adjusted through calibration. The posterior summary reports the median and 95% credible interval across the 50 best-fitting parameter sets identified from 1,500 calibration trials.*

## Table S2: cumulative intervention counts by scenario, 2030-2100

| scenario | cumulative_cancers | screens | ablations | LEEP_procedures | radiation_courses | therapeutic_vaccine_doses | prophylactic_vaccine_doses |
|---|---|---|---|---|---|---|---|
| No interventions | 304,320 | 0 | 0 | 0 | 0 | 0 | 0 |
| S&T&T 18% | 100,969 | 7,833,809 | 1,102,074 | 53,465 | 119,143 | 0 | 13,539,722 |
| S&T&T 35% | 94,345 | 14,323,237 | 1,880,674 | 88,911 | 214,607 | 0 | 13,558,636 |
| S&T&T 70% | 87,232 | 27,501,360 | 3,582,783 | 173,055 | 360,557 | 0 | 13,543,892 |
| S&T 18% | 94,755 | 7,698,433 | 2,079,792 | 94,570 | 118,398 | 0 | 13,557,296 |
| S&T 35% | 85,202 | 14,054,420 | 3,558,805 | 168,141 | 195,842 | 0 | 13,554,466 |
| S&T 70% | 62,753 | 25,477,863 | 5,538,219 | 263,008 | 304,709 | 0 | 13,541,807 |
| S&T 70%, 50% LTFU | 77,720 | 26,688,208 | 4,331,299 | 130,611 | 345,664 | 0 | 13,552,828 |
| S&TxV&T&T 18% | 99,632 | 7,721,517 | 1,063,054 | 50,338 | 109,761 | 2,819,969 | 13,546,424 |
| S&TxV&T&T 35% | 86,794 | 13,791,859 | 1,703,449 | 81,464 | 198,522 | 4,546,055 | 13,543,296 |
| S&TxV&T&T 70% | 71,161 | 25,458,502 | 2,825,628 | 130,164 | 313,346 | 7,481,890 | 13,555,508 |
| S&TxV 18% | 95,901 | 7,707,964 | 33,211 | 82,358 | 4,021 | 2,787,949 | 13,561,317 |
| S&TxV 35% | 86,422 | 14,011,082 | 67,912 | 131,355 | 7,893 | 4,745,917 | 13,538,680 |
| S&TxV 70% | 64,090 | 25,755,466 | 154,439 | 230,691 | 18,169 | 7,666,115 | 13,549,402 |
| HPV-Faster 18% | 91,185 | 7,654,797 | 969,080 | 47,210 | 104,250 | 0 | 13,570,252 |
| HPV-Faster 35% | 83,117 | 7,505,570 | 915,615 | 42,296 | 92,336 | 0 | 13,552,679 |
| HPV-Faster 70% | 60,879 | 7,200,266 | 718,433 | 34,998 | 81,911 | 0 | 13,541,658 |
| Mass TxV 90/0, 18% | 95,852 | 7,507,804 | 985,463 | 46,913 | 111,995 | 1,694,215 | 13,556,551 |
| Mass TxV 90/0, 35% | 91,083 | 7,620,692 | 1,000,802 | 48,998 | 111,399 | 3,298,626 | 13,546,573 |
| Mass TxV 90/0, 70% | 79,093 | 6,967,936 | 747,474 | 34,849 | 93,825 | 6,564,786 | 13,552,232 |
| Mass TxV 50/90, 18% | 93,466 | 7,733,729 | 1,019,567 | 47,657 | 113,782 | 1,694,215 | 13,552,083 |
| Mass TxV 50/90, 35% | 81,108 | 7,543,249 | 919,934 | 42,147 | 96,655 | 3,298,626 | 13,547,466 |
| Mass TxV 50/90, 70% | 61,685 | 6,851,325 | 615,970 | 28,594 | 70,741 | 6,564,786 | 13,542,254 |

*Caption: Cumulative counts of interventions delivered under each scenario over 2030-2100. Numbers show medians across the 50 best-fitting parameter sets. Columns cover: HPV DNA screens; ablative treatments; LEEP (large loop excision of the transformation zone) procedures; radiation courses for invasive cancer; therapeutic HPV vaccine doses; and prophylactic HPV vaccine doses (routine plus any campaign-linked doses).*
