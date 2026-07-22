# Dedicated archaic-introgression analysis: EV16A.SG

## Executive verdict

EV16A.SG is an **interesting but low-confidence Neanderthal-ancestry candidate**. The pipeline's standard estimate is **5.35%** (SE **2.11 percentage points**; 95% CI **1.22-9.49%**) from only **15,994 informative SNPs**. The signal does not collapse in the available genotype-level controls, but the interval is wide, the transversion-only estimate is much less precise, and the sample does not pass the pipeline's 200,000-informative-SNP high-confidence floor.

The defensible conclusion is therefore not that EV16A carried 5.35% Neanderthal ancestry as a settled biological fact, nor that EV16A had a recent Neanderthal ancestor. It is that EV16A has an elevated point estimate worth read-level and segment-based follow-up. The existing AADR genotype panel cannot supply that evidence.

There is **no Denisovan signal**: `D_Den = 0.00563`, SE `0.02599`, Z `0.22`. `D_Den` is a relative affinity statistic, not a percentage, so Denisovan and combined-archaic percentages are not estimable.

## Sample identity and provenance

| Field | Value |
|---|---|
| Genetic ID | EV16A.SG |
| Persistent AADR ID | 85430 |
| Individual / skeletal code | EV16A / EV16A |
| Material | Tooth |
| Group | Italy_Monteriggioni_IA_Etruscan |
| Locality | Monteriggioni, Tuscany (Siena), Italy |
| Date | 2549 +/- 58 BP; contextual range 700-500 BCE |
| Data type | SG; shotgun data, native pulldown on 3.2M SNP set |
| AADR 1240K SNPs hit | 35,204 (coverage field 0.0311) |
| Molecular sex | Unknown |
| AADR assessment | Pass |
| Contamination / damage covariates | Not available in the metadata used here |
| Publication | Ravasini & Trombetta et al., Genome Biology (2024), https://doi.org/10.1186/s13059-024-03430-4 |

## Core estimates

| Test | Estimate | Precision / range | Informative SNPs | Interpretation |
|---|---:|---:|---:|---|
| Standard Neanderthal f4-ratio | 5.3539% | SE 2.1094 pp; 95% CI 1.2195-9.4883% | 15,994 | Elevated point estimate, low power |
| 200-block bootstrap | 5.3583% mean | 2.5-97.5%: 1.8596-9.4790% | Block resampling | Consistent with the jackknife, still wide |
| Transversion-only | 6.9936% | SE 3.9214 pp; approximate 95% CI -0.692-14.679% | 4,941 | No damage-associated collapse, but very imprecise |
| Yoruba outgroup | 4.7346% | Point estimate | Panel recomputation | Same direction under alternate African baseline |
| Swapped Altai/Vindija roles | 6.5032% | Point estimate | Panel recomputation | Same direction; 1.77-pp reference range |
| Leave-one-chromosome-out | 4.3090-6.1315% | Range 1.8225 pp | 22 chromosome deletions | Elevated direction is not driven by one chromosome |
| Maximum single-block influence | 0.7211 pp | Dominant block 5 | 50 blocks | No one block explains the full estimate |
| Denisovan affinity | D = 0.005626 | SE 0.025991; Z = 0.216 | 17,776 | No evidence of unusual Denisovan affinity |

The transversion-only result is a damage proxy, not a read-level damage correction. Terminal-base trimming and higher mapping/base-quality filters require BAM/CRAM data.

## Etruscan and Monteriggioni context

The supplied values `+1.64` and `+1.44` are **standardized residuals**, not additional ancestry percentages:

- `+1.64` is EV16A's within-Etruscan residual z-score.
- `+1.44` is the residual after conditioning on genetic ancestry and geography.

Neither reaches the study's conventional absolute-z threshold of 2. Across the 75 Etruscan-context individuals, EV16A has the largest raw point estimate, but the Etruscan peer inverse-variance-weighted estimate excluding EV16A is **1.885 +/- 0.076%**. EV16A differs from that peer estimate by **3.469 percentage points**, corresponding to an approximate **z = 1.64** once EV16A's large SE is included. This is elevated but not a statistically exceptional Etruscan outlier.

Within the seven Monteriggioni samples, EV16A is also the largest raw estimate. The other six give an inverse-variance-weighted estimate of **2.675 +/- 0.528%**; EV16A's difference is approximately **z = 1.23**. All seven Monteriggioni results are below the pipeline's high-confidence SNP floor.

| Monteriggioni sample | Neanderthal estimate | SE (pp) | Informative SNPs | Ancestry-conditioned z |
|---|---:|---:|---:|---:|
| EV16A.SG | 5.354% | 2.109 | 15,994 | +1.444 |
| EV7A.SG | 3.528% | 0.899 | 106,360 | +1.389 |
| EV15A.SG | 3.477% | 1.876 | 24,567 | +0.625 |
| EV16D1.SG | 2.959% | 1.181 | 49,107 | +0.598 |
| EV19.SG | 2.061% | 1.627 | 31,201 | -0.153 |
| EV18.SG | 1.543% | 1.175 | 44,354 | -0.661 |
| EV16D2.SG | 0.894% | 2.023 | 13,986 | -0.690 |

In the full 17,143-person retained ancient-sample scan, EV16A ranks **74th by raw point estimate** (99.57th percentile). This raw rank is not a credibility-aware rank: low-information samples naturally occupy both tails more often because their estimates are noisy.

## Credibility assessment

Pipeline classification: **Broad analysis set / Low confidence**, artifact-risk score **63.62/100**.

The main penalties are:

- only 15,994 informative SNPs, 8.0% of the 200,000-SNP high-confidence floor;
- a 95% interval 8.27 percentage points wide;
- missing contamination and damage covariates;
- moderate reference and chromosome sensitivity in percentage-point terms;
- no validated segment call for EV16A.

The evidence in favor of a real elevation is narrower: standard, transversion-only, alternate-outgroup, swapped-reference, leave-one-chromosome-out, and block-bootstrap estimates all retain a positive/elevated direction. That makes EV16A a sensible follow-up target, but it does not overcome the low information content.

## Biological interpretation boundary

This result supports the statement: **"EV16A is a low-power candidate for above-average Neanderthal ancestry within the analyzed Etruscan dataset."**

It does not support any of the following:

- an exact claim that EV16A had 5.35% Neanderthal ancestry;
- a claim of a recent Neanderthal ancestor;
- a Denisovan ancestry percentage;
- a combined archaic ancestry percentage;
- a segment count or time-since-admixture estimate.

The most decisive next analysis would use the deposited read data (`ENA:PRJEB77116`) to apply terminal trimming, mapping/base-quality filters, contamination-aware checks, and a validated archaic-segment caller. Replication from a second extract/library would be especially valuable because the current genotype result is based on sparse coverage.

## Reproduction

The completed run used AADR v66.1 1240K and the repository's fixed seed `20260714`:

```powershell
.\.venv\Scripts\python.exe -m archaic.highest_archaic `
  --aadr-data PATH\TO\AADR `
  --metadata results\phase4_1240k_global_analysis.csv `
  --excluded results\phase2_1240k_global_excluded.csv `
  --config configs\highest_archaic.yaml `
  --output results\individual_EV16A `
  --subset EV16A.SG `
  --resume
```

Machine-readable results are in `all_sample_archaic_estimates.tsv` and `top_candidate_sensitivity_tests.tsv`; the run manifest records the configuration and input digest. The segment-follow-up table correctly reports that validated EV16A segment evidence is unavailable.
