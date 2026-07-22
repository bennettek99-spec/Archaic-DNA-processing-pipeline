# Highest archaic ancestry in the analyzed AADR subset

> Ancient individuals retained by the existing global Phase 2-4 scan. Neanderthal is an f4-ratio percentage; Denisovan is a relative D-statistic. Combined percentage is not estimable under the validated model.

## Conclusions

- Highest numerical Neanderthal estimate: **EV16A.SG**, 5.35% (95% CI 1.22-9.49%), Low confidence.
- No individual in this subset passes the pipeline's high/elite credibility thresholds; raw point estimates must therefore be treated as low-power screening results.
- Strongest Denisovan-related signal in this subset: **EV16A.SG**, D=0.00563, Z=0.22; no percentage is claimed.
- Highest combined archaic percentage: **not identifiable** because the validated Denisovan statistic is not on a percentage scale.
- Genotype sensitivity tests assess robustness of the genome-wide estimate but cannot establish recent admixture without validated segment or read-level evidence.

## Top credibility-aware Neanderthal results

| genetic_id | neanderthal_pct | neanderthal_ci_low_pct | alpha_nSNP | artifact_risk_score | credibility_class |
|---|---|---|---|---|---|

## Raw results most likely affected by artifacts

| genetic_id | neanderthal_pct | neanderthal_ci_low_pct | alpha_nSNP | artifact_risk_score | credibility_class |
|---|---|---|---|---|---|
| EV16A.SG | 5.3539 | 1.2195 | 15994 | 63.62 | Low confidence |

## Sensitivity-tested transversion ranking

| genetic_id | alpha_standard | alpha_transversion | alpha_alt_outgroup | loco_min | loco_max | bootstrap_q025 | artifact_risk_score | credibility_class |
|---|---|---|---|---|---|---|---|---|
| EV16A.SG | 0.0535 | 0.0699 | 0.0473 | 0.0431 | 0.0613 | 0.0186 | 63.62 | Low confidence |

## Context: Oase 1 and Bacho Kiro

| genetic_id | date_bp | neanderthal_pct | neanderthal_ci_low_pct | alpha_nSNP | alpha_transversion | loco_min | bootstrap_q025 | credibility_class |
|---|---|---|---|---|---|---|---|---|

No Oase 1 or Bacho Kiro comparison individual is present in this subset.

## Threshold provenance

- Broad informative-SNP floor: 10,000.
- High-confidence floor: 200,000 (existing pipeline threshold).
- Elite floor: 200,000, observed 75th percentile constrained by the high-confidence floor.
- Elite maximum SE: 2.109 percentage points (observed 25th percentile).

## Limitations

- EIGENSTRAT cannot support terminal-base trimming or higher base/mapping-quality thresholds; BAM/CRAM is required.
- One validated Denisovan reference cannot yield an absolute Denisovan f4-ratio or combined percentage.
- The metadata residual model is an outlier screen, not evidence of recent admixture.
- General segment detection is not validated. The separate Oase1 read-level workflow is contextual evidence only.
- Close-relative, batch, and genetic-ancestry-cluster controls need dedicated genotype/read analyses beyond duplicate-library and metadata residual checks.

## Reproduction

```bash
python -m archaic.highest_archaic --aadr-data PATH --metadata PATH --config configs/highest_archaic.yaml --output results/individual_EV16A --threads AUTO --resume
```
