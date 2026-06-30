# Genome-wide search for unexpected archaic introgression in ancient Eurasian genomes — findings

*Exploratory analysis on the AADR v66.p1 1240K panel. All results are hypotheses
subject to the stated limitations; no claim here is a validated biological
discovery.*

## Summary

Across **15,443 ancient Eurasian genomes**, after estimating each individual's
Neanderthal ancestry and asking whether it deviates from what is expected given
their **genetic ancestry, geography, and age**, **no individual shows a deviation
that survives multiple-testing correction.** The maximum standardized residual
(|z| = 3.90 over 9,252 high-confidence tests) is indistinguishable from the
expected extreme of that many draws from pure noise (3.87); zero samples pass
Bonferroni (z>4.55) or BH-FDR 5%. The nominal top candidates are explicable by
**known ancestry effects** (recent African/Levantine admixture → less Neanderthal),
not by novel archaic introgression.

This is the scientifically expected outcome: within-Eurasian variation in
Neanderthal ancestry is small (~0.1–0.4 percentage points between major groups)
and, at single-genome resolution on capture data, is dominated by measurement
noise. The value of the pipeline is that it **quantifies this rigorously** and
would surface a genuine outlier if one existed.

## Method (validated end-to-end)

| Quantity | Estimator |
|---|---|
| Neanderthal proportion | f4-ratio α = f4(Altai,Chimp; X,Mbuti) / f4(Altai,Chimp; Vindija,Mbuti) |
| Neanderthal affinity | D(X, Mbuti; Altai, Chimp) |
| Denisovan affinity | D(X, Mbuti; Denisova, Chimp) |
| Uncertainty | 50-block delete-one jackknife (per individual) |
| Ancestry | PCA (30k autosomal SNPs, randomized SVD); PC1 = West–East Eurasian (r=−0.67 with longitude) |
| Expectation | precision-weighted (1/SE²) mean of the K=80 nearest high-confidence neighbours in ancestry+geo+time space |
| Residual | z = (α_adj − E[α]) / √(SE_self² + σ_bio² + SE_E²), with σ_bio² = max(0, Var_neighbours − mean SE_n²) |

**Phase 1 validation (HO panel, 7/7 gates):** reproduced absolute Neanderthal
levels (~2–3% non-Africans; Yoruba −0.12%±0.28 ≈ 0), the East-Asian excess
(D(French,Han;Altai,Yoruba) Z=−2.3), the published rank order
Sardinian<French<Han<Papuan, the "100% Neanderthal" scale (Neanderthals-as-test
read 97–99%), and the Denisovan signal (Papuan Z=6.5, others ≈0).

## Dataset (Phase 2)

23,089 panel individuals → **15,443 retained** Eurasian ancients (excluded: 3,967
present-day, 1,950 non-Eurasian, 1,305 below the 30k-SNP floor, 402 contamination
CRITICAL/FAIL, 22 references). Median 618k SNPs/sample; **9,252 high-confidence**
(≥200k SNP, uncontaminated). Weighted-mean Neanderthal across all = **2.116%**.

## Bias control (Phase 4)

No gross technical confound: |weighted corr(α, covariate)| < 0.07 for usable-SNP
count, coverage, DNA damage, contamination, and age. The one real nuisance axis is
sequencing data type (shotgun 2.29% vs capture 2.00%, ~0.3pp), removed by an
additive offset before residual analysis.

## Result (Phases 5–6)

- Predictive model R² (ancestry PCs + geo + time) = **0.033** — within-ancestry
  Neanderthal variation is mostly measurement noise, as expected.
- Standardized residuals: 95th pct |z| = 1.58, 99th = 2.10, **max = 3.90**.
- Empirical tails are **narrower than N(0,1)** (1.7% beyond |z|>1.96 vs 5%
  expected) — the error model is conservative, so the null is not an artifact of
  inflated error bars.
- **0 / 9,252 pass Bonferroni or FDR.**
- Nominal extremes (e.g. `Georgia_Tkhina_19thCentury-oAfrican` at 0.41%,
  `Austria_ImperialRoman-oAfrica-oLevant`) are AADR-flagged ancestry outliers whose
  low Neanderthal is **expected from their actual admixed ancestry**, not anomalous.

## Robustness (Phase 9) and candidates (Phase 7)

The null is **robust**: 0 samples pass Bonferroni/FDR under every perturbation —
neighbour count K∈{40,80,160}, three 50% reference subsamples, and a tighter
≥400k-SNP floor — and the bootstrap (B=100) max|z| = 3.93 [3.79, 4.03] never
reaches z*=4.55. The residual *ranking* is itself stable (Spearman ρ≈0.99;
15–18/20 top candidates preserved), so the same individuals are consistently the
most deviant, but none rises above chance. The strongest "less-than-expected" hits
are AADR-flagged ancestry outliers (`Zana` = Georgia_Tkhina_19thC-oAfrican 0.41%,
Austria_ImperialRoman-oAfrica-oLevant, Hungary_EarlyAvar-oHighEastAsia) whose low
Neanderthal follows directly from recent African/Levantine/East-Asian admixture.
Per-candidate reports: `results/phase7_reports/`. Figures: `results/figures/`.

## Limitations

- **Single-genome capture data cannot resolve the small within-Eurasian signal.**
  Per-sample SE on α is ~0.3–0.6% (≥300k SNP) — comparable to or larger than the
  real between-group differences. Detecting a genuine individual outlier would
  require it to be *large* (>~2pp), which none are.
- Ancestry PCs are computed on a 30k-SNP subset (sufficient for major axes, not
  fine structure).
- Denisovan ancestry is relative-only (one high-coverage Denisovan; ≈0 in West
  Eurasia regardless).
- f4-ratio capped by Vindija coverage (528k SNP on 1240K).
- Duplicate individuals / close relatives are not de-duplicated; this is
  conservative for outlier calling (a duplicate neighbour pulls the expectation
  toward the sample, shrinking its residual).

## What a real discovery would look like here

An individual with high coverage (≥300k SNP), PASS contamination, α deviating from
its ancestry+geo+time expectation by ≫2pp at z>4.55, robust to neighbour-set and
SNP-subset perturbation (Phase 9), with no `-o` ancestry-outlier flag and no
matching low-Neanderthal admixture source. No such individual exists in this panel.

## Reproduce

```
python phase1_validate.py                 # estimator validation (HO)
python phase2_prepare.py   --panel 1240k  # dataset + QC + exclusion log
python phase3_estimate.py  --panel 1240k  # per-genome archaic estimates (~27 min)
python phase4_normalize.py --panel 1240k  # bias diagnostics + analysis table
python phase5_pca.py       --panel 1240k  # ancestry PCA
python phase6_outliers.py  --panel 1240k  # expectation + standardized residuals
python phase7_reports.py   1240k          # per-candidate investigation reports
python phase8_figures.py   1240k          # publication figures
python phase9_robustness.py --panel 1240k # robustness / sensitivity
```

Outputs in `results/`. See `README.md` for the full method and citations.
