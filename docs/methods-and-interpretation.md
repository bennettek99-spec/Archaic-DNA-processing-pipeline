# Methods and interpretation

## Core estimators

All core statistics use allele frequencies over autosomal AADR SNPs and
delete-one block-jackknife uncertainty over contiguous genomic blocks.

### Neanderthal proportion

The cross-sample Neanderthal estimate is an f4-ratio:

```text
alpha(X) = f4(Altai, Chimp; X, Mbuti)
           / f4(Altai, Chimp; Vindija, Mbuti)
```

Altai appears in the statistic and Vindija supplies the scale. The output is a
percentage with a block-jackknife standard error and confidence interval.

### Relative affinity channels

```text
D_Nea = D(X, Mbuti; Altai, Chimp)
D_Den = D(X, Mbuti; Denisova, Chimp)
```

These channels measure excess allele sharing relative to the African baseline.
`D_Den` supports ranking and significance testing, but a single clean Denisovan
calibration lineage does not identify an absolute Denisovan ancestry fraction.

### Differential population tests

Tests such as `D(Pop1, Pop2; Altai, Yoruba)` compare archaic sharing between
test populations. The African baseline is deliberate: substituting a deep
outgroup in this position can dilute the recent-introgression contrast with
ancestral variation.

## Non-negotiable interpretation boundary

- Neanderthal f4-ratio: percentage interpretation is supported within the
  documented calibration and panel limitations.
- Denisovan `D_Den`: relative affinity and Z-score only.
- Combined archaic percentage: not calculated by this pipeline.

Published Denisovan percentages can be cited as external context. If they are
combined with this pipeline's Neanderthal estimate in prose, the text must say
that the Denisovan component came from the literature and was not estimated by
this software.

## Validation

The estimator is checked through complementary evidence:

1. AADR controls, including African baselines, modern non-Africans, Papuan
   Denisovan affinity, and archaic genomes used as tests.
2. Comparison with published Neanderthal estimates.
3. Coalescent simulation with known introgression fractions.
4. Concordance checks against ADMIXTOOLS 2.
5. Synthetic packed-panel smoke tests for file-format and sign regressions.

See [VALIDATION.md](studies/VALIDATION.md) and
[SIMULATION_VALIDATION.md](studies/SIMULATION_VALIDATION.md).

## Credibility-aware ranking

The highest numerical point estimate is reported separately from the strongest
supported result. Credibility considers:

- informative SNP count and uncertainty;
- coverage, contamination, and damage metadata;
- duplicate-library and close-relative checks;
- transversion-only estimates;
- alternate African baseline and reference choices;
- leave-one-chromosome-out and block-influence behavior;
- bootstrap or subsampling stability;
- segment evidence when recent ancestry is claimed.

Contradictory sensitivity results remain visible. A sensitivity failure is not
silently averaged away.

## Panel limitations

AADR Human Origins and 1240K panels are ascertained SNP sets. They are powerful
for population comparisons but are not callable whole genomes. Consequences
include:

- panel-relative site distributions;
- limited individual resolution at low coverage;
- no defensible general haplotype caller from pseudo-haploid capture data;
- data-limited X-chromosome analysis when outgroup X genotypes are absent;
- the need for BAM/CRAM follow-up for read-level authentication and long
  introgressed segments.

## Evidence language

- Use **validated** for a method or control that passed explicit checks.
- Use **supported summary** for stable descriptive findings.
- Use **exploratory candidate** for outlier or locus-level follow-up.
- Use **inconclusive/data-limited** when the panel cannot identify the quantity.
- Do not use **discovery** for an unreplicated pipeline outlier.
