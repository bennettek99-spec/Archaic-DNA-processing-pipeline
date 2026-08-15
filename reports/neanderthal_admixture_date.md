# Exploratory date of shared Neanderthal admixture

## Result

The repaired pipeline independently places the main Neanderthal admixture event
carried by the Ranis/Zlaty kun population at approximately **48.5 thousand
years before present (ka BP)**.

This 48.5 ka value is the midpoint of the two primary point estimates, not a
formal pooled estimate:

| Ancient genome | Generations before sample | Calendar point estimate | Propagated 95% interval | Informative target SNPs | Fit R2 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Ranis13 (`RNI013.SG`) | 55.2 +/- 12.7 | 47.5 ka BP | 43.3-51.7 ka BP | 4,830 | 0.380 |
| Zlaty kun (`ZKU002.SG`) | 122.9 +/- 26.3 | 49.4 ka BP | 45.0-53.9 ka BP | 8,271 | 0.480 |

Because both individuals represent the same early population and share broad
age uncertainty, treating them as statistically independent would overstate
precision. A conservative union of the two propagated intervals is
approximately **43.3-53.9 ka BP**. The result is therefore reported as a
**supported exploratory estimate centered near 48.5 ka BP**, not as a new
high-precision date.

The calculation was not calibrated to reproduce a published date. Its result is
consistent with the 45-49 ka estimate from the high-coverage Ranis/Zlaty kun
study and the 50.5-43.5 ka extended gene-flow interval inferred from a broader
ancient/present-day segment analysis.

## What was calculated

At 8,419 AADR SNPs where Altai carries a chimp-polarized derived allele, at
least 80 of 142 African references are called, African derived frequency is at
most 1%, and the allele occurs in CEU/CHB controls, the pipeline measured
genotype covariance across genetic distance. It fit:

```text
C(d) = A exp(-lambda d) + c
```

from 0.02 to 10 cM. `lambda` estimates generations between Neanderthal gene
flow and the sampled human. Leave-one-autosome-out estimates supplied the
standard error. Calendar intervals propagated the genetic uncertainty, each
sample's AADR age uncertainty, and a uniform 25-33 year generation interval
(29-year point conversion).

## Replications and controls

| Genome | Generations +/- SE | Calendar point | Evidence assessment |
| --- | ---: | ---: | --- |
| Ust-Ishim | 88.7 +/- 56.1 | 46.9 ka BP | High uncertainty; generation CI reaches zero |
| Tianyuan | 178.1 +/- 154.8 | 44.7 ka BP | High uncertainty |
| Kostenki14 damage-restricted | 166.8 +/- 213.4 | 42.9 ka BP | High uncertainty |
| Bacho Kiro F6-620 | 69.6 +/- 31.5 | 45.3 ka BP | Separate recent-admixture control; not pooled |
| Oase1 damage-restricted | Not estimable | Not estimable | Only 481 target SNPs; no bins met the pair-count floor |

The later ancient samples are farther from the event and their 1240K curves are
less precise. Oase1's failure is scientifically informative: the array call set
cannot replace the repository's separate BAM/hmmix segment workflow for dating
its recent Neanderthal ancestor.

## Sensitivity results

The two primary point estimates remain near 48-50 ka under the informative
robustness settings:

- Strict African fixed-ancestral ascertainment: 47.55 and 49.12 ka BP.
- Pair-count-weighted fit: 47.51 and 49.19 ka BP.
- Coarser 0.02 cM bins: 47.51 and 49.44 ka BP.
- Excluding distances below 0.05 cM: 47.49 and 49.50 ka BP.
- Vindija-confirmed markers: Zlaty kun 48.87 ka BP; Ranis becomes data-limited.

Two sensitivity failures are retained:

- Transversion-only data leave just 1,187 and 1,723 target sites, so neither
  primary estimate is usable.
- Restricting the fit to 1 cM cannot identify the affine background and produces
  unstable chromosome-jackknife estimates. The 10 cM main fit is therefore
  required for this sparse panel.

See `results/admixture_dating/sensitivity_summary.tsv` for the complete audit.

## Interpretation boundary

The estimate dates ancestry-block breakdown, not the deeper divergence of
Neanderthal and modern-human lineages. It should be interpreted as the average
time of the shared introgression episode under a single-exponential model.
Extended or repeated gene flow can make this an intermediate effective date.

Remaining limitations include 1240K SNP ascertainment, no fine-scale map-error
correction, broad fossil-age uncertainty, non-independence of the primary
individuals, and low power in damage-resistant marker subsets. Whole-genome
segment analysis remains the preferred route for a definitive replication.

## Reproduction

```bash
python -m archaic.admixture_dating \
  --config configs/neanderthal_dating.yaml \
  --out results/admixture_dating
```

Primary references:

- Moorjani et al. 2016, https://doi.org/10.1073/pnas.1514696113
- Sumer et al. 2025, https://doi.org/10.1038/s41586-024-08420-x
- Iasi et al. 2024, https://doi.org/10.1126/science.adq3010
