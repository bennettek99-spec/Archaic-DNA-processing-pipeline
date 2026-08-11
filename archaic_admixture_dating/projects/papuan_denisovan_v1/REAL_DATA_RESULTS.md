# Real-data result: Skov et al. Papuan supplements

Status: **completed, with S5 event dating rejected as not estimable**

This run used the public supplementary files from Skov et al. (2018),
*Detecting archaic introgression using an unadmixed outgroup*:

- [Article](https://doi.org/10.1371/journal.pgen.1007641)
- S4 HMM parameter workbook, SHA-256
  `80b3ceb73d05ae3ffc3388380ab702213c4bfc011a79405edb3a0773ec4df8bc`
- S5 archaic-segment table, SHA-256
  `2a1e5a112208d97a29bd04d7b2ec04e3ce621103d519f383cf8b88143b63500f`

The files are public research supplements. They are not committed to Git.

## S5 interval-length reanalysis

The adapter retained Papuan HMM rows and treated `MeanProb` as a generic
archaic-state posterior, not a Denisovan-specific posterior. Physical lengths
were converted with the study's documented average recombination rate of
`1.2e-8` per base pair per generation. Denisovan affinity means greater sharing
with Denisova than Vindija; the strict sensitivity additionally requires
greater sharing with Denisova than Altai.

| Analysis | Retained tracts | Individuals | Median length | Verdict |
|---|---:|---:|---:|---|
| Broad Denisovan affinity | 78,937 | 89 | 0.2328 cM | Not estimable |
| Strict Denisovan affinity | 74,668 | 89 | 0.2274 cM | Not estimable |

Both analyses failed closed:

- the single-pulse fit reached its 100-generation lower bound and failed the
  exponential goodness-of-fit test;
- the younger component of both the two-pulse and continuous-flow fits reached
  the 50-generation lower bound;
- the best-BIC two-pulse model had poor parameter recovery and 0% model-family
  classification accuracy in the configured calibration.

The numerical lower-bound dates are diagnostic optimizer outputs, not
biological event estimates. In particular, they are not evidence for recent
Denisovan admixture or Denisovan survival.

## Matched S4 published-parameter benchmark

S4 contains HMM transition-parameter estimates for the exact same 89 Papuan
sample IDs found in the S5 HMM table. Filtering the `Human population
parameters` sheet to `Outgroup == Whole~world` gives:

| Quantity | Result |
|---|---:|
| Individuals | 89 |
| Median admixture-time parameter | 1,019.77 generations |
| Cross-individual range | 888.32-1,191.36 generations |
| Median at 29 years/generation | 29.57 kya |
| Median at 27-30 years/generation | 27.53-30.59 kya |
| 10,000-replicate bootstrap interval for the cross-individual median | 28.87-30.24 kya |

The bootstrap interval describes uncertainty in the median of the published
per-individual HMM estimates. It is not the full uncertainty interval for a
single admixture event and is not an independent raw-genome estimate.

## Interpretation

The defensible real-data result is:

1. the current independent-tract exponential models cannot date the exported
   S5 intervals;
2. the matched published S4 HMM parameters center near 29.6 kya;
3. neither result establishes when Denisovans went extinct or proves direct
   contact with a surviving Denisovan population at that date.

Independent event dating requires caller-aware transition likelihoods or
independently calibrated tract calls paired with a genomic recombination map.

## Validation

- 24 module tests passed.
- `pyflakes` passed.
- `git diff --check` passed.
- Both source files passed size and pinned SHA-256 verification.
- Broad and strict workflows reached `complete` checkpoint state.

Generated run directories are intentionally ignored:

- `archaic_admixture_dating/outputs/papuan_skov_real_v1_broad`
- `archaic_admixture_dating/outputs/papuan_skov_real_v1_strict`

## GRCh37 HapMap genetic-map sensitivity

A bounded Windows sensitivity run replaced the constant `1.2e-8` conversion
with linear interpolation on the chromosome-specific HapMap Phase II map
lifted to GRCh37 by Adam Auton. The 22 compressed autosomal maps total
36,704,757 bytes (35.0 MiB). Source commit:
`42e69e59f1d3a31725378e61d8f8e82a51e88d63`.

The run used the real S5 supplement, one thread, a 2 GiB configured ceiling,
20 bootstrap replicates, five complex-model bootstrap replicates, and two
calibration replicates for each M1-M10 scenario. It is a low-compute
sensitivity analysis; the small calibration is suitable for rejecting an
unstable result, not for positive model selection.

| Analysis | Retained tracts | Individuals | Median length | Verdict |
|---|---:|---:|---:|---|
| Broad Denisovan affinity | 75,360 | 89 | 0.3018 cM | Not estimable |
| Strict Denisovan affinity | 71,320 | 89 | 0.2911 cM | Not estimable |

Before source-class and confidence filtering, 1,528 autosomal Papuan HMM
intervals extended beyond the covered map ranges and were excluded rather
than extrapolated. Another 30 intervals spanned zero cM on the map and failed
the positive-length check.

In both analyses the single-pulse fit reached its 100-generation lower bound
and failed the exponential goodness-of-fit test (`KS p = 0`). The favored
two-pulse optimizer solutions contained a younger component at the
50-generation lower bound and had 0% model-family classification accuracy in
the bounded calibration. Their numerical dates are rejected diagnostics, not
admixture-event estimates.

This sensitivity changes the inferred genetic-length distribution but does
not rescue the independent-tract dating model. The defensible conclusion
remains `not-estimable`. Compact results are in
`MAP_AWARE_LIGHTWEIGHT_SUMMARY.tsv`; hashes for the uncommitted reference maps
are in `GRCH37_MAP_MANIFEST.tsv`.

Generated reports remain ignored:

- `archaic_admixture_dating/outputs/papuan_s5_map_grch37_light_broad/report/report.html`
- `archaic_admixture_dating/outputs/papuan_s5_map_grch37_light_strict/report/report.html`

## Failure anatomy and S4-S5 bridge

The full map-aware follow-up retained 75,360 broad Denisovan-affinity tracts
from the same 89 individuals. The distribution is not close to a stable single
exponential:

| Diagnostic | Result |
|---|---:|
| Effective decay at 0.02 cM | 83.90 generations |
| Effective decay at 0.50 cM | 42.55 generations |
| Effective decay at 1.00 cM | 36.63 generations |
| KS statistic at 0.02 cM | 0.2983 |
| KS p-value | 0 |
| Effective decay after removing longest 5% | 127.94 generations |
| Leave-one-chromosome-out range | 78.77-86.60 generations |
| Leave-one-individual-out range | 83.79-84.04 generations |

This pattern localizes the instability to the shape and long tail of the
distribution, not to one person or chromosome. These effective decay values
are unbounded diagnostics, not event dates.

The published Gower et al. (2021) GRCh37/hg19 Melanesian adaptive-
introgression candidate intervals overlapped 1,057 broad-affinity tracts.
Excluding them changed the effective decay only from 83.90 to 84.19
generations and left the KS p-value at zero. Candidate selected loci therefore
do not explain the model failure. The intervals are candidates, not a
definitive selected-locus catalogue.

All 89 `Whole~world` S4 individuals matched their map-aware S5 features. The
strongest cross-person association with the published S4 admixture-time
parameter was for Neanderthal-affinity total tract length (Spearman
`rho = 0.615`, 95% individual bootstrap interval `0.473-0.730`). The broad
Denisovan-affinity median-length association was weak within every source
cohort (`rho = 0.054`, `0.007`, and `0.179`). These are concordance checks on
related outputs from the same HMM study, not independent validation.

## Bounded observation-process calibration

The calibration generated 400 summary replicates: 20 replicates for each
M1-M10 approximation under the configured caller error and under a heavier
stress case (20% random false negatives and 25% length noise). Each replicate
requested 5,000 tracts.

Internal leave-one-out model classification averaged 92.5% under configured
error and 88.0% under stress. Single-pulse M1-M4 median relative recovery error
was 0.8-1.4% under configured error and 2.3-3.4% under stress. The calibration
can therefore recover and distinguish its own approximations.

Nevertheless, the real broad-affinity profile fell outside all nine 95%
feature envelopes for every model and error scenario. The closest tested
combination was stress-case M8, but it matched `0/9` features and was rejected.
No tested M1-M10 model is compatible, and the result remains **not estimable**.
This is posterior-predictive rejection, not a posterior probability over
demographic histories.

## External HMMix and raw chr22 gates

The official HMMix hg38 HGDP premade segment file is 399,257,886 bytes and
therefore below the 500 MB transfer cap. A 1 MiB HTTP range request verified
the schema. However, its Papuan samples use HGDP identifiers (for example,
`HGDP00540`), while none of the 89 S4/S5 individuals has an HGDP identifier.
The exact matched-person count is zero, so the remaining 398 MB was not
downloaded and no cross-build, unmatched population comparison was substituted.

A local 205,612,353-byte chromosome-22 1000 Genomes VCF was also inventoried.
It contains no Papuan population. WSL, bcftools, and vcftools are unavailable,
and the official hg19 HMMix ancestral, outgroup, and reference inputs are each
larger than the 500 MB cap. The raw caller pilot was therefore not run: there
is no authorized Papuan chr22 target in the local input, and the required
resource/tool gate does not pass.

Compact evidence is in `FAILURE_ANATOMY_SUMMARY.tsv`,
`OBSERVATION_CALIBRATION_SUMMARY.tsv`, and `EXTERNAL_PREFLIGHT.tsv`. Generated
HTML reports remain ignored under `archaic_admixture_dating/outputs/`.
