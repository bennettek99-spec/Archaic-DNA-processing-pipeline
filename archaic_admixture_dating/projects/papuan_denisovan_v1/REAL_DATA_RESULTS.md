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
