# Skov S5 segment-structure audit: the dating failure was an artifact

## Bottom line

The earlier `not estimable` verdict was **not** a property of the Papuan archaic
tract-length distribution. It was caused by a misreading of what the Skov et al.
2018 supplementary table S5 contains.

S5 is a **complete two-state segmentation of every decoded genome**, not a list
of archaic tracts. The previous analyses treated all 300,000+ exported intervals
as archaic tracts and filtered them on `MeanProb`, which is the posterior of
whichever state was decoded and therefore cannot fall below 0.5. Roughly 40% of
the "Denisovan tracts" that were being dated were **modern-human background
segments**, with a median length of 920 kb and a maximum of 20 Mb.

Once the hidden state is recovered and only archaic-state segments are kept, the
distribution stops being pathological:

| Property | Published pipeline (contaminated) | State-aware (corrected) |
|---|---:|---:|
| Denisovan-affinity tracts | 78,286 | 54,411 |
| Effective decay at 0.02 cM | 83.9 generations | 567.7 generations |
| Effective decay at 0.20 cM | 58.6 generations | 560.5 generations |
| Threshold spread over 0.02–0.20 cM | **1.43x** | **1.03x** |
| KS statistic at 0.05 cM | 0.2917 | 0.0204 |
| Subsampled rejection rate (nominal 0.05) | **1.000** | **0.025** |
| Per-person correlation with published S4 | rho = **−0.431** (wrong sign) | rho = **+0.479** (correct sign) |

The corrected distribution is statistically indistinguishable from a single
exponential. The contaminated one was rejected in every one of 200 subsamples.

After additionally correcting the length-dependent affinity labelling bias, the
estimate is **618.0 generations** (bootstrap 612.5–623.7), or 17.9 kya at 29
years per generation.

**That number is not an admixture date and must not be reported as one.** A
third bias remains: posterior decoding returns runs 1.56x longer than the HMM's
own fitted parameter implies, which accounts for essentially the whole remaining
gap. It cannot be removed from the exported segments without circularity. See
[The bias that remains](#the-bias-that-remains-and-why-the-answer-is-still-not-a-date).

## Evidence that S5 is a genome tiling

All figures are for the 89 Papuan `HMM` individuals, autosomes only
(304,080 segments).

| Test | Result | Expected if S5 listed archaic tracts |
|---|---:|---|
| Median decoded length per individual | 2.642 Gb | ~0.15 Gb |
| Segments per individual | 3,425 | similar |
| Neighbouring segments abutting exactly | 95.03% | near 0% |
| Median gap between neighbours | 0 bp | tens of kb |
| Coordinates off the 1-kb grid | 0 | n/a |
| Minimum `MeanProb` | 0.500000847 | can be any value |
| Genome length passing `MeanProb >= 0.8` | 99.15% | a few percent |

The `MeanProb` floor at 0.5 is the signature of a decoded-state posterior. The
`MeanProb >= 0.8` filter used previously removed 0.85% of decoded genome length,
so it never separated archaic from modern-human segments.

A second, smaller convention error was also found and fixed: the published `end`
column sits one 1-kb window past the segment, so `end - start` exceeds the
published `length` by 1,000 bp in 99.97% of rows. The audit rebuilds the
half-open interval as `[start, start + length)`.

## Recovering the hidden state

The Skov HMM emits the count of outgroup-private variants per 1-kb window with a
low rate in the modern-human state and a high rate in the archaic state. Fitting
`snps ~ Poisson(rate_state x length_kb)` as a two-component mixture recovers the
state with no hand-chosen cut-off:

- modern-human rate: **0.0256** private variants per kb
- archaic rate: **0.2245** private variants per kb
- contrast: **8.8-fold**; only **2.03%** of segments have a state posterior
  between 0.05 and 0.95
- EM converged in 7 iterations

Three independent checks confirm the assignment. None of them was used to fit it:

1. **HMM path alternation.** 98.61% of the 287,100 abutting neighbour pairs
   receive opposite states, as an alternating two-state path requires.
2. **Archaic ancestry fraction.** The archaic state covers **6.82%** of the
   genome (per-individual range 5.99%–7.56%), matching published Papuan
   Neanderthal-plus-Denisovan ancestry.
3. **Length scale.** Median archaic segment 59 kb versus 675 kb for the
   modern-human state.

## What the contamination did

The affinity rule compares raw counts of shared derived alleles. A 1 Mb
modern-human segment accumulates a handful of such alleles by incomplete lineage
sorting alone, and whenever the Denisova count happened to exceed the Vindija
count the segment was labelled Denisovan. That mislabelled **30,873**
modern-human segments as `denisovan_affinity_strict`, with a median length of
920 kb against 88 kb for genuine archaic-state segments in the same class.

Pooling two length distributions that differ tenfold in scale produces exactly
the pathology the earlier work documented:

- a non-exponential shape with a heavy tail;
- an effective decay rate that moves with the detection threshold;
- robustness to leave-one-chromosome-out and leave-one-individual-out, because
  the contamination is uniform across every individual and chromosome;
- indifference to removing published adaptive-introgression candidates;
- zero compatible models in the M1–M10 calibration, because the simulations
  generated pure archaic tracts while the real data was a mixture.

Every one of those observations is explained. The earlier work correctly ruled
out one chromosome, one individual, the map conversion, and selection; the cause
was upstream of all four.

## Corrected estimate and its status

State-aware Denisovan-affinity segments, GRCh37 map-aware, minimum 0.05 cM:

- n = 41,567 tracts across 89 individuals
- effective decay: **569.7 generations**
- individual-level bootstrap 95% interval (400 replicates): **563.3–576.2**
- KS statistic 0.0204; subsampled rejection rate 0.025 against a nominal 0.05
- a left-truncated two-exponential mixture puts **98.8%** of the weight on a
  single component, so there is no evidence of a second pulse in these segments

The estimate is interior to the optimizer bounds, threshold-stable, and
reproducible. Those were the three failures that produced the earlier
`not estimable` verdict, and all three are resolved.

## Residual bias

**569.7 generations is 16.5 kya at 29 years per generation. That is not a
credible Denisovan admixture date**, and it is roughly half the published S4 HMM
parameter (median 1,019.8 generations). The gap is not mysterious; it has a
measured cause.

Affinity labelling is strongly length-dependent, because short segments carry too
few private variants to show any archaic allele sharing:

| Segment-length decile | Median length (cM) | Median private SNPs | Fraction with any archaic sharing |
|---|---:|---:|---:|
| 1 (shortest) | 0.0026 | 5 | 21.3% |
| 5 | 0.0934 | 16 | 84.2% |
| 10 (longest) | 0.4456 | 42 | 97.4% |

Selecting Denisovan-affinity segments therefore discards short tracts
preferentially, lengthens the retained distribution, and biases the decay rate
towards the present. The unlabelled archaic remainder decays at **1,062
generations**, close to the published S4 parameter of 1,020, which is what a
short-segment-enriched subset should look like.

That 1,062 is **not** an alternative estimate. The unlabelled remainder fails the
same tests the corrected set passes (threshold spread 1.75x, subsampled rejection
rate 0.465), because it is enriched for HMM fragments and false positives. It
brackets the bias direction; it does not measure it.

Two controls confirm the diagnosis rather than a lucky subsetting:

| Analysis set | Threshold spread | Rejection rate |
|---|---:|---:|
| All archaic-state segments | 1.14x | 0.070 |
| State-aware Denisovan strict | 1.02x | 0.025 |
| State-aware Neanderthal affinity | 1.04x | 0.060 |
| State-aware unresolved affinity | 1.75x | 0.465 |
| Modern-human state, Denisovan-labelled | 1.00x | 0.765 |

The Neanderthal-affinity set behaves identically to the Denisovan one, as it
should if the fix is structural rather than tuned to one label. The
modern-human segments that the old pipeline labelled Denisovan are rejected on
their own, confirming they were a distinct population of intervals.

## Correcting the labelling bias

The bias is correctable, because the selection function is directly measurable on
the complete archaic state. For lengths above `T` drawn from an exponential
thinned by `c(l)`, the log-likelihood is

```
log L(lam) = -lam * sum(x) - n * log( integral_0^inf exp(-lam x) c(x + T) dx )
```

up to a constant, with `x = l - T` in Morgans. With `c` constant it collapses to
the naive estimator, so the correction is a strict generalisation of the
uncorrected fit.

**The estimator is validated against known truth.** On simulated data with a
saturating selection curve matching the real one, it recovers rates it was never
told, to within 5%:

| True rate | Uncorrected fit | Corrected fit |
|---:|---:|---:|
| 600 generations | 490 (−18%) | 600 (0%) |
| 1,000 generations | 800 (−20%) | 1,000 (0%) |
| 1,800 generations | 1,430 (−21%) | 1,790 (−0.6%) |

Applied to the real data:

| Minimum length | n | Uncorrected | Corrected | Unselected archaic | Ratio |
|---:|---:|---:|---:|---:|---:|
| 0.02 cM | 49,123 | 567.7 | 644.1 | 658.4 | 0.978 |
| 0.05 cM | 41,567 | 569.7 | **618.0** | 621.4 | 0.995 |
| 0.10 cM | 31,532 | 575.2 | 607.4 | 609.9 | 0.996 |
| 0.20 cM | 17,259 | 560.5 | 575.3 | 575.6 | 0.999 |

Headline: **618.0 generations**, individual-level bootstrap 612.5–623.7.

The decisive check is the last column. Reweighting the labelled subset
reproduces the unselected archaic distribution it was drawn from, to within
0.4–2.2%, at every threshold. Nothing was tuned to achieve that; it is what a
correct selection model must do, and it confirms the correction is removing the
labelling bias rather than manufacturing a number.

## The bias that remains, and why the answer is still not a date

618 generations is **17.9 kya** at 29 years per generation. That is still not a
credible Denisovan admixture date, and the reason is now measured rather than
suspected.

The Skov HMM fits a per-individual admixture-time parameter that sets the
geometric prior on archaic run length. Posterior decoding does not reproduce that
prior — it bridges weak evidence and merges runs. Comparing the two on the same
89 individuals measures the inflation directly, with no simulation:

| Quantity | Value |
|---|---:|
| Median fitted S4 parameter | 1,019.8 generations |
| Median decoded archaic decay | 655.3 generations |
| Ratio decoded / fitted | 0.640 (IQR 0.612–0.670) |
| Implied length inflation | **1.56x** |

Decoded runs are 1.56 times longer than the model's own admixture parameter
implies. That accounts for essentially the whole remaining gap: 618 x 1.56 = 966,
against a fitted 1,020.

**That multiplication must not be performed.** Rescaling by the S4 parameter
recovers the S4 parameter by construction and demonstrates nothing. The inflation
is a property of the decoder and can only be removed by simulating genotypes,
running the actual caller, and measuring what it returns for a known truth.

And even then the target would be 1,020 generations, or 29.6 kya — the value
Skov's own estimator produces, which is itself more recent than the 45–55 kya
consensus for Denisovan admixture. Closing that second gap is a question about
the estimator, not about this dataset.

## Summary of the three biases

| Bias | Size | Status |
|---|---:|---|
| Modern-human segments dated as archaic tracts | ~7x, and destroyed the fit entirely | **Fixed** |
| Length-dependent affinity labelling | 1.08x at 0.05 cM | **Fixed and validated** |
| Posterior-decoding run inflation | 1.56x | Measured, not correctable here |

## What still gates a published date

1. **Removed.** The length-dependent affinity ascertainment is corrected, with
   the estimator validated against known rates and against the unselected set.
2. Caller-aware simulation recovery: simulate genotypes, run the actual HMM, and
   confirm the observed distribution lies inside calibrated envelopes. This is
   now the sole remaining blocker and it is the only way to remove the 1.56x
   decoder inflation without circularity. The M1–M10 calibration must also be
   rerun against archaic-state segments, since it was previously compared
   against contaminated data.
3. Independent population-level replication (HGDP Papuans, hg38).

## Interpretation rules

Carried forward, with two corrections:

- **`MeanProb` is the posterior of the decoded state, not an archaic posterior.**
  It is bounded below by 0.5 and must never be used to select archaic segments.
  This supersedes the earlier rule that called it a "generic archaic-state
  posterior".
- **S5 rows are decoded segments, not tracts.** Roughly half are modern-human
  background.
- "Denisovan affinity" remains a relative sharing rule, not proof of a unique
  Denisovan donor. Its length bias is now corrected; its donor semantics are not.
- Effective decay values remain diagnostics, not event dates.
- Do not convert 618 generations into a biological admixture time, and do not
  multiply it by the 1.56x decoder inflation to reach one.
- Do not infer Denisovan survival or extinction timing.
- Raw and genotype data remain uncommitted.

## Reproduction

```bash
python -m archaic_admixture_dating.segment_structure \
  --s5 <path>/skov_2018_s5_segments.tsv \
  --s4 <path>/skov_2018_s4_parameters.xlsx \
  --genetic-map-dir <path>/grch37_hapmap_phase2 \
  --output archaic_admixture_dating/outputs/papuan_s5_segment_structure
```

No downloads. Runs in a few minutes on a laptop. Inputs are fingerprinted by
name, size, and SHA-256 in `provenance.json`; no absolute paths are recorded.

Source study: Skov et al. 2018, PLoS Genetics,
<https://doi.org/10.1371/journal.pgen.1007641>.
