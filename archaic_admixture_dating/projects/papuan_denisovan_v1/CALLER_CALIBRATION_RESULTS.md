# Caller-aware calibration: decoder inflation is not large enough to explain the gap

Status: **pilot complete, result is negative and fails closed**

## Bottom line

[SEGMENT_STRUCTURE_RESULTS.md](SEGMENT_STRUCTURE_RESULTS.md) ended with one
blocker: the corrected estimate of 618 generations could not be turned into a
date because posterior decoding returns runs 1.56x longer than the HMM's own
fitted parameter implies, and that inflation could only be removed by
simulating genotypes, running the caller, and measuring what it returns for a
known truth. That is now done.

The answer is that **decoder inflation, at the magnitude actually measured on
the real data, is not nearly large enough to close the gap.**

Reconciling the observed Papuan tract distribution with a 45 kya pulse would
require the caller to inflate runs by **2.72x**. The real data's own measured
inflation is **1.56x** — short by a factor of 1.74. Even reaching the Skov S4
value of 29.6 kya would require 2.12x.

The real decoded decay of 655.3 generations lies **below every point in the
calibration sweep**, including a pulse at 14.5 kya. No single Denisovan pulse
anywhere from 14.5 to 51.6 kya, under the published Jacobs et al. (2019)
demography and this observation model, produces a tract distribution that
decays as slowly as the real data does.

The formal inversion returns 315 generations (9.1 kya, 95% CI 2.6–13.3 kya).
**That is an extrapolation past the end of the curve and is not a date.** It is
reported because refusing to report it would hide how far outside the model the
real measurement sits.

## What was built

| Module | What it does |
|---|---|
| `skov_hmm.py` | Two-state Poisson HMM over outgroup-private window counts: Baum-Welch, posterior decoding, run extraction, truncated-exponential decay. numba-JIT'd. |
| `genotype_simulation.py` | Genotypes under the published Jacobs demography from stdpopsim, with the Denisovan pulse time freed. |
| `caller_calibration.py` | Sweep, curve fit, inversion with replicate bootstrap, and the required-inflation diagnostic. |
| `run_calibration.py` | The entry point that produced everything below. |

The published Jacobs model encodes a *young* Denisovan pulse (1027.6
generations, 29.8 kya) as a finding, so running it unmodified would assume the
thing under test. The two published Denisovan pulses are replaced by a single
pulse of the same total proportion (0.04) at the swept time. Every other
parameter — population sizes, split times, Neanderthal pulses, migration rates
— stays at its published value.

This replaces the previous tract-level observation model, which drew true
tracts from migration records and added false negatives and multiplicative
length noise. That error model is length-*preserving*; decoder inflation is
length-*creating*. No setting of those knobs could ever have reproduced it,
which is why the earlier M1–M10 calibration rejected every model at 0/9
features.

## Validation of the caller

On data simulated directly from a two-state Markov chain with known parameters:

| Quantity | Recovered | Truth |
|---|---:|---:|
| Modern-human rate | 0.0256 | 0.0256 |
| Archaic rate | 0.2233 | 0.2245 |
| Archaic fraction | 0.0677 | 0.0670 |
| Admixture parameter | 1537.8 | 1500.0 |

The decoder distortion is present and in the right direction: decoding returned
2,581 runs where the truth had 3,557, with median length 51 kb against a true
38 kb.

The JIT is what makes this feasible. A whole-genome forward-backward pass in
pure Python takes 57.7 s; JIT-compiled it takes 1.1 s, turning a 43-hour
per-replicate cost into about 5 s per individual.

## Observation-process fidelity

The simulation's nuisance parameters — effective variant density and outgroup
panel size — were fitted to the real analysis's *measurable* anchors before any
calibration was run. The pulse time was never tuned.

| Anchor | Simulated | Real | Ratio |
|---|---:|---:|---:|
| Modern-human rate (per kb) | 0.0286 | 0.0256 | 1.12 |
| Archaic rate (per kb) | 0.1990 | 0.2245 | 0.89 |
| Rate contrast | 6.97x | 8.8x | 0.79 |
| Decoded archaic fraction | 0.1029 | 0.0682 | 1.51 |
| Decoded / fitted ratio | 0.763 | 0.643 | 1.19 |

Rates match well; contrast and archaic fraction do not. Both residuals are
discussed under [Threats](#threats-to-the-conclusion).

## The calibration curve

Twelve pulse times, four replicates each, 10 Mb per replicate, 20 Papuan and
100 YRI diploids, variant density scaled to 0.40. Full table in
[CALLER_CALIBRATION_SUMMARY.tsv](CALLER_CALIBRATION_SUMMARY.tsv).

| True pulse (gen) | True (kya) | Fitted | Decoded decay | d/f |
|---:|---:|---:|---:|---:|
| 500 | 14.5 | 1166.7 | **784.9** | 0.677 |
| 700 | 20.3 | 1195.9 | 813.2 | 0.682 |
| 900 | 26.1 | 1309.7 | 937.4 | 0.732 |
| 1000 | 29.0 | 1283.7 | 1005.8 | 0.786 |
| 1250 | 36.2 | 1642.3 | 1331.1 | 0.812 |
| 1400 | 40.6 | 1627.2 | 1284.1 | 0.797 |
| 1550 | 45.0 | 1824.9 | 1470.6 | 0.802 |
| 1780 | 51.6 | 1790.7 | 1474.6 | 0.824 |

`decoded_decay = 0.5725 x true + 475.1`, R² = 0.916 on the point means.

The relationship is real and monotone: the caller does respond to the pulse
time. But the lowest decoded decay anywhere in the sweep is **784.9**, and the
real measurement is **655.3**.

## What would have to be true

| Candidate date | Predicted decoded decay | Extra compression needed | Required total inflation |
|---|---:|---:|---:|
| 29.6 kya (Skov S4 parameter) | 1059.0 | 0.619 | **2.12x** |
| 45 kya | 1362.5 | 0.481 | **2.72x** |
| 50 kya | 1462.1 | 0.448 | **2.92x** |
| 55 kya | 1561.1 | 0.420 | **3.12x** |

Against a simulated inflation of 1.31x and a **real, directly measured
inflation of 1.56x**.

This is the finding. The 1.56x decoder inflation documented in the previous
work is genuine, but it is roughly half of what would be needed to reconcile
the observed distribution with an old pulse. The remaining gap is not decoder
inflation.

## Over-calling does not explain it

The simulation calls more of the genome archaic than the real analysis does
(0.103 vs 0.068), and false positives are short, so they would bias simulated
decay upward and inflate every number in the table above. This was tested
directly by re-decoding the same replicates at increasing posterior thresholds
(pulse at 1550 generations, three seeds):

| Threshold | Archaic fraction | Decoded decay |
|---:|---:|---:|
| 0.50 | 0.0705 | 1707.9 |
| 0.70 | 0.0576 | 1678.4 |
| 0.80 | 0.0511 | 1737.7 |
| 0.90 | 0.0429 | 2034.4 |
| 0.99 | 0.0225 | 4256.6 |

At the threshold where the simulated archaic fraction matches the real 0.0682,
the simulated decay is 1707.9 against the real 655.3. Increasing specificity
moves the decay **further from** the real value, not towards it. The gap is not
an over-calling artifact.

## Threats to the conclusion

1. **Rate contrast is 6.97x against the real 8.8x.** Lower contrast means lower
   signal-to-noise, which should produce *more* decoder distortion, not less —
   so this residual works against the conclusion rather than for it. It cannot
   be tuned away, because contrast is set by the archaic divergence times in
   the published demography.
2. **Replicate scatter is large**: within-point SD of decoded decay is 167.8
   generations on ~130 tracts per replicate. The curve is well determined
   (R² = 0.916) but individual points are not.
3. **Pilot scale.** 10 Mb per replicate and 20 individuals, against 89
   individuals and whole genomes in the real analysis.
4. **Flat recombination rate** of 1.2e-8 in both simulation and length
   conversion. Self-consistent, but the real analysis used the GRCh37 HapMap
   map.
5. **This is a faithful reimplementation of the Skov model, not Skov's code.**
   Any behaviour specific to their implementation is not captured.
6. **Single pulses only.** Two-pulse and continuous-flow histories were not
   swept, and a mixture could in principle produce a slower-decaying
   distribution than any single pulse.
7. **The demography is held at published values.** Papuan effective size and
   bottleneck structure directly set tract-length structure, and a different
   demography would move the curve.
8. **The 655.3 anchor** is inherited from the state-aware reanalysis and
   carries that analysis's own caveats.

Of these, (6) is the most likely route to a reconciling explanation and is the
obvious next experiment: the required-inflation table says a single pulse
cannot do it, which is itself a reason to sweep mixtures.

## Interpretation rules

Carried forward from [SEGMENT_STRUCTURE_RESULTS.md](SEGMENT_STRUCTURE_RESULTS.md),
plus:

- **315 generations / 9.1 kya is not a date.** It is an extrapolation past the
  end of the calibration curve and is reported only to show how far outside the
  swept range the real measurement falls.
- Do not report the required-inflation figures as evidence for a young pulse.
  They measure a shortfall in the forward model, not a property of Denisovans.
- The decoder inflation measured in the earlier work (1.56x) remains correct.
  What is new is that it is insufficient, not that it is wrong.
- Effective decay values remain diagnostics, not event dates.
- Do not infer Denisovan survival or extinction timing.

## Reproduction

```bash
python archaic_admixture_dating/projects/papuan_denisovan_v1/run_calibration.py --output archaic_admixture_dating/outputs/papuan_caller_calibration_v1 --replicates 4
```

Requires the `calibration` extra (`msprime`, `stdpopsim`, `numba`). No
downloads and no controlled-access data. About 10 minutes on a laptop; the run
directory is gitignored, the summary table beside this file is tracked.
