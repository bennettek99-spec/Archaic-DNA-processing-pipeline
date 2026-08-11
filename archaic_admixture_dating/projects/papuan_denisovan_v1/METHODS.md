# Methods

## Tract input and QC

Imported tracts are normalized to the documented schema and validated before
use. Coordinates must be nonnegative with end greater than start. Genetic
length must be positive. Confidence, source-class probability, and callability
must fall between zero and one when present. Invalid rows are retained in an
exclusion table with reasons.

QC applies configurable minimum length, confidence, and callability thresholds.
Low-mappability, centromere, telomere, and selected-locus BED masks can be
applied without silently merging caller outputs. Overlapping tracts are
reported with provenance.

## Single pulse

For a pulse `t` generations ago, excess tract lengths above a detection limit
`Lmin`, expressed in Morgans, are modeled as exponential with rate `t`.
Conditioning on detection yields:

`log L(t) = n log(t) - t Σ(Li - Lmin)`.

Dates in years multiply generations by a configured generation time.

## Two pulses

The two-pulse density is an ordered mixture of two truncated exponentials. The
older date is constrained to exceed the younger date, preventing label
switching. The workflow warns when a weight approaches zero, dates collapse,
optimization fails, or tract count is inadequate.

## Prolonged flow

Prolonged flow is approximated by a uniform mixture of exponential rates over
an estimated younger-to-older interval. This is an interpretable V1
approximation rather than a claim of constant biological migration.

## Uncertainty

Chromosome-block and sample bootstrap are used for linked data. Generation
time, tract threshold, confidence, map, longest-tract exclusion, selected-locus
exclusion, caller error, and bottleneck assumptions are sensitivity axes.

## Simulations and model comparison

M1-M10 are configuration-defined. Deterministic tract-level simulations make
smoke and laptop calibration possible. Optional coalescent calibration is
required before strong real-data claims. Models are compared with likelihood,
AIC, BIC, stability, recovery, and posterior-predictive diagnostics where
available. A more complex model does not win solely from raw likelihood.
