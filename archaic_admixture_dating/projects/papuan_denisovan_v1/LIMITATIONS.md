# Limitations

- V1 imports established tract calls and does not independently validate every
  caller's posterior calibration.
- Short tracts from old pulses can fall below detection, biasing dates younger.
- Genetic map misspecification changes tract lengths and fitted generations.
- Phasing and genotype error can fragment or join tracts.
- Two pulses and prolonged flow can be statistically non-identifiable.
- Bottlenecks and modern-human population mixing can mimic young admixture.
- Selection can preserve a small number of long tracts that dominate a fit.
- A single Altai Denisovan is an incomplete proxy for donor diversity.
- Reference similarity does not uniquely identify geographic source.
- Tract-level simulation is a workflow/recovery approximation; strong claims
  require caller-aware coalescent calibration on appropriately governed data.
- Relative Denisovan affinity cannot be reported as an absolute ancestry
  percentage without a separately calibrated estimator.
- A young conditional tract date does not establish Denisovan survival to that
  date.

When models cannot be distinguished, sample size is inadequate, callability is
low, recovery is poor, or results depend on a few long tracts, the correct
result is `not distinguishable`, `not estimable`, or
`inconclusive/data-limited`.
