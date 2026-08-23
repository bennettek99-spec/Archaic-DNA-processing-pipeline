# Neanderthal admixture dating

## Question

This module estimates how many generations separate an ancient modern human
from a Neanderthal admixture event. It does not infer the human-Neanderthal
population split, and it does not date admixture from an f4-ratio alone.

The implementation follows the single-sample statistic of Moorjani et al.
(2016). At SNPs where the Altai Neanderthal carries a chimp-polarized derived
allele and an African reference panel is ancestral, the target genome's
genotype covariance is measured across pairs of SNPs at genetic distance `d`.
The curve

```text
C(d) = A exp(-lambda d) + c
```

is fit with `d` in Morgans. `lambda` is the number of generations from
admixture to the sampled person. Calendar time is calculated as:

```text
event age BP = sample age BP + lambda * generation interval
```

The default point conversion uses 29 years per generation and propagates a
25-33 year range together with chromosome-jackknife and AADR sample-age
uncertainty.

Primary method reference:

- Moorjani P. et al. 2016. *A genetic method for dating ancient genomes
  provides a direct estimate of human generation interval in the last 45,000
  years.* PNAS 113:5652-5657. https://doi.org/10.1073/pnas.1514696113

## AADR adaptation

The published ascertainment requires the African panel to be fixed ancestral.
In AADR's pseudo-haploid representation, one erroneous African call can remove
a SNP across a large panel. The default therefore allows at most 1% derived
frequency across 142 Mbuti/Yoruba/YRI references, requires at least 80 called
African samples, and requires the derived allele in at least one CEU or CHB
control. A strict zero-frequency run is retained as a sensitivity analysis.

The default fit spans 0.02-10 cM. The longer range is necessary for the sparse
1240K target calls to identify the affine background; the 0.02-1 cM sensitivity
is unstable and is reported rather than hidden.

## Prespecified targets

- Primary shared-event targets: `RNI013.SG` and `ZKU002.SG`.
- Early Eurasian replications: `Ust_Ishim.DG`, `Tianyuan.AG.BY.AA`, and
  `Kostenki14_d.AG.BY.AA`.
- Recent secondary-admixture controls: `F6-620.AG.BY.AA` and
  `Oase1_d.AG.BY.AA`.

Oase1 and Bacho Kiro are never pooled into the shared non-African event because
both have evidence of later Neanderthal ancestors. Oase1's sparse
damage-restricted 1240K call set is expected to fail closed; the separate
BAM/hmmix workflow is the appropriate high-resolution analysis.

## Run

Configure `aadr_dir` in ignored `config.local.yaml`, then run:

```bash
archaic-pipeline admixture-date
```

Equivalent module invocation:

```bash
python -m archaic.admixture_dating \
  --config configs/neanderthal_dating.yaml \
  --out results/admixture_dating
```

Useful sensitivity controls:

```bash
# Paper-strict African fixed-ancestral rule
python -m archaic.admixture_dating \
  --sample RNI013.SG --sample ZKU002.SG \
  --max-african-derived-frequency 0 \
  --out results/admixture_dating_sensitivity/strict_african

# Vindija-confirmed markers
python -m archaic.admixture_dating \
  --sample RNI013.SG --sample ZKU002.SG --require-vindija \
  --out results/admixture_dating_sensitivity/vindija_confirmed

# Damage-resistant but lower-power transversions
python -m archaic.admixture_dating \
  --sample RNI013.SG --sample ZKU002.SG --transversions-only \
  --out results/admixture_dating_sensitivity/transversions
```

## Outputs

- `estimates.tsv`: generations, jackknife uncertainty, calendar conversion,
  fit diagnostics, evidence status, and interpretation.
- `covariance_curves.tsv`: observed and fitted covariance by distance bin.
- `chromosome_jackknife.tsv`: leave-one-autosome-out generation estimates.
- `run_manifest.json`: configuration, Git identity, ascertainment counts, and
  a lightweight local-input fingerprint.
- `admixture_date_curves.png`: visual fit audit.

## Interpretation boundary

This is a supported exploratory 1240K estimate, not a definitive redating.
Important limitations are:

- sparse, ascertained 1240K rather than whole-genome data;
- no fine-scale genetic-map error correction;
- shared age uncertainty and relatedness among primary individuals;
- single-exponential pulse fitting, which returns an average date if gene flow
  lasted across multiple generations;
- low-power transversion-only and Vindija-confirmed subsets.

The single-sample statistic is downward biased when applied to present-day
genomes whose admixture event is more than roughly 2,000 generations old. This
module therefore does not present present-day single-genome dates as valid.

## SNP-density limitation (quantified)

`scripts/ns_dating_density.py` makes the 1240K-density limitation quantitative
rather than asserted. It simulates a chromosome with a known Neanderthal pulse
and recombination, ascertains Altai-derived/African-ancestral SNPs exactly as
the estimator does, then fits the covariance curve at whole-genome-like density
and at 1240K-like density (`results/dating_density_sim/density_recovery.csv`).

The short-distance (<2 cM) SNP pairs that pin down the decay rate `lambda`
collapse by roughly **625x** at 1240K density relative to whole-genome density.
At that pair count the recovered generations are unconstrained (jackknife SE of
the same order as the point estimate), while whole-genome density at least
bounds the rate. This is the quantitative form of the "sparse, ascertained 1240K"
limitation above, and the reason the checked-in run is labelled exploratory.

