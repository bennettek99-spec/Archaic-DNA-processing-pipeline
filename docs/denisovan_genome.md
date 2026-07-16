# Denisovan genome module

The main Phase 2–9 workflow excludes archaic reference genomes from its ancient
Eurasian target cohort. `archaic.denisovan_genome` is the reference-aware path
for running a Denisovan genome through the same validated AADR reader and
block-jackknife statistic engine.

## Default real-data run

```powershell
python -m archaic.denisovan_genome `
  --panel 1240k `
  --target Denisova3.DG `
  --output results/denisovan_genome

# Equivalent after an editable install:
archaic-pipeline denisovan-genome --panel 1240k
```

The default target is the high-coverage diploid `Denisova3.DG` call set. The
module also uses:

- `Denisova.SG`, a lower-density call set from the same specimen, strictly as a
  technical replicate;
- `Denisova3_snpAD.DG`, an alternate call set from the same specimen;
- `Denisova11.SG`, the known Neanderthal–Denisovan F1, as a biological control;
- `Denisova25.SG`, an older provisional Denisova Cave genome;
- Altai, Chagyrskaya, and Vindija Neanderthal comparators;
- Mbuti, Yoruba, Papuan, and French population controls.

## Outputs

The output directory contains a self-contained `report.html`, `RESULTS.md`, a
summary PNG, TSV tables for QC, pairwise concordance, f-statistics,
transversion sensitivity, chromosome robustness, and Denisovan-marker sharing,
plus a provenance-rich `run_manifest.json`.

## Interpretation boundary

The module does not report a Denisovan ancestry percentage for a Denisovan
genome. `Denisova3.DG` and `Denisova.SG` are different call sets from the same
biological specimen, so their agreement validates the data path but is not
independent ancestry evidence. The panel also lacks an independent second
high-quality Denisovan calibration genome appropriate for an absolute
Denisovan fraction. Reference-defined marker sharing is therefore labelled as
a lineage fingerprint, never an admixture percentage.
