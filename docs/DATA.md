# Data and artifact policy

The MIT license applies to this repository's source code and documentation. It
does not relicense AADR, ENA, or other upstream genomic data. Users are
responsible for the source dataset's access terms and citation requirements.

## AADR

The Allen Ancient DNA Resource genotype panels are large and are not tracked in
this repository. Obtain the relevant AADR release from the Reich Lab / Harvard
Dataverse, follow its terms, and cite Mallick et al. (2024) together with the
primary publications for samples used in an analysis.

The checked-in `config.yaml` contains a placeholder only. Put machine-local
paths in the ignored `config.local.yaml`, `ARCHAIC_CONFIG`, or
`ARCHAIC_AADR_DIR`.

Never commit AADR `.geno`, `.snp`, `.ind`, or `.anno` files.

## Intentionally tracked BAM exception

Five published, duplicate-removed, L35/MQ30 Sima de los Huesos BAMs are retained
under `data/sima_de_los_huesos_prjeb10597/`. They are inputs to the deliberately
exploratory adapter `tools/sima_de_los_huesos_scan.py` and are not part of the
main AADR Phase 2-9 workflow.

- ENA study: `PRJEB10597`
- Publication: Meyer et al. (2016), *Nature* 531:504-507
- Adapter restrictions: mapping quality >=30, read length >=35 bp, base quality
  >=30, deterministic pseudo-haploid calls at AADR 1240K positions
- Scientific boundary: the approximately 430 ka samples have extremely sparse
  nuclear coverage. Their outputs remain exploratory and below the normal AADR
  evidence floor.

| Sample | ENA run | File | Bytes | Publisher MD5 | Repository SHA-256 |
| --- | --- | --- | ---: | --- | --- |
| femurXIII | ERR995357 | `femurXIII.L35MQ30.bam` | 13,779,378 | `09d66d541a63b4e93bee59ad505770dd` | `bd5d4832c32fe4b649073186f490d437d37cd2a2ae07589091e9faecfd48be07` |
| femur fragment | ERR995361 | `femur_fragment.L35MQ30.bam` | 15,517,579 | `f21cd599ad206caeda0b24be5f180a8c` | `6ff13a3d9f35958d32557d8243d86e9377b70ea17f59a82d532ba56b178cd8ee` |
| incisor | ERR995358 | `incisor.L35MQ30.bam` | 55,948,406 | `660f1ab7f5ced32eb32169f0b71d6492` | `82f0e5c49e7c9b8123a644b6306df2211f5b4bbdf1f2dfbaa835ed3b322eb040` |
| molar | ERR995360 | `molar.L35MQ30.bam` | 5,155,276 | `9ecb7d4c628f812d62ba15a72df03f56` | `3f45904eaf1f4ef6bc8752c5ec668e1fe1866e662a631017a21c757b0d29a1f5` |
| scapula | ERR995359 | `scapula.L35MQ30.bam` | 8,325,626 | `0b70f2c1597e2d461094786e40b89c39` | `8a95d233b0e0002cd871c947242e0d3f6edd643a0b826dc4db8fbe16a63f74f4` |

The files stay in ordinary Git by project decision. Do not add further BAM/CRAM
inputs without documenting provenance, checksums, size, licensing, and why a
fetch recipe is insufficient.

## Generated results

The repository retains selected result tables, figures, HTML reports, and PDFs
when they provide compact, reviewable evidence for a documented analysis.

Every retained result should have:

1. a command or entry point that explains how it was produced;
2. input release/accession information;
3. a methods or provenance record;
4. an interpretation label: validated, supported summary, exploratory, or
   inconclusive/data-limited;
5. no personal genotype data or machine-local paths.

Temporary files, local PDF render workspaces, and partial reruns do not belong in
commits. Stage generated outputs explicitly rather than using `git add -A`.

### Which pipeline tables are tracked

The Phase 2-9 workflow writes a chain of per-sample tables, and only the end of
that chain is evidence. The **phase-4 analysis tables** stay tracked, because
every published study reads them:

| Tracked table | Read by |
| --- | --- |
| `results/phase4_1240k_global_analysis.csv` | global survey, Oceania transect, Neanderthal source contrast, highest-archaic scan |
| `results/phase4_1240k_analysis.csv` | Eurasian >5% survey |
| `results/phase2_1240k_global_excluded.csv` | highest-archaic scan |
| `results/highest_archaic/all_sample_archaic_estimates.tsv` | highest-archaic scan record |

The **upstream intermediates** that exist only to produce those tables are
ignored rather than tracked: phase-2 metadata and sample manifests, phase-3
estimate tables, the phase-5 PCA table, phase-6 residuals, and the
highest-archaic outlier dump. They are reproduced by
`python scripts/run_pipeline.py --panel 1240k` on a machine with the AADR
release available, and nothing downstream of phase 4 consumes them.
