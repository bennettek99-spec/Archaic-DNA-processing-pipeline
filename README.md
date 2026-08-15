# Archaic DNA Processing Pipeline

[![tests](https://github.com/bennettek99-spec/Archaic-DNA-processing-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/bennettek99-spec/Archaic-DNA-processing-pipeline/actions/workflows/ci.yml)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](pyproject.toml)
[![release](https://img.shields.io/github/v/release/bennettek99-spec/Archaic-DNA-processing-pipeline?include_prereleases)](https://github.com/bennettek99-spec/Archaic-DNA-processing-pipeline/releases)

A reproducible Python pipeline for studying Neanderthal and Denisovan-related
ancestry in ancient and present-day human genomes from the Allen Ancient DNA
Resource (AADR). It combines allele-frequency f-statistics, block-jackknife
uncertainty, technical-quality gates, ancestry-aware outlier models, sensitivity
analyses, and portable reports.

The pipeline is designed to separate a large numerical estimate from a credible
biological result. Low coverage, contamination, damage, relatedness, panel
ascertainment, and reference choice are surfaced rather than hidden.

> **Interpretation boundary:** the validated Neanderthal f4-ratio is reported as
> a percentage. `D_Den` is a relative Denisovan-affinity statistic, not a
> Denisovan percentage. The pipeline therefore does not calculate a combined
> Neanderthal-plus-Denisovan percentage.

## The two papers

These are the pipeline's headline outputs: one at the scale of the whole AADR,
one at the scale of a single individual.

### 📄 [Global AADR survey (PDF)](reports/global_archaic_survey/Global_archaic_survey.pdf)

The pipeline applied end to end to 21,109 quality-passing ancient and
present-day AADR genomes across Africa, the Americas, Eurasia, and Oceania:
global Neanderthal estimates, relative Denisovan-affinity controls,
coverage-aware filtering, artifact auditing, and continental comparisons, all
within the interpretation boundary above.

[Paper-style report](reports/global_archaic_survey/PAPER_global.md) ·
[supporting outputs](reports/global_archaic_survey/)

### 📄 [Iron Age Etruscan study (PDF)](reports/Etruscan_paper.pdf)

The same machinery aimed at one question, and declining to answer it. EV16A.SG
from Monteriggioni carries a raw 5.35% Neanderthal point estimate on only 15,994
informative SNPs, a 1.22-9.49% confidence interval. Transversion, alternate
outgroup, reference swap, per-chromosome, block-bootstrap, and local Etruscan
controls are all reported, and the individual still ends up classified **low
confidence** because coverage, read-level QC, and segment evidence do not support
it. That refusal is the result.

[Group-level findings](docs/studies/ETRUSCAN_FINDINGS.md) ·
[dedicated EV16A analysis](results/individual_EV16A/EV16A_dedicated_report.md) ·
[sensitivity table](results/individual_EV16A/top_candidate_sensitivity_tests.tsv)

### Also worth reading

- **[Remote Oceania transect](reports/oceania_transect/PAPER_oceania.md)** —
  Denisovan ancestry is normally measured on present-day genomes, so its arrival
  is inferred rather than observed. Vanuatu is the exception. Across 31 Vanuatu
  ancients in four dated horizons, pooled Denisovan affinity rises from
  indistinguishable from zero at founding to 75% of the present-day Papuan level.
  Ninety-six Guam/Marianas genomes, same expansion but without the Papuan-related
  influx, stay flat. The parallel Neanderthal rise appears in the control too, so
  it is reported as shared measurement drift rather than a result.
- **[Which Neanderthal? Altai vs Vindija](reports/neanderthal_source/PAPER_neanderthal_source.md)** —
  a source contrast across the AADR that returns a null with a stated detection
  limit rather than a finding.

## Quickstart

Five minutes, no data download. The synthetic smoke test drives the real
packed-genotype reader and the shared f-statistic engine, so a green run means
the install works end to end.

```bash
git clone https://github.com/bennettek99-spec/Archaic-DNA-processing-pipeline.git
cd Archaic-DNA-processing-pipeline
python -m venv .venv
```

Activate it — PowerShell `.\.venv\Scripts\Activate.ps1`, or bash
`source .venv/bin/activate` — then:

```bash
python -m pip install --upgrade pip
python -m pip install -e .
archaic-pipeline smoke-test
```

That prints a synthetic introgression estimate and exits 0. You now have five
console entry points:

| Command | What it does |
| --- | --- |
| `archaic-pipeline` | the Phase 2-9 workflow and every study subcommand |
| `archaic-highest` | highest-archaic credibility scan |
| `archaic-highest-segments` | segment support for scan candidates |
| `archaic-denisovan` | Denisovan reference-genome module |
| `archaic-admixture-dating` | Papuan Denisovan admixture dating V1 |

Run `archaic-pipeline --help` to see the subcommands.

**Next:** to work on the code, `python -m pip install -e ".[test,sim]"` then
`python -m pytest -q`. To run real analyses, you need an AADR panel — see
[Run with AADR](#run-with-aadr) below. To understand what the numbers mean before
trusting any of them, read
[methods and interpretation](docs/methods-and-interpretation.md).

## Status

- The core estimator passes seven AADR-based validation gates and is checked
  against published estimates, simulation, and ADMIXTOOLS 2.
- The Phase 2-9 AADR workflow, global survey, highest-archaic credibility scan,
  Etruscan study, Denisovan reference-genome module, focused Oase1 workflow,
  and exploratory single-genome Neanderthal admixture dating are implemented.
- Individual outlier findings remain hypotheses unless the documented coverage,
  damage, contamination, sensitivity, and—where appropriate—segment evidence
  support them.

See [methods and interpretation](docs/methods-and-interpretation.md) for the
statistical definitions and evidence limits.

## Run with AADR

The full workflow requires a locally obtained AADR v66.p1 Human Origins or
1240K panel. AADR genotype files are not redistributed by this repository.

Copy `config.yaml` to the ignored `config.local.yaml`, then set the directory
containing `v66.p1_1240K.{geno,snp,ind,anno}`:

```yaml
aadr_dir: "/path/to/aadr"
```

You can instead set `ARCHAIC_CONFIG=/path/to/config.yaml` or
`ARCHAIC_AADR_DIR=/path/to/aadr`.

```bash
archaic-pipeline validate --panel 1240k
archaic-pipeline all --panel 1240k
```

The orchestrator runs Phases 2-9 and generates
`reports/archaic_report_1240k.html`. The checked-in configuration contains no
developer-specific AADR path and fails with a setup message when none is
provided.

Detailed setup and command examples are in the
[getting-started guide](docs/getting-started.md).

## Scientific scope

The core statistics are:

1. **Neanderthal proportion:** an Altai/Vindija-scaled f4-ratio with a
   block-jackknife confidence interval.
2. **Neanderthal affinity:** `D(X, Mbuti; Altai, Chimp)`, used as a relative
   significance channel.
3. **Differential sharing:** population contrasts such as
   `D(Pop1, Pop2; Altai, Yoruba)`.
4. **Denisovan affinity:** `D(X, Mbuti; Denisova, Chimp)`, reported only as a
   relative statistic and Z-score.

Important limitations:

- A single low-coverage ancient genome cannot support fine-scale percentage
  comparisons at population-level resolution.
- Human Origins and 1240K are ascertained panels, not callable whole genomes.
- A high point estimate is not evidence of recent admixture without technical
  stability and, where relevant, long-segment evidence.
- Published Denisovan percentages may be discussed as external context, but
  they are never presented as outputs calculated by this pipeline.

## Analyses and reports

| Analysis | Entry point | Main documentation or report |
| --- | --- | --- |
| Core validation | `archaic-pipeline validate --panel 1240k` | [VALIDATION.md](docs/studies/VALIDATION.md) |
| Full Phase 2-9 workflow | `archaic-pipeline all --panel 1240k` | `reports/archaic_report_1240k.html` |
| Highest-archaic credibility scan | `archaic-pipeline highest-archaic` | [module guide](docs/highest_archaic.md) |
| Denisovan reference genome | `archaic-pipeline denisovan-genome --panel 1240k` | [module guide](docs/denisovan_genome.md) |
| Papuan Denisovan admixture dating V1 | `archaic-admixture-dating run-all --profile smoke --resume` | [module guide](archaic_admixture_dating/README.md) |
| Remote Oceania archaic time transect | `archaic-pipeline oceania-transect --panel 1240k` | [report](reports/oceania_transect/PAPER_oceania.md) |
| Which Neanderthal? Altai vs Vindija source contrast | `archaic-pipeline neanderthal-source --panel 1240k` | [report](reports/neanderthal_source/PAPER_neanderthal_source.md) |
| Global and Eurasian >5% survey | `python scripts/global_archaic_survey.py` | [study index](docs/studies-and-reports.md) |
| Neanderthal admixture dating | `archaic-pipeline admixture-date` | [module guide](docs/neanderthal_admixture_dating.md) |
| Oase1 segment analysis | `python scripts/oase1_haplotype.py` | [study index](docs/studies-and-reports.md) |
| West-Eurasian source ancestry | `archaic-pipeline ancestry` | [study index](docs/studies-and-reports.md) |
| Etruscan case study | `python scripts/etruscan_study.py` | [study index](docs/studies-and-reports.md) |

The [study and report index](docs/studies-and-reports.md) labels which outputs
are validation artifacts, supported summaries, exploratory results, or
data-limited negative/inconclusive findings.

## Data and repository policy

Large AADR panels remain outside Git. Five small, published Sima de los Huesos
BAMs from ENA study PRJEB10597 are intentionally retained as a documented
exception for the exploratory adapter in `tools/sima_de_los_huesos_scan.py`.
Their accessions, publisher MD5 values, repository SHA-256 values, and scientific
limits are recorded in [docs/DATA.md](docs/DATA.md).

See [docs/DATA.md](docs/DATA.md) before adding data or generated results — it
also records which pipeline tables are tracked and which are regenerated. The
source code is MIT licensed; upstream datasets retain their own terms and
citation requirements.

## Documentation

- [Documentation index](docs/index.md)
- [Getting started](docs/getting-started.md)
- [Methods and interpretation](docs/methods-and-interpretation.md)
- [Studies and reports](docs/studies-and-reports.md)
- [Data and artifact policy](docs/DATA.md)
- [Simulation validation](docs/studies/SIMULATION_VALIDATION.md)
- [Papuan Denisovan admixture dating V1](archaic_admixture_dating/README.md)
- [Roadmap](docs/ROADMAP_GENETICS.md)
- [Changelog](CHANGELOG.md)
- [Contributing](CONTRIBUTING.md)
- [Release process](RELEASING.md)

## Citation and license

Use the repository's [CITATION.cff](CITATION.cff) metadata when citing the
software, and cite Mallick et al. (2024) plus the primary publications for AADR
samples used in an analysis.

Code is released under the [MIT License](LICENSE). Genomic data are not
relicensed by this repository.
