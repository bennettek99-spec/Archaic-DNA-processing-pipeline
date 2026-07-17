# Studies and reports

This page maps each analysis to its entry point, retained output, and intended
evidence label. It distinguishes generated artifacts from independently
reviewed scientific conclusions.

| Analysis | Command | Canonical output | Status |
| --- | --- | --- | --- |
| Core AADR validation | `archaic-pipeline validate --panel 1240k` | [VALIDATION.md](../VALIDATION.md) | Validated method/control suite |
| Coalescent validation | `python validate_simulation.py` | [SIMULATION_VALIDATION.md](../SIMULATION_VALIDATION.md) | Validated simulation benchmark |
| Phase 2-9 AADR workflow | `archaic-pipeline all --panel 1240k` | `reports/archaic_report_1240k.html` | Supported summaries plus explicitly labelled candidates |
| Highest-archaic scan | `archaic-pipeline highest-archaic` | `reports/highest_archaic_ancestry_report.md` | Credibility-aware ranking; candidate interpretation |
| Denisovan reference genome | `archaic-pipeline denisovan-genome --panel 1240k` | `results/denisovan_genome/report.html` | Reference/lineage characterization; no Denisovan percentage |
| Global >5% survey | `python global_archaic_survey.py` | `reports/global_archaic_survey/PAPER_global.md` | Neanderthal threshold result supported; Denisovan context is relative or literature-derived |
| Eurasian >5% survey | `python high_archaic_survey.py` | `reports/high_archaic_survey/PAPER.md` | Supported high-confidence threshold result; raw crossings audited as artifacts |
| Oase1 array segments | `python oase1_haplotype.py` | `reports/oase1_haplotype/PAPER_oase1.md` | Focused contextual evidence; read-level workflow preferred for confirmation |
| West-Eurasian ancestry | `archaic-pipeline ancestry` | `reports/ancestry/PAPER_ancestry.md` | Supported descriptive model comparison |
| Etruscan case study | `python etruscan_study.py` | [ETRUSCAN_FINDINGS.md](../ETRUSCAN_FINDINGS.md) | Group-level supported summary; individual/locus claims exploratory |
| FADS analysis | `python fads_report.py` | [FADS_REPORT.md](../FADS_REPORT.md) | Selection context; not proof of adaptive introgression |
| Local archaic windows | `python local_archaic_scan.py --panel 1240k` | [LOCAL_ARCHAIC_REPORT.md](../LOCAL_ARCHAIC_REPORT.md) | Peaks are exploratory; desert recovery failed |
| X-chromosome depletion | `python xchrom_depletion.py --panel ho` | [XCHROM_REPORT.md](../XCHROM_REPORT.md) | Inconclusive/data-limited |
| Denisovan population survey | `python denisovan_survey.py --panel 1240k` | [DENISOVAN_REPORT.md](../DENISOVAN_REPORT.md) | Validated relative-affinity positive control; no percentage |
| Sima de los Huesos BAM adapter | `python tools/sima_de_los_huesos_scan.py` | `output/sima_de_los_huesos/` when run | Exploratory and below normal AADR evidence floor |

## Highest-archaic reporting

The highest-archaic workflow deliberately reports:

- the raw numerical Neanderthal maximum;
- the credibility-aware supported maximum;
- the strongest Denisovan-affinity statistic without converting it to a
  percentage;
- technical-risk flags and contradictory sensitivity tests;
- candidate-specific reports and a provenance manifest.

See the [module guide](highest_archaic.md).

## Global Denisovan context

The global survey recovers high relative Denisovan affinity in Oceania. Any
statement that living Oceanians exceed 5% combined archaic ancestry depends on
published Denisovan percentage estimates, not a Denisovan percentage calculated
by this repository. Reports and summaries must keep that source distinction
explicit.

## Generated artifacts

HTML reports and PDFs are retained for portability, while CSV/TSV tables are
retained when they provide inspectable evidence. A file's presence in `results/`
or `reports/` does not elevate it from exploratory to validated. The status in
this index and the limitations in the report govern interpretation.
