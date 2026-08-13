# Documentation

This documentation separates operating instructions, statistical definitions,
and study-specific findings so the README can remain a reliable entry point.

## Start here

- [Getting started](getting-started.md): installation, synthetic smoke testing,
  AADR configuration, and the full workflow.
- [Methods and interpretation](methods-and-interpretation.md): estimators,
  validation, uncertainty, and claims the data cannot support.
- [Studies and reports](studies-and-reports.md): canonical entry points and the
  status of each retained study output.
- [Data and artifact policy](DATA.md): AADR handling, the documented
  PRJEB10597 BAM exception, checksums, and result-retention rules.

## Module guides

- [Highest-archaic AADR scan](highest_archaic.md)
- [Denisovan reference-genome module](denisovan_genome.md)

## Validation and project records

- [External validation](studies/VALIDATION.md)
- [Simulation validation](studies/SIMULATION_VALIDATION.md)
- [Changelog](../CHANGELOG.md)
- [Roadmap](ROADMAP_GENETICS.md)
- [Release process](../RELEASING.md)

## Evidence labels

Documentation uses four labels consistently:

- **Validated:** a method or expected control has passed prespecified checks.
- **Supported summary:** a result is stable enough for descriptive reporting but
  is not by itself a causal or recent-admixture claim.
- **Exploratory:** a candidate or analysis requires independent follow-up.
- **Inconclusive/data-limited:** the available panel cannot identify the
  requested quantity with defensible uncertainty.
