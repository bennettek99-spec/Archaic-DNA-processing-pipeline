# Candidate data sources

This is a discovery and governance checklist, not authorization to download or
an assertion that every source can be combined.

| Candidate | Official source | Access | Approximate size | Format/phasing | Coverage | Automated V1 download | V1 suitability |
|---|---|---|---:|---|---|---|---|
| SGDP modern genomes | [Simons Foundation SGDP](https://www.simonsfoundation.org/simons-genome-diversity-project/) and [IGSR SGDP collection](https://www.internationalgenome.org/data-portal/data-collection/SGDP/) | public project, but verify current human-subject and reuse terms | official SGDP page describes roughly 10 TB for the full historical distribution | genome VCF/BAM; release-dependent phasing | 279 genomes from 130 populations in the public project, including Oceania | no by default; select authorized files/components only | strong candidate if current reuse and community-governance review permit |
| HGDP/1000 Genomes high-coverage panel | [IGSR HGDP collection](https://www.internationalgenome.org/data-portal/data-collection/HGDP) | portal publishes a collection-specific reuse policy | many thousands of files; whole collection is not laptop scale | PCR-free high-coverage alignments and release-specific callsets | African and East Asian comparisons; sample coverage is collection-specific | component-only after manifest sizing | comparison references |
| Published Papuan Denisovan tract tables | publication repository or author archive | publication-specific | MB to low GB | tabular caller output | study-specific | only when the official source permits | preferred laptop V1 input |
| High-coverage Altai Denisovan | [Max Planck EVA Denisovan genome project](https://www.eva.mpg.de/genetics/genome-projects/denisova/) | official page states public release without passwords; still record exact file terms | alignments and raw reads are large; genotype derivatives are smaller | hg19/GRCh37 alignments and genotype calls | approximately 30x single Denisovan individual | explicit component only, with size/checksum preflight | archaic reference, not a unique donor proxy |
| SHAPEIT/1000 Genomes genetic maps | versioned official caller/project distribution; record exact URL in run manifest | public release-specific | usually hundreds of MB compressed | chromosome map text with cM positions | genome-wide, assembly-specific | component-wise after checksum/size preflight | dating map candidate; compare at least one alternate map/scale |
| Genome accessibility/problem-region masks | ENCODE/Genome in a Bottle or study-specific mask | public with versioned terms | MB to GB | BED | assembly-specific | yes only with pinned version/checksum | QC masking |

For every actual run, add exact URLs, assembly, release, sample count,
population labels, file sizes, checksums, citations, phasing state, reuse
terms, and automated-download permission to the run input manifest.

## Tract-caller decision

IBDmix is the leading V1 adapter because it is established and does not require
phased target genomes. admixfrog is the preferred secondary sensitivity
candidate where suitable inputs and compute are available. ArchaicSeeker 2.0
and hmmix remain candidates. V1 does not create a new caller.

## Citations to pin during a real-data run

The final report must cite the exact dataset release, tract caller and version,
Altai Denisovan source, recombination map, masks, and primary method papers.
These cannot be finalized before an authorized input dataset is selected.

## Live verification note

The portal links above were verified on 2026-07-25. File counts, sizes, access
terms, and endpoints can change; rerun source discovery and governance review
before any real download.
