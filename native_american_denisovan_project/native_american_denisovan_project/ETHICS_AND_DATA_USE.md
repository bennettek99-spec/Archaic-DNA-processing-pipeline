# Ethics and data use — Native American Denisovan ancestry study

This project analyses genetic data from ancient Indigenous individuals of the
Americas and from present-day populations. It does so under the constraints
below. These are binding project rules, not aspirations.

## Data sources and access

- All genotype data are from the **Allen Ancient DNA Resource (AADR) v66.p1**,
  a publicly released secondary-analysis resource (Mallick et al. 2024, Sci. Data).
  AADR aggregates published, consented, publicly deposited genotype calls; the
  project uses the local copy at `C:/Users/benne/aadr_v66/` for analysis only.
- The project **redistributes no individual-level genotype data**. Outputs are
  aggregate statistics, metadata manifests, diagnostic-site lists, and figures.
- The project does **not** re-consent, re-access, or pull raw reads (BAM/FASTQ);
  terminal-base / read-level damage filters cannot be reconstructed from AADR
  calls and are stated as a limitation, not worked around.

## Indigenous data and identity

- The project makes **no claim** about tribal identity, cultural affiliation,
  biological descent, or migration history of any Indigenous group from
  Denisovan ancestry or any other statistic.
- Ancient genomes are treated as **data points from published studies**, not as
  property of this project; sample IDs are retained as publication identifiers
  (required for reproducibility) but are not relabelled as "ours."
- Geographic precision in published manifests is **reduced** (coordinates
  rounded) to limit site re-identification where ethically appropriate, while
  retaining enough resolution for regional grouping.
- The project recognises that AADR metadata does **not** document
  tribal/community permissions for every included individual; secondary
  analysis of public data is permitted by the data generators' release terms,
  but the absence of community-level consent is a real limitation and is noted
  rather than ignored.
- Results are framed at the **population/regional** level; Native American
  populations are **not** treated as genetically uniform, and no result is
  presented as characterising a specific tribe or nation.

## What the project will not do

- Infer tribal identity from genetics.
- Present ancient genomes as the property of the project or its authors.
- Redistribute controlled-access or raw data.
- Make claims about Indigenous identity based on Denisovan ancestry.
- Frame Native American populations as genetically uniform.
- Interpret genetic affinity as proof of direct cultural or biological descent.
- Use "Papuan-like ancestry" where the result only indicates allele-sharing
  similarity (see `docs/statistic_interpretation.md` language discipline).
- Force a positive result; "no replication" or "unresolvable" is an acceptable
  and scientifically useful outcome.

## Modern admixed American data

- Present-day admixed American populations (e.g. 1000 Genomes PUR/CLM/MXL/PEL)
  are used only as flagged **admixed** controls and are never used to infer
  unadmixed Native American Denisovan ancestry without local-ancestry
  restriction (out of scope on AADR pseudo-haploid data; deferred to a
  whole-genome tier).
- Pre-contact **ancient** genomes are prioritised for historical inference
  because modern admixture can seriously distort low-level Denisovan signals
  (Browning et al. 2018 explicitly flag admixture/LD false positives).

## Reproducibility vs. redistribution

- All scripts are committed to the repository; the AADR data are **not** (they
  are machine-local, gitignored, and large). A user reproduces the analysis by
  pointing `config/paths.yaml` at their own AADR checkout — no data are bundled.
- Provenance (panel, SNP counts, software versions, commands) is recorded in
  every output manifest, mirroring the production pipeline's convention.

## Limitations of this ethics statement

- It governs this project's use of public AADR data; it does not substitute for
  the original data generators' consent and community-engagement processes,
  which remain the responsibility of the primary studies.
- It does not constitute legal advice; the AADR terms of use are the binding
  instrument and may evolve with future releases.
