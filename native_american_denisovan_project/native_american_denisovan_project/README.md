# Native American Denisovan ancestry study — research module

A **separate, experimental** research module (not part of the production
archaic-introgression pipeline) that investigates whether Native Americans
carry a Denisovan ancestry component more closely related to the
**Papuan-associated** Denisovan component than to the principal **East-Asian**
component — the hypothesis suggested (but **not** demonstrated) by Figure 4 of
Browning et al. 2018.

**Current stage: feasibility.** Per the project brief (Section 35), the full
workflow is NOT implemented until the feasibility report is reviewed. This
repository contains the feasibility-stage deliverables plus the scaffolding and
empirical probes that ground the feasibility assessment.

> **Honest framing.** The feasibility probe already shows total Denisovan
> affinity in Native Americans is null (replicating the production survey), a
> "more-Denisovan-than-Han" signal is suggestive but unstable in a few ancient
> Siberians, and one plausible statistic failed validation. The expected
> outcomes are C/D/E, not a forced positive. **No claim of "Papuan-like
> Denisovan ancestry in Native Americans" is made.**

## What is in this module

```
native_american_denisovan_project/
├── FEASIBILITY_REPORT.md            <-- START HERE (the deliverable)
├── ETHICS_AND_DATA_USE.md
├── README.md                        (this file)
├── .gitignore
├── config/  default.yaml, paths.example.yaml, populations.yaml
├── data/manifests/                  aadr_inventory_full.tsv, aadr_inventory_relevant.tsv
│   └── diagnostic_sites/            denisovan_diagnostic_sites.tsv (12,180 Set-B sites)
├── scripts/
│   ├── 01_inventory_aadr.py         scan .anno -> sample/population inventory (Tables 1,2)
│   ├── 02_diagnostic_snp_power.py   Denisovan-diagnostic SNP counts + callability (power)
│   └── 03_dstat_probe.py            Tier-1/2/3 block-jackknife D-stat feasibility probe
├── tests/  test_dstat_sign.py       sign/orientation gate (5/5 pass) — MANDATORY before results
├── results/tables/                  table1..table4 + diagnostic SNP tables (TSV)
├── results/logs/                    inventory_summary.txt, diagnostic_snp_power.txt, dstat_probe.txt
└── docs/
    ├── browning_figure4_interpretation.md   Section 4 critical interpretation
    ├── statistic_interpretation.md          3 candidate statistics + signs + the failed one
    ├── synthetic_validation_design.md       sign/artifact/power test design
    └── literature_matrix.tsv                16-paper verified literature matrix
```

## Key empirical findings (feasibility stage)

1. **Data exist.** AADR v66.p1 1240K has **951 usable ancient American**
   individuals plus all key anchors (USR1, Kolyma1, MA1, Yana1, Sumidouro6) and
   strong present-day controls (Papuan, Han, Karitiana, French, Mbuti, Yoruba).
2. **Power is not the binding constraint.** **12,180** Denisovan-diagnostic
   Set-B SNPs on 1240K (10,942 strict Set-A; ~2,000 transversions); each key
   ancient individual carries ~9–12k of them.
3. **Total Denisovan affinity is null in Native Americans.**
   `D(X, Mbuti; Denisova, Chimp)` all-sites: Papuan Z=5.92 (control), Karitiana
   Z=0.18, USR1 Z=−0.48, Kolyma1 Z=1.33, MA1 Z=1.06, Yana1 Z=0.44, Sumidouro6
   Z=0.24. No ancient American is significant.
4. **Transversion-only has a baseline shift** (French +1.55%, Z=2.45) — only
   TV *contrasts* are interpretable.
5. **A "more-Denisovan-than-Han" signal is suggestive but unstable**:
   `D(X, Han; Denisova, Chimp)` — Kolyma1 Z=2.53 (all-sites); Yana1 Z=3.77,
   Karitiana Z=2.05 (TV). Hypothesis-generating, not a result.
6. **One statistic failed validation and is excluded.**
   `D(Altai, Denisova; X, Mbuti)` scores French +3.57% but Papuan +0.35% — a
   Neanderthal indicator anti-correlated with Denisovan ancestry. Locked in
   the passing sign-gate.

## Running (uses the sibling pipeline's validated venv — read-only reuse)

```powershell
# copy the machine-local path file (gitignored) and edit if needed
cp config/paths.example.yaml config/paths.yaml

# 1. inventory (Tables 1,2) — reads only the .anno file
..\archaic-introgression\.venv\Scripts\python.exe scripts\01_inventory_aadr.py

# 2. diagnostic-SNP power estimate (reads a handful of archaic/African refs)
..\archaic-introgression\.venv\Scripts\python.exe scripts\02_diagnostic_snp_power.py

# 3. D-stat feasibility probe (Tier 1/2/3)
..\archaic-introgression\.venv\Scripts\python.exe scripts\03_dstat_probe.py

# sign-gate (mandatory before any real-data S3 result)
..\archaic-introgression\.venv\Scripts\python.exe tests\test_dstat_sign.py
```

`PYTHONIOENCODING=utf-8` is recommended (Windows console is cp1252).

## Relationship to the production pipeline

This module is **independent**: it imports the validated, read-only
`archaic.*` reader/statistics engine (TGENO reader, block-jackknife f4/D,
`.anno` parser) by adding `../archaic-introgression` to `sys.path`, but it does
**not** modify, overwrite, or restructure the production pipeline. Integration
into the main pipeline is deferred until manual review, per the brief.

## Primary statistic (preregistered)

`S3`: *After conditioning on East-Asian-related and Ancient-North-Eurasian-related
ancestry, Native American populations show greater sharing with Papuan-associated
Denisovan markers than with East-Asian-associated Denisovan markers.* — an
allele-sharing-similarity statement, not a descent or percentage claim. Full
derivation and the excluded statistic in `docs/statistic_interpretation.md`.

## References (verified)

- Browning et al. 2018, Cell 173:53-61.e9. DOI 10.1016/j.cell.2018.02.031 (PMC5866234).
- Qin & Stoneking 2015, Mol Biol Evol 32:2665-2674. DOI 10.1093/molbev/msv141.
- Jacobs et al. 2019, Cell 177:1010-1021. DOI 10.1016/j.cell.2019.02.035.
- Full 16-entry matrix: `docs/literature_matrix.tsv`.
