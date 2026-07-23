# Feasibility report — Native American Denisovan ancestry study

**Stage:** feasibility / pre-implementation (project brief Section 35, item 10 — stop here for review).
**Module:** `native_american_denisovan_project/` (renamed from `GLMNATIVEAMERICANDENISOVAN/`).
**Status of the hypothesis:** **not yet tested; the feasibility probe already leans toward a weak/absent or unresolvable outcome.** No claim of "Papuan-like Denisovan ancestry in Native Americans" is made or implied.

This report is grounded in (a) verified primary literature fetched this session, (b) a read-only scan of the local AADR v66.p1 1240K annotation, and (c) block-jackknife D-statistics and diagnostic-SNP counts computed directly on the AADR 1240K panel reusing the validated production `archaic.*` engine. All numbers below are reproducible from the scripts in `scripts/`.

---

## 1. Exact interpretation of the hypothesis

The central question is whether Native Americans carry a Denisovan ancestry component more closely related to the **Papuan-associated** Denisovan component than to the principal **East-Asian** Denisovan component (Browning et al. 2018). Per the brief and the Section-4 analysis (`docs/browning_figure4_interpretation.md`), this is **not** what Browning et al. 2018 Figure 4 demonstrated. Browning's American panels were the **admixed** 1000 Genomes populations (PUR/CLM/MXL/PEL); the authors explicitly attributed their lower match rates to **admixture/LD false positives** and any Denisovan segments to **East-Asian-related Native ancestors**, and the formal two-component test was **not significant in any American population**. The project hypothesis is therefore an **unreplicated over-reading of Figure 4's layout**, and the project's job is to test whether it survives controlled statistics — with the real possibility that the correct outcome is "no replication" (D) or "unresolvable" (E).

The hypothesis is reframed as a falsifiable, preregistered contrast (Section 25): *after conditioning on East-Asian-related and Ancient-North-Eurasian-related ancestry, Native American populations show greater sharing with Papuan-associated Denisovan markers than with East-Asian-associated Denisovan markers.* The statistic is **allele-sharing similarity**, not "descent from Papuans" and not a percentage.

## 2. Relevant literature (full matrix in `docs/literature_matrix.tsv`)

Three core papers verified by fetching primary records this session:
- **Browning et al. 2018** (Cell 173:53-61.e9, DOI 10.1016/j.cell.2018.02.031, PMC5866234): reference-free Sprime; two Denisovan components (high-altair-affinity ~1/3 of East Asian segments; moderate-affinity dominant in Papuans/South Asians); two-component test significant **only** in CHS/CHB/CDX/JPT (+FIN/Punjabi), **not Americans**; American signal attributed to admixture/LD + East-Asian-related Native ancestors.
- **Qin & Stoneking 2015** (Mol Biol Evol 32:2665-2674, DOI 10.1093/molbev/msv141, PMID 26104010): "very low level of Denisovan ancestry across Eastern Eurasian and Native American (EE/NA)" **equally correlated** with New Guinea and Australian ancestry → common source for EE/NA and Oceanian Denisovan ancestry. This is already the hypothesis, at SNP resolution.
- **Jacobs et al. 2019** (Cell 177:1010-1021, DOI 10.1016/j.cell.2019.02.035): ≥2 deeply divergent Denisovan ancestries in Papuans → the "structured Denisovan source" model (Model 4).

Supporting/contextual papers (matrix has 16 entries): Reich 2010/2011, Meyer 2012, Huerta-Sanchez 2014 (EPAS1), Skoglund & Jakobsson 2011, Patterson 2012 (antisymmetry), Mallick 2024 (AADR), Raghavan 2014/2015 (ANE/Siberian), Moreno-Mayar 2018 (First Americans), Lazaridis 2014 (ANE), Skoglund & Mathieson 2018, Browning & Browning 2015 (tract dating needs phased WGS).

## 3. Available datasets (read-only, this machine)

- **AADR v66.p1 1240K**: 23,089 individuals × 1,233,013 SNPs, TGENO packed format, at `C:/Users/benne/aadr_v66/`. Pseudo-haploid ancient calls; present-day diploid. `.anno` metadata parsed by the validated `archaic.anno` loader.
- **AADR v66.p1 Human Origins (HO)**: 27,594 × 584,131 (more present-day diversity; available for ascertainment/positive-control checks).
- **Archaic references present in AADR 1240K**: Denisova.SG (574k autosomal SNPs), Denisova3.DG, Denisova3_snpAD.DG, Denisova11.SG, Denisova25.SG, AltaiNeanderthal.DG (1.15M SNPs), VindijaG1_final.SG (528k), Chagyrskaya8.DG (1.06M), Chimp.REF (1.10M). No second independent high-coverage Denisovan (Denisova 11 is an F1 Neanderthal×Denisovan; Denisova 25 is provisional) → an **absolute Denisovan-ancestry percentage is not identifiable**, matching the existing pipeline's stated limitation.
- **African outgroups**: Mbuti (n=15), Yoruba/YRI (n=127) — plus Mende, Dinka, Ju_hoan, San in HO.

## 4. Number of suitable ancient Native American and related individuals

From `scripts/01_inventory_aadr.py` (SNP floor 30,000 1240K; the pipeline's high-confidence threshold):

| group | n usable ancient | date range (BP) | median 1240K SNPs | median cov |
|------|----:|----|----:|----:|
| Ancient South America | 374 | 500–11,885 | 481,950 | 0.71× |
| Ancient Caribbean | 245 | 510–3,188 | 706,670 | 2.00× |
| Ancient North America | 155 | 500–12,712 | 549,351 | 0.81× |
| Ancient Mesoamerica | 154 | 500–9,534 | 335,835 | 0.52× |
| Ancient Arctic/Beringian | 23 | 505–3,885 | 649,435 | 0.98× |
| **Total usable ancient American** | **951** | | | |
| Ancient Paleo-Siberian | 47 | — | 761,870 | 1.50× |
| Ancient Siberian | 272 | — | 666,845 | 1.42× |
| Ancient North Eurasian | 8 | — | 370,190 | 0.41× |
| Ancient East Asian | 1,253 | — | 512,514 | 0.89× |
| Ancient Jomon Japan | 16 | — | 703,503 | 1.29× |

Key named individuals **confirmed present**: USR1.SG (Ancient Beringian, 11,425 BP, 1.15M SNPs), Kolyma1.SG (Ancient Paleo-Siberian, 9,775 BP), MA1.SG (Mal'ta / ANE, 24,320 BP), Yana1.SG (Upper Paleolithic Siberian, 31,850 BP, 1.15M SNPs), Sumidouro6.SG (Lagoa Santa, ~10k BP). **Not found** under the guessed IDs Anzick-1.SG / SpiritCave.SG (real IDs differ — to be resolved in the manifest build).

Present-day controls (1240K): Papuan (n≈32), Han (46), Japanese (31), Dai (14), Karitiana (16), French (31), Yoruba (127), Mbuti (15). HO adds more Australasian/Andamanese/Negrito diversity.

## 5. Denisovan-informative SNP overlap with 1240K (empirical)

From `scripts/02_diagnostic_snp_power.py` (frequency-based, polarisation-invariant, identical rule to the validated `denisovan_genome._marker_sharing`):

- **Set A** (strict: Denisova high ≥0.90, African mean ≤0.10, Neanderthal mean ≤0.10): **10,942** autosomal SNPs on 1240K.
- **Set B** (shared: Neanderthal ≤0.50): **12,180** SNPs.
- **Set F** (shared-archaic negative control: Nea mean ≥0.90): 13,852 SNPs.
- **Transversion-only**: Set A **1,738**, Set B **1,998**, Set F 2,480 (across a 228,560-SNP TV sub-panel).

Per-individual callability of Set-B markers (the power number):
| sample | role | callable Set-B markers |
|------|------|----:|
| USR1.SG | Ancient Beringian | 12,180 |
| Kolyma1.SG | Ancient Paleo-Siberian | 12,177 |
| MA1.SG | ANE (Mal'ta) | 8,949 |
| Yana1.SG | UP Siberian (Yana) | 12,178 |
| Sumidouro6.SG | Lagoa Santa | 10,857 |
| [pool] Papuan | positive control | 12,126 |

**Power is NOT the binding constraint.** There are ~12k diagnostic markers and every key ancient individual carries ~9–12k of them. The constraint is **signal size**: Native American Denisovan ancestry is small, and raw marker frequencies are dominated by ancestral polymorphism (Mbuti baseline reads 1.90% and Yoruba 3.68% at the "Denisovan" marker despite zero Denisovan ancestry; French 5.61%), which is precisely why naive marker-counting is misleading and D-statistics (which normalize out the baseline) are required.

## 6. Power concerns — what the feasibility probe actually showed

`scripts/03_dstat_probe.py` ran block-jackknife (50-block) D-statistics on the key individuals. Three findings that determine feasibility:

**(a) Total Denisovan affinity is null in Native Americans/Beringians (replicates the production survey).** `D(X, Mbuti; Denisova, Chimp)` all-sites: Papuan Z=5.92 (control), Karitiana Z=0.18, USR1 Z=−0.48, Kolyma1 Z=1.33, MA1 Z=1.06, Yana1 Z=0.44, Sumidouro6 Z=0.24. **No ancient American/Beringian is significant.** A total-affinity test alone cannot support the hypothesis.

**(b) Transversion-only analysis has a baseline shift.** In TV-only, French becomes +1.55% (Z=2.45) and Han/Japanese/Dai turn small-positive — so TV absolute values are **not comparable** to all-sites values and only **TV contrasts** (Tier-2) are interpretable. This is a methodological caveat any result must state.

**(c) The contrast `D(X, Han; Denisova, Chimp)` shows a suggestive but unstable "more Denisovan than Han" signal.** All-sites: Kolyma1 Z=2.53 (significant). TV-only: Yana1 Z=3.77, Karitiana Z=2.05; MA1 Z=1.77, Kolyma1 Z=1.86 (trend). This is the only place a "different-from-Han" signal peeks through, and it is **unstable across all-sites vs TV and across individuals**. It is hypothesis-generating, not a result, and exactly the kind of thing the full Tier-5 diagnostic-marker + ancestry-stratified analysis must adjudicate.

**(d) One candidate statistic failed validation and is excluded.** `D(Altai, Denisova; X, Mbuti)` (the "obvious" Neanderthal-conditioned form) scores French +3.57% (Z=6.33) but **Papuan +0.35% (Z=0.44)** — it is a **Neanderthal indicator, anti-correlated with Denisovan ancestry**, and would have produced a wrong conclusion if used uncritically. This is locked in `tests/test_dstat_sign.py` (5/5 pass) and documented in `docs/statistic_interpretation.md`. It validates the project's "derive sign first" rule and is the single most important methodological outcome of the feasibility stage.

## 7. Methods possible with AADR (Tier 1–5, 14, 16–18)

- **Tier 1–3** D-statistics with block jackknife (all-sites, transversions-only, multiple outgroups): directly on AADR using the validated engine. ✓
- **Tier 4–5** diagnostic-SNP sets (Set A/B + Papuan/Han-associated subsets from published tracts, training/validation partitioned) + allele-sharing profile correlations, conditional marker enrichment, f4 contrasts, residual/regression, mixture fitting with block jackknife: all feasible on AADR pseudo-haploid frequencies. ✓
- **Tier 14** ancestry-stratified regression (Denoise using published qpAdm/ADMIXTURE East-Asian + ANE fractions): feasible; qpAdm via the pipeline's Python qpAdm or ADMIXTOOLS 2 (installed, R 4.6.1). ✓
- **Tier 16** downsampling curves: feasible (high-coverage genomes projected to 1240K and subsampled). ✓
- **Tier 17** transversion, TV-only, mappability mask, repeat exclusion, CpG: feasible to the extent derivable from the .snp file + the 1240K capture design (terminal-base and read-level filters **cannot** be reconstructed from AADR genotype calls — stated limitation).
- **Tier 18** per-autosome + leave-one-chromosome-out: feasible. ✓

## 8. Methods requiring whole genomes / phasing (Tier 21) — NOT attempted on AADR

Sprime/S*, IBDmix, hmmix, ArchaicSeeker, **tract-length dating**, haplotype clustering, archaic-tract phylogenies, and modern local-ancestry restriction of admixed Americans all require **phased whole genomes**. AADR 1240K pseudo-haploid calls cannot support them, and the brief explicitly forbids haplotype-length dating on pseudo-haploid 1240K. The whole-genome tier is therefore a **future-data** item (Section "Future data needed"), not a feasibility-stage deliverable.

## 9. Ethical / access restrictions (full statement in `ETHICS_AND_DATA_USE.md`)

- AADR v66.p1 is **publicly released** for secondary analysis; individual-level genotypes are redistributable under the AADR terms by the data generators, but this project **redistributes no genotype data** — only aggregate statistics, metadata manifests with reduced geographic precision, and diagnostic-site lists.
- Ancient Native American individuals require care: no tribal/community-permission documentation is carried in AADR metadata; the project makes **no claim of identity, cultural affiliation, or descent** from Denisovan ancestry, and treats archaeological dates as ranges, not points.
- Sample names are retained as publication IDs (needed for reproducibility) but coordinates are rounded; no modern admixed American individual is used without flagging the admixture.

## 10. Recommended minimum viable analysis (Section 33), scoped to feasibility

1. Build the curated manifest (already started: `data/manifests/aadr_inventory_relevant.tsv`) and resolve real IDs for Anzick-1 / Spirit Cave / other named samples. ✓-ish
2. Run the **S1** and **S2** D-statistics across the full time transect (951 ancient Americans + Siberians + present-day controls), all-sites and transversions-only. *(S1 probe done on a subset; S2 not yet run.)*
3. Construct the high-confidence Set-A diagnostic set (10,942 sites already enumerated) and the Papuan-associated / East-Asian-associated subsets from **published** tract lists (training/validation split).
4. Compute the **S3** conditional contrast across the transect, with block-bootstrap CIs, leave-one-chromosome-out, and selected-locus (EPAS1 etc.) exclusions.
5. Run coverage-matched downsampling to set reliable/provisional/insufficient thresholds.
6. Regression of the Denisovan signal on East-Asian + ANE ancestry fractions.
7. Cautious written interpretation against Models 0–6 and Outcomes A–E.

The sign gate (`tests/test_dstat_sign.py`, passing) is a **mandatory precondition** before step 4 produces any reported number.

## 11. Estimated computational requirements

- The inventory + diagnostic-SNP + probe scripts ran in **minutes** on this machine (15.6 GB RAM); reading a handful of individuals from the 7 GB TGENO via the individual-major reader is cheap.
- Full Tier-1/2 over ~1,200 ancient individuals × ~560k shared SNPs: a few hours (the production Phase 3 over all ~15k Eurasians took ~27 min; this is a strict subset).
- Tier-5 diagnostic-marker profiles over the transect: minutes (only ~12k markers per individual, vectorised).
- Simulations (msprime WGS projected to 1240K for Models 0–6 power/FPR): the largest item, ~hours; deferred to the full workflow.

## 12. Reasons the project could fail (honest prior)

1. **Signal too small**: total Denisovan affinity is already null in Native Americans (probe finding a); the composition contrast may be below detection even with 12k markers.
2. **Ancestral polymorphism / Neanderthal confounding**: raw marker frequencies are dominated by it (Yoruba 3.68% baseline); if the Set-A filtering does not fully remove it, S3 is biased.
3. **Transversion baseline shift** (finding b): TV-only results are only interpretable as contrasts; if a result rests on TV absolute values it is unsound.
4. **Unstable single-individual signals** (finding c): Kolyma1/Yana1/Karitiana each show one significant contrast in one mode; if the signal does not generalize across individuals and chromosomes it is an artifact.
5. **Wrong-statistic trap** (finding d): already cost one statistic; further naive forms (e.g. `D(NativeAmerican, Han; Papuan, African)` labelled "Denisovan") measure broader affinity, not Denisovan composition.
6. **Ascertainment circularity**: defining "Papuan-associated" markers with Papuans then "finding" Papuan enrichment is circular; must use published tracts or partitioned training sets.
7. **Modern admixture** (Browning's own caveat): present-day admixed Americans must be local-ancestry-restricted or excluded; ancient genomes are prioritised but lower-coverage and pseudo-haploid.
8. **Single Denisovan reference**: only the Altai Denisovan is high-coverage; "Papuan-like vs Han-like" is inferred from match-rate differences to one genome, not a true phylogeny — the whole-genome tract tier would be needed to resolve, and it is out of reach on AADR.
9. **Outcome E is scientifically acceptable**: if S1 is null, S2 is null, and S3 is null after conditioning, the correct conclusion is "no replication / unresolvable with AADR," not a forced positive.

## Bottom line / recommendation

The data and methods exist to **test** the hypothesis rigorously on AADR 1240K. The feasibility probe already shows that (i) **total** Denisovan ancestry in Native Americans is null, (ii) a "more-Denisovan-than-Han" signal is **suggestive but unstable** in a few ancient Siberians/Karitiana, and (iii) one plausible statistic is **wrong** and was caught only by the sign gate. The recommended next step is the **minimum viable analysis** (Section 10 above), gated on the passing sign tests, with the explicit expectation that the most likely outcomes are C (weak/suggestive), D (no replication), or E (unresolvable), and that **Outcome A is unlikely given the probe**. The project will not, at this stage, assert that Native Americans carry a Papuan-like Denisovan component.

*Stop point: per project brief Section 35, the full workflow is NOT implemented until this feasibility report is reviewed.*
