# Genome-wide detection of unexpected archaic introgression in ancient Eurasian genomes

A rigorous, transparent, pure-Python pipeline that estimates Neanderthal (and,
where feasible, Denisovan) ancestry for ancient Eurasian genomes in a local AADR
database, and flags individuals whose archaic ancestry is statistically unusual
*given their ancestry, geography, age, and context* — not simply the highest.

**Status: Phase 1 (validation) COMPLETE — estimator validated (7/7 gates pass).**
Findings from later phases are exploratory hypotheses, never discoveries, until
technical explanations are ruled out and compared with the literature.

---

## Data

* **Panel:** AADR v66.p1 **Human Origins** (`C:\Users\benne\aadr_v66\v66.p1_HO.*`),
  27,594 individuals × 584,131 SNPs (579,720 autosomal), packed **TGENO**
  (individual-major) layout — ideal for our "few individuals, all SNPs" access.
* **Archaic references present in-panel:** `AltaiNeanderthal.DG` (high-cov),
  `VindijaG1_final.SG` (high-cov, second Neanderthal used as the "100%" scale),
  `Denisova.SG` (high-cov Denisovan), late Neanderthals (Mezmaiskaya1/2, Spy,
  Goyet, Les_Cottès), and `Denisova11.SG` (a Nea×Den F1 — *not* used as a clean
  Denisovan).
* **Outgroup / baseline:** `Chimp_HO.HO` (rooting), Mbuti & Yoruba (African,
  ~no Neanderthal — the introgression baseline).

## Method

All statistics are allele-frequency f-statistics (Patterson et al. 2012) on the
autosomes, with standard errors from a **50-block delete-one jackknife** over
contiguous equal-SNP-count genomic blocks (Busing et al. 1999). f4 and the
normalised D are invariant to per-SNP allele polarisation, so a single consistent
allele coding is used and no explicit "derived allele" call is needed.

1. **Neanderthal proportion (absolute, cross-sample comparable)** — direct
   f4-ratio:
   ```
   alpha(X) = f4(Altai, Chimp ; X, Mbuti) / f4(Altai, Chimp ; Vindija, Mbuti)
   ```
   Two *different* high-coverage Neanderthals — Altai in the statistic, Vindija
   as the "100% Neanderthal" yardstick — so the estimate is not biased by using
   one genome as both source and scale (Reich et al. 2009; Patterson et al. 2012;
   Petr et al. 2019). The ratio form cancels the test sample's private drift,
   making `alpha` comparable across individuals of differing coverage (a raw f4
   or D would not).

2. **Neanderthal affinity (relative, for significance)** —
   `D(X, Mbuti ; Altai, Chimp)` (>0 ⇒ excess Neanderthal-derived sharing vs an
   African baseline).

3. **Differential Neanderthal between two test pops** —
   `D(Pop1, Pop2 ; Altai, Yoruba)`. Note the **African outgroup (Yoruba)**: it
   isolates *recent* Neanderthal-derived sharing because Africans carry the
   modern-human ancestral background with ~no Neanderthal. Chimp in that position
   dilutes the signal with deep ancestral variation (verified empirically — see
   `probe_eastasian.py`).

4. **Denisovan affinity (relative only)** — `D(X, Mbuti ; Denisova, Chimp)`.
   With a single high-coverage Denisovan there is no second-archaic scale, so we
   report relative affinity / significance, **not** an absolute Denisovan
   fraction.

## Assumptions (and why they hold here)

* Allele coding is consistent per SNP across all individuals (guaranteed by the
  packed `.geno`); f4/D are polarisation-invariant, so the EIGENSTRAT allele1/2
  convention need not be resolved.
* Mbuti/Yoruba ≈ zero Neanderthal ancestry (the standard baseline; African
  back-to-Africa Neanderthal is ~0.1–0.3% and treated as the null floor — Yoruba
  reads −0.12% ± 0.28% here, consistent with ~0).
* Altai and Vindija bracket the true introgressing Neanderthal lineage; the
  f4-ratio yields a slightly scaled but well-calibrated proportion (validated to
  ~2–3% for non-Africans, matching the literature).
* 50 contiguous blocks (~11.6k SNPs each) are large enough to absorb LD.

## Validation results (Phase 1 gate — `phase1_validate.py`)

| Test population | category | α (Neanderthal) | D_Nea (Z) | D_Den (Z) | SNPs |
|---|---|---:|---:|---:|---:|
| French | modern W-Eur | 2.75% ± 0.43 | 0.0219 (5.8) | −0.003 (−1.0) | 266,712 |
| Sardinian | modern W-Eur | 2.47% ± 0.43 | 0.0204 (5.4) | −0.005 (−1.3) | 266,712 |
| Han | modern E-Asia | 2.81% ± 0.46 | 0.0229 (5.2) | −0.004 (−1.0) | 266,712 |
| Papuan | modern Oceania | 3.23% ± 0.50 | 0.0334 (6.7) | **+0.0355 (6.5)** | 266,712 |
| Karitiana | modern America | 2.72% ± 0.54 | 0.0241 (5.1) | 0.000 (0.0) | 266,712 |
| **Yoruba** | African control | **−0.12% ± 0.28** | −0.004 (−1.8) | −0.005 (−1.9) | 266,712 |
| Loschbour (WHG) | ancient Eur | 3.85% ± 0.59 | 0.0263 (5.0) | 0.003 (0.4) | 168,754 |
| Stuttgart (LBK) | ancient Eur | 2.55% ± 0.65 | 0.0212 (3.5) | −0.005 (−0.8) | 239,404 |
| Kostenki14 (UP) | ancient Eur | 2.71% ± 0.69 | 0.0218 (3.6) | −0.003 (−0.5) | 248,342 |
| Sunghir (UP) | ancient Eur | 2.79% ± 0.48 | 0.0222 (3.9) | −0.002 (−0.4) | 266,711 |
| MA1 / Mal'ta (UP) | ancient Eur | 2.79% ± 0.73 | 0.0245 (4.1) | 0.005 (0.8) | 195,264 |
| **Mezmaiskaya1** (Neanderthal as test) | control | **99.19%** | 0.861 (179) | 0.347 (53) | 127,255 |
| **Spy** (Neanderthal as test) | control | **97.10%** | 0.852 (189) | 0.355 (41) | 100,431 |

**Gates (all PASS):**
- **G1** non-African α in [1.2%, 3.5%] — 7/8 (Loschbour 3.85% just over, low coverage).
- **G2** East-Asian excess `D(French, Han; Altai, Yoruba) = −0.0071, Z = −2.3` (Han more Neanderthal); Sardinian–Han Z = −2.8.
- **G3** African control α(Yoruba) = −0.12% ≈ 0.
- **G4** Neanderthal-as-test reads **99% / 97%** (the "100% Neanderthal" scale is correct).
- **G5** D_Nea Z > 3 for all non-Africans; |Z(Yoruba)| = 1.8.
- **G6** Denisovan channel: Papuan +0.0355 (Z 6.5) ≫ everyone else ≈ 0.
- **G7** rank order **Sardinian < French < Han < Papuan** — the canonical published gradient.

These reproduce: Green et al. 2010 / Wall et al. 2013 (East-Asian excess), Lazaridis
et al. 2014 (Loschbour/Stuttgart ~1.8–2.2%), Seguin-Orlando 2014 (Kostenki14),
Sikora 2017 (Sunghir), Raghavan 2014 (Mal'ta), Reich 2010 / Meyer 2012 (Papuan
Denisovan).

## Quantified limitations (read before trusting any individual)

* **Per-individual standard error is ~0.4–0.7% on α.** The known East-Asian
  *population-level* excess (~0.4 percentage points) is only ~Z 2–3 even with
  hundreds of HGDP samples per population. **A single low-coverage genome cannot
  support individual-level Neanderthal-% claims at that resolution.** Outliers
  must be interpreted at the group level, or as hypotheses for shotgun/1240k
  follow-up.
* **HO ascertainment** (SNPs polymorphic in modern humans, not enriched for
  archaic alleles) limits archaic-informative sites; the f4-ratio is further
  **capped at ~266k SNPs by Vindija coverage**. The **1240k** panel (~1.15M SNPs)
  would roughly double usable sites and materially improve individual-level
  power — recommended before Phase 6 outlier work is taken seriously.
* **Denisovan is relative-only** within West Eurasia (signal ≈ 0 there anyway).
* Low-coverage and contaminated samples will dominate naive outlier lists; Phase 4
  (normalisation on coverage/SNPs/contamination/damage from the `.anno`) and
  Phase 9 (robustness) exist precisely to remove these.

## Run

```bash
python phase1_validate.py     # validation gate (this README's table)
python probe_eastasian.py     # methods check: outgroup choice for the E-Asian excess
python validate_published.py  # EXTERNAL validation vs published papers -> VALIDATION.md
python generate_report.py --panel 1240k   # self-contained HTML executive summary -> reports/
```

After any run, `generate_report.py` reads whatever `results/` artifacts exist and
writes a portable, figure-embedded HTML highlights page to **`reports/archaic_report_<panel>.html`**
(every section is optional, so it works for partial runs too).

External validation against the literature (Oase1, Ust'-Ishim, Fu 2016 UP series,
modern populations): see **`VALIDATION.md`** — r=0.87 (0.96 excl. Oase1) vs published,
16/17 within range, Oase1 recovered at 6–9%, Papuan Denisovan confirmed.

## Full pipeline (one command)

```bash
python run_pipeline.py --panel 1240k        # Phases 2→9 + HTML report, in order
python -m pytest tests/ -q                  # validate the core f-statistics math
```

## Sub-study: Etruscan archaic introgression (+ research paper)

```bash
python etruscan_study.py     # core analysis -> ETRUSCAN_FINDINGS.md + reports/etruscan_report.html
python etruscan_paper.py     # expanded paper -> PAPER.md + reports/etruscan_paper.html
```

Situates the 75 Etruscans in the Italian time transect and a regional comparison panel
(Anatolia/Steppe/Latins/Romans/Greeks/Sicily/Sardinia + moderns) and asks: did archaic
ancestry change over time, do Etruscans differ from neighbours, are any individuals outliers
*after conditioning on ancestry+geography*, do outliers track genetic ancestry, and **which
archaic genes** shifted in frequency over time. Findings: Etruscan Neanderthal ~2.0%,
indistinguishable from Latins (D-stat Z=0.2); flat over time; the steppe/Levantine/East-Med
ancestry outliers are *not* archaic outliers (decoupled); **FADS1-2** (dietary selection) is
the one introgression locus surviving Bonferroni — a candidate, not proof. The full write-up
is **`PAPER.md`** (manuscript) / `reports/etruscan_paper.html`.

The paper uses **population mean-genome profiles** (`archaic/profiles.py`): per-cohort mean
allele-frequency vectors that average out single-genome noise, giving tight group-level
archaic estimates, cohort-vs-cohort differential tests, and genetic-distance maps.

## Layout

```
archaic/
  lib_eigenstrat.py   memory-mapped AADR EIGENSTRAT/TGENO reader (self-contained copy)
  stats.py            D, f4, f4-ratio + block-jackknife (scalar & vectorised batch)
  panel.py            load panel, select samples, per-pop autosomal allele freqs
  anno.py             parse the AADR .anno metadata
  refs.py             per-panel config (HO / 1240K reference samples, QC thresholds)
  loci.py             locus/gene-level archaic-allele analysis (adaptive-introgression panel)
  profiles.py         population "mean-genome" profiles + group-level archaic + distances
phase1_validate.py … phase9_robustness.py   the nine pipeline stages
etruscan_study.py / etruscan_paper.py        Etruscan sub-study + manuscript (PAPER.md)
validate_published.py external validation vs the literature  -> VALIDATION.md
generate_report.py    self-contained HTML executive summary  -> reports/
etruscan_study.py     focused Etruscan sub-study              -> reports/etruscan_report.html
run_pipeline.py       orchestrator (Phases 2→9 + report)
probe_eastasian.py    methods probe (differential-statistic outgroup choice)
tests/test_stats.py   pytest unit tests for the f-statistics
results/   figures/   reports/             outputs
```

Provenance: reader adapted from the validated `ancient-pca/lib_eigenstrat.py`.

## Key citations

Patterson et al. 2012 *Genetics* 192:1065 (f-statistics/ABBA-BABA) · Green et al.
2010 *Science* 328:710 (Neanderthal draft, D) · Reich et al. 2010 *Nature* 468:1053
(Denisova) · Meyer et al. 2012 *Science* 338:222 (high-cov Denisovan, f4-ratio) ·
Prüfer et al. 2014 *Nature* 505:43 (Altai) · Prüfer et al. 2017 *Science* 358:655
(Vindija 33.19) · Wall et al. 2013 *Genetics* 194:199 (East-Asian excess) ·
Sankararaman et al. 2014 *Nature* 507:354 · Petr et al. 2019 *PNAS* 116:1639 ·
Lazaridis et al. 2014 *Nature* 513:409 (Loschbour/Stuttgart) · Seguin-Orlando et
al. 2014 *Science* 346:1113 (Kostenki14) · Sikora et al. 2017 *Science* 358:659
(Sunghir) · Raghavan et al. 2014 *Nature* 505:87 (Mal'ta) · Busing et al. 1999
*Stat. Comput.* 9:3 (delete-m jackknife) · Mallick et al. 2024 *Sci. Data* (AADR).
```
