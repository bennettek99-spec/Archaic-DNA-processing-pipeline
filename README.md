# Genome-wide detection of unexpected archaic introgression in ancient Eurasian genomes

[![tests](https://github.com/bennettek99-spec/Archaic-DNA-processing-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/bennettek99-spec/Archaic-DNA-processing-pipeline/actions/workflows/ci.yml)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](pyproject.toml)

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

## Configure (no hard-coded paths)

Point the pipeline at your AADR download by editing **`config.yaml`** (or set
`ARCHAIC_CONFIG=/path/to/your.yaml`):

```yaml
aadr_dir: "/path/to/aadr"      # folder with v66.p1_1240K.{geno,snp,ind,anno}
```

Everything else (panel prefixes, QC thresholds) is in the same file; nothing is
hard-coded in the code.

## Full pipeline (one command)

```bash
pip install -e .                            # installable package (pyproject.toml)
python run_pipeline.py --panel 1240k        # Phases 2→9 + HTML report, in order
python -m pytest tests/ -q                  # unit tests (math, simulation, kinship, config)
```

## Validation & robustness (for reuse / publication)

```bash
python validate_simulation.py               # ground-truth recovery -> SIMULATION_VALIDATION.md
Rscript tools/admixtools_concordance.R <prefix> results/admixtools_results.csv
python tools/compare_admixtools.py          # concordance vs ADMIXTOOLS 2
```

- **Simulation (msprime)** — recovers a *known* Neanderthal fraction (calibration
  0.97·true + 0.01pp, r=0.99); confirms SE ∝ 1/√SNPs, the detector's null
  false-positive rate, and its power. The ground-truth backbone for a methods paper.
- **ADMIXTOOLS 2 concordance** — the same f-statistics computed by the field-standard
  tool on the same files; `compare_admixtools.py` reports the differences.
- **Kinship** (`archaic.kinship.prune`) — READ-style relatedness/duplicate removal
  before group estimates (e.g. flags the Tarquinia necropolis family clusters among
  the Etruscans: 75 → 69 independent, 1 duplicate + 6 first-degree).
- **Release** — `pip`-installable, `CITATION.cff`, `.zenodo.json`, CI; see `RELEASING.md`
  for minting a Zenodo DOI.

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

```bash
python make_pdf.py PAPER.md reports/Etruscan_paper.pdf   # manuscript -> styled PDF
python fads_report.py        # FADS deep-dive   -> FADS_REPORT.md + reports/fads_report.html
python mutation_diversity.py # worldwide diversity -> DIVERSITY_REPORT.md + reports/diversity_report.html
```

- **FADS deep-dive** — the candidate selection locus, traced SNP-by-SNP (rs174537/rs174550) and
  shown to be the pan-European post-Neolithic FADS dietary sweep (archaic allele ~90% in
  Mesolithic HGs, declining; p≈0 across 5,936 Europeans), not adaptive introgression.
- **Worldwide diversity** — heterozygosity across 20 populations (the measurable proxy for
  "mutation rate"; rate per se needs trios). Africans highest, Native Americans lowest, the
  serial-founder gradient (r=−0.74), and a 1240K-vs-Human-Origins cross-check that *demonstrates*
  SNP-ascertainment bias (the San jump from #17 to #4).

## Sub-study: >5% archaic-ancestry survey (Eurasia + whole AADR)

```bash
python high_archaic_survey.py     # Eurasia (15,443 ancients) -> reports/high_archaic_survey/

python phase2_prepare.py --panel 1240k --scope global       # keep present-day + all continents
python phase3_estimate.py --panel 1240k \
       --meta results/phase2_1240k_global_metadata.csv \
       --out  results/phase3_1240k_global_estimates.csv     # only the new non-Eurasian genomes
python global_archaic_survey.py   # whole AADR (21,109 genomes) -> reports/global_archaic_survey/
```

Asks a simple absolute-threshold question distinct from Phase 6's conditioned outlier
search: who, across the AADR, carries **>5% archaic ancestry**? **Eurasia** (15,443
ancients): 110 raw crossings but **0 of 9,252 high-confidence genomes** — every crossing
is a low-coverage artifact (median 21,604 SNPs vs 272,044 panel-wide); the one individual
whose 95% CI excludes 5% is **Oase1** (9.8%, a documented recent-Neanderthal-ancestor
case, see below). Neanderthal ancestry is elevated (~2.5–3.2%) in the earliest Upper
Paleolithic (Bacho Kiro, Goyet, Muierii2, Ust'-Ishim), settling to the ~2.1% Holocene
baseline. A relative Denisovan D-statistic recovers a coherent Island Southeast Asian
cluster (3 Indonesian genomes).

Extending to the **whole AADR** (21,109 genomes, all continents, present-day included) via
a new `phase2_prepare.py --scope global` mode (keeps present-day + every continent, tags
`continent`/`is_modern`; `phase3_estimate.py --meta` lets Phase 3 run against that
metadata, reusing already-computed Eurasian estimates and only estimating the newcomers)
adds two things: an out-of-Africa **negative control** — mean α=0.29%, the most-negative
Denisovan affinity of any continent, on data never used to calibrate the estimator — and
the global **archaic maximum**: **Oceania** (highest Neanderthal at 2.3% *and* 31.8% of
genomes flagged Denisovan, up to Z=7.2 in Papuans). Still 0 high-confidence genomes exceed
5% Neanderthal on any continent, but Oceanians are the only humans whose *combined*
Neanderthal+Denisovan ancestry robustly clears 5%.

Full write-ups (paper-style: Abstract → Conclusion, reproducible Methods, figures) and
rendered PDFs:
- **`reports/high_archaic_survey/PAPER.md`** / `Eurasia_high_archaic_survey.pdf`
- **`reports/global_archaic_survey/PAPER_global.md`** / `Global_archaic_survey.pdf`

## Sub-study: Oase1 haplotype segment analysis (+ BAM pipeline)

```bash
python oase1_haplotype.py                    # array-resolution segments -> reports/oase1_haplotype/
bash oase1_bam_pipeline/run_oase1_hmmix.sh   # full read-level hmmix pipeline (Linux/macOS/WSL)
```

Oase1's genome-wide Neanderthal proportion (9.8% ± 2.2%) cannot by itself distinguish a
recent Neanderthal ancestor from an unusual background — the discriminator is *segment
length*. Using Oase1's 1240K genotypes (both AADR libraries merged), the in-panel
Altai/Vindija/African references, and the panel's own genetic map, a 2-state Viterbi HMM
finds Oase1's Neanderthal ancestry sits in **9 segments totalling 122 cM, longest 36.5 cM
(chr5) and 30.4 cM (chr9)** — no penecontemporaneous ~2% Eurasian (Ust'-Ishim, Kostenki14,
Loschbour) has a single segment over 18.5 cM. The longest-segment length implies a
Neanderthal ancestor within ~3 generations (an upper bound, since array density fragments
true blocks), consistent with Fu et al.'s (2015) published estimate of 4–6 generations
from the original shotgun genome.

A complete, runnable **read-level pipeline** (`oase1_bam_pipeline/`, using `hmmix` — Skov
et al. 2018) is provided for full-resolution confirmation: it downloads the Oase1 reads
from ENA **PRJEB8987**, calls variants against hg19, and runs `create_ingroup` → `train`
(haploid) → `decode -admixpop` to segment and annotate each block's archaic source
(Neanderthal vs Denisovan). Requires `samtools`/`bcftools`/`hmmix` on Linux/macOS/WSL —
not available on this project's Windows dev box, hence the separate script rather than a
direct run.

Full write-up: **`reports/oase1_haplotype/PAPER_oase1.md`** / `Oase1_haplotype_analysis.pdf`.

## Rendering a paper to PDF

`make_pdf.py` is a generic Markdown→PDF renderer (reportlab Platypus): headings, tables,
bullet lists, fenced code blocks, horizontal rules, markdown/bare hyperlinks, and inline
`![caption](figure.png)` images (resolved relative to the source `.md`'s own directory).

```bash
python make_pdf.py reports/high_archaic_survey/PAPER.md reports/high_archaic_survey/Eurasia_high_archaic_survey.pdf
python make_pdf.py reports/global_archaic_survey/PAPER_global.md reports/global_archaic_survey/Global_archaic_survey.pdf
python make_pdf.py reports/oase1_haplotype/PAPER_oase1.md reports/oase1_haplotype/Oase1_haplotype_analysis.pdf
python make_pdf.py                                        # legacy default: PAPER.md -> reports/Etruscan_paper.pdf
```

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
  simulate.py         msprime coalescent simulation with known archaic introgression
  kinship.py          READ-style relatedness / duplicate detection and pruning
  config.py           load paths + thresholds from config.yaml (nothing hard-coded)
phase1_validate.py … phase9_robustness.py   the nine pipeline stages
etruscan_study.py / etruscan_paper.py        Etruscan sub-study + manuscript (PAPER.md)
validate_simulation.py                       ground-truth simulation validation
tools/admixtools_concordance.R + compare_admixtools.py   ADMIXTOOLS 2 concordance
config.yaml · pyproject.toml · CITATION.cff · .zenodo.json · RELEASING.md   reuse/release
validate_published.py external validation vs the literature  -> VALIDATION.md
generate_report.py    self-contained HTML executive summary  -> reports/
etruscan_study.py     focused Etruscan sub-study              -> reports/etruscan_report.html
run_pipeline.py       orchestrator (Phases 2→9 + report)
probe_eastasian.py    methods probe (differential-statistic outgroup choice)
high_archaic_survey.py    >5% archaic survey, Eurasia    -> reports/high_archaic_survey/
global_archaic_survey.py  >5% archaic survey, whole AADR -> reports/global_archaic_survey/
  (needs phase2_prepare.py --scope global + phase3_estimate.py --meta, see above)
oase1_haplotype.py         Oase1 array-resolution segment/karyogram analysis -> reports/oase1_haplotype/
oase1_bam_pipeline/        runnable read-level hmmix pipeline (Linux/macOS/WSL; see its README.md)
make_pdf.py                generic Markdown -> PDF renderer (any PAPER*.md, not just the Etruscan one)
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
