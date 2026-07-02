# Archaic (Neanderthal and Denisovan) introgression in a single present-day genome: a validated-pipeline case study

*Bennett Kuhn — Archaic-DNA Introgression Pipeline (v0.3.0). Analysis run 2026-07-01 on the AADR v66.p1 1240K reference panel.*

## Abstract

Every person of non-African descent carries a small residue of DNA inherited from Neanderthals, and some populations additionally carry Denisovan DNA — the genetic legacy of interbreeding events as our ancestors expanded out of Africa ~50,000 years ago. This study estimates that legacy in a single present-day personal genome (a MyHeritage low-pass whole-genome-sequencing dataset, 560,853 autosomal sites, build37), using the same archaic reference genomes (Altai and Vindija Neanderthals; the Denisova genome; a chimpanzee outgroup) and the same allele-frequency machinery that the pipeline applies to ancient genomes from the Allen Ancient DNA Resource (AADR). We first show that the pipeline's standard ancient-genome estimator — the f4-ratio against the archaic genomes — **fails on consumer-array data**: on the sites a consumer array shares with the archaic references, populations *known* to be ~2% Neanderthal (French, English) read approximately zero, because the array under-samples archaic-informative sites. We therefore use the statistic that is valid for sparse array data — an **archaic-allele match rate** over curated Neanderthal-marker alleles (derived, carried by both Neanderthals, and essentially absent in Africans), internally calibrated to a percentage against reference populations of known archaic ancestry measured on the identical SNP set. The result: **Neanderthal ancestry = 2.25% ± 0.17%**, statistically indistinguishable from the West-Eurasian reference mean (2.08%) and stable (2.0–2.6%) across eight independent marker definitions. Denisovan ancestry is **negligible** (the only reference population showing a Denisovan excess is Papuan), consistent with the subject's West-Eurasian ancestry. **Total archaic introgression ≈ 2.3%.** The estimate is a hypothesis about one genome, not a discovery; its main value is as a worked, honestly-bounded demonstration of how much — and how little — a single consumer genome can say about archaic ancestry.

## 1. Introduction

Anatomically modern humans expanding out of Africa met and interbred with at least two archaic human groups: Neanderthals, across western Eurasia, and Denisovans, further east. As a result, present-day people outside sub-Saharan Africa carry roughly 1.5–2.5% Neanderthal ancestry (Green et al. 2010; Prüfer et al. 2014, 2017), and people from Oceania and parts of East and South Asia additionally carry up to ~3–5% Denisovan ancestry (Reich et al. 2010; Meyer et al. 2012). Sub-Saharan Africans carry little or none, providing the natural zero-point against which the introgression signal is measured.

The quantity has become familiar to the public because direct-to-consumer testing companies report a "Neanderthal" number. But estimating archaic ancestry rigorously from a *single* consumer genome is harder than the marketing implies, for two reasons. First, the archaic reference genomes are covered at only a subset of sites, and consumer arrays (or low-pass imputation to a modern reference panel) are deliberately ascertained toward common, well-behaved human polymorphisms — precisely *not* the sites that are most informative about archaic ancestry. Second, a single diploid genome yields per-site allele "frequencies" of only 0, ½ or 1, which are far noisier than the smooth frequencies of a reference population.

This study asks a concrete question — *how much Neanderthal and Denisovan DNA does one specific person carry?* — and answers it with the pipeline's validated toolkit, while being explicit about which statistics work on consumer data and which do not. The subject (hereafter "Bennett") is the pipeline's author; his MyHeritage-estimated ancestry is ~62% British Isles, ~17% Italian, ~14% Germanic/Dutch/Danish and ~7% French — entirely West-Eurasian, which sets a clear prior expectation (~2% Neanderthal, ~0 Denisovan) against which the analysis can be checked.

## 2. Materials and methods

### 2.1 Data

**Personal genome.** MyHeritage raw DNA export (format MHv1.0; method: low-pass whole-genome sequencing; reference build37/hg19), 560,853 genotyped autosomal sites reported on the forward strand.

**Reference panel.** AADR v66.p1 "1240K" (23,089 individuals × 1,233,013 SNPs, build37). In-panel archaic references: `AltaiNeanderthal.DG`, `VindijaG1_final.SG`, `Denisova.SG`; outgroup `Chimp.REF`; African baselines `Mbuti`, `Yoruba`. Present-day context populations: Orcadian, English, French, Basque, Spanish, Italian_North, Sardinian (West-Eurasian); Han (East Asian); Papuan (Oceanian, the positive control for Denisovan); Yoruba and Mbuti (African, the zero-point).

### 2.2 Aligning a consumer genome to the panel

Each panel SNP is defined by two alleles (`a1`, `a2`); the panel counts the `a1` allele, and all f-statistics used here are invariant to which allele is counted as long as the coding is *consistent across populations at each site* (Patterson et al. 2012). We therefore matched the personal genome to the panel by **(chromosome, position)** — robust to differing SNP-name schemes — and converted each diploid call into the dosage (0, 1 or 2 copies) of that site's `a1` allele, i.e. a frequency in {0, ½, 1}, exactly as the panel treats a single-genome reference. Calls whose alleles were consistent with the panel site only after reverse-complementing were flipped; **strand-ambiguous sites (a1/a2 = A/T or C/G) were dropped** because their strand cannot be resolved by allele identity; triallelic mismatches and no-calls were dropped.

Of 560,853 consumer sites, **164,427 matched a 1240K autosomal SNP by position** (the 1240K captures a specific SNP set that only partly overlaps consumer arrays), 503 were dropped as strand-ambiguous and none as mismatches, leaving **163,924 usable sites**. As an internal check that the allele coding was consistent with the panel, the personal genome's per-site `a1` dosage correlates **r = +0.639** with French `a1` frequency across shared sites — strongly positive, as required for a West-Eurasian genome coded the same way as the panel (a coding error would invert this to a large negative correlation).

### 2.3 Why the ancient-genome estimator fails here, and what replaces it

For ancient shotgun genomes the pipeline estimates Neanderthal ancestry with the **f4-ratio**

  α(X) = f4(Altai, Chimp; X, Mbuti) / f4(Altai, Chimp; Vindija, Mbuti),

which uses a *second* Neanderthal (Vindija) to set the "100% Neanderthal" scale (Petr et al. 2019). This statistic requires Vindija to be covered at the same sites and draws its power from tens of thousands of archaic-informative SNPs. A consumer array shares too few: only **14.5%** of the strongly archaic-informative sites (where Altai differs sharply from Chimp) survive the overlap. The consequence is decisive — restricted to the personal genome's SNP set, the f4-ratio reads **−0.29% for French and −0.41% for English**, populations known to carry ~2% Neanderthal. The f4-ratio is not *wrong*, it is simply *underpowered and biased toward zero* on array data, and any single-genome number it produces (here 0.03% ± 0.78%) is uninformative.

We therefore use the statistic that is valid for sparse array data — an **archaic-allele match rate**. This is the same idea underlying tag-SNP counts used by consumer companies and by Vernot & Akey (2014), but built from the in-panel archaic genomes rather than a published introgression map:

1. **Define Neanderthal-marker alleles** on the panel: sites where the *derived* allele (differing from the chimpanzee ancestral state) is carried by **both** the Altai and Vindija Neanderthals and is **essentially absent in Africans** (frequency < 5% in Yoruba and Mbuti). Such alleles are the classic signature of introgression: archaic-derived, African-rare, and therefore present in living non-Africans predominantly because of admixture. The panel yields 5,303 such sites; **394 fall on the personal array** (primary threshold).

2. **Measure the match rate** = the target genome's mean frequency of the archaic-derived allele over these marker sites, with a 50-block delete-one genome jackknife for the standard error.

3. **Calibrate to a percentage.** Because the marker set is a fixed ascertainment, the match rate is linearly related to the introgression fraction. We fit `known α ~ match rate` across the reference populations (whose true Neanderthal % is the full-panel f4-ratio) on the *identical* marker sites, and read the personal genome off that line. Papuan is excluded from the calibration (its Denisovan ancestry perturbs Neanderthal-marker sharing); the African references are retained because they anchor the zero-point. The **entire calibrated prediction is block-jackknifed** (each block deletion re-fits the calibration and re-predicts), so the reported error propagates both the genome's sampling noise and the calibration uncertainty.

The same construction with Denisova-specific markers (derived, carried by the Denisovan but *not* by the Neanderthals, African-rare) gives a Denisovan match rate. Denisovan ancestry is reported as a *relative* signal (there is only one high-coverage Denisovan genome, so it is not calibrated to an absolute percentage).

## 3. Results

### 3.1 Neanderthal ancestry ≈ 2.25%

The personal genome carries the Neanderthal-marker alleles at a **match rate of 18.1% ± 1.5%**, within the West-Eurasian reference band (16.3–17.9%) and far above the African baseline (Yoruba 2.1%, Mbuti 0.9%). The internal calibration is excellent (r = 0.990; α = 0.138 × match − 0.246), and placing the personal genome on it gives:

> **Neanderthal ancestry = 2.25% ± 0.17%.**

This is statistically indistinguishable from the West-Eurasian reference mean of **2.08%** (Figure 1) — the subject sits at the upper edge of the reference cloud, but well within one standard error of it, as expected for a single noisy genome. It is also consistent, within error, with the ~1.5–2.1% reported for West Europeans in the literature (Prüfer et al. 2017).

![Figure 1. Neanderthal ancestry (calibrated archaic-allele match rate, with jackknife standard errors) for the personal genome (red) against present-day reference populations; grey = West-Eurasian references, dark = cross-continental references. The dashed line is the West-Eurasian mean (2.08%). Bennett (2.25%) is statistically indistinguishable from the West-Eurasian cloud and far above the African baseline (near zero).](fig_pg1_neanderthal.png)

Figure 2 shows the calibration itself: the reference populations fall on a tight line from the African zero-point (~2% match → ~0% ancestry) up through the West-Eurasian and East-Asian references (~16–17% match → ~2% ancestry), and the personal genome is placed on that line with its jackknife error bar.

![Figure 2. Internal calibration of the match-rate estimator. Each point is a reference population: x = archaic-allele match rate on the personal genome's marker sites, y = its Neanderthal ancestry (full-panel f4-ratio "truth"; Papuan, in purple, is shown but excluded from the fit). The dashed line is the calibration; the red star is the personal genome placed on it (2.25% ± 0.17%).](fig_pg2_calibration.png)

### 3.2 Reference comparison

| Population | n | match rate (%) | calibrated α (%) | full-panel α "truth" (%) |
|---|---|---|---|---|
| **Bennett (this study)** | 1 | **18.1 ± 1.5** | **2.25 ± 0.17** | — |
| Orcadian | 17 | 17.1 | 2.12 | 2.11 |
| English | 12 | 17.4 | 2.15 | 2.02 |
| French | 31 | 16.6 | 2.04 | 2.09 |
| Basque | 25 | 17.9 | 2.22 | 2.04 |
| Spanish | 55 | 16.3 | 2.00 | 2.05 |
| Italian (North) | 24 | 16.4 | 2.01 | 2.04 |
| Sardinian | 31 | 16.5 | 2.03 | 2.03 |
| Han | 46 | 16.8 | 2.07 | 2.29 |
| Papuan | 32 | 15.3 | (1.86) | 2.94 |
| Yoruba | 24 | 2.1 | 0.04 | −0.11 |
| Mbuti | 15 | 0.9 | −0.12 | 0.00 |

The calibrated column recovers each reference's known Neanderthal ancestry to within its error (Africans → ~0; West-Eurasians → ~2.0–2.2%), confirming the estimator is well-behaved. The one deliberate exception is Papuan, whose Neanderthal-marker match rate is depressed relative to its true ancestry (hence its exclusion from the calibration and the parenthetical value) — a known consequence of its additional Denisovan ancestry and deeper drift.

### 3.3 Denisovan ancestry is negligible

On Denisova-specific markers, the only reference population with an elevated match rate is **Papuan (19.4%)**, the textbook Denisovan carrier; all West-Eurasian references sit at a common baseline (~16.8%) that reflects background archaic sharing rather than Denisovan admixture. The personal genome's Denisovan match rate (14.1% ± 4.6%) is **at or below** that West-Eurasian baseline — i.e. **no Denisovan excess**. The direct Denisovan D-statistic on the personal genome agrees (D(Bennett, Yoruba; Denisova, Chimp) Z = 0.3, indistinguishable from zero). This is exactly the expectation for a West-Eurasian genome: Denisovan ancestry is essentially confined to Oceanian, and to a lesser extent East/South Asian, populations. **Total archaic introgression is therefore essentially all Neanderthal: ≈ 2.3%.**

### 3.4 Robustness

The estimate does not depend on the precise marker definition. Sweeping the African-rarity cutoff (2%, 5%, 10%, 15%) and the archaic-fixation threshold (0.95, 0.99) gives eight independent marker sets (241 to 1,094 sites) and eight estimates spanning **2.00%–2.63%** (mean 2.26%), each with calibration r ≥ 0.986:

| African-rare cutoff | archaic threshold | marker SNPs | calibration r | Neanderthal % |
|---|---|---|---|---|
| 2% | 0.95 / 0.99 | 241 | 0.992 | 2.63 |
| 5% (primary) | 0.95 / 0.99 | 394 | 0.990 | 2.25 |
| 10% | 0.95 / 0.99 | 697 | 0.989 | 2.14 |
| 15% | 0.95 / 0.99 | 1,094 | 0.986 | 2.00 |

The stricter, cleaner marker sets give slightly higher values (~2.6%) and the broadest, noisier sets slightly lower (~2.0%), bracketing the primary estimate. A per-chromosome breakdown (Figure 3) scatters around the genome-wide value as expected given only ~10–50 marker sites per chromosome; no single chromosome drives the result. Together these place the subject's Neanderthal ancestry robustly in the **~2.0–2.3%** range.

![Figure 3. Per-chromosome Neanderthal estimate for the personal genome (point size proportional to the number of marker SNPs on each chromosome). The dashed line and band are the genome-wide estimate ± SE. The per-chromosome scatter is sampling noise from the small marker counts, not real chromosomal structure.](fig_pg3_perchrom.png)

### 3.5 Methods comparison

For completeness, the ancient-genome statistics were also computed directly on the personal genome. The f4-ratio gives **0.03% ± 0.78%** (70,495 SNP) and the Neanderthal D-statistic D(Bennett, Yoruba; Altai, Chimp) gives **Z = 2.7** — a detectable but weak affinity, and a point estimate that is uninformative because, as shown in §2.3, these statistics collapse toward zero on the array's SNP ascertainment. Their disagreement with the match-rate estimate is not a contradiction: it is the quantitative demonstration that the f4-ratio is the wrong tool for a consumer genome, and that the calibrated match rate is the right one.

## 4. Discussion

The headline number — **~2.25% Neanderthal, ~0% Denisovan, ~2.3% archaic in total** — is unremarkable, and that is the point. It is exactly what a West-Eurasian genome should carry, it matches the reference populations that share the subject's ancestry, and it is stable to how the estimate is constructed. In absolute terms, roughly 2% of a ~3.1-gigabase autosomal genome corresponds to on the order of ~60 megabases of sequence tracing to Neanderthal ancestors — distributed as short introgressed haplotypes scattered across the chromosomes, the surviving fragments of admixture that happened some 50,000–60,000 years ago. Consistent with this, the subject's independent chromosome-painting project recovered ~1.2% Neanderthal on chromosome 1 alone; segment-painting and allele-frequency methods measure related but not identical quantities, and both land in the same low-single-digit-percent regime.

The scientifically interesting content of this study is methodological. It makes concrete a fact that is easy to overlook: **the statistic that is gold-standard for a high-coverage ancient genome is the wrong statistic for a consumer array.** The f4-ratio's power comes from the second archaic genome and from archaic-informative sites, both of which are largely absent from consumer data; forcing it onto that data returns ~0% even for people who are demonstrably ~2% Neanderthal. The match-rate estimator sidesteps both problems by (i) using only curated, high-confidence marker alleles and (ii) calibrating internally against populations of known ancestry on the identical sites — turning a biased measurement into an unbiased one. This is the same principle consumer companies use when they count "Neanderthal variants," made explicit and validated here against the AADR panel.

The estimate should be read as a **hypothesis about one genome**, in keeping with this pipeline's standing rule that findings are hypotheses until technical causes are excluded. Here the technical causes largely *have* been examined — array ascertainment (quantified and corrected), allele coding (checked), marker definition (swept), calibration (cross-validated against nine reference populations) — and the result survives them. What it cannot claim is fine resolution: a single consumer genome can confidently place a person in the "West-Eurasian, ~2% Neanderthal, no Denisovan" bin, but it cannot resolve whether that person carries 2.0% or 2.3% versus a neighbour, because the between-individual differences within a continental group (tenths of a percent) are smaller than a single genome's error bar.

## 5. Limitations

- **Single diploid genome.** Per-site "frequencies" are 0/½/1, so the sampling noise is intrinsically larger than for a reference population; the ±0.17% error bar reflects this and the ~2.0–2.6% robustness range is the honest envelope.
- **Array/imputation ascertainment.** Consumer data under-samples archaic-informative sites and (for low-pass imputation) is imputed against a modern reference panel, which can introduce mild reference bias toward common human alleles; the calibration corrects for the ascertainment but cannot fully undo imputation bias.
- **Calibrated, not absolute.** The Neanderthal percentage is anchored to the reference populations' f4-ratio values, which themselves carry a small (~+0.2 percentage-point) absolute-scale offset quantified elsewhere by simulation; relative comparisons are unbiased, the absolute decimal is approximate.
- **Denisovan is relative only.** With a single high-coverage Denisovan genome, the Denisovan signal is reported as an affinity, not a calibrated fraction; the safe conclusion is "no detectable Denisovan excess," not a hard 0.0%.
- **No haplotype phasing or segment maps.** This study measures a genome-wide fraction, not the physical locations, lengths, or functional content of introgressed segments; those require phased sequence data and a segment caller (e.g. hmmix / IBDmix), which is the natural follow-up.

## 6. Conclusion

Applying the pipeline's validated archaic references to a single present-day personal genome, and using the estimator appropriate to sparse consumer data, the subject carries **~2.25% ± 0.17% Neanderthal ancestry and no detectable Denisovan ancestry — a total archaic introgression of about 2.3%.** The value is typical for West-Eurasian ancestry, robust to methodological choices, and internally calibrated against reference populations of known archaic content. Beyond the personal number, the study documents a transferable lesson for anyone estimating archaic ancestry from consumer genotypes: the f4-ratio that works for ancient genomes silently fails on arrays, and a calibrated archaic-allele match rate is the statistic that recovers the right answer.

## Data and code availability

Pipeline and this analysis: `personal_genome_study.py` and `archaic/consumer_dna.py` in the Archaic-DNA Introgression Pipeline (github.com/bennettek99-spec/Archaic-DNA-processing-pipeline). Reference data: Allen Ancient DNA Resource v66.p1 (1240K), Harvard Dataverse. Personal genotype data are private and not redistributed. Results, tables and figures: `results/personal/` and `reports/personal_genome/`.

## Selected references

- Green R.E. et al. (2010) A draft sequence of the Neandertal genome. *Science* 328:710.
- Reich D. et al. (2010) Genetic history of an archaic hominin group from Denisova Cave. *Nature* 468:1053.
- Meyer M. et al. (2012) A high-coverage genome sequence from an archaic Denisovan individual. *Science* 338:222.
- Patterson N. et al. (2012) Ancient admixture in human history. *Genetics* 192:1065.
- Prüfer K. et al. (2014) The complete genome sequence of a Neanderthal from the Altai Mountains. *Nature* 505:43.
- Vernot B., Akey J.M. (2014) Resurrecting surviving Neandertal lineages from modern human genomes. *Science* 343:1017.
- Prüfer K. et al. (2017) A high-coverage Neandertal genome from Vindija Cave, Croatia. *Science* 358:655.
- Petr M. et al. (2019) Limits of long-term selection against Neandertal introgression. *PNAS* 116:1639.
- Mallick S. et al. (2024) The Allen Ancient DNA Resource (AADR). *Scientific Data* 11:182.
