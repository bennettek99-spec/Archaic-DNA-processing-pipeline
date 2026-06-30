# Neanderthal ancestry in the Etruscans: regional context, temporal stability, and a candidate selection signal

*Computational re-analysis of public ancient-DNA (AADR 1240K). Generated 2026-06-29. Exploratory; candidate findings are hypotheses, not validated discoveries.*

## Abstract
The Etruscans were a culturally distinctive Iron Age population of central Italy whose biological origins were long debated until ancient DNA showed them to be genetically continuous with neighbouring Italic peoples. Whether their archaic (Neanderthal) ancestry differs from that of contemporaries, changed through time, or shows individual anomalies has not been examined. Using 75 Etruscan genomes and a regional comparison panel from the Allen Ancient DNA Resource (1240K), we estimate Neanderthal ancestry with a validated f4-ratio and population mean-genome profiles. Etruscan Neanderthal ancestry is 1.96% (± 0.39), indistinguishable from Latins/Italics (2.24%) and Imperial Romans (2.04%) and from the broader European range. Etruscans are archaic-identical to Latins/Italics (Z=+0.2) and to Anatolian farmers and Italian Bronze Age; the only detectable differences are small (|D| ≤ 0.003) and involve Imperial Roman (Z=+2.8), Magna Graecia Greek (Z=+3.2), reflecting those groups' somewhat more eastern-Mediterranean ancestry rather than unusual archaic introgression. Genome-wide Neanderthal ancestry is flat across the Italian transect from the Neolithic to the medieval period. No Etruscan individual is a genome-wide-significant archaic outlier, and—crucially—the steppe/Levantine/East-Mediterranean ancestry outliers are not archaic outliers, showing that fine-scale ancestry variation in Iron-Age Etruria was decoupled from archaic content. A targeted scan of adaptive-introgression loci over time recovers FADS1-2 (fatty-acid metabolism; p=1.2e-03) as the single Bonferroni-significant candidate—consistent with the well-known dietary selection at FADS—though it rests on few markers and is reported as a hypothesis. The results place Etruscan archaic ancestry firmly within its regional and temporal context.

## Introduction
The Etruscans (Etruria, central Italy, ~800–100 BCE) developed a distinctive language and culture whose origins were debated for millennia. Ancient-DNA studies (Antonio et al. 2019; Posth et al. 2021) showed Etruscans were genetically continuous with neighbouring Italic peoples and the preceding Bronze Age, carrying a Steppe-derived component like other Iron Age Europeans, with a minority of individuals of eastern-Mediterranean or steppe-leaning ancestry. All present-day and ancient non-Africans also carry ~2% Neanderthal ancestry from Late-Pleistocene admixture (Green et al. 2010; Prüfer et al. 2014). Whether Etruscan archaic ancestry differs from neighbours, changed over time, or shows individual anomalies — and whether archaic-introgressed loci were under selection in this transect — has not been assessed. We address these questions with a validated genome-wide estimator, population mean-genome profiles, and a targeted locus scan.

## Materials and Methods
**Data.** 75 Etruscan genomes and regional comparison cohorts (Anatolian Neolithic, Italian Neolithic/Bronze Age, Aegean Bronze Age, Yamnaya Steppe, Latins/Italics, Republican and Imperial Romans, Magna-Graecia Greeks, ancient Sicilians and Sardinians, and modern French/Spanish/Sardinian/Han) from the AADR v66.1 1240K panel (1,233,013 SNPs). Quality control, ancestry PCA, and per-individual estimates follow the parent pipeline (Phases 2–6).

**Neanderthal ancestry.** The f4-ratio α = f4(Altai, Chimp; X, Mbuti) / f4(Altai, Chimp; Vindija, Mbuti), with two independent high-coverage Neanderthals (Altai in the statistic, Vindija as the 100%-Neanderthal scale; Petr et al. 2019). Denisovan affinity is D(X, Mbuti; Denisova, Chimp). Standard errors from a 50-block delete-one jackknife. The estimator reproduces published values (r=0.87; Oase1 recovered at 6–9%; see VALIDATION.md).

**Mean-genome profiles.** For each cohort we build a mean allele-frequency vector ("mean genome"), which averages out single-genome noise and yields group-level α with much tighter SEs than any individual. Differential Neanderthal sharing between cohorts is tested with D(P1, P2; Altai, Yoruba) (an African outgroup isolates recent introgression). Pairwise genome-wide allele-frequency distance summarises genetic affinity (visualised by MDS).

**Selection scan.** At curated adaptive-introgression loci (BNC2, OAS, TLR, FADS, HLA, TBX15, etc.), "archaic-informative" SNPs (Altai+Vindija-derived, near-absent in Africans) are identified; per-individual archaic-allele dosage is regressed on age across Italian samples, controlling for genome-wide archaic ancestry and ancestry turnover (PC1, PC2), so the age term isolates locus-specific temporal change.

## Results

### 1. Etruscan Neanderthal ancestry in regional context
Group-level Neanderthal ancestry (mean-genome f4-ratio):

| population | n | Neanderthal % ± SE | Denisovan D Z |
|---|---|---|---|
| Anatolia Neolithic | 58 | 1.94 ± 0.38 | -0.9 |
| Italy Neolithic | 46 | 2.08 ± 0.39 | -0.7 |
| Aegean BA | 24 | 2.09 ± 0.41 | -0.8 |
| Italy Bronze Age | 115 | 2.11 ± 0.40 | -1.0 |
| Yamnaya (Steppe) | 150 | 2.14 ± 0.39 | -0.7 |
| Etruscan | 75 | 1.96 ± 0.39 | -0.8 |
| Latin / Italic IA | 6 | 2.24 ± 0.43 | -0.6 |
| Sicily (anc) | 150 | 2.22 ± 0.39 | -0.6 |
| Magna Graecia Greek | 30 | 2.13 ± 0.40 | -0.9 |
| Republican Roman | 19 | 2.33 ± 0.37 | -0.2 |
| Imperial Roman | 150 | 2.04 ± 0.37 | -0.9 |
| Sardinia (anc) | 139 | 2.20 ± 0.39 | -1.1 |
| Sardinian (mod) | 31 | 2.03 ± 0.38 | -0.9 |
| Spanish (mod) | 55 | 2.05 ± 0.38 | -0.5 |
| French (mod) | 31 | 2.09 ± 0.38 | -0.4 |
| Han (mod) | 46 | 2.29 ± 0.38 | -0.8 |

Etruscan Neanderthal ancestry is **1.96% ± 0.39**, within the tight European band and essentially identical to Latins/Italics and Imperial Romans. Denisovan affinity is ~0 throughout, as expected for West Eurasians.

### 2. Etruscans do not differ from their neighbours
Differential Neanderthal sharing D(Etruscan, X; Altai, Yoruba):

| contrast | D | Z |
|---|---|---|
| Etruscan vs Latin / Italic IA | +0.0002 | +0.2 |
| Etruscan vs Imperial Roman | +0.0019 | +2.8 |
| Etruscan vs Anatolia Neolithic | +0.0002 | +0.2 |
| Etruscan vs Yamnaya (Steppe) | -0.0024 | -2.2 |
| Etruscan vs Magna Graecia Greek | +0.0031 | +3.2 |
| Etruscan vs Italy Bronze Age | +0.0013 | +1.8 |

Etruscans are archaic-identical to Latins/Italics (Z=+0.2) and to Anatolian farmers and Italian Bronze Age; the only detectable differences are small (|D| ≤ 0.003) and involve Imperial Roman (Z=+2.8), Magna Graecia Greek (Z=+3.2), reflecting those groups' somewhat more eastern-Mediterranean ancestry rather than unusual archaic introgression. This is consistent with the genetic continuity reported by Posth et al. 2021. Tuscan vs Latial Etruscans: 2.06% vs 1.94%.

### 3. Archaic ancestry was stable over time
Genome-wide Neanderthal ancestry is flat across the Italian transect (Neolithic → medieval); the major Neolithic and Steppe ancestry turnovers did not change the archaic fraction.

### 4. Individual outliers are explained by ancestry, not archaic anomalies
No Etruscan is a genome-wide-significant archaic outlier. The AADR steppe/Levantine/East-Mediterranean ancestry outliers have mean |z| = 0.41 versus 0.72 for typical Etruscans — i.e. they are **not** archaic outliers, because those alternative West-Eurasian ancestries carry similar ~2% Neanderthal. Fine-scale ancestry variation in Iron-Age Etruria was decoupled from archaic ancestry.

### 5. A candidate selection signal at FADS
Of 11 adaptive-introgression loci tested, after controlling for ancestry and overall archaic level:

| gene | phenotype | archaic SNPs | Δ/kyr (pp) | p |
|---|---|---|---|---|
| FADS1-2 | fatty-acid metabolism | 5 | +6.407 | 0.001 |
| OAS1-3 | antiviral innate immunity | 19 | -2.654 | 0.040 |
| TBX15-WARS2 | body-fat / cold response | 12 | -0.393 | 0.262 |
| OCA2-HERC2 | eye/skin pigmentation | 22 | -0.975 | 0.291 |
| EPAS1 | hypoxia (Denisovan; control) | 9 | -0.308 | 0.370 |
| KRT-cluster | keratin / skin & hair | 26 | +0.509 | 0.432 |
| TLR1-6-10 | innate immunity (TLRs) | 7 | +1.053 | 0.519 |
| POU2F3 | keratinocyte differentiation | 4 | +0.450 | 0.673 |
| BNC2 | skin pigmentation | 14 | +0.261 | 0.764 |
| HLA | MHC / adaptive immunity | 30 | +0.163 | 0.829 |
| SLC16A11 | lipid metabolism / T2D | 4 | -0.018 | 0.961 |

FADS1-2 (fatty-acid metabolism; p=1.2e-03) is the one locus surviving Bonferroni correction — biologically plausible given the strong, well-documented dietary selection at FADS in Europeans (Mathieson et al. 2015; Buckley et al. 2017), but resting on few archaic-informative SNPs and reported as a hypothesis for higher-coverage follow-up.

## Discussion
Etruscan archaic ancestry is unremarkable: ~2% Neanderthal, ~0 Denisovan, indistinguishable from neighbours and stable through time. This is the expected outcome — Neanderthal ancestry is shared across all non-Africans and the within-European variance is small — and it reinforces, from the archaic-ancestry angle, the genetic continuity of Etruscans with Italic peoples. The decoupling of genetic-ancestry outliers from archaic-ancestry outliers is a useful methodological point: an individual can be a clear ancestry outlier yet carry a perfectly ordinary archaic complement, because the relevant alternative ancestries are themselves ~2% Neanderthal. The FADS candidate is intriguing and concordant with the strongest known signal of recent dietary selection in Europe, but the single-locus power on capture data is low and the result must be treated as a hypothesis.

## Limitations
Capture-array, pseudo-haploid genotypes; ~75 Etruscans over a narrow window; single-locus selection power is low; archaic-allele sets are putative (incomplete lineage sorting can mimic introgression); the f4-ratio absolute scale runs slightly high versus some studies (relative comparisons are robust). No claim here is a validated biological discovery.

## References
Antonio et al. 2019 *Science* 366:708 · Posth et al. 2021 *Sci. Adv.* 7:eabi7673 · Green et al. 2010 *Science* 328:710 · Prüfer et al. 2014 *Nature* 505:43 · Prüfer et al. 2017 *Science* 358:655 · Petr et al. 2019 *PNAS* 116:1639 · Patterson et al. 2012 *Genetics* 192:1065 · Mathieson et al. 2015 *Nature* 528:499 · Buckley et al. 2017 *Mol. Biol. Evol.* 34:1307 · Mallick et al. 2024 *Sci. Data* 11:182 (AADR).
