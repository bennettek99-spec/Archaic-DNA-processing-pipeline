# Etruscan archaic-introgression study

*Exploratory; results are candidates, not proof of selection.*

## A. Genome-wide Neanderthal ancestry over time (Italy)
Weighted trend vs age: **+0.0081 pp / 1000 yr** (p=0.79). Per-bin means:

| era | n | Neanderthal % ± SE |
|---|---|---|
| Neolithic/Copper | 42 | 2.08 ± 0.11 |
| Bronze Age | 82 | 2.00 ± 0.07 |
| Iron Age (Etruscan/Latin) | 229 | 2.03 ± 0.05 |
| Roman | 144 | 1.97 ± 0.05 |
| Late Antique/Medieval | 116 | 1.99 ± 0.06 |

## B. Etruscan outliers & ancestry
Mean |within-Etruscan z|: AADR '-o' ancestry outliers **0.45** vs typical Etruscans **0.68**. Of 1 Etruscans with |z|>2 within the group, **0** fall below |z|=2 once genetic ancestry + geography are conditioned on — i.e. their apparent deviation is **explained by ancestry**, not anomalous archaic introgression.

Top within-Etruscan residuals:

| id | group | BP | Nea% | z(within) | z(ancestry-cond) | -o |
|---|---|---|---|---|---|---|
| ITTQ14.SG | Italy_Lazio_BA_IA_a_Etruscan | 2930 | 0.56 | -2.42 | -2.81 |  |
| EV7A.SG | Italy_Monteriggioni_IA_Etruscan | 2312 | 3.53 | +1.82 | +1.39 |  |
| EV16A.SG | Italy_Monteriggioni_IA_Etruscan | 2549 | 5.35 | +1.64 | +1.44 |  |
| TAQ016.AG | Italy_Lazio_IA_c_IA_d_Etruscan | 2151 | 2.70 | +1.46 | +1.44 |  |
| MAG001.AG | Italy_Tuscany_Etruscan | 2617 | 0.93 | -1.45 | -1.69 |  |
| TAQ018.AG | Italy_Lazio_IA_c_IA_d_Etruscan | 2151 | 1.17 | -1.39 | -1.61 |  |
| CSN002.AG | Italy_Tuscany_IA_b_IA_c_Etruscan-o | 2342 | 4.69 | +1.39 | +1.37 | YES |
| R473.SG | Italy_Lazio_IA_b_Etruscan | 2600 | 1.19 | -1.39 | -1.92 |  |
| R10338.SG | Italy_Lazio_IA_Etruscan | 2328 | 2.59 | +1.33 | +0.58 |  |
| TAQ004.AG | Italy_Lazio_IA_c_IA_d_Etruscan | 2151 | 2.64 | +1.31 | +1.12 |  |
| ITTQ10.SG | Italy_Lazio_TarquiniaCivita_IA_a_E | 2746 | 1.13 | -1.27 | -1.92 |  |
| TAQ007.AG | Italy_Lazio_IA_c_Etruscan | 2240 | 1.18 | -1.11 | -0.93 |  |

## C. Archaic alleles / genes over time
14 adaptive-introgression loci tested; **2** with nominal p<0.05 (Bonferroni threshold 0.0036). Age coefficient controls for genome-wide archaic ancestry, so it isolates locus-specific change.

| gene | phenotype | archaic SNPs | Δ/kyr | p | direction |
|---|---|---|---|---|---|
| FADS1-2 | fatty-acid metabolism | 5 | +7.409pp | 0.000 | falling toward present |
| OAS1-3 | antiviral innate immunity | 19 | -3.040pp | 0.021 | rising toward present |
| TBX15-WARS2 | body-fat / cold response | 12 | -0.378pp | 0.292 | rising toward present |
| OCA2-HERC2 | eye/skin pigmentation | 22 | -0.921pp | 0.328 | rising toward present |
| KRT-cluster | keratin / skin & hair | 26 | +0.634pp | 0.334 | falling toward present |
| TLR1-6-10 | innate immunity (TLRs) | 7 | +1.341pp | 0.421 | falling toward present |
| POU2F3 | keratinocyte differentiation | 4 | +0.847pp | 0.444 | falling toward present |
| EPAS1 | hypoxia (Denisovan; control) | 9 | -0.261pp | 0.459 | rising toward present |
| BNC2 | skin pigmentation | 14 | +0.171pp | 0.848 | falling toward present |
| HLA | MHC / adaptive immunity | 30 | +0.065pp | 0.932 | falling toward present |
| SLC16A11 | lipid metabolism / T2D | 4 | +0.017pp | 0.964 | falling toward present |

## D. Verdict
Genome-wide Neanderthal ancestry is essentially flat across the Italian transect (+0.0081 pp/kyr, p=0.79). No Etruscan individual is a significant archaic-ancestry outlier; crucially, the AADR steppe/Levantine/East-Mediterranean *genetic-ancestry* outliers are **not** archaic-ancestry outliers (mean |z| 0.45 vs 0.68 for typical Etruscans) — those alternative West-Eurasian ancestries carry similar ~2% Neanderthal, so the fine-scale ancestry variation in Iron-Age Etruria did **not** translate into archaic-ancestry differences. 1 of 11 introgression loci survive Bonferroni after controlling for ancestry (PCs) and overall archaic level: **FADS1-2** (p=2.6e-04, falling toward present, 5 archaic SNPs). These rest on few archaic-informative SNPs and are flagged as selection *candidates/hypotheses* for higher-coverage follow-up, not proof. FADS1-2 is a well-known target of strong dietary selection in Europeans (Mathieson 2015; Buckley 2017), so a temporal shift there is biologically plausible.

## E. ADMIXTOOLS 2 ancestry-model cross-check

*Authoritative tool cross-check; exploratory. The pure-Python qpAdm in this
pipeline is a simplified rotating-outgroup form and is not the publication
instrument — the ADMIXTOOLS 2 run below is.*

Four Italian target cohorts (Etruscan, Latin, ImperialRoman, ItalyBA) were
modelled with the canonical West-Eurasian sources (Anatolia_N, Yamnaya, WHG, and
with/without Iran_N) against seven distal outgroups (Mbuti, Han, Papuan,
Karitiana, Natufian, Ust_Ishim, MA1), using the same PLINK export and cached f2
statistics as the concordance validation (`tools/admixtools_concordance.R`).
Results are in `results/etruscan_admixtools/{qpwave,qpadm_popdrop}.csv`, produced
by `tools/etruscan_qpwave_qpadm.R`.

The outcome is a formal rejection at full SNP density, consistent with the
pipeline's documented handling of >1M-SNP qpAdm:

- **qpWave** rejects every tested rank (p ≤ 5e-05 for all four cohorts and both
  source sets): the left/right sets are connected by more ancestry streams than
  any 3- or 4-source model admits.
- **qpAdm popdrop ladder** rejects the full 3-way model and every reduced
  (drop-one-source) model for all four cohorts (all p < 0.05; most < 1e-6), and
  rejects the 4-way (Iran_N added) model and its drops as well.

The interpretation boundary therefore stands: the mixture *proportions* reported
by the pipeline's ancestry decomposition are descriptive summaries of a model
that is formally rejected at >1M SNPs, not an accepted ancestry model. They are
kept because they reproduce the well-established Steppe-migration signal with no
manual tuning, but they are not presented as a validated qpAdm fit.
