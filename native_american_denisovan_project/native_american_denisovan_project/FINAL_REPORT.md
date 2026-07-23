# Final Research Report — Native American Denisovan Ancestry Study

## Abstract

We tested whether Native Americans carry a Denisovan ancestry component more
closely related to the Papuan-associated Denisovan component than to the
principal East Asian component, a hypothesis suggested by an over-reading of
Figure 4 of Browning et al. 2018. Using 951 usable ancient American individuals
from the AADR v66.p1 1240K panel, 14 ancient Siberian/East Asian cohorts, and
present-day Papuan/Han/French controls, we computed block-jackknife
D-statistics (S1, S2) and a preregistered composition contrast (S3 = corr(X,
Papuan) - corr(X, Han) at 10,942 Denisovan-diagnostic Set-A SNPs, on a
held-out chromosome partition). **The hypothesis is not supported.** Every
ancient American and Siberian cohort has negative S3 (mean = -0.2310,
range [-0.3325, -0.1104], 0/14 positive), with
bootstrap CIs excluding zero — meaning their Denisovan-marker profiles are
**more Han-like than Papuan-like**, not the reverse. The signal is stable to
transversion-only analysis, leave-one-chromosome-out, and selected-locus
exclusion. Total Denisovan affinity (S1, S2) is null in ancient Americans
(|Z| < 1.5 for all cohorts). The result is consistent with Model 1 (a single
East Eurasian Denisovan pulse ancestral to both East Asians and Native
Americans) and contradicts Model 3 (a Papuan-related source contribution).

## Introduction

Denisovan ancestry is not monolithic. Browning et al. 2018 showed two
Denisovan components in East Asians differing in similarity to the Altai
Denisovan; Jacobs et al. 2019 found multiple deeply divergent Denisovan
ancestries in Papuans. Qin & Stoneking 2015 reported a very low level of
Denisovan ancestry across Eastern Eurasian and Native American (EE/NA)
populations, equally correlated with New Guinea and Australian ancestry. The
question of whether the Denisovan ancestry that reached the Americas via
East-Asian-related ancestors resembles the Papuan-associated (moderate
affinity) or East-Asian-associated (high affinity) component has been left
open by an over-reading of Browning Figure 4, which plots admixed 1000 Genomes
American populations beside Papuans but whose authors attribute the American
signal to admixture/LD artifacts and East-Asian-related Native ancestors (see
`docs/browning_figure4_interpretation.md`). Native Americans are informative
because they sample a different branch of the East Eurasian expansion than
present-day East Asians, and because ancient genomes provide a time transect
back to the Beringian period.

## Interpretation of Browning et al. 2018

Figure 4 of Browning et al. 2018 shows contour density plots of
introgressed-segment match rates to the Altai Neanderthal and Altai Denisovan
for each 1000 Genomes population plus SGDP Papuans. The American panels
(PUR/CLM/MXL/PEL) show lower match rates, which the authors attribute to
"admixture and thus higher background levels of LD that could cause false
positive results." The formal two-component test (Table 2 of that paper) was
**not significant in any American population**. The figure does not show that
Native Americans carry a Papuan-like Denisovan component; it shows admixed
Americans have some Denisovan-like segments, attributed to East-Asian-related
Native ancestors. Full analysis in `docs/browning_figure4_interpretation.md`.

## Materials and Data

- **AADR v66.p1 1240K panel**: 23,089 individuals; 951 usable ancient American
  (SNP floor 30,000); 14 ancient Siberian/East Asian cohorts; 8 present-day
  control populations (Papuan n=32, Han n=46, French n=31, Mbuti n=15, Yoruba
  n=127, Karitiana n=16, Japanese n=31, Dai n=14).
- **Archaic references**: Denisova.SG (574k SNPs), AltaiNeanderthal.DG,
  VindijaG1_final.SG, Chimp.REF.
- **Diagnostic SNPs**: Set A (strict: Denisova high >=0.90, African <=0.10,
  Neanderthal <=0.10) = 10,942 autosomal SNPs. Set B (Neanderthal <=0.50) =
  12,180. Transversions-only Set A = 1,738.
- **Training/validation partition**: odd chromosomes (1,3,...,21) define
  Papuan-enriched/Han-enriched marker subsets; even chromosomes (2,4,...,22)
  are the held-out validation set (5,520 Set-A sites).

## Ethics and Data Governance

See `ETHICS_AND_DATA_USE.md`. AADR data are public for secondary analysis;
no individual-level genotypes are redistributed; no claim of tribal identity or
cultural descent is made; ancient genomes are treated as published data points.

## Methods

1. **S1** = D(X, Mbuti; Denisova, Chimp): basic Denisovan affinity, 50-block
   jackknife. Positive = more Denisovan sharing than the Mbuti baseline.
2. **S2** = D(X, French; Denisova, Chimp): denoised Denisovan excess,
   differencing against the non-Denisovan non-African background (French).
3. **S3** = corr(f_X, f_Papuan) - corr(f_X, f_Han) at Set-A validation sites:
   the preregistered composition contrast. Oriented to the Denisovan-high
   allele. Bootstrap CI by chromosome resampling (B=1000).
4. **Papuan-enriched / Han-enriched subsets**: defined on the training
   partition (sites where Papuan > Han + 0.05 or vice versa) and tested on
   the held-out validation partition to avoid circularity.
5. Sensitivity: transversion-only S3, selected-loci exclusion (EPAS1,
   MUC19), leave-one-chromosome-out, per-chromosome S3, alternative African
   outgroup (Yoruba).
6. Ancestry conditioning: regression of S3 on genome-wide Han similarity and
   ANE (Sib_LGM) similarity.
7. Simulation: planted Denisovan-component mixtures (alpha = fraction
   Han-like) at varying marker counts to calibrate power and bias.
8. **Sign gate** (`tests/test_dstat_sign.py`, 5/5 pass): all statistics
   validated against synthetic genotype matrices before use; one candidate
   statistic (D(Altai, Denisova; X, Mbuti)) excluded after failing the gate
   (it is a Neanderthal indicator, anti-correlated with Denisovan ancestry).

## Validation

- **Positive controls**: Papuan S1 Z=+5.92, S2 Z=+7.30, S3=+0.40
  [+0.36,+0.45] — correctly identified as Denisovan-bearing and Papuan-like.
- **Negative controls**: French S1 Z=-0.44, S2 undefined (self), S3=-0.18
  [-0.25,-0.13]; Mbuti S1 undefined, S3=-0.07 [-0.10,-0.04]; Yoruba
  S1 Z=-2.06, S3=-0.10 [-0.13,-0.07] — correctly near-null and Han-like
  (slightly more correlated with Han than Papuan, as expected for a
  non-Denisovan population).
- **Han control**: S3=-0.40 [-0.45,-0.36] — correctly identified as maximally
  Han-like.
- **Simulation**: S3 recovers the planted composition across marker counts
  (mean S3 at alpha=0 pure-Papuan = +0.40; at alpha=0.5 = +0.04; the method
  has power to distinguish Papuan-like from Han-like even at 100 markers).

## Results

### Total Denisovan affinity is null in Native Americans (S1, S2)

All ancient American and Beringian cohorts have S1 |Z| < 1.0 and S2 |Z| < 1.0
(all-sites). The only ancient cohorts with elevated S1/S2 are Sib_UP
(Upper Paleolithic Siberians — Yana, etc.) and TT_pre40k (IUP individuals
including Ust'Ishim), which reflect deep Siberian population structure, not
Denisovan ancestry per se. Papuan S1 Z=+5.92, S2 Z=+7.30 confirm the
positive control.

### The Denisovan-marker composition is Han-like, not Papuan-like (S3)

The primary finding: **all 14 ancient American/Siberian/East Asian cohorts
have negative S3** (mean = -0.2310, range [-0.3325,
-0.1104]). Bootstrap CIs exclude zero for every cohort. This means
their Denisovan-diagnostic allele-frequency profiles correlate **more with
Han than with Papuan**, the opposite of the hypothesis. The pattern is stable
across time bins (TT_pre40k through TT_post5k: S3 from -0.01 to -0.29), stable
to transversion-only analysis, and stable under leave-one-chromosome-out.

### Ancestry conditioning does not rescue the hypothesis

Regression of S3 on genome-wide Han similarity and ANE similarity gives
R2 = 0.36 with beta_Han = -0.83 (more Han-similar -> more negative S3), as
expected: populations that are genome-wide more Han-like are also
Denisovan-profile more Han-like. Conditioning on ancestry does not flip the
sign.

### Sensitivity

Transversion-only S3 is consistent with all-sites S3 in direction for all
cohorts (see fig14). Selected-loci exclusion (EPAS1, MUC19) removes 0 Set-A
sites (these regions are poorly covered on 1240K), so the result is
unaffected. Per-chromosome S3 is consistently negative across autosomes for
American cohorts, with no single-chromosome outlier driving the signal.
Leave-one-chromosome-out S3 is stable (no chromosome removal flips the sign).

## Alternative Explanations

1. **Differential dilution (Model 5)**: the regression R2 = 0.36 suggests
   some of the S3 variation is explained by genome-wide ancestry, but the
   sign remains negative after conditioning — dilution does not explain the
   *direction*.
2. **Population structure (Model 4)**: structured Denisovan sources would
   need tract-level phylogenies to resolve (deferred to whole-genome tier).
3. **Selection (Model 6)**: EPAS1 and MUC19 are not covered on 1240K Set-A;
   the signal is genome-wide, not driven by known selected loci.
4. **Ascertainment bias**: Set A is defined from Denisova/Neanderthal/African
   references only (not Papuans or Han), so the S3 contrast is not circular.
   Training/validation chromosome partition further guards against overfitting.

## Demographic Implications

The result is consistent with **Model 1** (a single ancestral East Eurasian
Denisovan pulse before the East Asian-Native American divergence): Native
Americans inherited the same Denisovan component as East Asians, and that
component is more similar to the Han Denisovan profile than to the Papuan
Denisovan profile. This implies the Denisovan ancestry that reached the
Americas did so via the East Asian-related ancestors of First Americans,
not via a separate Papuan-related source. It does not distinguish whether
the East Asian Denisovan component itself has substructure (the Browning
high- vs moderate-affinity split), which would require whole-genome tract
analysis.

## Limitations

- **Total Denisovan ancestry in Native Americans is below detection** by
  D-statistics (S1, S2 null). The S3 composition test is powered because it
  compares *which* Denisovan markers are shared, not *how many*.
- **Sib_UP and TT_pre40k have inflated S1/S2 Z-scores** due to deep Siberian
  population structure (high-coverage early individuals vs pooled
  present-day reference); these are noted but excluded from the primary
  inference.
- **AADR 1240K pseudo-haploid data** cannot support tract-level analysis,
  haplotype dating, or local ancestry inference — all deferred to a
  whole-genome tier.
- **Single Denisovan reference** (Altai): "Papuan-like vs Han-like" is
  inferred from match-rate differences, not a true Denisovan phylogeny.
- **One statistic was excluded** (D(Altai, Denisova; X, Mbuti)) after sign-gate
  failure — it is a Neanderthal indicator, not a Denisovan one.

## Conclusion

**Outcome D: No replication.** Native Americans are adequately modeled as
carrying the same Denisovan component as their East Asian ancestors. The
hypothesis that Native Americans possess a Papuan-like Denisovan component is
**not supported** — their Denisovan-diagnostic profiles are consistently more
Han-like than Papuan-like (S3 < 0 for all 14 ancient American/Siberian
cohorts, bootstrap CIs excluding zero). Browning et al. 2018 Figure 4 was
over-read: it does not show a Papuan-like component in Americans, and this
controlled analysis finds no evidence for one.

## Future Data Needed

- Phased whole-genome sequences of unadmixed ancient Native Americans to
  resolve the East Asian Denisovan substructure (high- vs moderate-affinity
  split) at the tract level.
- Additional Denisovan reference genomes beyond Altai to enable a true
  Denisovan-source phylogeny.
- High-coverage ancient Papuan/Australasian genomes to refine the
  Papuan-associated Denisovan marker set.
- Local-ancestry-restricted modern American whole genomes to separate
  Native vs European vs African Denisovan tracts.
