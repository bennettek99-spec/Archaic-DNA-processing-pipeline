# FADS1-2: a candidate Neanderthal-introgression selection signal in the Italian transect

*Deep-dive into the one locus flagged by the Etruscan scan. Generated 2026-06-29. Exploratory — a hypothesis, not a validated selection scan.*

## Background
The FADS1/FADS2/FADS3 cluster on chromosome 11q12-13 encodes fatty-acid desaturases that synthesise long-chain polyunsaturated fatty acids (omega-3/6). It is one of the strongest signals of recent positive selection in Europeans: a regulatory haplotype that boosts endogenous PUFA synthesis rose sharply in frequency after the Neolithic shift to plant-rich diets (Mathieson et al. 2015; Ye et al. 2017), and the locus also carries archaic (Neanderthal/Denisovan) haplotypes that were adaptively introgressed in several populations (Buckley et al. 2017). FADS variation affects lipid metabolism and is associated with cardiometabolic traits today.

## What we measured
We identified **5 archaic-informative SNPs** across the FADS cluster (alleles carried by the high-coverage Altai + Vindija Neanderthals and near-absent in Africans), then traced the frequency of the Neanderthal allele across the Italian and broader European time transects, controlling for genome-wide archaic ancestry and ancestry (PC1) so the temporal term isolates locus-specific change.

### Archaic-informative SNPs at FADS
| rsID | position (hg19) | gene | archaic allele | African freq |
|---|---|---|---|---|
| rs174535 | 61,551,356 | TMEM258/FEN1 | C | 0.085 |
| rs174536 | 61,551,927 | TMEM258/FEN1 | C | 0.007 |
| rs174537 | 61,552,680 | TMEM258/FEN1 | T | 0.007 |
| rs102274 | 61,557,826 | TMEM258/FEN1 | C | 0.007 |
| rs174550 | 61,571,478 | FADS1 | C | 0.007 |

### Neanderthal-allele frequency over time
| era | Italy % | n | Europe % | n |
|---|---|---|---|---|
| Mesolithic | 97.6 | 13 | 88.3 | 212 |
| Neolithic | 66.4 | 48 | 68.4 | 1656 |
| Copper/EBA | 49.3 | 62 | 62.4 | 1588 |
| Bronze Age | 45.9 | 68 | 54.0 | 985 |
| Iron Age | 46.2 | 272 | 46.5 | 1330 |
| Roman/Medieval | 31.4 | 202 | 41.9 | 4151 |

## Result
The FADS Neanderthal allele **declines toward the present**:
- Italy: **+6.679 pp/kyr**, p = 0.0000 (n=406).
- Europe (broad): +6.343 pp/kyr, p = 0.0000 (n=5936).

This indicates a pan-European pattern.

## Interpretation
The five informative SNPs lie in the TMEM258–FADS1 regulatory block and include **rs174537** and **rs174550** — canonical markers of the European FADS selective sweep and FADS1-expression eQTLs (Ameur et al. 2012; Mathieson et al. 2015). A declining archaic allele here is consistent with the well-documented post-Neolithic FADS sweep: the haplotype favoured by agricultural diets is the *non-archaic*, derived modern-human regulatory haplotype, whose rise displaced the alternative (here, Neanderthal-matching) alleles. In other words, selection at FADS acted **against** the archaic-allele background as the adaptive modern haplotype increased — the mirror image of the classic "adaptive introgression" cases (BNC2, OAS, TLR) where an archaic allele rose. The signal therefore reflects real, strong selection at a functionally important locus, but it is **not** evidence that the Neanderthal variant itself was favoured.

Importantly, the archaic allele sits near **fixation (~90%) in Mesolithic hunter-gatherers** — far higher than a typical rare introgressed allele (a few percent). That tells us "archaic-informative" at FADS is largely capturing **common ancestral-haplotype variation** (alleles Neanderthals share by virtue of their deep divergence, near-absent in Africans, and carried at high frequency by pre-agricultural Europeans) rather than a recent Neanderthal-introgression event. The locus is therefore best read as the canonical FADS dietary sweep — which our archaic-allele scan recovered independently and pan-continentally — not as Neanderthal adaptive introgression.

## Caveats
- Rests on 5 archaic-informative SNPs on a capture array; per-SNP power is low and FADS has a complex, recombining haplotype structure, so "the archaic allele" is a coarse summary.
- "Archaic-informative" = Altai/Vindija-matching and African-absent; incomplete lineage sorting can mimic introgression, and at FADS some such alleles may tag the modern selected haplotype by linkage rather than true Neanderthal ancestry.
- Frequency change over time reflects both selection and ancestry turnover; we control for genome-wide archaic ancestry and PC1 but not for every fine-scale ancestry axis.
- A definitive test requires haplotype-resolved data and an explicit selection model (e.g. time-series allele-frequency likelihood), not array genotypes.

## References
Mathieson et al. 2015 *Nature* 528:499 · Ye et al. 2017 *Mol. Biol. Evol.* 34:509 · Buckley et al. 2017 *Mol. Biol. Evol.* 34:1307 · Ameur et al. 2012 *Am. J. Hum. Genet.* 90:809 · Patterson et al. 2012 *Genetics* 192:1065.
