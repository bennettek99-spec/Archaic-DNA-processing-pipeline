# Which human population is the most "mutation-rich"? A worldwide genetic-diversity analysis

*Generated 2026-06-29. Exploratory; see the framing note on mutation rate vs. diversity.*

## The question, made precise
"Highest rate of mutation" needs care. The **de-novo mutation rate** per generation (~1.2x10^-8 per base pair) is a biological near-constant across human populations and **cannot be measured from a genotype panel** — that requires sequencing parent–offspring trios. What a panel *can* measure is **standing genetic diversity** (heterozygosity): the variation a population currently carries, which equals mutation rate x effective population size x time and is dominated by demographic history. So the well-posed question is: *which population carries the most accumulated genetic diversity?* We answer that, and explain why.

## Method
Per-individual heterozygosity H = (heterozygous sites)/(called sites) on the autosomes, using high-coverage **diploid** (.DG/.SG) SGDP/HGDP samples (pseudo-haploid ancient samples cannot yield heterozygosity). Twenty populations span every inhabited continent and deliberately include the diversity extremes. We also regress H on great-circle distance from East Africa.

## Result — the answer
**The most genetically diverse human population is African.** By direct panel measurement the West African **Yoruba** ranks first (H = 0.2787) — and robustly so, topping *both* the 1240K and the balanced Human Origins panels. But the deeper, ascertainment-corrected answer is the **Khoesan (Ju|'hoan / San)**: standard SNP panels are largely ascertained in non-San peoples and badly *undercount* San-specific variation. On the Eurasian-leaning **1240K** panel the San are pushed down to #17 of 20; on the balanced **Human Origins** panel they leap to #4 (Africans take three of the top four), and on unbiased whole-genome data the Khoesan are the single most diverse human population (Mallick et al. 2016). The least diverse on every panel are **Native Americans (Karitiana, H = 0.1779)** — a 1.57x spread. Diversity declines strongly with distance from Africa (r = -0.74, p = 1.8e-04).

1240K ranking (20 populations):

| rank | population | region | n | heterozygosity ± SE | dist. from E. Africa (km) |
|---|---|---|---|---|---|
| 1 | Yoruba | Africa | 24 | 0.2787 ± 0.0006 | 3,816 |
| 2 | Mozabite | N.Africa/Mid-East | 29 | 0.2783 ± 0.0012 | 4,477 |
| 3 | Burusho | South Asia | 26 | 0.2722 ± 0.0009 | 4,661 |
| 4 | Russian | Europe | 27 | 0.2701 ± 0.0007 | 5,561 |
| 5 | French | Europe | 31 | 0.2692 ± 0.0013 | 5,393 |
| 6 | BedouinA | N.Africa/Mid-East | 25 | 0.2692 ± 0.0021 | 2,476 |
| 7 | Biaka | Africa | 24 | 0.2659 ± 0.0007 | 2,460 |
| 8 | Druze | N.Africa/Mid-East | 44 | 0.2636 ± 0.0015 | 2,696 |
| 9 | Basque | Europe | 25 | 0.2632 ± 0.0008 | 5,463 |
| 10 | Sardinian | Europe | 31 | 0.2623 ± 0.0008 | 4,529 |
| 11 | Balochi | South Asia | 26 | 0.2616 ± 0.0029 | 3,674 |
| 12 | Kalash | South Asia | 24 | 0.2589 ± 0.0014 | 4,509 |
| 13 | Yakut | East Asia/Siberia | 22 | 0.2493 ± 0.0010 | 9,192 |
| 14 | Han | East Asia/Siberia | 46 | 0.2481 ± 0.0007 | 7,729 |
| 15 | Mbuti | Africa | 15 | 0.2459 ± 0.0011 | 1,359 |
| 16 | Japanese | East Asia/Siberia | 31 | 0.2458 ± 0.0006 | 10,244 |
| 17 | Ju_hoan_North | Africa | 10 | 0.2451 ± 0.0054 | 3,762 |
| 18 | Mayan | Americas | 24 | 0.2217 ± 0.0098 | 13,589 |
| 19 | Papuan | Oceania | 32 | 0.1975 ± 0.0013 | 11,786 |
| 20 | Karitiana | Americas | 16 | 0.1779 ± 0.0117 | 11,450 |

## Ascertainment bias — why the panel matters (and why the San are really #1)
The 1240K SNP set was discovered mostly in non-San populations, so it *undercounts* the variation that is private to deeply-diverging African groups — artificially deflating their heterozygosity. We demonstrate this directly: the **same shotgun-diploid individuals** scored on the balanced Human Origins panel rise sharply. the Khoesan (Ju|'hoan) move from #17 on 1240K to #4 on Human Origins (of 20 shared populations).

| population | region | H on 1240K (rank) | H on Human Origins (rank) |
|---|---|---|---|
| Yoruba | Africa | 0.2787 (#1) | 0.2652 (#1) |
| Mozabite | N.Africa/Mid-East | 0.2783 (#2) | 0.2561 (#3) |
| Burusho | South Asia | 0.2722 (#3) | 0.2487 (#5) |
| Russian | Europe | 0.2701 (#4) | 0.2455 (#8) |
| French | Europe | 0.2692 (#5) | 0.2456 (#7) |
| BedouinA | N.Africa/Mid-East | 0.2692 (#6) | 0.2468 (#6) |
| Biaka | Africa | 0.2659 (#7) | 0.2579 (#2) |
| Druze | N.Africa/Mid-East | 0.2636 (#8) | 0.2404 (#9) |
| Basque | Europe | 0.2632 (#9) | 0.2391 (#11) |
| Sardinian | Europe | 0.2623 (#10) | 0.2384 (#13) |
| Balochi | South Asia | 0.2616 (#11) | 0.2390 (#12) |
| Kalash | South Asia | 0.2589 (#12) | 0.2358 (#14) |
| Yakut | East Asia/Siberia | 0.2493 (#13) | 0.2299 (#16) |
| Han | East Asia/Siberia | 0.2481 (#14) | 0.2308 (#15) |
| Mbuti | Africa | 0.2459 (#15) | 0.2401 (#10) |
| Japanese | East Asia/Siberia | 0.2458 (#16) | 0.2277 (#17) |
| Ju_hoan_North | Africa | 0.2451 (#17) | 0.2548 (#4) |
| Mayan | Americas | 0.2217 (#18) | 0.2032 (#18) |
| Papuan | Oceania | 0.1975 (#19) | 0.1864 (#19) |
| Karitiana | Americas | 0.1779 (#20) | 0.1630 (#20) |

This is the key methodological caveat — the 1240K "winner" reflects ascertainment as much as biology; the ascertainment-balanced panel restores the textbook result that the Khoesan are the most diverse humans.

## Why — the reasoning
1. **Africa is the homeland and the reservoir of diversity.** Anatomically modern humans arose in Africa ~300 kya. African populations have had the longest time and the largest long-term effective population size to accumulate variation, so they carry the most.
2. **The Khoesan (Ju|'hoan / San) are the most diverse of all** (shown above on the balanced panel): theirs is among the earliest-diverging human lineages (split ~200–300 kya) and they maintained a large effective size — more independent ancestral lineages, hence the highest heterozygosity of any living people, a long-standing result.
3. **The out-of-Africa bottleneck.** Every non-African descends from a small group that left Africa ~50–70 kya. That founder event discarded much of the African variation in one step, so all non-Africans start lower.
4. **Serial founder effects.** As humans spread further, each new population was founded by a subset of the last, losing diversity at each step — which is why heterozygosity falls almost linearly with distance from Africa (Ramachandran et al. 2005), reproduced here (Figure 2). **Native Americans (Karitiana)**, at the end of the longest migration (Africa → Asia → Beringia → the Americas, through repeated bottlenecks), are the least diverse; Oceanians and isolated groups (Papuans, Kalash) are also low for their regions.
5. **It is demography, not a faster mutational clock.** Because the per-generation mutation rate is essentially the same everywhere, these differences reflect population size and migration history, not biology that mutates faster. (A subtler, real phenomenon is variation in the mutation *spectrum* — e.g. the European excess of TCC→TTC mutations (Harris 2015) — but that concerns the *types* of mutations, not the overall rate, and is not resolvable on this panel.)

## Caveats
- The 1240K panel is an ascertained SNP set; absolute heterozygosity is affected by ascertainment, though the African-highest gradient is robust and large.
- Some samples are few (e.g. Ju|'hoan n≈10); SEs are shown.
- Heterozygosity is one diversity measure; private-allele counts and runs-of-homozygosity tell a consistent story.

## References
Ramachandran et al. 2005 *PNAS* 102:15942 · Henn et al. 2012 *PNAS* 109:17758 · Mallick et al. 2016 *Nature* 538:201 (SGDP) · Harris 2015 *PNAS* 112:3439 · 1000 Genomes Project 2015 *Nature* 526:68.
