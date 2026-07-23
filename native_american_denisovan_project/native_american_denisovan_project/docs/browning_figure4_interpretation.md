# Critical interpretation of Browning et al. 2018, Figure 4

*Source verified this session: PMC5866234 full text (PMID 29551270, DOI 10.1016/j.cell.2018.02.031), Cell 173(1):53-61.e9. Quotation is kept to brief necessary phrases; methods and findings are otherwise paraphrased.*

## What Figure 4 actually is

Figure 4 is a grid of **two-way contour density plots**. For each population panel, every axis is the **match proportion** of an introgressed segment to one archaic reference genome:

- horizontal axis: proportion of the segment's putative archaic-specific alleles that match the **Altai Neanderthal**;
- vertical axis: proportion that match the **Altai Denisovan**.

"Match proportion" = (alleles in the segment matching the archaic genome) / (alleles comparable, i.e. not masked in that archaic genome). Only segments with **at least 10 variants not masked in the Neanderthal genome AND at least 10 not masked in the Denisovan genome** are plotted. The segments themselves were called by **Sprime**, a **reference-free** S*-like method that infers archaic-specific alleles from high LD between alleles that are rare/absent in the Yoruba outgroup — it does **not** use the Altai genome to find them; the reference is only used afterwards to label the segments.

The grid rows are: (1) European, (2) East Asian, (3) South Asian, (4) **"American and SGDP Papuan populations."** The American populations are the four admixed 1000 Genomes American groups — PUR (Puerto Rican), CLM (Colombian), MXL (Mexican-American), PEL (Peruvian) — **not** unadmixed or ancient Native Americans. There are **no ancient genomes** in Figure 4.

## What the paper says Figure 4 shows

In every population there is a large cluster at high Neanderthal-match / low Denisovan-match (Neanderthal-introgressed segments), a small near-zero cluster (false positives), and in Asian/Papuan panels a **third cluster** with high Denisovan-match / low Neanderthal-match (Denisovan-introgressed segments). In Japanese and the three Chinese populations the Denisovan cluster is **wide and bimodal**, and a formal two-component Gaussian-mixture test is significant (Bonferroni p<0.0026) **only** in CHS, CHB, CDX, JPT (plus FIN and Punjabi at uncorrected p<0.05) — **Table 2 lists no American population**. Roughly one-third of East-Asian Denisovan segments fall in the **high-altair-affinity** component (match ~0.80 to Altai Denisovan); the rest are **moderate-affinity** (match ~0.46-0.52).

The authors' interpretation: the high-affinity component is "primarily present in East Asians"; the moderate-affinity component is the major part of Denisovan ancestry in **Papuans and South Asians**. A consistent scenario would have the high-affinity component introgressing into East Asia after the East/South Asia split, and the moderate-affinity component entering Asia via migration from Papuan-related ancestors.

## What Figure 4 says about Americans specifically

Two sentences in the paper directly bear on the American panels, and both **deflate** the "Native Americans carry a Papuan-like Denisovan component" reading:

1. "The match rate to the Altai Neanderthal and Altai Denisovan genomes is **lower in the American populations** than in the other 1000 Genomes populations (Figure 3). This is likely due to the fact that the American populations are **admixed and thus have higher background levels of LD that could cause false positive results**." The authors themselves attribute the American signal to an **admixture/LD artifact**, not to a real Denisovan source.

2. "Figures 4 and 5 also indicate that several other populations may carry a small proportion of segments introgressed from Denisovans. These include the Finns… and **admixed American populations whose Native American ancestors are related to East Asians**." Any American Denisovan segments are explicitly attributed to **Native American ancestors that are East-Asian-related** — i.e. the *same* source as the East Asian Denisovan component, **not** a distinct Papuan-like source.

The two-component test was **not** significant in any American population, and no American population appears in Table 2. Figure 4 therefore **does not decompose** American Denisovan segments into Papuan-like vs Han-like components at all.

---

## Strongest defensible interpretation of Browning et al. Figure 4

- Figure 4 establishes, with a reference-free method, that **Denisovan introgression is not monolithic**: East Asians carry a detectable Denisovan component that is itself bimodal in similarity to the Altai Denisovan (a high-altair-affinity ~1/3 and a moderate-affinity ~2/3), while Papuans/South Asians are dominated by the moderate-affinity component.
- It shows that **admixed 1000 Genomes American populations contain a small number of Denisovan-like segments** (visible in the row-4 panels), which is consistent with their Native American fraction being East-Asian-related and therefore carrying some East-Asian-type Denisovan ancestry.
- It motivates — but does **not** resolve — the question of whether the Denisovan ancestry that reached the Americas via East-Asian-related ancestors is the high-altair-affinity component, the moderate-affinity (Papuan-like) component, or both.
- The only defensible *population-history* statement Figure 4 supports about the Americas is: **admixed Americans show some Denisovan-like segments, plausibly inherited from East-Asian-related Native American ancestors; their match rates are depressed, most likely by admixture-driven LD, and no two-component structure was demonstrated in Americans.**

## Claims that Figure 4 does not justify

- It does **not** show that **unadmixed or ancient Native Americans** carry Denisovan ancestry — the plotted American samples are admixed 1000 Genomes populations (PUR/CLM/MXL/PEL), and no ancient genome is shown.
- It does **not** show that Native American Denisovan ancestry is **Papuan-like** (moderate-affinity). Placing Americans in the same Figure-4 row as Papuans is a layout choice, not evidence of shared component; the authors do not claim Americans carry the moderate-affinity component, and the two-component test was negative in Americans.
- It does **not** show that Americans carry the **high-altair-affinity** component either — no American two-component decomposition is reported.
- It does **not** establish that a **distinct Denisovan pulse** reached the Americas. The authors' own wording ties American segments to East-Asian-related ancestors and to admixture artifacts.
- It does **not** provide a **genome-wide ancestry proportion** for Americans — it shows segment-level match rates, not a calibrated percentage.
- It does **not** give **confidence intervals / bootstrap distributions** for the American panels' cluster positions; the contours are descriptive density plots, and the only formal significance test (Table 2) excludes Americans.
- It does **not** rule out that the few American Denisovan-like segments are **shared ancestral polymorphism, Neanderthal–Denisovan shared derived alleles, or LD false positives** (the authors flag LD false positives themselves).
- It does **not** speak to the **time** the signal entered the Americas (no ancient individuals, no tract-length dating in Americans, and the paper states haplotype lengths did not significantly differ between components).

## Implication for this project

The project's central hypothesis — *"Do Native Americans carry a Denisovan ancestry component more closely related to the Papuan-associated component than to the principal East Asian component?"* — is **not** a claim Browning et al. made or demonstrated. It is a hypothesis **generated by an over-reading of Figure 4's layout** (Americans plotted beside Papuans) and by the broader two-component Denisovan model. The project must therefore treat the hypothesis as **unreplicated** and assess whether the data can test it at all, with the prior expectation that Browning's own evidence points toward an East-Asian-related (or artifactual) source for any American Denisovan signal, not a Papuan-like one.
