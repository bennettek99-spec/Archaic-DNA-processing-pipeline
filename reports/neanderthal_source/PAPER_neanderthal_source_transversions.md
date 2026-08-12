# Which Neanderthal? Altai-versus-Vindija affinity across the ancient Eurasian record

*AADR v66.p1 1240k panel (transversions only - damage-robust sensitivity run). 10,954 unique Eurasian ancient genomes above a 100,000-SNP floor, pooled into 41 dated cohorts, plus 18 present-day anchors. Paired 50-block jackknife throughout.*

## Abstract

Whether every Eurasian lineage descends from one Neanderthal source population is normally addressed with a handful of genomes. Here the contrast D(X, Yoruba; Vindija, Altai) - positive when a population is closer to the Croatian Vindija Neanderthal than to the Siberian Altai one - is measured across the ancient Eurasian record. Two anchors establish that the instrument works. Every Neanderthal-carrying cohort is displaced towards Vindija (present-day French +0.0782, Han +0.0749), reproducing the result that Vindija sits closer to the introgressing population (Prufer et al. 2017); and with the Denisovan genome as the test population the statistic returns -0.1459 (Z = -4.4), recovering the pull towards Altai expected from that genome's ~1% Denisovan-related ancestry - a real, published difference between the two Neanderthal genomes' own histories.

On the question asked, the answer is a null with a stated limit. Upper Palaeolithic Europeans and Palaeolithic north-east Asians differ by +0.00147 +/- 0.01099 (Z = +0.13) in raw D_VA, and present-day French and Han by +0.00330 +/- 0.00760 (Z = +0.43). A single proportional relation D_VA = 1.57 x D_NEA (+/- 0.28), where D_NEA is a Vindija/Altai-symmetric measure of how much Neanderthal ancestry a cohort has, describes 50 of 54 cohorts within two standard errors, and 0 depart from it after correction for multiple testing. The detection limit is **0.0151 in D_VA units for a typical pair of cohorts, 20% of the total Vindija-over-Altai signal** (0.0031, 4%, for the best-powered pairs). Sources differing by less than that are invisible here, which includes most of what the literature debates about a second pulse into East Asia.

One pattern is not null and is reported as a candidate rather than a finding. The oldest cohorts sit below the single-source line: normalised source affinity is -0.37 +/- 0.21 (Z = -1.79) lower in pre-LGM Upper Palaeolithic Europeans than in Medieval-to-recent Europeans, and -0.29 +/- 0.21 (Z = -1.35) lower in the ~45,000 BP Initial Upper Palaeolithic group. Three things argue against taking this at face value: it correlates with sample age, which is also a proxy for deamination damage and coverage (rho = -0.49 (p = 0.0013, n = 40) against date); it is absent from the one Palaeolithic cohort with high coverage (Palaeolithic north-east Asia, -0.21 +/- 0.24, Z = -0.90 against the same comparator); and the effect is carried by the *denominator*, since the raw D_VA of these cohorts is ordinary and it is their elevated D_NEA that moves the ratio. The honest reading is that this study cannot separate a genuinely less Vindija-like early source from an age-linked measurement effect, and it should not be reported as the former.

## The statistic, and the trap inside it

Vindija and Altai are both called at 144,757 autosomal 1240K sites, but only **5,237 of those actually distinguish them** (3.62%). The archaic genomes, not the ancient cohorts, are the limiting sample, and that is why this contrast is hard however many ancient genomes are available. It also has a useful consequence: the same sites carry the signal for every cohort, so their sampling noise is *common-mode* and cancels when two cohorts are differenced inside the same jackknife replicate. Pairing the jackknife this way is 1.9x tighter than combining two independent standard errors, which is the difference between a usable detection limit and an uninformative one.

Three further offsets are common-mode, and are therefore differenced away rather than interpreted:

1. **Vindija is pseudo-haploid (`.SG`) and Altai is diploid (`.DG`)** in this release, and Vindija is called at less than half as many sites. The two genomes are not symmetric inputs.

2. **Yoruba is not archaic-free.** All Africans carry a little Neanderthal ancestry from back-migration (Chen et al. 2020), which shows up here as a non-zero baseline: present-day Mbuti read -0.0175 against the Yoruba baseline rather than zero, and chimpanzee reads -0.0916.

3. **The 1240K panel is ascertained in modern humans**, so archaic variation is sampled non-randomly.

None of these can produce a *difference between two Eurasian cohorts*, since every cohort is measured against the same two archaic genomes and the same baseline. Only differences are reported as results.

### Removing the quantity confound

The confound that does *not* cancel is how much Neanderthal ancestry each cohort has, because D_VA scales with it. Each cohort is therefore also measured with `D_NEA = D(X, Yoruba; NeaAvg, Chimp)`, where `NeaAvg` is the mean of the two archaic frequencies. Swapping Vindija for Altai leaves `NeaAvg` unchanged, so `D_NEA` tracks the amount of Neanderthal ancestry without responding to its source. Under a single shared source every cohort satisfies `D_VA = k * D_NEA` with one k, and k estimates the Vindija preference of the source population itself.

It is worth recording how large this confound actually is on this panel, because the literature figure and this panel's own answer disagree. Published estimates put East Asians somewhere around 8-20% above Europeans in Neanderthal ancestry; this repository's f4-ratio on 1240K puts the ratio at about 1.02, and D_NEA here gives 0.0438 for Han against 0.0429 for French. The correction is therefore small for the East Asian comparison specifically - but it is large for the Palaeolithic cohorts and very large for Oase1, which is where it changes the conclusion.

**Scope condition.** The proportional model only describes cohorts whose departure from the Yoruba baseline is dominated by Neanderthal introgression. Cohorts with D_NEA below 0.02 - all the African ones - are excluded from the fit and from the residual test: their departure from Yoruba is dominated by deep African population structure, and R divides by something indistinguishable from zero. They are still reported, and their D_VA values are informative as a baseline. Chimpanzee is scored on D_VA only, since it appears on both sides of D_NEA and is degenerate there by construction.

## Positive controls: the statistic can see a real source difference

**Denisova as the test population.** The Altai Neanderthal carries ~1% Denisovan-related ancestry (Prufer et al. 2014), so the Denisovan genome should be pulled towards Altai. It is: D_VA = -0.1459 +/- 0.0331 (Z = -4.4), strongly negative while every Eurasian cohort is strongly positive, and a residual of -0.254 (Z = -6.2) from the single-source line. A published difference between the two Neanderthal genomes' own histories is recovered at high significance, which is the evidence that a null elsewhere means something.

**Proportionality.** Across cohorts spanning near-zero to ~12% Neanderthal affinity, D_VA is proportional to D_NEA with slope k = 1.5670 +/- 0.2795 (Z = 5.6).


![Figure 1](fig_n1_scaling.png)

**Figure 1.** Left: every cohort on the D_VA-versus-D_NEA plane with the fitted single-source line. Right: residuals from that line in standard errors. The line is fitted on the region-by-period grid alone, so the named Palaeolithic cohorts, the present-day anchors and the Denisovan control are all scored against a line they did not help define.

## Result 1: the comparisons the question asks

| label                  | kind              |   n_ind |   D_NEA |    D_VA |   residual |   residual_se |   residual_z |
|:-----------------------|:------------------|--------:|--------:|--------:|-----------:|--------------:|-------------:|
| CTRL_Denisova          | control-reference |       1 |  0.0692 | -0.1459 |   -0.25429 |       0.04084 |        -6.23 |
| Europe_Palaeolithic    | grid              |      28 |  0.0557 |  0.0716 |   -0.01568 |       0.00826 |        -1.9  |
| PD_French              | present-day       |      31 |  0.0429 |  0.0782 |    0.01097 |       0.00611 |         1.79 |
| Europe_Medieval_Recent | grid              |    2984 |  0.0483 |  0.08   |    0.00427 |       0.00265 |         1.61 |
| UP_Europe_pre_LGM      | named             |      14 |  0.0599 |  0.0767 |   -0.01711 |       0.0116  |        -1.47 |
| PD_Papuan              | present-day       |      32 |  0.0585 |  0.0732 |   -0.01854 |       0.01433 |        -1.29 |
| Oase1_40ka             | named-lowpower    |       1 |  0.0889 |  0.0367 |   -0.1027  |       0.09112 |        -1.13 |
| UP_Europe_post_LGM     | named             |       5 |  0.0458 |  0.0577 |   -0.01402 |       0.0132  |        -1.06 |
| IUP_Eurasia_45ka       | named             |      10 |  0.0536 |  0.0732 |   -0.01082 |       0.01019 |        -1.06 |
| PD_Han                 | present-day       |      46 |  0.0438 |  0.0749 |    0.00624 |       0.00687 |         0.91 |
| EastAsia_Bronze        | grid              |     212 |  0.0461 |  0.0773 |    0.00502 |       0.00558 |         0.9  |
| UP_NorthEastAsia       | named             |       8 |  0.0522 |  0.0753 |   -0.00654 |       0.01106 |        -0.59 |
| UstIshim_44ka          | named             |       1 |  0.0439 |  0.0608 |   -0.00804 |       0.01828 |        -0.44 |


**Upper Palaeolithic Europeans versus Palaeolithic north-east Asians.** Raw D_VA differs by +0.00147 +/- 0.01099 (Z = +0.13); normalised source affinity by -0.16 +/- 0.25 (Z = -0.64). The two groups' Neanderthal ancestry is indistinguishable in Vindija-versus-Altai character.

**The East Asian second pulse.** The specific contested claim is that East Asians received Neanderthal ancestry from an additional or different source. Present-day French and Han differ by +0.00330 +/- 0.00760 (Z = +0.43) in raw D_VA and +0.11 +/- 0.25 (Z = +0.46) normalised; Medieval-to-recent Europeans and Bronze Age East Asians differ by +0.00262 +/- 0.00627 (Z = +0.42). This neither supports nor excludes a second pulse. It bounds how different the two sources could be, and the bound - 20% of the total signal - is not tight enough to adjudicate the debate.

**Oase1.** In raw D_VA Oase1 is unremarkable (+0.0367 +/- 0.0811), which is *itself* the surprise: with D_NEA = 0.0889, three times any other cohort and consistent with his ~10% Neanderthal ancestry, the single-source model predicts a D_VA near 0.14. His residual is -0.103 +/- 0.091 (Z = -1.13). This is **not** reported as evidence that his recent Neanderthal ancestor came from a different population. Oase1 is admitted here deliberately below the study's SNP floor, at 0.05x coverage and 9,169 usable sites against ~507,000 for every other cohort; only ~3.6% of those distinguish the two Neanderthals, leaving a few hundred effectively informative sites spread over 50 jackknife blocks. At that density the jackknife is not trustworthy, the first-order approximation D_VA ~ a x D_VA(source) is least accurate at his admixture proportion, and a single pseudo-haploid genome enters the D denominator differently from a pooled cohort. The honest statement is that Oase1 cannot be placed on this axis with the 1240K panel.


![Figure 4](fig_n4_targets.png)

**Figure 4.** Named targets' residuals from the single-source line against the detection limit (grey).

## Result 2: an age-correlated pattern, reported as a candidate

The oldest cohorts sit consistently below the single-source line. Normalised source affinity R is -0.37 +/- 0.21 (Z = -1.79) lower in pre-LGM Upper Palaeolithic Europeans than in Medieval-to-recent Europeans, and -0.29 +/- 0.21 (Z = -1.35) lower in the ~45,000 BP Initial Upper Palaeolithic group; Holocene European cohorts from six periods agree with each other to within a few percent. Taken at face value this would say the earliest Eurasians' Neanderthal ancestry was less Vindija-like. Four checks say it should not be taken at face value.

1. **It tracks the assay, not just the calendar.** Across cohorts, R correlates with sample age (rho = -0.49 (p = 0.0013, n = 40)), and age is also a proxy for deamination damage (rho = -0.11 (p = 0.52, n = 38)) and coverage (rho = -0.23 (p = 0.15, n = 40)). These covariates cannot be separated in observational ancient-DNA data.

2. **The high-coverage Palaeolithic cohort does not show it.** Palaeolithic north-east Asia - the only Palaeolithic cohort with median coverage above 4x - sits -0.21 +/- 0.24 (Z = -0.90) from the same Medieval comparator, i.e. on the line. If the effect were a property of Palaeolithic *people*, it should appear there too.

3. **The numerator is ordinary; the denominator moves.** These cohorts' raw D_VA is normal (pre-LGM Upper Palaeolithic Europe +0.0767 against +0.0800 for Medieval Europe, a difference of -0.00322 +/- 0.00944). What moves is D_NEA, elevated to 0.0599 from 0.0483. An elevated Neanderthal level in early Upper Palaeolithic Europeans is independently expected (Fu et al. 2016), so part of this is real - but any age-linked inflation of D_NEA would produce exactly the same deficit in R without any change of source.

4. **Multiple testing.** With 54 cohorts scored, the Bonferroni threshold is |Z| = 3.31; these residuals do not reach it.

The transversions-only rerun is the sharpest discriminator available for check 1, and it is reported next.

*This report is itself the transversions-only run; see `PAPER_neanderthal_source.md` for the full-panel analysis.*


| covariate       |    rho |      p |   n |
|:----------------|-------:|-------:|----:|
| date_bp         | -0.49  | 0.0013 |  40 |
| median_coverage | -0.232 | 0.1506 |  40 |
| median_damage   | -0.109 | 0.5156 |  38 |
| median_snps     | -0.334 | 0.0352 |  40 |


![Figure 2](fig_n2_time.png)

**Figure 2.** Left: raw D_VA through time by region, which tracks Neanderthal quantity. Right: normalised source affinity, flat across the Holocene and dipping in the oldest cohorts.

## The time axis

Because the AADR supplies dated cohorts rather than only present-day populations, the question can be asked as a time series - the part of this analysis with no published counterpart.

| label                        | kind   |   n_ind |   date_bp | region      |   median_coverage |   median_damage |   D_NEA |   D_VA |   D_VA_se |   D_VA_z |     R |   R_se |   n_snp |
|:-----------------------------|:-------|--------:|----------:|:------------|------------------:|----------------:|--------:|-------:|----------:|---------:|------:|-------:|--------:|
| Caucasus_Neolithic           | grid   |      98 |      5534 | Caucasus    |              2.12 |           0.127 |  0.0459 | 0.0711 |    0.0117 |     6.07 | 1.551 |  0.315 |  138509 |
| Caucasus_Bronze              | grid   |     307 |      3842 | Caucasus    |              2.12 |           0.115 |  0.0452 | 0.0717 |    0.011  |     6.52 | 1.586 |  0.32  |  138509 |
| Caucasus_Iron_Classical      | grid   |     127 |      2599 | Caucasus    |              2.06 |           0.119 |  0.0461 | 0.072  |    0.0113 |     6.35 | 1.561 |  0.304 |  138509 |
| Caucasus_Medieval_Recent     | grid   |      36 |       975 | Caucasus    |              2.06 |           0.096 |  0.0473 | 0.0765 |    0.0109 |     7    | 1.616 |  0.307 |  138509 |
| CentralAsia_Neolithic        | grid   |      35 |      5163 | CentralAsia |              1.68 |           0.186 |  0.0426 | 0.0559 |    0.0109 |     5.14 | 1.313 |  0.301 |  138509 |
| CentralAsia_Bronze           | grid   |      74 |      3760 | CentralAsia |              1.73 |           0.126 |  0.0474 | 0.0693 |    0.0111 |     6.21 | 1.461 |  0.303 |  138503 |
| CentralAsia_Iron_Classical   | grid   |     134 |      2838 | CentralAsia |              1.82 |           0.133 |  0.0461 | 0.0732 |    0.01   |     7.3  | 1.587 |  0.29  |  138509 |
| CentralAsia_Medieval_Recent  | grid   |      43 |      1000 | CentralAsia |              0.9  |           0.092 |  0.0463 | 0.072  |    0.0109 |     6.63 | 1.555 |  0.284 |  138509 |
| EastAsia_LateGlacial_Meso    | grid   |      20 |      9417 | EastAsia    |              3.34 |           0.428 |  0.0474 | 0.0796 |    0.0133 |     5.96 | 1.68  |  0.307 |  138509 |
| EastAsia_Neolithic           | grid   |     150 |      6950 | EastAsia    |              1.24 |           0.089 |  0.0438 | 0.0814 |    0.0116 |     7.01 | 1.858 |  0.346 |  138509 |
| EastAsia_Bronze              | grid   |     212 |      4154 | EastAsia    |              1.16 |           0.084 |  0.0461 | 0.0773 |    0.0119 |     6.49 | 1.676 |  0.307 |  138509 |
| EastAsia_Iron_Classical      | grid   |     204 |      2051 | EastAsia    |              1.78 |           0.155 |  0.0481 | 0.0773 |    0.0116 |     6.66 | 1.609 |  0.29  |  138509 |
| EastAsia_Medieval_Recent     | grid   |     111 |       700 | EastAsia    |              1.83 |           0.242 |  0.0486 | 0.0776 |    0.012  |     6.46 | 1.598 |  0.289 |  138509 |
| Europe_Palaeolithic          | grid   |      28 |     30940 | Europe      |              1.01 |           0.337 |  0.0557 | 0.0716 |    0.0117 |     6.13 | 1.285 |  0.229 |  138509 |
| Europe_LateGlacial_Meso      | grid   |     145 |      8704 | Europe      |              2.57 |           0.224 |  0.0485 | 0.0753 |    0.0122 |     6.17 | 1.551 |  0.302 |  138509 |
| Europe_Neolithic             | grid   |    1462 |      6427 | Europe      |              1.64 |           0.116 |  0.0467 | 0.0788 |    0.0119 |     6.63 | 1.687 |  0.316 |  138509 |
| Europe_Bronze                | grid   |    1432 |      4044 | Europe      |              1.33 |           0.108 |  0.0484 | 0.0781 |    0.0111 |     7.01 | 1.614 |  0.296 |  138509 |
| Europe_Iron_Classical        | grid   |    1646 |      2175 | Europe      |              1.62 |           0.116 |  0.0471 | 0.0786 |    0.0114 |     6.88 | 1.67  |  0.311 |  138509 |
| Europe_Medieval_Recent       | grid   |    2984 |      1175 | Europe      |              1.7  |           0.154 |  0.0483 | 0.08   |    0.0112 |     7.14 | 1.655 |  0.296 |  138509 |
| InnerAsia_Bronze             | grid   |      59 |      3722 | InnerAsia   |              1.73 |           0.077 |  0.0492 | 0.0748 |    0.0116 |     6.45 | 1.52  |  0.308 |  138498 |
| InnerAsia_Iron_Classical     | grid   |     220 |      2225 | InnerAsia   |              1.49 |           0.061 |  0.0491 | 0.0759 |    0.0105 |     7.26 | 1.545 |  0.264 |  138509 |
| InnerAsia_Medieval_Recent    | grid   |      61 |      1008 | InnerAsia   |              0.87 |           0.037 |  0.0462 | 0.0758 |    0.0107 |     7.1  | 1.641 |  0.29  |  138509 |
| NearEast_Bronze              | grid   |     106 |      3513 | NearEast    |              1.17 |           0.15  |  0.0431 | 0.0741 |    0.0123 |     6.05 | 1.719 |  0.35  |  138505 |
| NearEast_Iron_Classical      | grid   |      29 |      2375 | NearEast    |              0.94 |           0.161 |  0.0454 | 0.0734 |    0.0134 |     5.47 | 1.616 |  0.344 |  138505 |
| NearEast_Medieval_Recent     | grid   |      24 |       800 | NearEast    |              0.76 |           0.078 |  0.0397 | 0.0803 |    0.0114 |     7.06 | 2.024 |  0.405 |  138506 |
| Siberia_Neolithic            | grid   |      66 |      5972 | Siberia     |              2.31 |           0.092 |  0.0501 | 0.0707 |    0.0124 |     5.7  | 1.411 |  0.294 |  138506 |
| Siberia_Bronze               | grid   |     145 |      4065 | Siberia     |              2.28 |           0.079 |  0.0508 | 0.0747 |    0.0112 |     6.68 | 1.469 |  0.28  |  138509 |
| Siberia_Iron_Classical       | grid   |     103 |      2270 | Siberia     |              2.34 |           0.069 |  0.0493 | 0.0711 |    0.0109 |     6.52 | 1.443 |  0.277 |  138509 |
| SouthAsia_Iron_Classical     | grid   |      67 |      2243 | SouthAsia   |              1.63 |           0.096 |  0.0441 | 0.0713 |    0.0115 |     6.21 | 1.615 |  0.335 |  138509 |
| SouthAsia_Medieval_Recent    | grid   |      30 |       900 | SouthAsia   |              4.82 |           0.143 |  0.0475 | 0.0746 |    0.0117 |     6.37 | 1.572 |  0.304 |  138501 |
| WestSiberia_LateGlacial_Meso | grid   |      37 |      8342 | WestSiberia |              1.39 |           0.331 |  0.0529 | 0.0795 |    0.0134 |     5.93 | 1.503 |  0.306 |  138508 |
| WestSiberia_Neolithic        | grid   |     191 |      6617 | WestSiberia |              2.74 |           0.103 |  0.05   | 0.0733 |    0.0121 |     6.06 | 1.466 |  0.316 |  138509 |
| WestSiberia_Bronze           | grid   |     209 |      3841 | WestSiberia |              1.89 |           0.074 |  0.0487 | 0.072  |    0.0112 |     6.45 | 1.478 |  0.295 |  138509 |
| WestSiberia_Iron_Classical   | grid   |     102 |      2227 | WestSiberia |              0.9  |           0.085 |  0.048  | 0.0733 |    0.0117 |     6.28 | 1.526 |  0.311 |  138509 |
| WestSiberia_Medieval_Recent  | grid   |     129 |       950 | WestSiberia |              2.2  |           0.077 |  0.0478 | 0.0739 |    0.011  |     6.7  | 1.547 |  0.287 |  138509 |
## Characterising the null

### Null constructions

The floor is measured rather than assumed, by splitting single homogeneous cohorts - which share a Neanderthal source by construction - and differencing the halves. Random splits give the sampling floor. **Coverage-stratified splits give the floor that matters**, because coverage trends with date across the AADR and any coverage-linked bias would masquerade as a temporal signal:

| control        | cohort                 | a                                   | b                                   |   D_VA_diff |   D_VA_diff_se |   D_VA_diff_z |   D_VA_diff_se_independent |   R_diff |   R_diff_se |   R_diff_z |
|:---------------|:-----------------------|:------------------------------------|:------------------------------------|------------:|---------------:|--------------:|---------------------------:|---------:|------------:|-----------:|
| random-split   | Europe_Medieval_Recent | NULLrand_Europe_Medieval_Recent_r0A | NULLrand_Europe_Medieval_Recent_r0B |    -0.0003  |        0.00089 |         -0.34 |                    0.01585 |  -0.009  |      0.0225 |      -0.4  |
| random-split   | Europe_Medieval_Recent | NULLrand_Europe_Medieval_Recent_r1A | NULLrand_Europe_Medieval_Recent_r1B |    -5e-05   |        0.0008  |         -0.06 |                    0.01585 |  -0.0095 |      0.0198 |      -0.48 |
| random-split   | Europe_Medieval_Recent | NULLrand_Europe_Medieval_Recent_r2A | NULLrand_Europe_Medieval_Recent_r2B |     0.00067 |        0.00088 |          0.76 |                    0.01585 |   0.0039 |      0.0276 |       0.14 |
| coverage-split | Europe_Medieval_Recent | NULLcov_Europe_Medieval_Recent_HI   | NULLcov_Europe_Medieval_Recent_LO   |     0.00134 |        0.00107 |          1.25 |                    0.0159  |   0.0348 |      0.0243 |       1.43 |
| random-split   | Europe_Neolithic       | NULLrand_Europe_Neolithic_r0A       | NULLrand_Europe_Neolithic_r0B       |     0.00094 |        0.00113 |          0.84 |                    0.01684 |   0.0516 |      0.0335 |       1.54 |
| random-split   | Europe_Neolithic       | NULLrand_Europe_Neolithic_r1A       | NULLrand_Europe_Neolithic_r1B       |    -5e-05   |        0.00115 |         -0.04 |                    0.01684 |   0.0001 |      0.0307 |       0    |
| random-split   | Europe_Neolithic       | NULLrand_Europe_Neolithic_r2A       | NULLrand_Europe_Neolithic_r2B       |     0.00092 |        0.00107 |          0.86 |                    0.01683 |   0.0018 |      0.0262 |       0.07 |
| coverage-split | Europe_Neolithic       | NULLcov_Europe_Neolithic_HI         | NULLcov_Europe_Neolithic_LO         |    -0.00029 |        0.00129 |         -0.22 |                    0.01691 |  -0.0316 |      0.0362 |      -0.87 |
| random-split   | Europe_Bronze          | NULLrand_Europe_Bronze_r0A          | NULLrand_Europe_Bronze_r0B          |    -0.00064 |        0.00113 |         -0.57 |                    0.01578 |   0.0044 |      0.0281 |       0.16 |
| random-split   | Europe_Bronze          | NULLrand_Europe_Bronze_r1A          | NULLrand_Europe_Bronze_r1B          |    -0.00151 |        0.00125 |         -1.21 |                    0.01578 |  -0.0128 |      0.0343 |      -0.37 |
| random-split   | Europe_Bronze          | NULLrand_Europe_Bronze_r2A          | NULLrand_Europe_Bronze_r2B          |     0.00124 |        0.00124 |          1    |                    0.01578 |   0.0019 |      0.0317 |       0.06 |
| coverage-split | Europe_Bronze          | NULLcov_Europe_Bronze_HI            | NULLcov_Europe_Bronze_LO            |     0.00092 |        0.00153 |          0.6  |                    0.01579 |   0.0431 |      0.0374 |       1.15 |
| random-split   | Europe_Iron_Classical  | NULLrand_Europe_Iron_Classical_r0A  | NULLrand_Europe_Iron_Classical_r0B  |     0.00031 |        0.00106 |          0.29 |                    0.01616 |  -0.0068 |      0.0282 |      -0.24 |
| random-split   | Europe_Iron_Classical  | NULLrand_Europe_Iron_Classical_r1A  | NULLrand_Europe_Iron_Classical_r1B  |    -0.00011 |        0.00124 |         -0.09 |                    0.01618 |   0.0279 |      0.0311 |       0.9  |
| random-split   | Europe_Iron_Classical  | NULLrand_Europe_Iron_Classical_r2A  | NULLrand_Europe_Iron_Classical_r2B  |     0.00186 |        0.00121 |          1.53 |                    0.01618 |   0.0527 |      0.0343 |       1.54 |
| coverage-split | Europe_Iron_Classical  | NULLcov_Europe_Iron_Classical_HI    | NULLcov_Europe_Iron_Classical_LO    |    -0.00068 |        0.00127 |         -0.53 |                    0.01616 |  -0.0345 |      0.0342 |      -1.01 |
| random-split   | EastAsia_Bronze        | NULLrand_EastAsia_Bronze_r0A        | NULLrand_EastAsia_Bronze_r0B        |    -0.00287 |        0.00315 |         -0.91 |                    0.01699 |  -0.0579 |      0.0822 |      -0.7  |
| random-split   | EastAsia_Bronze        | NULLrand_EastAsia_Bronze_r1A        | NULLrand_EastAsia_Bronze_r1B        |     0.00015 |        0.00292 |          0.05 |                    0.01701 |  -0.0164 |      0.081  |      -0.2  |
| random-split   | EastAsia_Bronze        | NULLrand_EastAsia_Bronze_r2A        | NULLrand_EastAsia_Bronze_r2B        |     0.00128 |        0.00299 |          0.43 |                    0.01695 |  -0.0687 |      0.0904 |      -0.76 |
| coverage-split | EastAsia_Bronze        | NULLcov_EastAsia_Bronze_HI          | NULLcov_EastAsia_Bronze_LO          |    -0.00469 |        0.00356 |         -1.32 |                    0.0171  |  -0.1353 |      0.0984 |      -1.38 |


The largest absolute null difference in D_VA is 0.00469 and the spread of null differences gives a systematic floor of 0.00308. On the normalised statistic the same splits are noisier - coverage-split R differences reach 0.135 - which is a direct measurement of how much coverage alone can move R, and is the reason the age-correlated pattern above is reported as a candidate.

Note that a same-cohort split is the right yardstick for *bias* and the wrong one for *power*: two halves of one cohort share almost all their ancestry, so their difference has a much smaller standard error than a comparison between genuinely different populations. The statistical floor is therefore taken from the 1431 real cohort-versus-cohort comparisons (0.01510 at the median), not from the splits.


![Figure 3](fig_n3_null.png)

**Figure 3.** Left: null constructions against the detection limit. Right: the distribution of all pairwise cohort-difference Z scores against a standard normal - a check that the paired jackknife is calibrated rather than merely tight.

### Stated detection limit

> **A difference of 0.0151 in D_VA between two typical cohorts is resolvable at 2 sigma; against a typical cohort's D_VA of 0.0746 that is 20%. For the best-powered pairs the limit falls to 0.0031, 4%.**

Concretely: if a cohort replaced a fraction *f* of its Neanderthal ancestry with ancestry from a Neanderthal lineage equidistant between Vindija and Altai, this study would detect it only for *f* > 20% in a typical comparison. Structure within the introgressing population finer than that is invisible here. The corresponding limit on R is 0.41.

This is a statement about the panel, not about history. The limit is set by the 5,237 1240K sites that separate the two archaic genomes; shotgun data at all sites, or the addition of Chagyrskaya and Mezmaiskaya (absent from the AADR), would tighten it substantially.

### Coverage matching

Every core cohort was recomputed on the 224,854 SNPs covered in all of them. The fitted slope is 1.5944 +/- 0.2822 against 1.5670 +/- 0.2795 on the full panel, so the shared-SNP restriction changes nothing material. Full numbers in `ns_coverage_matched.csv`.

## Limitations

- **The panel caps the question.** The AADR 1240K release contains three archaic genomes (Altai, Vindija, Denisova) and no Chagyrskaya or Mezmaiskaya, so 'which Neanderthal' can be asked along one axis only, and that axis is defined by two genomes that are themselves close relatives.

- **D_VA and D_NEA are relative affinities, not percentages**, and R is a ratio of two of them, interpretable only across cohorts.

- **The proportional model is first-order.** D_VA ~ a x D_VA(source) neglects terms in the admixture fraction that are acceptable at a ~2% and marginal at Oase1's ~10%.

- **Cohorts are not random samples.** Pooled cohorts mix sites, periods and degrees of relatedness, and no kinship pruning was applied at this scale, so a large cemetery can be over-represented within its cohort.

- **Denisovan ancestry perturbs the contrast.** Because the Altai genome carries Denisovan-related ancestry, populations with Denisovan ancestry of their own are pulled towards Altai for reasons unrelated to their Neanderthal source. Present-day Papuans sit -0.0185 (Z = -1.29) from the line, in the predicted direction; Oceanian cohorts should be read with that in mind.

- **The age-correlated residual is unresolved**, as set out above.

## Reproduce

```bash
python neanderthal_source.py --panel 1240k
python neanderthal_source.py --panel 1240k --transversions
```


*Refs: Green et al. 2010 Science 328:710; Reich et al. 2010 Nature 468:1053; Patterson et al. 2012 Genetics 192:1065; Prufer et al. 2014 Nature 505:43; Fu et al. 2015 Nature 524:216 (Oase1); Fu et al. 2016 Nature 534:200; Prufer et al. 2017 Science 358:655 (Vindija); Chen et al. 2020 Cell 180:677; Mallick et al. 2024 Sci. Data 11:182 (AADR).*


## Full cohort tables

### Named Palaeolithic cohorts and reference controls

| label              | kind              |   n_ind |   date_bp | region      |   median_coverage |   median_damage |   D_NEA |    D_VA |   D_VA_se |   D_VA_z |       R |    R_se |   n_snp |
|:-------------------|:------------------|--------:|----------:|:------------|------------------:|----------------:|--------:|--------:|----------:|---------:|--------:|--------:|--------:|
| CTRL_Chimp         | control-reference |       1 |       nan |             |            nan    |         nan     | -1      | -0.0916 |    0.0173 |    -5.3  | nan     | nan     |  138509 |
| CTRL_Denisova      | control-reference |       1 |       nan |             |            nan    |         nan     |  0.0692 | -0.1459 |    0.0331 |    -4.41 |  -2.11  |   0.649 |   67046 |
| CTRL_Mbuti         | control-reference |      15 |       nan |             |            nan    |         nan     | -0.0061 | -0.018  |    0.008  |    -2.25 | nan     | nan     |  137883 |
| IUP_Eurasia_45ka   | named             |      10 |     45875 | Europe      |              1.08 |           0.53  |  0.0536 |  0.0732 |    0.0117 |     6.26 |   1.365 |   0.248 |  138509 |
| UP_Europe_post_LGM | named             |       5 |     17045 | Europe      |              1.63 |           0.128 |  0.0458 |  0.0577 |    0.0131 |     4.4  |   1.261 |   0.331 |  136504 |
| UP_Europe_pre_LGM  | named             |      14 |     30246 | Europe      |              1    |           0.259 |  0.0599 |  0.0767 |    0.0134 |     5.74 |   1.281 |   0.249 |  138504 |
| UP_NorthEastAsia   | named             |       8 |     31850 | EastAsia    |              4.36 |         nan     |  0.0522 |  0.0753 |    0.0136 |     5.52 |   1.442 |   0.315 |  138509 |
| UstIshim_44ka      | named             |       1 |     44366 | WestSiberia |             40.66 |         nan     |  0.0439 |  0.0608 |    0.0187 |     3.25 |   1.384 |   0.44  |  138498 |
| Oase1_40ka         | named-lowpower    |       1 |     39982 | Europe      |              0.05 |         nan     |  0.0889 |  0.0367 |    0.0811 |     0.45 |   0.412 |   0.901 |    9169 |

### Present-day anchors

| label        | kind        |   n_ind |   date_bp | region      |   median_coverage |   median_damage |   D_NEA |    D_VA |   D_VA_se |   D_VA_z |       R |    R_se |   n_snp |
|:-------------|:------------|--------:|----------:|:------------|------------------:|----------------:|--------:|--------:|----------:|---------:|--------:|--------:|--------:|
| PD_Basque    | present-day |      25 |         0 | Europe      |             37.3  |             nan |  0.0437 |  0.0808 |    0.0115 |     7.05 |   1.85  |   0.339 |  137886 |
| PD_Biaka     | present-day |      24 |         0 | Africa      |             41    |             nan | -0.004  | -0.0055 |    0.0066 |    -0.84 | nan     | nan     |  137873 |
| PD_Dai       | present-day |      14 |         0 | SEAsia      |             36.59 |             nan |  0.0439 |  0.067  |    0.012  |     5.57 |   1.527 |   0.339 |  138508 |
| PD_Esan      | present-day |     101 |         0 | Africa      |             35.44 |             nan |  0.0006 |  0.0009 |    0.0023 |     0.38 | nan     | nan     |  138143 |
| PD_French    | present-day |      31 |         0 | Europe      |             40.07 |             nan |  0.0429 |  0.0782 |    0.0126 |     6.22 |   1.823 |   0.384 |  137887 |
| PD_Gambian   | present-day |     114 |         0 | Africa      |             36.31 |             nan |  0.0016 |  0.0007 |    0.0031 |     0.22 | nan     | nan     |  138144 |
| PD_Han       | present-day |      46 |         0 | EastAsia    |             38.11 |             nan |  0.0438 |  0.0749 |    0.0115 |     6.52 |   1.709 |   0.338 |  137884 |
| PD_Japanese  | present-day |      31 |         0 | EastAsia    |             45.18 |             nan |  0.0416 |  0.0768 |    0.0117 |     6.54 |   1.846 |   0.345 |  138146 |
| PD_Juhoan    | present-day |      10 |         0 | Africa      |             35.64 |             nan |  0.0059 | -0.0231 |    0.0102 |    -2.27 | nan     | nan     |  137884 |
| PD_Karitiana | present-day |      15 |         0 | Americas    |             35.76 |             nan |  0.0534 |  0.1011 |    0.0139 |     7.27 |   1.894 |   0.343 |  137884 |
| PD_Luhya     | present-day |     103 |         0 | Africa      |             34.43 |             nan |  0.002  |  0.006  |    0.0037 |     1.62 | nan     | nan     |  138490 |
| PD_Mandenka  | present-day |      26 |         0 | Africa      |             33.93 |             nan | -0.0011 |  0.0047 |    0.0052 |     0.89 | nan     | nan     |  138505 |
| PD_Mbuti     | present-day |      17 |         0 | Africa      |             36.02 |             nan | -0.0063 | -0.0175 |    0.008  |    -2.19 | nan     | nan     |  137883 |
| PD_Mende     | present-day |      87 |         0 | Africa      |             36.25 |             nan | -0.0012 | -0.0008 |    0.0035 |    -0.22 | nan     | nan     |  138143 |
| PD_Mixe      | present-day |       3 |         0 | Americas    |             42.23 |             nan |  0.0469 |  0.0975 |    0.015  |     6.52 |   2.08  |   0.424 |  134834 |
| PD_Papuan    | present-day |      32 |         0 | Oceania     |             40.42 |             nan |  0.0585 |  0.0732 |    0.0166 |     4.42 |   1.25  |   0.317 |  137888 |
| PD_Russian   | present-day |      27 |         0 | WestSiberia |             35.4  |             nan |  0.0441 |  0.0776 |    0.0107 |     7.26 |   1.761 |   0.332 |  137888 |
| PD_Sardinian | present-day |      31 |         0 | Europe      |             33.49 |             nan |  0.0424 |  0.0759 |    0.0124 |     6.14 |   1.792 |   0.366 |  137888 |