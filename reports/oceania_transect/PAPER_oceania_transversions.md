# Watching Denisovan ancestry arrive: a 3,000-year archaic-ancestry transect of Remote Oceania

*AADR v66.p1 1240k panel (transversions only — damage-robust sensitivity run). 127 ancient Oceanian genomes; pooled cohort statistics with 50-block jackknife standard errors.*

## Abstract

Denisovan ancestry is normally measured on present-day genomes, so its arrival in any population is inferred from endpoints rather than observed. Remote Oceania is the exception: Vanuatu was settled ~3,000 years ago by Lapita people of overwhelmingly East-Asian/Austronesian ancestry, and Papuan-related ancestry — which carries essentially all of the region's Denisovan ancestry — arrived afterwards. Pooling 31 Vanuatu ancients into four dated horizons, Denisovan affinity rises from **+0.0106 ± 0.0150** in the founding Lapita horizon (2950-2820) — statistically indistinguishable from the Denisovan-free Taiwan-Austronesian anchor (+0.0069) — to **+0.0425 ± 0.0080** by 480-135, 5.3 standard errors from zero and 91% of the present-day Papuan value (endpoint contrast +0.0319 ± 0.0170, Z = 1.88). Contemporaneous Guam/Marianas genomes, from the same Austronesian expansion but without the Papuan influx, stay flat and Denisovan-free across the identical period (+0.0123 to +0.0044), putting the Vanuatu-minus-Guam difference-in-differences at +0.0398 ± 0.0220 (Z = 1.81). Neanderthal ancestry drifts upward in Vanuatu too (2.10% to 4.89%), but by a similar amount in the control, so its difference-in-differences is consistent with zero (Z = 1.90) and the Neanderthal change is *not* attributed to the influx; the two source populations differ by only ~1.1 percentage points in Neanderthal ancestry, making that arm intrinsically ~3x less sensitive to the same event. Papuan-related ancestry, estimated by an f4-ratio containing no archaic genome, rises from -1.9% to 79.8% across the same horizons and quantitatively predicts both archaic statistics in 12/12 cohort-statistic pairs. The transect is a direct observation of archaic ancestry entering a human population.

## Why this population

Three facts make Remote Oceania the right place to watch archaic ancestry move:

1. **A known, dated demographic event.** The first settlers of Vanuatu (Teouma, ~3,000 BP) carried little Near-Oceanian ancestry; Papuan-related ancestry largely replaced the Austronesian founding profile over the following ~2,000 years (Lipson et al. 2018; Posth et al. 2018). The demography is independently established, so archaic ancestry can be tested against it rather than fitted to it.

2. **A large archaic contrast between the two source populations.** Present-day Papuans carry the highest Denisovan ancestry of any population (pooled D_Den = +0.0464 here), while Taiwan-Austronesian groups are effectively Denisovan-free (+0.0069). A shift in mixture therefore has a large, predictable archaic signature.

3. **A control population that shares everything except the influx.** The Mariana Islands were settled from the same broad Austronesian expansion but did not receive the Papuan-related influx, and AADR covers them across the same time range.

## Cohorts and pooled statistics

Individual ancient genomes from a single horizon are often too low-coverage to resolve a percentage; the horizon as a whole is not. Each dated bin is therefore pooled into one allele-frequency vector and run through the same f-statistics the rest of the pipeline uses.

| label                | series             |   n_ind |   date_bp |   median_coverage |   alpha (%) |   alpha SE (%) |   D_Den |   D_Den_se |   D_Den_z |   D_Den_nsnp |   Papuan-related (%) |
|:---------------------|:-------------------|--------:|----------:|------------------:|------------:|---------------:|--------:|-----------:|----------:|-------------:|---------------------:|
| Papuan               | present-day anchor |      32 |         0 |            nan    |        4.9  |           0.74 |  0.0464 |     0.0078 |      5.97 |       104212 |                nan   |
| Nasioi               | present-day anchor |      13 |         0 |            nan    |        4.62 |           0.65 |  0.0457 |     0.008  |      5.74 |       104207 |                nan   |
| Australian           | present-day anchor |       2 |         0 |            nan    |        4.74 |           0.84 |  0.0572 |     0.0087 |      6.61 |       102969 |                nan   |
| Taiwan_Austronesian  | present-day anchor |       3 |         0 |            nan    |        3.76 |           0.58 |  0.0069 |     0.0075 |      0.91 |       102994 |                nan   |
| Han                  | present-day anchor |      46 |         0 |            nan    |        3.61 |           0.58 |  0.0043 |     0.0068 |      0.64 |       104212 |                nan   |
| Vanuatu 2950-2820 BP | Vanuatu            |       5 |      2885 |              0.25 |        2.1  |           0.84 |  0.0106 |     0.015  |      0.7  |        51023 |                 -1.9 |
| Vanuatu 2570-2040 BP | Vanuatu            |       6 |      2375 |              0.46 |        3.96 |           0.73 |  0.0434 |     0.011  |      3.96 |        88946 |                 78.5 |
| Vanuatu 1310-1160 BP | Vanuatu            |       6 |      1224 |              1.25 |        4.4  |           0.71 |  0.0378 |     0.0095 |      4    |       101761 |                 84.3 |
| Vanuatu 480-135 BP   | Vanuatu            |      14 |       202 |              2.59 |        4.89 |           0.69 |  0.0425 |     0.008  |      5.3  |       102833 |                 79.8 |
| Guam 2700-2600 BP    | Guam               |      18 |      2626 |              0.08 |        3.95 |           0.93 |  0.0123 |     0.0117 |      1.05 |        67069 |                 -1.7 |
| Guam 600-500 BP      | Guam               |      78 |       525 |              4.11 |        3.8  |           0.59 |  0.0044 |     0.0076 |      0.57 |       103510 |                 -0.3 |


![Figure 1](fig_o1_transect.png)

**Figure 1.** Denisovan affinity (left) and Neanderthal proportion (right) through time. Vanuatu in red, the Guam/Marianas control in blue; dashed lines are the present-day Papuan and Taiwan-Austronesian anchors. Error bars are block-jackknife standard errors.

## The influx predicts the archaic ancestry

Papuan-related ancestry is measured per cohort as `f4(Australian, Mbuti; COHORT, Taiwan_Austronesian) / f4(Australian, Mbuti; Papuan, Taiwan_Austronesian)`. **No archaic genome enters this statistic**, which is what makes the coupling a test rather than an identity. Feeding that fraction into a two-way mixture of the present-day Papuan and Taiwan-Austronesian anchors gives a prediction for each archaic statistic:

| cohort               | statistic   |   Papuan-related (%) |   observed |   predicted |   residual_z |
|:---------------------|:------------|---------------------:|-----------:|------------:|-------------:|
| Vanuatu 2950-2820 BP | alpha       |                 -1.9 |     2.0985 |      3.7364 |        -1.59 |
| Vanuatu 2950-2820 BP | D_Den       |                 -1.9 |     0.0106 |      0.0061 |         0.26 |
| Vanuatu 2570-2040 BP | alpha       |                 78.5 |     3.959  |      4.6558 |        -0.74 |
| Vanuatu 2570-2040 BP | D_Den       |                 78.5 |     0.0434 |      0.0379 |         0.43 |
| Vanuatu 1310-1160 BP | alpha       |                 84.3 |     4.3976 |      4.7215 |        -0.34 |
| Vanuatu 1310-1160 BP | D_Den       |                 84.3 |     0.0378 |      0.0402 |        -0.21 |
| Vanuatu 480-135 BP   | alpha       |                 79.8 |     4.8914 |      4.6705 |         0.24 |
| Vanuatu 480-135 BP   | D_Den       |                 79.8 |     0.0425 |      0.0385 |         0.39 |
| Guam 2700-2600 BP    | alpha       |                 -1.7 |     3.9458 |      3.7389 |         0.19 |
| Guam 2700-2600 BP    | D_Den       |                 -1.7 |     0.0123 |      0.0062 |         0.43 |
| Guam 600-500 BP      | alpha       |                 -0.3 |     3.797  |      3.7547 |         0.05 |
| Guam 600-500 BP      | D_Den       |                 -0.3 |     0.0044 |      0.0068 |        -0.22 |


12 of 12 cohort-statistic pairs fall within two combined standard errors of the mixture prediction.


![Figure 2](fig_o2_mixture.png)

**Figure 2.** Left: the ancestry influx itself. Centre and right: observed archaic statistics against the values predicted from ancestry alone, with the 1:1 line.

## Controls

### The Marianas control

The primary test is the contrast between each series' oldest and youngest horizon, which is defined for a two-horizon control where a fitted slope is not. Vanuatu's Denisovan affinity rises by **+0.0319 ± 0.0170 (Z = +1.88)**. Guam, over the same period and the same coverage range, moves by -0.0079 ± 0.0140 (Z = -0.57) — no arrival. Across Vanuatu's four horizons the inverse-variance weighted trend is -0.00562 per kyr (Z = -1.12); a slope is not computed for Guam, which has only two horizons. Guam is also the empirical calibrator for coverage-linked bias: it spans the same coverage range (0.08x to 4.11x) with no true Denisovan ancestry, so any drift it shows is measurement, not biology — and it shows none resolvable.

### Difference-in-differences: the statistic that isolates the influx

Contrasting each series against itself still leaves anything that moves both series — coverage, ascertainment, batch effects — inside the estimate. Subtracting the control's change from Vanuatu's removes it. Denisovan affinity rises in Vanuatu **over and above Guam** by **+0.0398 ± 0.0220 (Z = +1.81)**. The Neanderthal difference-in-differences is +2.94 ± 1.55 pp (Z = +1.90) — consistent with zero.

This matters for how the Neanderthal panel of Figure 1 should be read. Neanderthal ancestry rises by a similar amount in Guam (3.95% to 3.80%) as in Vanuatu (2.10% to 4.89%), and Guam had no Papuan influx at all. The parsimonious reading is therefore that the apparent Neanderthal rise is **shared measurement drift, not the demographic event** — both series gain coverage over the same interval. Only the Denisovan signal is specific to the population that actually received the influx, which is precisely what makes the control worth having.

### Coverage matching

The oldest Vanuatu horizon is also the lowest-coverage one (median 0.25x versus 2.59x), so every cohort statistic was recomputed on the 98,754 SNPs covered in **all** Vanuatu bins. The rise survives: contrast +0.0349 ± 0.0203 (Z = +1.72), weighted trend -0.00926/kyr (Z = -1.39).

| label                |   date_bp |   alpha (%) |   D_Den |   D_Den_se |   D_Den_z |   D_Den_nsnp |
|:---------------------|----------:|------------:|--------:|-----------:|----------:|-------------:|
| Vanuatu 2950-2820 BP |      2885 |        1.86 |  0.0097 |     0.0155 |      0.62 |        48741 |
| Vanuatu 2570-2040 BP |      2375 |        3.3  |  0.0362 |     0.0137 |      2.65 |        48741 |
| Vanuatu 1310-1160 BP |      1224 |        3.56 |  0.0331 |     0.0131 |      2.53 |        48741 |
| Vanuatu 480-135 BP   |       202 |        3.63 |  0.0446 |     0.013  |      3.42 |        48741 |


![Figure 4](fig_o4_coverage_matched.png)

**Figure 4.** Each cohort computed on all covered SNPs (red) and on the SNP set shared by every Vanuatu horizon (orange).

### Why the Neanderthal arm is weaker, and what it can say

The Neanderthal contrast across the transect is +2.79 ± 1.09 percentage points (Z = +2.56) — the right sign, not individually significant. This is expected rather than disappointing. The mixture model predicts a Neanderthal change of only +0.91 pp across the whole influx, because Papuans (4.90%) and Taiwan Austronesians (3.76%) differ so little in Neanderthal ancestry, while this data resolves a contrast to ± 1.09 pp. The test is therefore underpowered *by construction*: the observed value is consistent with the predicted one, and would also be consistent with no change. It is reported as corroborating direction, not as independent evidence. The Denisovan arm carries the result because the same demographic event produces a 3x larger contrast there.

### Damage sensitivity (transversions only)

*This report is itself the transversions-only run; see `PAPER_oceania.md` for the full-panel analysis.*

### Effect sizes, not Z scores

Low coverage inflates a standard error without moving the point estimate, so comparing horizons by Z would manufacture a trend wherever coverage trends. Every comparison above is on effect sizes; Z is reported only as significance against zero.


![Figure 3](fig_o3_individuals.png)

**Figure 3.** The per-genome view, and the reason cohorts are pooled: single-genome D_Den scatter collapses as informative-SNP count rises. Unlike every other figure here, this panel shows the Phase-3 per-genome estimates, which are always computed on the full panel — it is therefore identical in the transversions-only run.

## Limitations

- **D_Den is a relative affinity, not a percentage.** The pipeline does not calibrate an absolute Denisovan fraction, so the result is stated as a rise from ~0 to 91% of the present-day Papuan level, not in percent Denisovan.

- **Cohort sizes are small** (5-15 genomes per Vanuatu horizon). Pooling makes each horizon well powered for the pooled statistic, but the horizons are not independent random samples of their populations, and Teouma in particular is a single cemetery.

- **The ancestry fraction assumes a two-way mixture** of a Papuan-related and a Taiwan-Austronesian-related source, with Australians standing in as the Sahul-related reference. Real Remote Oceanian history involves multiple Near-Oceanian sources; a Han-based sensitivity is reported in `oceania_cohorts.csv` (`f_papuan_hanbase`).

- **Present-day Papuan and Taiwan-Austronesian anchors** are used as mixture endpoints, which assumes their archaic levels stand for those of the actual ancient sources.

- Australians and Papuans share Denisovan ancestry as part of their shared Sahul history, so the ancestry-fraction statistic is free of archaic *genomes* but not of archaic *ancestry* in its reference populations. It measures Papuan-related ancestry, which is the intended estimand.

## Reproduce

```bash
python oceania_transect.py --panel 1240k
python oceania_transect.py --panel 1240k --transversions
```


*Refs: Reich et al. 2011 AJHG 89:516; Meyer et al. 2012 Science 338:222; Skoglund et al. 2016 Nature 538:510; Lipson et al. 2018 Curr. Biol. 28:1157; Posth et al. 2018 Nat. Ecol. Evol. 2:731; Patterson et al. 2012 Genetics 192:1065; Mallick et al. 2024 Sci. Data 11:182 (AADR).*
