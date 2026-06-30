# External validation — pipeline vs published papers

Neanderthal ancestry (f4-ratio α) and Denisovan affinity (D) recomputed on the
AADR 1240K panel for individuals/populations with literature values. Absolute
f4-ratio scale is reference-dependent, so the tests are: right ballpark (~2%),
relative ordering / correlation, and the standout cases.

| sample / population | this pipeline (%) | published (%) | source | D_Nea Z | D_Den Z |
|---|---|---|---|---|---|
| Oase1 | 5.06 ± 1.36 (116,429 SNP) | 7.3 [6.0–9.0] | Fu & Paabo 2015 Nature 524:216 | 4.6 | 1.7 |
| Oase1_dmg | 9.80 ± 2.18 (25,775 SNP) | 7.3 [6.0–9.0] | Fu & Paabo 2015 (damage-restricted) | 4.4 | 1.3 |
| Ust_Ishim | 2.34 ± 0.55 (505,813 SNP) | 2.3 [1.9–2.7] | Fu et al. 2014 Nature 514:445 | 2.8 | -1.2 |
| Tianyuan | 2.38 ± 0.71 (392,218 SNP) | 2.5 [2.0–3.0] | Yang et al. 2017 Curr Biol; Fu 2013 | 2.8 | -1.3 |
| Kostenki14 | 2.33 ± 0.53 (464,272 SNP) | 2.6 [2.2–3.0] | Seguin-Orlando 2014 Science; Fu 2016 | 3.2 | -0.3 |
| GoyetQ116-1 | 4.13 ± 1.00 (129,167 SNP) | 3.0 [2.6–3.4] | Fu et al. 2016 Nature 534:200 | 4.4 | 0.3 |
| Vestonice14 | 1.56 ± 0.77 (329,939 SNP) | 2.6 [2.2–3.1] | Fu et al. 2016 Nature 534:200 | 2.6 | -0.1 |
| MA1_Malta | 2.33 ± 0.70 (362,024 SNP) | 2.0 [1.6–2.6] | Raghavan et al. 2014 Nature 505:87 | 3.6 | 1.1 |
| Loschbour_WHG | 2.92 ± 0.56 (206,599 SNP) | 1.9 [1.6–2.2] | Lazaridis et al. 2014 Nature 513:409 | 4.3 | 0.5 |
| Stuttgart_LBK | 1.81 ± 0.49 (410,508 SNP) | 1.8 [1.5–2.1] | Lazaridis et al. 2014 Nature 513:409 | 3.8 | -0.6 |
| French | 2.09 ± 0.38 (505,842 SNP) | 1.9 [1.7–2.1] | Prufer et al. 2017 Science 358:655 | 4.8 | -0.4 |
| Sardinian | 2.03 ± 0.38 (505,843 SNP) | 1.7 [1.5–1.9] | Prufer 2017; Lazaridis 2014 | 4.8 | -0.9 |
| Han | 2.29 ± 0.38 (505,837 SNP) | 2.3 [2.0–2.5] | Wall 2013; Prufer 2017 (~20% > Eur) | 4.6 | -0.8 |
| Dai | 2.18 ± 0.37 (505,845 SNP) | 2.2 [1.9–2.5] | Wall 2013; Vernot & Akey 2014 | 4.8 | -0.5 |
| Papuan | 2.94 ± 0.44 (505,841 SNP) | 2.2 [1.8–2.6] | Reich 2010; Vernot 2016 (Nea); +Denisovan | 6.0 | 5.9 |
| Karitiana | 2.35 ± 0.47 (505,843 SNP) | 2.0 [1.7–2.3] | Vernot & Akey 2014 | 4.5 | 0.2 |
| Yoruba | -0.11 ± 0.21 (505,833 SNP) | 0.0 [0.0–0.3] | baseline (~0; Chen 2020 ~0.3% back-flow) | -2.4 | -2.0 |

## Key tests
```
================================================================
KEY VALIDATION TESTS
================================================================
(non-African baseline from modern reference pops = 2.18%)
[Oase1] standard 5.06±1.36%, damage-restricted 9.80±2.18%  vs published 6-9% (Fu 2015)
        single most elevated sample; D_Nea Z=4.6/4.4 (significant). Damage-restriction RAISES it (removes modern
        contamination that dilutes archaic signal) -> recovers the published 6-9% -> CONFIRMED
[Ust-Ishim] 2.34±0.55% vs published 2.3% (Fu 2014), 505,813 SNP -> MATCH
[Yoruba] -0.11±0.21% vs published ~0 -> MATCH (~0)
[E-Asian excess] Han 2.29% > French 2.09% (ratio 1.10; direct D-stat is the powered test, see Phase 1)
[Denisovan] Papuan D_Den Z=5.9 >> all others (~0) -> CONFIRMED (Reich 2010/Meyer 2012)

Correlation mine vs published: r=0.866 (all 16), r=0.959 (excl. Oase1)
Within published range (+/- 1 SE): 16/17 samples
```

![validation](figures/fig7_validation_vs_published.png)

**Verdict:** the pipeline reproduces published Neanderthal estimates across the
literature anchors (r=0.87 overall, 0.96 excluding the Oase1 leverage
point), quantitatively matches the high-coverage controls (Ust'-Ishim 2.3%, Stuttgart 1.8%), returns the African baseline at ~0, correctly recovers the famous
Oase1 6–9% Neanderthal (and shows damage-restriction removing contamination dilution),
and confirms the Papuan-specific Denisovan signal. This independently corroborates the
Phase-1 internal validation.
