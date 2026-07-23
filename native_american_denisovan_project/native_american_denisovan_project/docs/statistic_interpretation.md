# Candidate primary statistics — expected sign and interpretation

This document derives the expected sign of each proposed statistic **before** it is used on real data, per the project rule that no statistic enters results until its interpretation survives synthetic tests. It also records one statistic that **failed** empirical sign-checking in the feasibility probe and is therefore excluded.

## Sign conventions (the antisymmetry trap)

D(W,X;Y,Z) is antisymmetric: `D(W,X;Y,Z) = -D(X,W;Y,Z) = -D(W,X;Z,Y)` (Patterson et al. 2012). The project spec writes the basic Denisovan statistic as `D(African, X; Denisovan, Chimp)`, i.e. `D(Mbuti, X; Denisova, Chimp)`. The validated production pipeline (`denisovan_survey.py`) uses the opposite first-pair order, `D(X, Mbuti; Denisova, Chimp)`, because that makes the Papuan positive control **positive** (Papuan ≈ +3.45% in the existing survey). The two are exact negatives:

```
D_pipeline(X) = D(X, Mbuti; Denisova, Chimp) = - D_spec(X) = - D(Mbuti, X; Denisova, Chimp)
```

Throughout this project we use the **pipeline convention** (positive = more Denisovan sharing than the Mbuti baseline) and state the conversion explicitly on every table. The feasibility probe confirmed Papuan = +0.0345 (Z=5.92) on all sites and +0.0464 (Z=5.97) on transversions only, reproducing the published positive control.

Per-SNP algebra: `D(W,X;Y,Z)` numerator = `(pW-pX)(pY-pZ)`, denominator = `(pW+pX-2pWpX)(pY+pZ-2pYpZ)`. The D-statistic is **polarisation-invariant** (flipping the counted allele at a SNP sends every p→1−p, leaving (pW−pX)(pY−pZ) and the denominator unchanged), so we never need an ancestral-allele call for D; we only need it for the diagnostic-marker orientation.

## Statistic 1 (Tier 1) — basic Denisovan affinity

`S1(X) = D(X, Mbuti; Denisova, Chimp)`

- **Per-SNP**: `(pX - pMbuti)(pDenisova - pChimp)`.
- **Expected sign**: positive where X carries alleles the Denisovan shares against chimp that Mbuti lacks — i.e. populations with Denisovan ancestry. Papuan strongly positive (verified). French ≈ 0 (verified −0.0016 all sites). East Asians weakly positive in principle (~0.2% Denisovan) but in practice ~0 because the statistic conflates Denisovan ancestry with **shared ancestral polymorphism + Neanderthal sharing** that all non-Africans carry.
- **Null hypothesis**: X and Mbuti share the Denisovan-derived allele equally (no excess Denisovan affinity).
- **What a non-zero value means**: X has *some* excess allele-sharing with Denisova relative to Mbuti — not necessarily Denisovan introgression (could be ancestral polymorphism, Neanderthal–Denisovan shared derived alleles, reference bias, or drift).
- **Probe result (all sites)**: Papuan Z=5.92 (control); Karitiana Z=0.18; USR1 Z=−0.48; Kolyma1 Z=1.33; MA1 Z=1.06; Yana1 Z=0.44; Sumidouro6 Z=0.24. **No ancient American/Beringian is significant.** This replicates the production survey's controlled near-null for *total* Denisovan affinity.
- **Weakness**: it measures total Denisovan-related sharing, not the *composition* (Papuan-like vs Han-like). It cannot, by itself, answer the project's central question.
- **Use in project**: baseline / positive-control scale only; never as the deciding statistic.

## Statistic 2 (Tier 2/3) — denoised Denisovan excess vs a non-Denisovan non-African baseline

`S2(X) = D(X, French; Denisova, Chimp)`

- **Per-SNP**: `(pX - pFrench)(pDenisova - pChimp)`.
- **Rationale**: French carries the shared non-African background (Neanderthal introgression + ancestral polymorphism common to all out-of-Africa populations) but ~0 Denisovan ancestry. Differencing X against French therefore removes the shared non-African baseline that contaminates S1, leaving X's **Denisovan excess** as the dominant signal. This is the right "is there any Denisovan ancestry in X at all" test on AADR data.
- **Expected sign**: positive where X has Denisovan ancestry beyond the shared non-African background. Papuan strongly positive; East Asians small-positive (their ~0.2%); French itself = 0 by construction (self-comparison); a Native American with no Denisovan ancestry → ~0; a Native American carrying *any* Denisovan component → positive.
- **Null**: X and French share Denisovan-derived alleles equally.
- **Caveat**: still conflates the two Denisovan components (it cannot tell Papuan-like from Han-like); it is a *total-excess* test, conditioned on the non-African background. It is the appropriate replacement for the naive Statistic 3 (below), which failed.
- **Use in project**: the principal "is there detectable Denisovan ancestry in Native Americans at all" statistic, run all-sites and transversions-only, across the time transect.

## Statistic 3 (Tier 5, primary) — composition: Papuan-associated vs East-Asian-associated diagnostic-marker contrast

`S3(X) = corr( f_X[Set_A] , f_Papuan[Set_A] ) - corr( f_X[Set_A] , f_Han[Set_A] )`

restricted to the **Set A** diagnostic sites (Denisova high, Neanderthal low, African low), where `f_X[Set A]` is X's oriented Denisovan-allele frequency vector. Equivalently a regression of X's Set-A frequencies onto Papuan and Han profiles with ANE and East-Asian ancestry fractions as covariates (Section 14). The primary, preregistered contrast (Section 25) is:

> After conditioning on East-Asian-related and Ancient-North-Eurasian-related ancestry, Native American populations show greater sharing with Papuan-associated Denisovan markers than with East-Asian-associated Denisovan markers.

- **Expected sign under the hypothesis (Papuan-like component in Native Americans)**: positive (X's profile is more Papuan-correlated than Han-correlated at Denisovan-diagnostic sites, after conditioning).
- **Expected sign under the null (single East-Asian-like component)**: ≈ 0 (X's profile tracks Han, conditional on ancestry).
- **Expected sign under "dilution/structure"**: depends on the covariates; this is why conditioning on East-Asian + ANE ancestry is mandatory, not optional.
- **Null**: the Papuan-association and Han-association contrasts are equal after conditioning.
- **What a positive value means**: X's Denisovan-marker profile more closely resembles Papuans' than Han's — **allele-sharing similarity**, NOT "X descends from Papuans" and NOT "X has p% Papuan ancestry." Per the language discipline (Section 30): "greater affinity to markers enriched in Papuan Denisovan-derived segments."
- **Circularity guard**: Set A is defined from Denisova/Neanderthal/African references only — **not** from Papuans or Han — so using it to compare X's Papuan- vs Han-similarity is not circular. The Papuan-associated / East-Asian-associated *subsets* of Set A (Sets D and E, Tier 5) must be defined on a **training partition** and tested on a held-out partition, or taken from published tract calls, to avoid defining markers with Papuans and then "finding" Papuan enrichment.
- **Use in project**: the deciding statistic. Run with block-bootstrap CIs, leave-one-chromosome-out, transversions-only, and selected-locus exclusions before any interpretation.

## Statistic that FAILED validation — excluded from results

`S_bad(X) = D(Altai, Denisova; X, Mbuti)` — the "obvious" Neanderthal-conditioned form one might reach for to isolate Denisovan-specific sharing.

- **Per-SNP**: `(pAltai - pDenisova)(pX - pMbuti)`.
- **Expected under the naive intuition** (what one might *hope* it measures): positive where X shares alleles with Denisova that it does not share with Neanderthal — i.e. Denisovan-specific ancestry.
- **Empirical probe result (all sites)** — this is why it is excluded:
  - French  = +0.0357, Z = 6.33
  - Han     = +0.0389, Z = 6.51
  - Karitiana = +0.0345, Z = 5.67
  - Papuan  = +0.0035, Z = 0.44   ← **the Papuan positive control is NULL**
  - Yoruba  = +0.0017, Z = 0.47
- **What this means**: the statistic is dominated by **Neanderthal ancestry in X**, not Denisovan ancestry. At sites where Altai (Neanderthal) differs from Denisova, a population with Neanderthal ancestry (every non-African, including French) has `pX - pMbuti > 0`, giving a large positive value; a population with **Denisovan** ancestry (Papuan) is pulled *toward Denisova*, **reducing** `pAltai - pDenisova`-weighted sharing, so the statistic is **suppressed** in Papuans — the opposite of "isolating Denisovan ancestry." French (0% Denisovan) scores higher than Papuan (5% Denisovan). This is a decisive sign failure.
- **Conclusion**: `D(Altai, Denisova; X, Mbuti)` is a **Neanderthal indicator** under the project's sign convention, **anti-correlated with Denisovan ancestry**, and is excluded from all results. This is exactly the failure mode the project spec warns about (Section 10: "Do not include a statistic in final results unless its interpretation passes these tests"). It is retained here as a documented negative result and as a synthetic-test case (see `synthetic_validation_design.md`).

## Summary table

| ID | statistic | expected + | probe check | role |
|----|-----------|-----------|-------------|------|
| S1 | D(X, Mbuti; Denisova, Chimp) | Papuan large +; French ~0; EAs ~0 | Papuan Z=5.9, French Z=−0.4, Karitiana Z=0.2 | positive-control scale; total affinity only |
| S2 | D(X, French; Denisova, Chimp) | Papuan large +; EAs small +; French 0 | to be run | denoised Denisovan-excess detection |
| S3 | corr(f_X, f_Papuan) − corr(f_X, f_Han) on Set A, conditional | + if Papuan-like; 0 if Han-like | to be run on partitioned Set A | **primary** composition test |
| S_bad | D(Altai, Denisova; X, Mbuti) | (intuition: + for Denisovan) — **wrong** | French Z=6.3, Papuan Z=0.4 | **excluded** (Neanderthal indicator) |

All four will be locked in `tests/test_dstat_sign.py` against synthetic genotype matrices with known introgression before any real-data result is reported.
