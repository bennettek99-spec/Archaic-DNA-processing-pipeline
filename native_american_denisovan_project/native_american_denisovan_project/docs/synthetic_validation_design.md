# Synthetic validation design

Before any real-data statistic is reported, the implementation must prove on **synthetic genotype matrices with known ground truth** that each statistic (a) has the expected sign, (b) is robust to the artifacts AADR data actually carries, and (c) has controlled false-positive and bias behaviour at the site counts real ancient genomes contribute. This document specifies the tests; the executable skeleton is `tests/test_dstat_sign.py`.

The rationale is concrete: the feasibility probe already found a statistic (`D(Altai, Denisova; X, Mbuti)`, the "obvious" Neanderthal-conditioned form) whose **intuitive** sign is wrong — it scores French (0% Denisovan) higher than Papuan (5% Denisovan). No statistic whose sign has not been confirmed on a synthetic matrix with a planted Denisovan pulse may appear in the results.

## A. Sign/orientation tests (mandatory gate)

Build a small synthetic genotype matrix (e.g. 2000 SNPs × 8 populations: Chimp, Denisovan-reference, Neanderthal-reference, African-outgroup, X_no_introgression, X_neanderthal_only, X_denisovan_only, X_both) with controlled allele frequencies. Plant introgression by copying a fraction α of the Denisovan-reference alleles into X_denisovan_only (and analogously for Neanderthal), so the ground-truth ancestry of every X is known by construction. Then assert:

| test | planted state | expected statistic outcome |
|------|--------------|------------------------------|
| T1 null | no archaic introgression in X | S1(X)≈0, S2(X)≈0, S3(X)≈0 within MC noise |
| T2 Neanderthal only | X has Neanderthal ancestry, 0 Denisovan | S1(X)≈0 (Neanderthal cancels), S2(X)≈0, S3(X)≈0; **S_bad(X) is large +** (proving it is a Neanderthal indicator) |
| T3 Denisovan only | X has Denisovan ancestry | S1(X)>0, S2(X)>0, S3(X) sign tracks which reference the planted Denisovan resembles; **S_bad(X) suppressed** (Papuan-style anti-correlation reproduced) |
| T4 both | X has both | S1,S2 positive; S3 dominated by the Denisovan component |
| T5 allele flip | count allele flipped on a subset of SNPs | all statistics unchanged (polarisation-invariance) — a regression test for the reader |

A statistic that fails its row is **excluded from the project** (as `S_bad` already is).

## B. Artifact-robustness tests

Using the same matrix generator, perturb the X populations and confirm the statistic is stable or fails predictably:

| artifact | how simulated | required behaviour |
|----------|---------------|------------------------|
| unequal missing data | randomly drop 10–60% of X genotypes | S1,S2 bias within ±1 SE; site count reported |
| pseudo-haploid sampling | force X to 0/1 (drop hets) | S1,S2 within ±1 SE of diploid truth; pseudo-haploid inflates S3 variance but not its mean |
| sequencing error | flip X alleles at rate 1e-3, 1e-2 | S1,S2 shift < 1e-3 at 1e-3; document breakdown at 1e-2 |
| reference bias | bias X toward the reference (counted) allele by 1–5% | S1 shifts positive; transversion-only reduces it; report as a known confound |
| ancestral-allele misidentification | deliberately mis-orient 5% of Set-A markers | S3 mean shifts by a bounded amount; flagged in limitations |
| ascertainment to 1240K | restrict synthetic WGS to the 1240K SNP list | S3 variance inflates as informative-site count drops; feeds the downsampling tests |

## C. Power / downsampling tests (Tier 16)

From a high-coverage synthetic genome with a *known* Denisovan fraction, project onto the real 1240K Set-A marker list, then **subsample** the diagnostic markers to {100, 500, 1000, 5000, 10000} and recompute S2 and S3. Record:

- bias(D) and SE at each count;
- false-positive rate (fraction of null-simulation runs with |Z|>3);
- power (fraction of α=0.2% Denisovan-simulation runs with |Z|>3) at each count;
- the minimum diagnostic-marker count at which S3 retains the planted direction.

These curves define the **reliable / provisional / insufficient** thresholds (Section 16) that gate whether a given ancient individual is even reportable. The feasibility probe already shows the *upper* end: ~12,180 Set-B markers exist on 1240K and each representative ancient individual carries ~9,000–12,180 of them, so power is limited by **signal size** (Native American Denisovan ancestry is small), not by marker count. The downsampling tests will quantify how small is still detectable.

## D. Population-model simulations (Tier 23, deferred to full workflow)

Use `archaic.simulate` (msprime) to generate whole genomes under Models 0–6 (Section 22) with a planted Denisovan pulse of known divergence, then project to 1240K and run S2/S3. Record power, false-positive rate, bias and CI coverage for **distinguishing Model 1 (single pulse) from Model 3 (Papuan-related source) from Model 5 (dilution)**. This is the step that converts "suggestive" into "supported" or "rejected"; it is part of the full workflow, **not** the feasibility stage, and is listed here only to fix the design contract.

## E. What the feasibility stage delivers vs defers

- **Feasibility stage (this stage):** the sign/orientation gate (Section A) and the artifact-robustness design (Section B) are **specified** here and a runnable skeleton (`tests/test_dstat_sign.py`) is provided; the `S_bad` failure is documented from real data as proof the gate is necessary.
- **Deferred to full workflow:** executing B–D at scale, the downsampling curves, and the msprime model grid. Per Section 35 of the project brief, the workflow stops at the feasibility report until the sign-gate and a minimum-viable S1+S2 run are reviewed.
