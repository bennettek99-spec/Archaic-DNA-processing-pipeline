#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_dstat_sign.py — sign/orientation gate for the project's candidate statistics.

This is the executable skeleton for Section A of docs/synthetic_validation_design.md.
It builds a SYNTHETIC allele-frequency matrix with KNOWN planted introgression
(no AADR data required) and asserts that each statistic has the expected sign,
*before* any real-data result may be reported. It encodes the failure already
found in the feasibility probe: D(Altai, Denisova; X, Mbuti) is a Neanderthal
indicator (French > Papuan), not a Denisovan-isolating statistic.

Run (from the module root, using the pipeline venv):
  ..\\archaic-introgression\\.venv\\Scripts\\python.exe -m pytest tests\\test_dstat_sign.py -q
or directly:
  ..\\archaic-introgression\\.venv\\Scripts\\python.exe tests\\test_dstat_sign.py
"""
from __future__ import annotations
import os, sys
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
for _c in (os.path.join(os.path.dirname(_REPO), "archaic-introgression"), os.path.dirname(_REPO)):
    if os.path.isdir(os.path.join(_c, "archaic")):
        sys.path.insert(0, _c)
        break
from archaic import stats as st


def _d(freq, W, X, Y, Z, n_blocks=20):
    """Block-jackknife D on a synthetic freq dict (polarisation-invariant)."""
    block = st.assign_blocks(len(freq[W]), n_blocks)
    return st.dstat(freq, W, X, Y, Z, block, n_blocks)["theta"]


def synthetic_freqs(n_snp=4000, alpha_den=0.05, alpha_nea=0.02, seed=1):
    """Build freq arrays for Chimp, Den, Nea, Afr, X_null, X_nea, X_den, X_both.

    Introgression is planted by copying the archaic allele into X at fraction alpha.
    'Derived' is the archaic allele where it differs from chimp/african. The matrix
    is constructed so that the D-statistic sign is determined by construction."""
    rng = np.random.default_rng(seed)
    # ancestral (chimp/african) allele frequency baseline
    p_afr = rng.random(n_snp) * 0.5
    p_chimp = p_afr.copy()
    # archaic alleles: Denisovan and Neanderthal diverged from african at a subset
    archaic_site = rng.random(n_snp) < 0.5
    p_den = np.where(archaic_site, rng.random(n_snp) * 0.5 + 0.5, p_afr)
    p_nea = np.where(archaic_site, rng.random(n_snp) * 0.5 + 0.5, p_afr)
    # ensure den and nea differ on a subset (the 'denisovan-specific' sites)
    den_specific = archaic_site & (rng.random(n_snp) < 0.5)
    p_nea = np.where(den_specific, p_afr, p_nea)   # nea stays ancestral at den-specific
    nea_specific = archaic_site & (rng.random(n_snp) < 0.5) & ~den_specific
    p_den = np.where(nea_specific, p_afr, p_den)   # den stays ancestral at nea-specific

    def admix(p_base, p_arch, alpha):
        return (1 - alpha) * p_base + alpha * p_arch

    pX_null = p_afr.copy()
    pX_nea = admix(p_afr, p_nea, alpha_nea)
    pX_den = admix(p_afr, p_den, alpha_den)
    pX_both = admix(admix(p_afr, p_nea, alpha_nea), p_den, alpha_den)

    return {
        "Chimp": p_chimp, "Den": p_den, "Nea": p_nea, "Afr": p_afr,
        "X_null": pX_null, "X_nea": pX_nea, "X_den": pX_den, "X_both": pX_both,
    }


def test_S1_basic_denisovan_affinity_signs():
    """S1 = D(X, Afr; Den, Chimp): + for Denisovan-bearing, ~0 for null/Neanderthal."""
    f = synthetic_freqs()
    assert _d(f, "X_den", "Afr", "Den", "Chimp") > 0.01, "X_den should be +"
    assert abs(_d(f, "X_null", "Afr", "Den", "Chimp")) < 0.02, "X_null ~0"
    # Neanderthal-only should not read strongly + (the key non-confound)
    assert _d(f, "X_nea", "Afr", "Den", "Chimp") < 0.02, "X_nea must not inflate S1"


def test_S2_denoised_excess_signs():
    """S2 = D(X, WestEur-like; Den, Chimp). Use X_nea as the West-Eurasian baseline
    (Neanderthal, no Denisovan). X_den should be +; X_nea self ~0."""
    f = synthetic_freqs()
    s2_den = _d(f, "X_den", "X_nea", "Den", "Chimp")
    s2_null = _d(f, "X_null", "X_nea", "Den", "Chimp")
    assert s2_den > 0.01, "X_den should exceed the non-Denisovan non-African baseline"
    assert abs(s2_null) < 0.02, "X_null should be indistinguishable from baseline"


def test_S_bad_is_neanderthal_indicator_not_denisovan():
    """S_bad = D(Altai=Nea, Den; X, Afr) must NOT be used as a Denisovan statistic.
    Encodes the feasibility-probe failure: it is large for Neanderthal-bearing X and
    SUPPRESSED for Denisovan-bearing X. We assert the qualitative anti-pattern that
    disqualified it: X_nea >= X_den (Neanderthal scores higher than Denisovan)."""
    f = synthetic_freqs()
    sbad_nea = _d(f, "Nea", "Den", "X_nea", "Afr")
    sbad_den = _d(f, "Nea", "Den", "X_den", "Afr")
    assert sbad_nea > 0, "Neanderthal-bearing X gives large + on S_bad (it IS a Nea indicator)"
    assert sbad_den <= sbad_nea + 1e-9, (
        "S_bad is suppressed in Denisovan-bearing X relative to Neanderthal-bearing X; "
        "it cannot isolate Denisovan ancestry")


def test_polarisation_invariance():
    """Flipping the counted allele at every SNP must leave D unchanged (regression
    test for the reader/orientation logic)."""
    f = synthetic_freqs()
    D_before = _d(f, "X_den", "Afr", "Den", "Chimp")
    f2 = {k: 1.0 - v for k, v in f.items()}
    D_after = _d(f2, "X_den", "Afr", "Den", "Chimp")
    assert abs(D_before - D_after) < 1e-9, "D must be polarisation-invariant"


def test_sign_convention_antisymmetry():
    """D(X, Mbuti; Den, Chimp) == - D(Mbuti, X; Den, Chimp). Locks the convention
    conversion documented in statistic_interpretation.md."""
    f = synthetic_freqs()
    a = _d(f, "X_den", "Afr", "Den", "Chimp")
    b = _d(f, "Afr", "X_den", "Den", "Chimp")
    assert abs(a + b) < 1e-9


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failures = 0
    for fn in fns:
        try:
            fn(); print(f"PASS  {fn.__name__}")
        except AssertionError as e:
            failures += 1; print(f"FAIL  {fn.__name__}: {e}")
    print(f"\n{len(fns)-failures}/{len(fns)} sign-gate tests passed")
    sys.exit(1 if failures else 0)
