"""
Unit tests for archaic.qpadm (constrained + unconstrained) and archaic.ancestry
(source library / model competition) — synthetic data only, no AADR required.
Run: pytest -q  (or run this file directly; see tests/test_modules.py convention)
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from archaic import qpadm as qp, ancestry as anc


def _synthetic_freq(seed=0, n_snp=20000, missing=0.05):
    """Random allele frequencies for a small qpAdm system, with a Target that is
    an EXACT known linear mixture of S1/S2/S3 (so recovery can be checked)."""
    rng = np.random.default_rng(seed)
    names = ["S1", "S2", "S3", "R0", "R1", "R2", "R3"]
    freq = {n: np.clip(rng.beta(2, 2, n_snp), 0.001, 0.999) for n in names}
    true_w = np.array([0.5, 0.3, 0.2])
    freq["Target"] = true_w[0] * freq["S1"] + true_w[1] * freq["S2"] + true_w[2] * freq["S3"]
    for n in list(freq):
        miss = rng.random(n_snp) < missing
        freq[n] = freq[n].copy()
        freq[n][miss] = np.nan
    block = (np.arange(n_snp) * 50 // n_snp).astype(np.int32)
    return freq, block, true_w


def test_qpadm_recovers_exact_mixture():
    freq, block, true_w = _synthetic_freq(seed=1)
    r = qp.qpadm(freq, "Target", ["S1", "S2", "S3"], ["R0", "R1", "R2", "R3"], block, 50)
    assert np.allclose(r["weights"], true_w, atol=1e-6)
    assert r["feasible"]


def test_qpadm_constrained_matches_unconstrained_when_feasible():
    # when the unconstrained solution is already on the simplex, the constrained
    # fit should recover essentially the same weights (SLSQP finds the same optimum)
    freq, block, true_w = _synthetic_freq(seed=2)
    r_free = qp.qpadm(freq, "Target", ["S1", "S2", "S3"], ["R0", "R1", "R2", "R3"], block, 50)
    r_con = qp.qpadm_constrained(freq, "Target", ["S1", "S2", "S3"], ["R0", "R1", "R2", "R3"], block, 50)
    assert np.allclose(r_free["weights"], r_con["weights"], atol=1e-3)
    assert abs(r_con["weights"].sum() - 1.0) < 1e-6
    assert (r_con["weights"] >= -1e-9).all()


def test_qpadm_constrained_clamps_negative_weight_onto_simplex():
    # directly unit-test the solver on a system whose OLS solution is infeasible
    A = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    true = np.array([-0.4, 0.6])           # infeasible target mixture
    b = A @ true
    x = qp._solve_constrained(A, b)
    assert (x >= -1e-9).all()
    assert x.sum() <= 1 + 1e-9


def test_build_system_blockwise_matches_naive_bruteforce():
    """Regression guard: the vectorised block-sum system builder must reproduce
    the naive per-block-recompute f4 exactly (this is the perf-critical path)."""
    freq, block, _ = _synthetic_freq(seed=4, n_snp=8000)
    sources = ["S1", "S2", "S3"]; outgroups = ["R0", "R1", "R2", "R3"]
    sys_ = qp._build_system(freq, "Target", sources, outgroups, block, 50)

    S1, others = sources[0], sources[1:]
    R0, Rj = outgroups[0], outgroups[1:]
    need = ["Target"] + sources + outgroups
    mask = np.all([np.isfinite(freq[p]) for p in need], axis=0)

    def f4_naive(W, X, Y, Z, sel):
        a = (freq[W] - freq[X]) * (freq[Y] - freq[Z])
        return float(np.nanmean(a[sel]))

    b_naive = np.array([f4_naive("Target", S1, R0, r, mask) for r in Rj])
    A_naive = np.array([[f4_naive(s, S1, R0, r, mask) for s in others] for r in Rj])
    assert np.allclose(sys_["b"], b_naive, atol=1e-10)
    assert np.allclose(sys_["A"], A_naive, atol=1e-10)

    # spot-check one leave-one-block-out replicate
    bl = 5
    sel = mask & (block != bl)
    b_loo_naive = np.array([f4_naive("Target", S1, R0, r, sel) for r in Rj])
    assert np.allclose(sys_["b_loo"][bl], b_loo_naive, atol=1e-9)


def test_default_outgroups_excludes_model_sources():
    outs = anc.default_outgroups(["WHG", "EHG", "CHG", "Anatolia_N"])
    for s in ["WHG", "EHG", "CHG", "Anatolia_N"]:
        assert s not in outs
    # always-safe distal outgroups still present
    for base in anc.BASE_RIGHT:
        assert base in outs


def test_decompose_best_prefers_the_true_generating_model():
    """decompose_best() should rank the model that actually generated Target
    (S1/S2/S3, an exact mixture) above an unrelated 2-source model built from
    independent random populations. decompose_best() looks up the canonical
    BASE_RIGHT outgroup names internally, so the synthetic freq dict must
    include them (as independent random populations) for it to run at all."""
    freq, block, true_w = _synthetic_freq(seed=5)
    n_snp = len(freq["S1"])
    rng = np.random.default_rng(5)
    for name in list(anc.BASE_RIGHT) + ["Unrelated1", "Unrelated2"]:
        freq[name] = np.clip(rng.beta(2, 2, n_snp), 0.001, 0.999)
    models = {
        "true_model": ["S1", "S2", "S3"],
        "wrong_model": ["Unrelated1", "Unrelated2"],
    }
    ranked = anc.decompose_best(freq, "Target", models, block, 50)
    assert ranked[0]["model"] == "true_model"
    assert ranked[0]["ok"] and ranked[0]["free"]["feasible"]


def test_weights_dict_flattens_constrained_result():
    freq, block, true_w = _synthetic_freq(seed=6)
    result = anc.decompose(freq, "Target", ["S1", "S2", "S3"], ["R0", "R1", "R2", "R3"], block, 50)
    wd = anc.weights_dict(result, which="constrained")
    assert set(wd) == {"S1", "S2", "S3"}
    total = sum(w for w, _ in wd.values())
    assert abs(total - 1.0) < 1e-6


def test_qpwave_returns_rank_tests():
    freq, block, _ = _synthetic_freq(seed=7)
    rows = qp.qpwave(freq, ["Target", "S1", "S2"], ["R0", "R1", "R2", "R3"], block, 50)
    assert rows
    assert {r["rank"] for r in rows} == {0, 1}
    assert all(r["n_snp"] > 0 for r in rows)
    assert all(r["dof"] > 0 for r in rows)


def test_model_rejection_table_recovers_plausible_vs_rejected():
    """The true generating model (S1/S2/S3) must be labelled plausible and a
    wrong model (unrelated sources) rejected, in a per-target x model table."""
    freq, block, true_w = _synthetic_freq(seed=8)
    n_snp = len(freq["S1"])
    rng = np.random.default_rng(8)
    for name in list(anc.BASE_RIGHT) + ["Unrelated1", "Unrelated2"]:
        freq[name] = np.clip(rng.beta(2, 2, n_snp), 0.001, 0.999)
    models = {
        "true_model": ["S1", "S2", "S3"],
        "wrong_model": ["Unrelated1", "Unrelated2"],
    }
    rows = anc.model_rejection_table(freq, ["Target"], models, block, 50)
    by = {r["model"]: r for r in rows}
    assert by["true_model"]["status"] == "plausible"
    assert by["wrong_model"]["status"] == "rejected"
    assert all(r["n_snp"] > 0 for r in rows)


def test_model_rejection_table_skips_missing_targets():
    freq, block, _ = _synthetic_freq(seed=9)
    rows = anc.model_rejection_table(freq, ["NotPresent"], {"west3": ["S1", "S2", "S3"]},
                                     block, 50)
    assert rows == []


if __name__ == "__main__":
    tests = [v for k, v in list(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"OK {t.__name__}")
    print(f"\n{len(tests)}/{len(tests)} tests passed")
