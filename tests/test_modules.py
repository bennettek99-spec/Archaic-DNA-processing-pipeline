"""
Unit tests for the new modules (simulation, kinship, config) — runnable in CI
without the AADR data. Run: pytest -q
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from archaic import stats as st, kinship as kin, config as cfg


def test_simulate_recovers_alpha():
    sim = __import__("archaic.simulate", fromlist=["simulate"])
    freq, block, _, info = sim.simulate(0.04, n_test=15, n_afr=12, n_blocks=25, seed=3)
    a = st.f4_ratio(freq, "Altai", "Chimp", "X", "Mbuti", "Vindija", block, info["n_blocks"])
    assert abs(a["theta"] - 0.04) < 0.02          # recovers known truth within tolerance
    assert a["theta"] > 0.015                       # clearly non-zero


def test_simulate_null_is_zero():
    sim = __import__("archaic.simulate", fromlist=["simulate"])
    freq, block, _, info = sim.simulate(0.0, n_test=15, n_afr=12, n_blocks=25, seed=4)
    a = st.f4_ratio(freq, "Altai", "Chimp", "X", "Mbuti", "Vindija", block, info["n_blocks"])
    assert abs(a["theta"]) < 0.012                  # ~0 under no introgression


def test_kinship_haploid_resolution():
    G = np.array([[0, 2, -1], [1, 1, 0]])           # (snp, ind): 0/2 -> 0/1, het random, -1 kept
    A = kin._haploid_alleles(G, seed=0)
    assert A[0, 0] == 0 and A[0, 1] == 1 and A[0, 2] == -1
    assert A[1, 0] in (0, 1) and A[1, 2] == 0


def test_kinship_classify_thresholds():
    # synthetic P0: pair (0,1) identical, (0,2) unrelated; median sets the scale
    P0 = np.full((3, 3), np.nan)
    P0[0, 1] = P0[1, 0] = 0.10                       # very low -> identical
    P0[0, 2] = P0[2, 0] = 0.25
    P0[1, 2] = P0[2, 1] = 0.25                       # median ~0.25 -> norm(0,1)=0.4
    _, pairs = kin.classify(P0)
    degs = {(i, j): d for i, j, _, d in pairs}
    assert degs.get((0, 1)) == "identical/duplicate"
    assert (0, 2) not in degs and (1, 2) not in degs  # unrelated


def test_config_loads_panels():
    c = cfg.load_config()
    assert "1240k" in c["panels"] and "ho" in c["panels"]
    assert c["panels"]["1240k"]["snp_floor"] == 30000
