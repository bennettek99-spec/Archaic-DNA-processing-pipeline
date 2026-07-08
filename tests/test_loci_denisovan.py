"""
Unit test for archaic.loci.denisovan_informative — checks the selection criteria
hold for every SNP it returns, on a real on-disk synthetic panel. Run: pytest -q
"""
import os, sys, shutil, tempfile
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from archaic.synthetic import write_synthetic_panel
from archaic.panel import Panel
from archaic import loci as loci_mod

REFS = dict(Denisova=dict(ids=["Denisova"]), Altai=dict(ids=["Altai"]),
            Vindija=dict(ids=["Vindija"]), Mbuti=dict(pops=["Mbuti"]),
            Yoruba=dict(pops=["Yoruba"]))          # Yoruba absent in synthetic -> empty


def test_denisovan_informative_criteria():
    tmp = tempfile.mkdtemp(prefix="archaic_deni_")
    try:
        prefix, _ = write_synthetic_panel(tmp, n_snp=8000, seed=7)
        panel = Panel(prefix, autosomes_only=True)
        di = loci_mod.denisovan_informative(panel, panel.snp_rows, REFS,
                                            den_thresh=0.9, afr_thresh=0.1, nea_max=0.5)
        # consistent shapes
        k = len(di["rows"])
        for key in ("den_is_a1", "p_den", "p_afr", "p_nea"):
            assert len(di[key]) == k
        assert np.all(np.isin(di["rows"], panel.snp_rows))
        if k:
            a1 = di["den_is_a1"]
            den = np.where(a1, di["p_den"], 1 - di["p_den"])
            afr = np.where(a1, di["p_afr"], 1 - di["p_afr"])
            nea = np.where(a1, di["p_nea"], 1 - di["p_nea"])
            assert np.all(den >= 0.9)          # Denisova ~fixed for the archaic allele
            assert np.all(afr <= 0.1)          # ~absent in Africans
            assert np.all(nea <= 0.5)          # not carried by Neanderthals (Deni-specific)
    finally:
        del panel; shutil.rmtree(tmp, ignore_errors=True)
