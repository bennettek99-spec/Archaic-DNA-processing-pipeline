import os
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from archaic import cohort_rules as cr
from archaic import validation as val
from archaic.manifest import infer_aadr_release, manifest_record, write_frozen_manifest
from archaic.snp_filters import transversion_mask


def test_manifest_record_and_hash_are_stable():
    row = {
        "genetic_id": "I1.AG",
        "master_id": "S1",
        "group_id": "Italy_Etruscan",
        "locality": "Tarquinia",
        "country": "Italy",
        "lat": 42.25,
        "lon": 11.75,
        "date_bp": 2600,
        "date_sd": 50,
        "full_date": "2600 +/- 50 BP",
        "snps_1240k": 250000,
        "coverage": 0.11,
        "damage": 0.04,
        "angsd": "[0.001,0.006]",
        "hapconx": "",
        "mol_sex": "M",
        "assessment": "PASS",
    }
    rec = manifest_record(
        row,
        panel="1240k",
        panel_prefix=r"C:\aadr\v66.p1_1240K",
        snps_col="snps_1240k",
        status="retained",
    )
    assert rec["aadr_release"] == "v66.p1"
    assert rec["snp_count"] == 250000
    with tempfile.TemporaryDirectory(dir=os.getcwd()) as td:
        out = os.path.join(td, "manifest.csv")
        digest1 = write_frozen_manifest([rec], out)
        digest2 = write_frozen_manifest([rec], out)
        assert digest1 == digest2
        assert os.path.exists(out + ".sha256")


def test_infer_aadr_release_unknown_fallback():
    assert infer_aadr_release("v54.1.p1_HO") == "v54"
    assert infer_aadr_release("custom_panel") == "unknown"


def test_cohort_rules_and_duplicate_pruning():
    df = pd.DataFrame(
        [
            dict(genetic_id="I1.AG", group_id="Italy_Etruscan", country="Italy", date_bp=2600, alpha_nSNP=200000),
            dict(genetic_id="I1.SG", group_id="Italy_Etruscan", country="Italy", date_bp=2600, alpha_nSNP=300000),
            dict(genetic_id="L1.AG", group_id="Latini_IA", country="Italy", date_bp=2700, alpha_nSNP=150000),
            dict(genetic_id="B1.AG", group_id="Italy_EBA", country="Italy", date_bp=4100, alpha_nSNP=120000),
        ]
    )
    out = cr.add_population_test_keep(cr.apply_cohort_rules(df))
    assert out.loc[out.genetic_id == "I1.AG", "archaeological_cohort"].iloc[0] == "Etruscan_context"
    assert out.loc[out.genetic_id == "L1.AG", "archaeological_cohort"].iloc[0] == "Latin_context"
    assert out.loc[out.genetic_id == "B1.AG", "archaeological_cohort"].iloc[0] == "Preceding_Bronze_Age_Italy"
    assert not bool(out.loc[out.genetic_id == "I1.AG", "population_test_keep"].iloc[0])
    assert bool(out.loc[out.genetic_id == "I1.SG", "population_test_keep"].iloc[0])
    assert cr.duplicate_root("Oase1.AG.BY.AA") == "Oase1"


def test_validation_metrics_and_threshold_sensitivity():
    df = pd.DataFrame(
        [
            dict(name="A", category="x", my_Nea=2.1, my_SE=0.2, pub_Nea=2.0, pub_lo=1.8, pub_hi=2.2, source="s"),
            dict(name="B", category="x", my_Nea=3.0, my_SE=0.3, pub_Nea=2.5, pub_lo=2.2, pub_hi=2.8, source="s"),
            dict(name="C", category="x", my_Nea=0.1, my_SE=0.2, pub_Nea=0.0, pub_lo=0.0, pub_hi=0.3, source="s"),
        ]
    )
    summary = val.validation_summary(df, label="toy")
    assert summary["n"] == 3
    assert np.isclose(summary["mae_pct_points"], (0.1 + 0.5 + 0.1) / 3)
    assert summary["bias_pct_points"] > 0
    ba = val.bland_altman_table(df)
    assert set(["method_mean", "method_diff"]).issubset(ba.columns)
    sens = val.threshold_sensitivity(df, thresholds=(2.0, 3.0))
    assert sens.loc[sens.threshold_pct == 2.0, "agreement"].iloc[0] == 1.0


def test_transversion_mask():
    snp = pd.DataFrame({"a1": ["A", "A", "C", "G"], "a2": ["G", "C", "T", "T"]})
    assert transversion_mask(snp) == [False, True, False, True]
