from pathlib import Path

import numpy as np
import pandas as pd

from archaic.denisovan_genome import (
    _marker_sharing,
    _portable_path,
    _transversion_sensitivity,
    genotype_pair_metrics,
)


class _DummyPanel:
    n_snp = 5


def test_public_provenance_paths_do_not_leak_external_directories(tmp_path):
    external = tmp_path.parent / "private" / "v66.p1_1240K.geno"
    assert _portable_path(external) == "v66.p1_1240K.geno"
    local = Path.cwd() / "results" / "denisovan_genome"
    assert _portable_path(local) == "results/denisovan_genome"


def test_genotype_pair_metrics_handles_missing_and_dosage_distance():
    a = np.array([0, 1, 2, -1, 0])
    b = np.array([0, 2, 0, 1, -1])
    result = genotype_pair_metrics(a, b)
    assert result["joint_calls"] == 3
    assert result["exact_concordance"] == 1 / 3
    assert result["mean_allele_distance"] == 0.5
    assert result["opposite_homozygote_rate"] == 1 / 3


def test_genotype_pair_metrics_empty_overlap_is_explicit():
    result = genotype_pair_metrics([-1, -1], [0, 2])
    assert result["joint_calls"] == 0
    assert np.isnan(result["exact_concordance"])


def test_marker_sharing_orients_the_denisovan_allele():
    freq = {
        "Denisova_replicate": np.array([1.0, 0.0, 1.0, 0.0, 1.0]),
        "Altai_Neanderthal": np.array([0.0, 1.0, 0.4, 0.6, 0.0]),
        "Vindija_Neanderthal": np.array([0.0, 1.0, 0.4, 0.6, 0.0]),
        "Mbuti": np.array([0.0, 1.0, 0.0, 1.0, 0.4]),
        "Yoruba": np.array([0.0, 1.0, 0.0, 1.0, 0.4]),
        "Target": np.array([1.0, 0.0, 0.5, 0.0, 1.0]),
    }
    table = _marker_sharing(freq, _DummyPanel(), "all_snps")
    target = table[table["sample"] == "Target"].iloc[0]
    # The last site fails the African <=10% criterion, leaving four markers.
    assert target["reference_defined_markers"] == 4
    assert target["n_callable_markers"] == 4
    assert target["mean_denisovan_marker_allele"] == 0.875
    assert not bool(target["is_admixture_percentage"])


def test_transversion_sensitivity_tracks_direction_and_retention():
    full = pd.DataFrame(
        [{"test_id": "a", "estimate": 0.2, "se": 0.02, "z": 10.0, "n_snp": 100}]
    )
    tv = pd.DataFrame(
        [{"test_id": "a", "estimate": 0.18, "se": 0.04, "z": 4.5, "n_snp": 25}]
    )
    out = _transversion_sensitivity(full, tv).iloc[0]
    assert bool(out["same_direction"])
    assert out["tv_site_retention"] == 0.25
    assert np.isclose(out["estimate_delta"], -0.02)
