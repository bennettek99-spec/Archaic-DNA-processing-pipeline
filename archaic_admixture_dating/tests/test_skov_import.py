from __future__ import annotations

import pandas as pd

from archaic_admixture_dating.tract_filtering import filter_tracts
from archaic_admixture_dating.tract_import import import_tracts


def _skov_row(**changes):
    value = {
        "name": "PAPUAN_1",
        "chrom": 1,
        "start": 1_000_000,
        "end": 1_100_000,
        "length": 100_000,
        "snps": 30,
        "pop": "Papuans",
        "country": "Papua New Guinea",
        "region": "Melanesia",
        "MeanProb": 0.91,
        "Shared_with_Altai": 4,
        "Shared_with_Denisova": 12,
        "Shared_with_Vindija": 5,
        "outgroup": "SubAfricans",
        "method": "HMM",
    }
    value.update(changes)
    return value


def test_skov_hmm_import_is_papuan_individual_level_and_map_explicit(tmp_path):
    source = tmp_path / "skov.tsv"
    pd.DataFrame(
        [
            _skov_row(),
            _skov_row(
                name="PAPUAN_2",
                Shared_with_Altai=15,
                Shared_with_Denisova=12,
                Shared_with_Vindija=5,
            ),
            _skov_row(name="NON_HMM", method="Sstar"),
            _skov_row(name="NON_PAPUAN", pop="Han", region="EastAsia"),
        ]
    ).to_csv(source, sep="\t", index=False)

    valid, excluded = import_tracts(source, caller="skov_hmm")

    assert excluded.empty
    assert valid["sample_id"].tolist() == ["PAPUAN_1", "PAPUAN_2"]
    assert valid["posterior_archaic"].tolist() == [0.91, 0.91]
    assert valid["posterior_denisovan"].isna().all()
    assert abs(valid.loc[0, "length_cm"] - 0.12) < 1e-12
    assert valid.loc[0, "source_class"] == "denisovan_affinity_strict"
    assert valid.loc[1, "source_class"] == "denisovan_affinity"
    assert "1.2e-8" in valid.loc[0, "genetic_length_method"]


def test_qc_uses_generic_archaic_posterior_without_calling_it_denisovan(tmp_path):
    source = tmp_path / "skov.tsv"
    pd.DataFrame(
        [
            _skov_row(name="PASS", MeanProb=0.91),
            _skov_row(name="FAIL", MeanProb=0.70),
        ]
    ).to_csv(source, sep="\t", index=False)
    valid, _ = import_tracts(source, caller="skov_hmm")

    retained, excluded = filter_tracts(
        valid,
        minimum_length_cm=0.02,
        minimum_confidence=0.8,
        minimum_callable_fraction=0.8,
    )

    assert retained["sample_id"].tolist() == ["PASS"]
    assert excluded["sample_id"].tolist() == ["FAIL"]
    assert excluded["_exclusion_reason"].tolist() == ["below_confidence"]
