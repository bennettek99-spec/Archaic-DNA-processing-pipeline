from __future__ import annotations

import pandas as pd

from archaic_admixture_dating.tract_filtering import overlapping_pairs
from archaic_admixture_dating.tract_schema import STANDARD_COLUMNS, validate_tracts


def _row(**changes):
    value = {
        "sample_id": "P1",
        "population": "Papuan_fixture",
        "chromosome": "chr21",
        "start_bp": 100,
        "end_bp": 200,
        "start_cm": 1.0,
        "end_cm": 1.1,
        "posterior_denisovan": 0.9,
        "caller": "fixture",
        "callable_fraction": 0.95,
    }
    value.update(changes)
    return value


def test_schema_calculates_lengths_and_rejects_invalid_coordinates():
    frame = pd.DataFrame([_row(), _row(start_bp=300, end_bp=200), _row(chromosome="chr99")])
    valid, excluded = validate_tracts(frame)
    assert len(valid) == 1
    assert valid.loc[0, "chromosome"] == "21"
    assert valid.loc[0, "length_bp"] == 100
    assert abs(valid.loc[0, "length_cm"] - 0.1) < 1e-12
    assert len(excluded) == 2
    reasons = excluded["_exclusion_reason"].tolist()
    assert any("end_not_greater_than_start" in reason for reason in reasons)
    assert any("invalid_chromosome" in reason for reason in reasons)
    assert set(STANDARD_COLUMNS).issubset(valid.columns)


def test_overlaps_are_reported_without_merging_callers():
    tracts, _ = validate_tracts(
        pd.DataFrame(
            [
                _row(start_bp=100, end_bp=200, caller="ibdmix"),
                _row(start_bp=150, end_bp=250, caller="admixfrog"),
            ]
        )
    )
    overlaps = overlapping_pairs(tracts)
    assert len(tracts) == 2
    assert len(overlaps) == 1
    assert not bool(overlaps.loc[0, "same_caller"])
