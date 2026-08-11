from __future__ import annotations

import gzip
from argparse import Namespace

import pandas as pd
import pytest

from archaic_admixture_dating.genetic_map import apply_genetic_map, load_genetic_map
from archaic_admixture_dating.cli import _input_manifest
from archaic_admixture_dating.tract_import import import_tracts


def _write_map(path, chromosome="22"):
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)\n")
        handle.write(f"chr{chromosome}\t100\t10\t0.0\n")
        handle.write(f"chr{chromosome}\t200\t10\t1.0\n")
        handle.write(f"chr{chromosome}\t300\t10\t2.0\n")


def test_map_interpolation_and_outside_range_fail_closed(tmp_path):
    map_path = tmp_path / "genetic_map_GRCh37_chr22.txt.gz"
    _write_map(map_path)
    tracts = pd.DataFrame(
        {
            "chromosome": ["22", "22"],
            "start_bp": [150, 50],
            "end_bp": [250, 150],
        }
    )

    mapped = apply_genetic_map(tracts, tmp_path)

    assert mapped.loc[0, "start_cm"] == pytest.approx(0.5)
    assert mapped.loc[0, "end_cm"] == pytest.approx(1.5)
    assert mapped.loc[0, "length_cm"] == pytest.approx(1.0)
    assert mapped.loc[0, "genetic_map_status"] == "interpolated"
    assert pd.isna(mapped.loc[1, "length_cm"])
    assert mapped.loc[1, "genetic_map_status"] == "outside_map_range"


def test_map_loader_rejects_nonmonotonic_cm(tmp_path):
    path = tmp_path / "bad.txt.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("Position(bp)\tMap(cM)\n100\t1.0\n200\t0.5\n")
    with pytest.raises(ValueError, match="not monotonic"):
        load_genetic_map(path)


def test_skov_import_can_replace_constant_rate_with_genetic_map(tmp_path):
    _write_map(tmp_path / "genetic_map_GRCh37_chr22.txt.gz")
    source = tmp_path / "skov.tsv"
    pd.DataFrame(
        [
            {
                "name": "PAPUAN_1",
                "pop": "Papuans",
                "chrom": 22,
                "start": 150,
                "end": 250,
                "length": 100,
                "MeanProb": 0.91,
                "Shared_with_Altai": 4,
                "Shared_with_Denisova": 12,
                "Shared_with_Vindija": 5,
                "method": "HMM",
            },
            {
                "name": "PAPUAN_X",
                "pop": "Papuans",
                "chrom": "X",
                "start": 150,
                "end": 250,
                "length": 100,
                "MeanProb": 0.91,
                "Shared_with_Altai": 4,
                "Shared_with_Denisova": 12,
                "Shared_with_Vindija": 5,
                "method": "HMM",
            },
        ]
    ).to_csv(source, sep="\t", index=False)

    valid, excluded = import_tracts(
        source,
        caller="skov_hmm",
        chromosomes=["22"],
        genetic_map_directory=tmp_path,
    )

    assert excluded.empty
    assert valid["sample_id"].tolist() == ["PAPUAN_1"]
    assert valid.loc[0, "length_cm"] == pytest.approx(1.0)
    assert valid.loc[0, "genetic_map_status"] == "interpolated"
    assert "GRCh37" in valid.loc[0, "genetic_length_method"]


def test_map_manifest_is_portable_and_hashed(tmp_path):
    source = tmp_path / "tracts.tsv"
    source.write_text("placeholder\n", encoding="utf-8")
    map_path = tmp_path / "genetic_map_GRCh37_chr22.txt.gz"
    _write_map(map_path)
    args = Namespace(
        input=str(source),
        caller="generic",
        population=None,
        chromosomes=["22"],
        genetic_map_dir=str(tmp_path),
        genetic_map_pattern="genetic_map_GRCh37_chr{chromosome}.txt.gz",
        genetic_map_build="GRCh37",
    )

    manifest = _input_manifest(args)

    assert manifest["input"]["filename"] == "tracts.tsv"
    assert "path" not in manifest["input"]
    assert manifest["genetic_map"]["files"][0]["filename"] == map_path.name
    assert len(manifest["genetic_map"]["files"][0]["sha256"]) == 64
