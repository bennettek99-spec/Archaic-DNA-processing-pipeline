from pathlib import Path

import numpy as np
import pandas as pd

from archaic.highest_archaic import (
    load_analysis, load_settings, parse_contamination, prepare,
    rank_denisovan, rank_neanderthal,
)
from archaic.highest_archaic_segments import canonical_id, summarize_segments


FIXTURE = Path(__file__).parent / "fixtures" / "highest_archaic_small.tsv"


def test_parse_contamination_uses_lower_bound():
    assert parse_contamination("[0.005,0.012]") == 0.005
    assert np.isnan(parse_contamination("n/a"))


def test_analysis_sets_are_distribution_derived_and_modern_is_removed(tmp_path):
    df = load_analysis(FIXTURE)
    assert "MODERN.HO" not in set(df.genetic_id)
    cfg = load_settings()
    out, thresholds = prepare(df, cfg, str(tmp_path / "missing_panel"))
    assert thresholds["elite_min_snps"] >= thresholds["high_min_snps"]
    assert out.loc[out.genetic_id == "RAW1.AG", "analysis_set"].iloc[0] == "broad"
    assert out.loc[out.genetic_id == "HC1.SG", "analysis_set"].iloc[0] in {
        "high-confidence", "elite-confidence"
    }
    assert "contamination" in out.loc[
        out.genetic_id == "CONTAM.AG", "analysis_set_fail_reasons"
    ].iloc[0]


def test_rankings_separate_raw_from_supported_lower_bound(tmp_path):
    df = load_analysis(FIXTURE)
    out, _ = prepare(df, load_settings(), str(tmp_path / "missing_panel"))
    assert rank_neanderthal(out, 1).iloc[0].genetic_id == "RAW1.AG"
    assert rank_neanderthal(out, 1, lcb=True).iloc[0].genetic_id == "HC1.SG"
    assert rank_denisovan(out, 1).iloc[0].genetic_id == "DEN1.SG"


def test_combined_percentage_is_not_fabricated(tmp_path):
    df = load_analysis(FIXTURE)
    out, _ = prepare(df, load_settings(), str(tmp_path / "missing_panel"))
    assert out.combined_archaic_pct.isna().all()
    assert out.denisovan_pct.isna().all()
    assert out.combined_status.str.contains("not_estimable").all()


def test_segment_followup_summarizes_only_supplied_calls():
    seg=pd.DataFrame({"sample":["Oase1","Oase1"],"chrom":[5,9],"span_cM":[36.5,30.4]})
    out=summarize_segments(seg).iloc[0]
    assert out.n_segments == 2
    assert np.isclose(out.total_segment_length,66.9)
    assert canonical_id("Oase1_d.AG.BY.AA") == "Oase1"
