"""Optional, honest segment follow-up for highest-archaic candidates.

This module summarizes a supplied segment table; it does not call segments from
EIGENSTRAT genotypes. The repository's existing Oase1 array/HMM output can be
used as contextual evidence. Other samples require a validated caller and, for
read-quality controls, BAM/CRAM input.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


def canonical_id(value):
    s=re.sub(r"\.BY\.AA$","",str(value)); s=s.rsplit(".",1)[0] if "." in s else s
    return re.sub(r"_d$","",s,flags=re.IGNORECASE)


def summarize_segments(frame):
    rows=[]
    for sample,g in frame.groupby("sample"):
        if "span_cM" in g:
            length=pd.to_numeric(g["span_cM"],errors="coerce")
            unit="cM"
        elif {"start_cM","end_cM"}.issubset(g.columns):
            length=pd.to_numeric(g.end_cM,errors="coerce")-pd.to_numeric(g.start_cM,errors="coerce"); unit="cM"
        elif {"start","end"}.issubset(g.columns):
            length=(pd.to_numeric(g.end,errors="coerce")-pd.to_numeric(g.start,errors="coerce"))/1e6; unit="Mb"
        else:
            raise ValueError("Segment table needs span_cM, start_cM/end_cM, or start/end columns.")
        assignment=(g["assignment"].fillna("unassigned").value_counts().idxmax()
                    if "assignment" in g and len(g) else "unassigned")
        rows.append({"segment_sample":sample,"canonical_id":canonical_id(sample),"n_segments":int(length.notna().sum()),
                     "total_segment_length":float(length.sum()),"mean_segment_length":float(length.mean()),
                     "max_segment_length":float(length.max()),"length_unit":unit,
                     "chromosome_distribution":";".join(f"{k}:{v}" for k,v in g["chrom"].astype(str).value_counts().sort_index().items()) if "chrom" in g else "",
                     "majority_assignment":assignment})
    return pd.DataFrame(rows)


def write_segment_followup(candidate_ids,segments,output):
    path=Path(segments); summaries=summarize_segments(pd.read_csv(path)) if path.exists() else pd.DataFrame()
    by_id={r.canonical_id:r for r in summaries.itertuples(index=False)} if len(summaries) else {}
    rows=[]
    for sid in candidate_ids:
        hit=by_id.get(canonical_id(sid))
        if hit:
            row=hit._asdict(); row.update(genetic_id=sid,segment_status="existing_array_proxy",source=str(path))
        else:
            row={"genetic_id":sid,"canonical_id":canonical_id(sid),"segment_status":"not_available: validated general caller/read data required","source":""}
        rows.append(row)
    out=pd.DataFrame(rows); Path(output).parent.mkdir(parents=True,exist_ok=True); out.to_csv(output,sep="\t",index=False); return out


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument("--candidates",required=True,help="newline list or comma-separated IDs")
    ap.add_argument("--segments",default=str(Path(__file__).resolve().parent.parent/"reports"/"oase1_haplotype"/"oase1_segments.csv")); ap.add_argument("--output",required=True)
    args=ap.parse_args(argv); p=Path(args.candidates)
    ids=[x.strip() for x in (p.read_text(encoding="utf-8").splitlines() if p.exists() else args.candidates.split(",")) if x.strip()]
    write_segment_followup(ids,args.segments,args.output); return 0


if __name__ == "__main__":
    raise SystemExit(main())
