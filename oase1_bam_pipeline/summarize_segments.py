#!/usr/bin/env python
"""
Summarise an hmmix `decode` segment table for Oase1: how much of the genome is
archaic, how long the segments are (cM), how many exceed Fu et al. 2015's 50 cM
recent-ancestry criterion, which archaic source each segment matches, and the
implied number of generations to the Neanderthal ancestor.

hmmix decode (with -extrainfo -admixpop) writes a whitespace table whose columns
vary by version but always include chrom/start/end/length/state and, when a
genetic map or admixpop is supplied, a genetic length and per-archaic SNP-match
columns. This parser is tolerant: it locates columns by name.

Usage: python summarize_segments.py work/segments.Oase1.txt
Segment length -> generations: mean introgressed length L(cM) ~= 100/g.
"""
import sys
import numpy as np
import pandas as pd

CM_PER_MB = 1.3          # autosomal average, used only if no genetic map column

path = sys.argv[1] if len(sys.argv) > 1 else "work/segments.Oase1.txt"
df = pd.read_csv(path, sep=r"\s+", engine="python")
cols = {c.lower(): c for c in df.columns}


def col(*names):
    for n in names:
        if n in cols:
            return cols[n]
    return None


state_c = col("state")
arch = df[df[state_c].astype(str).str.lower().str.contains("archaic")] if state_c else df

# segment length in cM: prefer an explicit genetic-length column, else bp*rate
gcol = col("genetic_length", "length_cm", "cm")
if gcol:
    length_cm = arch[gcol].astype(float).values
else:
    lc = col("length")
    scol, ecol = col("start"), col("end")
    if lc:
        length_bp = arch[lc].astype(float).values
    else:
        length_bp = (arch[ecol] - arch[scol]).astype(float).values
    length_cm = length_bp / 1e6 * CM_PER_MB

length_cm = length_cm[np.isfinite(length_cm)]
n = len(length_cm)
tot = length_cm.sum()
longest = length_cm.max() if n else np.nan
top3 = np.sort(length_cm)[::-1][:3].mean() if n else np.nan

# archaic source: tally best-matching archaic if per-archaic SNP columns exist
src = []
for a in ["altai", "vindija", "chagyrskaya", "denisova"]:
    c = col(a, f"{a}_snps", f"shared_with_{a}")
    if c:
        src.append((a, arch[c].astype(float).sum()))

print("Oase1 hmmix segment summary")
print("=" * 40)
print(f"archaic segments        : {n}")
print(f"total archaic length    : {tot:.1f} cM")
print(f"longest segment         : {longest:.1f} cM")
print(f"mean of top-3 segments  : {top3:.1f} cM")
print(f"segments > 50 cM        : {(length_cm > 50).sum()}   (Fu et al. recent-ancestry criterion)")
print(f"segments > 30 cM        : {(length_cm > 30).sum()}")
if n:
    print(f"\nGenerations to Neanderthal ancestor (L ~= 100/g):")
    print(f"  from longest  : g ~= {100/longest:.1f}")
    print(f"  from top-3    : g ~= {100/top3:.1f}")
    print(f"  (Fu et al. 2015 estimate: 4-6 generations)")
if src:
    tot_src = sum(v for _, v in src) or 1
    print("\nArchaic-source composition of segments (SNP matches):")
    for a, v in src:
        print(f"  {a.capitalize():12s} {v:8.0f}  ({100*v/tot_src:5.1f}%)")
