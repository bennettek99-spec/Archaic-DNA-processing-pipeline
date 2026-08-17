#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fetch the four high-coverage archaic genomes as hg19 BCFs.

Zenodo record 7246376, "Archaic variants from Altai, Vindija, Chagyrskaya and
Denisova (hg19)": one BCF per chromosome carrying all four samples, ~1.31 GB for
the 22 autosomes. This is the small distribution. The primary releases at
ftp.eva.mpg.de are all-sites VCFs — Chagyrskaya alone is 43 GB of autosomes at
about 14 bytes per base — and fetching only the panel's sites out of them by
range request does not help, because at one target site per ~2.4 kb almost every
compressed block is touched anyway.

hg19 matters: the AADR 1240K panel is on hg19, so positions join directly with
no liftover. Zenodo record 13368126 is the same data on hg38 and is the wrong
one for this pipeline.

Why this exists at all: the download is deliberately not committed, so the
`.gitignore` entry that excludes it points here. Files are size-checked and
skipped when already present, so re-running costs nothing and an interrupted run
resumes.

Attribution: the Chagyrskaya genome is Mafessoni et al. 2020 (PNAS 117:15132);
Altai is Prufer et al. 2014 (Nature 505:43); Vindija 33.19 is Prufer et al. 2017
(Science 358:655); Denisova is Meyer et al. 2012 (Science 338:222). The EVA
copies carry a Ft. Lauderdale request reserving first genome-wide analysis to
the data producers; for all four that first analysis is the published paper
above, so the reservation is discharged and ordinary citation applies.

Run: PYTHONIOENCODING=utf-8 python scripts/fetch_archaic_bcf.py
"""
import argparse
import json
import os
import sys
import time
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEST = os.path.join(ROOT, "data", "archaic_hg19")
RECORD = 7246376
SAMPLES = ("AltaiNeandertal", "Vindija33.19", "Denisova", "Chagyrskaya-Phalanx")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", default="1-22",
                    help="e.g. '1-22', '22', '1,2,22'")
    ap.add_argument("--dest", default=DEST)
    args = ap.parse_args()

    want = []
    for part in args.chroms.split(","):
        if "-" in part:
            a, b = part.split("-")
            want.extend(range(int(a), int(b) + 1))
        else:
            want.append(int(part))

    os.makedirs(args.dest, exist_ok=True)
    url = f"https://zenodo.org/api/records/{RECORD}"
    with urllib.request.urlopen(url, timeout=60) as r:
        meta = json.load(r)
    files = {f["key"]: f for f in meta["files"]}
    print(f"Zenodo {RECORD}: {meta['metadata']['title']}")
    print(f"samples expected: {', '.join(SAMPLES)}")

    total = fetched = 0
    for c in want:
        for key in (f"highcov_ind_{c}.bcf", f"highcov_ind_{c}.bcf.csi"):
            f = files.get(key)
            if f is None:
                print(f"  !! {key} not in the record")
                continue
            out = os.path.join(args.dest, key)
            total += f["size"]
            if os.path.exists(out) and os.path.getsize(out) == f["size"]:
                continue
            t0 = time.time()
            urllib.request.urlretrieve(f["links"]["self"], out)
            got = os.path.getsize(out)
            if got != f["size"]:
                os.remove(out)
                sys.exit(f"size mismatch for {key}: {got} != {f['size']}")
            fetched += got
            print(f"  {key:28s} {got/2**20:7.1f} MB  {time.time()-t0:5.1f}s",
                  flush=True)
    print(f"complete: {total/2**30:.2f} GB on disk, "
          f"{fetched/2**30:.2f} GB fetched this run")


if __name__ == "__main__":
    main()
