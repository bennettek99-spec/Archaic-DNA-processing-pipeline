#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 3 — estimate archaic introgression for every retained genome (vectorised
+ parallelised).

For each Phase-2 sample:
  alpha_Nea  Neanderthal proportion  = f4(Altai,Chimp; X,Mbuti)/f4(Altai,Chimp; Vindija,Mbuti)
  D_Nea      Neanderthal affinity    = D(X, Mbuti; Altai, Chimp)     (+ jackknife Z)
  D_Den      Denisovan affinity      = D(X, Mbuti; Denisova, Chimp)  (+ jackknife Z)

Reference allele frequencies and all reference-only per-SNP constants are computed
ONCE; test individuals are processed in chunks with a fully vectorised block
jackknife (archaic.stats.batch_jackknife_ratio). Chunks are distributed across
CPU cores via ThreadPoolExecutor (numpy releases the GIL for the heavy
operations, and the memmap-backed TGENO reads are also GIL-free). Output:

  results/phase3_<panel>_estimates.csv      (resumable; skips ids already written)

Usage:
  python phase3_estimate.py --panel 1240k
  python phase3_estimate.py --panel 1240k --limit 200 --out results/_tmp.csv   # test
"""
import os, sys, time, argparse, csv
import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archaic.panel import Panel
from archaic import stats as st
from archaic.refs import PANELS
from archaic.log_utils import get_logger

log = get_logger("archaic.phase3")

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
N_BLOCKS = 50
F = np.float32


def _process_chunk(chunk, panel_prefix, starts, pMb, axN, axD, denYN, denYD, a_vmb,
                   snp_rows, transversions_only, F):
    """Process one chunk of individuals and return list of result dicts."""
    panel = Panel(panel_prefix, transversions_only=transversions_only)
    cols = np.array([c for _, c in chunk], dtype=np.int64)
    G = panel.pg.read(snp_rows, cols)
    pX = G.astype(F); pX[G < 0] = np.nan; pX *= F(0.5)
    diffXM = pX - pMb[:, None]
    denX = pX + pMb[:, None] - 2.0 * pX * pMb[:, None]

    num = axN[:, None] * diffXM
    den = np.broadcast_to(a_vmb[:, None], num.shape)
    a_t, a_se, _, a_n = st.batch_jackknife_ratio(num, den, starts)

    num = diffXM * axN[:, None]
    den = denX * denYN[:, None]
    dn_t, dn_se, dn_z, dn_n = st.batch_jackknife_ratio(num, den, starts)

    num = diffXM * axD[:, None]
    den = denX * denYD[:, None]
    dd_t, dd_se, dd_z, dd_n = st.batch_jackknife_ratio(num, den, starts)

    return [dict(
        genetic_id=gid,
        alpha_Nea=round(float(a_t[k]), 6), alpha_SE=round(float(a_se[k]), 6),
        alpha_nSNP=int(a_n[k]),
        D_Nea=round(float(dn_t[k]), 6), D_Nea_SE=round(float(dn_se[k]), 6),
        D_Nea_Z=round(float(dn_z[k]), 3), D_Nea_nSNP=int(dn_n[k]),
        D_Den=round(float(dd_t[k]), 6), D_Den_SE=round(float(dd_se[k]), 6),
        D_Den_Z=round(float(dd_z[k]), 3), D_Den_nSNP=int(dd_n[k]),
    ) for k, (gid, _) in enumerate(chunk)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", choices=list(PANELS), default="1240k")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--chunk", type=int, default=128)
    ap.add_argument("--workers", type=int, default=0,
                    help="parallel worker count (0 = os.cpu_count())")
    ap.add_argument("--out", default="")
    ap.add_argument("--meta", default="",
                    help="override Phase-2 metadata CSV (e.g. the global scope one)")
    ap.add_argument("--transversions-only", action="store_true",
                    help="use only transversion SNPs for a damage-robust comparison run")
    args = ap.parse_args()
    cfg = PANELS[args.panel]
    n_workers = args.workers or os.cpu_count() or 4

    meta_path = args.meta or os.path.join(RESULTS, f"phase2_{args.panel}_metadata.csv")
    meta = pd.read_csv(meta_path)
    ids = meta["genetic_id"].tolist()
    if args.limit:
        ids = ids[:args.limit]
    log.info(f"Phase 3 — panel={args.panel}  samples={len(ids):,}  "
             f"chunk={args.chunk}  workers={n_workers}")

    panel = Panel(cfg["prefix"], transversions_only=args.transversions_only)
    starts = st.block_starts(panel.n_snp, N_BLOCKS)

    log.info("Reference allele frequencies (computed once)...")
    rf, ri = panel.frequencies({k: cfg["refs"][k] for k in
                                ["Altai", "Vindija", "Denisova", "Chimp", "Mbuti"]})
    for k in ["Altai", "Vindija", "Denisova", "Chimp", "Mbuti"]:
        log.info(f"  {k:9s} SNPs={ri[k]['n_snp_covered']:,}")

    pAlt = rf["Altai"].astype(F); pChi = rf["Chimp"].astype(F)
    pVin = rf["Vindija"].astype(F); pDen = rf["Denisova"].astype(F)
    pMb = rf["Mbuti"].astype(F)
    axN = pAlt - pChi
    axD = pDen - pChi
    denYN = pAlt + pChi - 2.0 * pAlt * pChi
    denYD = pDen + pChi - 2.0 * pDen * pChi
    a_vmb = axN * (pVin - pMb)

    col_of = panel._id_to_col
    use = [(i, col_of[i]) for i in ids if i in col_of]
    missing = [i for i in ids if i not in col_of]
    if missing:
        log.warning(f"{len(missing)} ids not in .ind (skipped)")

    default_name = (f"phase3_{args.panel}_transversions_estimates.csv"
                    if args.transversions_only else f"phase3_{args.panel}_estimates.csv")
    out_path = args.out or os.path.join(RESULTS, default_name)
    fields = ["genetic_id", "alpha_Nea", "alpha_SE", "alpha_nSNP",
              "D_Nea", "D_Nea_SE", "D_Nea_Z", "D_Nea_nSNP",
              "D_Den", "D_Den_SE", "D_Den_Z", "D_Den_nSNP"]
    done_ids = set()
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        try:
            done_ids = set(pd.read_csv(out_path, usecols=["genetic_id"])["genetic_id"])
        except Exception:
            done_ids = set()
    if done_ids:
        before = len(use)
        use = [(i, c) for (i, c) in use if i not in done_ids]
        log.info(f"resume: {len(done_ids):,} done, {len(use):,}/{before:,} remaining")
    if not use:
        log.info("Nothing to do.")
        return

    fh = open(out_path, "a", newline="")
    writer = csv.DictWriter(fh, fieldnames=fields)
    if not done_ids:
        writer.writeheader()

    # precompute shared constants for worker threads
    pMb_arr = pMb[:, None]
    a_vmb_arr = a_vmb[:, None]

    t0 = time.time()
    done = 0
    chunks = [use[c0:c0 + args.chunk] for c0 in range(0, len(use), args.chunk)]
    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        futs = {}
        for chunk in chunks:
            fut = pool.submit(_process_chunk, chunk, cfg["prefix"], starts,
                              pMb_arr, axN[:, None], axD[:, None],
                              denYN[:, None], denYD[:, None], a_vmb_arr,
                              panel.snp_rows, args.transversions_only, F)
            futs[fut] = len(chunk)
        for fut in as_completed(futs):
            results = fut.result()
            for row in results:
                writer.writerow(row)
            done += len(results)
            fh.flush()
            rate = done / (time.time() - t0)
            remaining = sum(s for f, s in futs.items() if not f.done())
            if rate > 0 and remaining > 0:
                log.info(f"{done:,}/{len(use):,}  ({rate:.0f}/s, "
                         f"ETA {remaining/rate/60:.1f} min)")
    fh.close()
    log.info(f"Done in {(time.time()-t0)/60:.1f} min -> {out_path}")


if __name__ == "__main__":
    main()
