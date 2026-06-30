#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
export_plink.py — export a cohort subset of the AADR (TGENO) panel to PLINK
binary (.bed/.bim/.fam). PLINK is the universal input for ADMIXTURE, LEA, and
ADMIXTOOLS 2, so this bridges the pipeline to those external tools.

By default it exports the Etruscan-context analysis set (Etruscans + neighbours +
archaic references + outgroups), autosomes only, optionally SNP-thinned for
clustering tools.

    python export_plink.py --out results/plink/etruscan_set [--snp-step 10]

The .fam family-ID column carries a population label (group_id / pop), so the
downstream tools can group individuals.
"""
import os, sys, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archaic.panel import Panel
from archaic.refs import PANELS

# curated export set: clean LABEL -> selector (kind, value). First match wins, so
# specific groups (Etruscan, Latin) take precedence over broad ones.
_whg = lambda g: any(s in g for s in ("loschbour", "villabruna", "bichon", "labrana",
                                      "iberia_mesolithic", "france_mesolithic", "england_mesolithic"))
_latin = lambda g: ("latini" in g) or (("lazio_ia" in g) and ("etruscan" not in g))
DEFAULT = [
    # archaic refs + outgroups (single samples / pops) — assign first
    ("Altai", "id", "AltaiNeanderthal.DG"), ("Vindija", "id", "VindijaG1_final.SG"),
    ("Denisova", "id", "Denisova.SG"), ("Chimp", "id", "Chimp.REF"),
    ("Ust_Ishim", "id", "Ust_Ishim.DG"), ("MA1", "id", "MA1.SG"),
    ("Mbuti", "pop", "Mbuti"), ("Yoruba", "pop", "Yoruba"), ("Han", "pop", "Han"),
    ("Papuan", "pop", "Papuan"), ("Karitiana", "pop", "Karitiana"),
    ("French", "pop", "French"), ("Sardinian", "pop", "Sardinian"),
    # ancient sources / targets — specific before broad
    ("Etruscan", "grp", "etruscan"), ("Latin", "grpfn", _latin),
    ("ImperialRoman", "grp", "imperialroman"), ("RepublicanRoman", "grp", "republic"),
    ("Anatolia_N", "grp", "turkey_n"), ("Yamnaya", "grp", "yamnaya"),
    ("WHG", "grpfn", _whg), ("Iran_N", "grp", "iran_ganjdareh_n"),
    ("Natufian", "grp", "israel_natufian"), ("ItalyBA", "grpfn",
        lambda g: "italy" in g and any(b in g for b in ("_ba", "_eba", "_mba", "_lba"))),
]
MAXN = 80  # cap individuals per population
# dosage (copies of allele1) -> PLINK 2-bit code (A1=allele1):
#   dosage 2 (hom A1) -> 0b00=0 ; dosage 1 (het) -> 0b10=2 ; dosage 0 (hom A2) -> 0b11=3 ; missing -> 0b01=1
_CODE = np.array([3, 2, 0], dtype=np.uint8)  # index by dosage 0,1,2 ; missing handled separately


def write_bed(path, panel, snp_rows, cols, chunk=20000):
    n_ind = len(cols)
    nb = (n_ind + 3) // 4
    with open(path, "wb") as fh:
        fh.write(bytes([0x6c, 0x1b, 0x01]))            # magic, SNP-major
        for s in range(0, len(snp_rows), chunk):
            G = panel.pg.read(snp_rows[s:s + chunk], cols)        # (c, n_ind) int8
            P = np.where(G < 0, np.uint8(1), _CODE[np.clip(G, 0, 2)]).astype(np.uint8)
            if nb * 4 != n_ind:
                P = np.concatenate([P, np.zeros((P.shape[0], nb * 4 - n_ind), np.uint8)], 1)
            P = P.reshape(P.shape[0], nb, 4)
            packed = (P[:, :, 0] | (P[:, :, 1] << 2) | (P[:, :, 2] << 4) | (P[:, :, 3] << 6))
            fh.write(packed.astype(np.uint8).tobytes())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", default="1240k")
    ap.add_argument("--out", default="results/plink/etruscan_set")
    ap.add_argument("--snp-step", type=int, default=1, help="keep every Nth SNP (thinning)")
    args = ap.parse_args()
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    panel = Panel(PANELS[args.panel]["prefix"])
    meta_path = os.path.join("results", f"phase2_{args.panel}_metadata.csv")

    # gather columns + clean labels (first curated match wins; cap per population)
    gl = panel.ind["pop"].str.lower()
    pops = panel.ind["pop"].values
    rng = np.random.default_rng(0)
    idset, labels = [], {}
    for label, kind, val in DEFAULT:
        if kind == "pop":
            idx = np.where(pops == val)[0]
        elif kind == "id":
            c = panel._id_to_col.get(val)
            idx = np.array([c]) if c is not None else np.array([], int)
        elif kind == "grpfn":
            idx = np.where(gl.map(val).fillna(False).values)[0]
        else:  # grp substring
            idx = np.where(gl.str.contains(val.lower(), na=False).values)[0]
        idx = [int(i) for i in idx if int(i) not in labels]
        if len(idx) > MAXN:
            idx = sorted(rng.choice(idx, MAXN, replace=False).tolist())
        for i in idx:
            labels[i] = label
            idset.append(i)
    cols = np.array(sorted(idset), dtype=np.int64)
    print(f"exporting {len(cols)} individuals across {len(set(labels.values()))} populations")

    # SNP set (autosomes, optionally thinned)
    snp_rows = panel.snp_rows[::args.snp_step]
    snp = panel.snp.loc[snp_rows]
    print(f"SNPs: {len(snp_rows):,} (step {args.snp_step})")

    # .bim
    with open(args.out + ".bim", "w") as fh:
        for r, (_, s) in zip(snp_rows, snp.iterrows()):
            fh.write(f"{s['chrom']}\t{s['name']}\t{s['gpos']}\t{s['pos']}\t{s['a1']}\t{s['a2']}\n")
    # .fam (FID = population label, IID = genetic id)
    sexmap = {"M": "1", "F": "2"}
    with open(args.out + ".fam", "w") as fh:
        for i in cols:
            iid = panel.ind["id"].values[i]
            fam = labels[i].replace(" ", "_")
            sx = sexmap.get(panel.ind["sex"].values[i], "0")
            fh.write(f"{fam}\t{iid}\t0\t0\t{sx}\t-9\n")
    # .bed
    write_bed(args.out + ".bed", panel, snp_rows, cols)
    sz = os.path.getsize(args.out + ".bed") / 1e6
    print(f"Wrote {args.out}.bed/.bim/.fam  ({sz:.0f} MB .bed)")
    print("Use with:  admixtools2 (extract_f2), ADMIXTURE, or LEA (after ped2geno).")


if __name__ == "__main__":
    main()
