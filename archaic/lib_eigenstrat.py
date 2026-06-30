"""
lib_eigenstrat.py
=================
Minimal, transparent reader for AADR EIGENSTRAT files (.ind / .snp / .geno),
including BOTH packed .geno binary layouts that EIGENSOFT/AADR emit.

(Copied verbatim from the validated reader in ../ancient-pca/lib_eigenstrat.py so
this package is self-contained and reproducible. Do not diverge the two without
reason.)

Why we read the packed format directly in Python instead of using PLINK /
EIGENSOFT: it keeps the whole pipeline inspectable (no opaque format
conversions) and lets us pull ONLY the SNPs and individuals we actually need,
so a 3.8 GB / 27,594-individual file never has to be fully loaded.

Two packed layouts (auto-detected from the 5-byte magic):
  * b"GENO"  (packedancestrymap) -- SNP-major: one record per SNP, each record
    rlen = max(48, ceil(nind/4)) bytes; record 0 is the header.
  * b"TGENO" (transposed packed)  -- individual-major: a 48-byte header, then one
    record per INDIVIDUAL of rlen = ceil(nsnp/4) bytes, SNPs packed 2-bits each.
    AADR v66.p1 Human Origins ships in THIS layout (ideal for "few individuals,
    all SNPs", which is exactly our access pattern).

In both layouts each genotype is 2 bits, 4 per byte, the FIRST element in the
HIGH bits, value 0/1/2 = copies of the counted allele and 3 = missing.
"""

import os
import numpy as np
import pandas as pd


# ---------------------------------------------------------------- .ind --------
def read_ind(path):
    """Return a DataFrame with columns: id, sex, pop (one row per individual)."""
    ids, sexes, pops = [], [], []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 3:
                continue
            ids.append(parts[0])
            sexes.append(parts[1])
            pops.append(parts[2])
    return pd.DataFrame({"id": ids, "sex": sexes, "pop": pops})


# ---------------------------------------------------------------- .snp --------
def read_snp(path):
    """
    Read the .snp file.  Columns: name, chrom, gpos, pos, a1, a2.
    chrom is kept as a string; pos is an int.  The row index is the SNP's
    position in the .geno file (0-based), which is what the packed reader needs.
    """
    df = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        names=["name", "chrom", "gpos", "pos", "a1", "a2"],
        dtype={"name": str, "chrom": str, "gpos": float,
               "pos": np.int64, "a1": str, "a2": str},
    )
    df.index.name = "geno_row"          # 0..nsnp-1, matches packed .geno order
    return df


# ------------------------------------------------------ packed .geno reader ---
class PackedGeno:
    """
    Memory-mapped reader supporting both 'GENO' and 'TGENO' packed layouts.

        pg = PackedGeno(geno_path, nind, nsnp)
        block = pg.read(snp_rows, ind_cols)  # int8 [len(snp), len(ind)],
                                             # 0/1/2 = dosage, -1 = missing
    """

    def __init__(self, geno_path, nind, nsnp):
        self.path = geno_path
        self.nind = int(nind)
        self.nsnp = int(nsnp)
        fsize = os.path.getsize(geno_path)

        self.mm = np.memmap(geno_path, dtype=np.uint8, mode="r")
        magic = self.mm[:5].tobytes()

        # parse nind/nsnp from the header text for a cross-check
        htxt = self.mm[:48].tobytes().split(b"\x00")[0].decode("ascii", "ignore")
        try:
            h = htxt.split()
            self.header_nind, self.header_nsnp = int(h[1]), int(h[2])
        except Exception:
            self.header_nind = self.header_nsnp = None

        if magic.startswith(b"TGENO"):
            self.transposed = True
            self.rlen = (self.nsnp + 3) // 4                # ceil(nsnp/4)
            self.hlen = fsize - self.nind * self.rlen        # header bytes (==48)
            if self.hlen < 0 or fsize != self.hlen + self.nind * self.rlen:
                raise ValueError(
                    f"{geno_path}: TGENO size {fsize} inconsistent with "
                    f"nind={nind}, nsnp={nsnp}, rlen={self.rlen}.")
            self.data = self.mm[self.hlen: self.hlen + self.nind * self.rlen
                                ].reshape(self.nind, self.rlen)
        elif magic.startswith(b"GENO"):
            self.transposed = False
            self.hlen = self.rlen = max(48, (self.nind + 3) // 4)  # ceil(nind/4)
            if fsize != self.rlen * (self.nsnp + 1):
                raise ValueError(
                    f"{geno_path}: GENO size {fsize} != rlen*(nsnp+1)="
                    f"{self.rlen * (self.nsnp + 1)}.")
            self.data = self.mm[self.rlen: self.rlen + self.nsnp * self.rlen
                                ].reshape(self.nsnp, self.rlen)
        else:
            raise ValueError(f"{geno_path}: unknown magic {magic!r} "
                             f"(expected b'GENO' or b'TGENO').")

    def read(self, snp_rows, ind_cols, ind_chunk=256, snp_chunk=20000):
        """
        snp_rows : 1-D int array of .geno SNP-row indices to fetch.
        ind_cols : 1-D int array of individual indices to fetch.
        Returns int8 array (len(snp_rows), len(ind_cols)); 0/1/2, -1 = missing.
        """
        snp_rows = np.asarray(snp_rows, dtype=np.int64)
        ind_cols = np.asarray(ind_cols, dtype=np.int64)
        out = np.empty((len(snp_rows), len(ind_cols)), dtype=np.int8)

        if self.transposed:
            # data is individual-major -> extract SNP bits, looping individuals
            sbyte = snp_rows // 4
            sshift = ((3 - (snp_rows % 4)) * 2).astype(np.uint8)    # (n_snp,)
            for c in range(0, len(ind_cols), ind_chunk):
                rows = self.data[ind_cols[c:c + ind_chunk]]          # (k, rlen)
                sub = rows[:, sbyte]                                 # (k, n_snp)
                vals = ((sub >> sshift) & 3).astype(np.int8)         # per-SNP shift
                vals[vals == 3] = -1
                out[:, c:c + rows.shape[0]] = vals.T
        else:
            # data is SNP-major -> extract individual bits, looping SNPs
            ibyte = ind_cols // 4
            ishift = ((3 - (ind_cols % 4)) * 2).astype(np.uint8)     # (n_ind,)
            for s in range(0, len(snp_rows), snp_chunk):
                rows = self.data[snp_rows[s:s + snp_chunk]]          # (c, rlen)
                sub = rows[:, ibyte]                                 # (c, n_ind)
                vals = ((sub >> ishift) & 3).astype(np.int8)         # per-ind shift
                vals[vals == 3] = -1
                out[s:s + rows.shape[0]] = vals
        return out
