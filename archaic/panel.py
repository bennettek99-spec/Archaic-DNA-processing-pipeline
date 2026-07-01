"""
panel.py — load the AADR panel and compute per-population allele frequencies
on the autosomes, pulling only the individuals we ask for.

A "population" here is just a named set of individuals, selected by exact .ind
individual-ID (for single reference genomes like AltaiNeanderthal.DG) and/or by
exact group label (for present-day HGDP pops like 'Han'). Frequencies are the
mean dosage / 2 over non-missing genotypes; a SNP with no data in a population is
NaN there and is dropped from any statistic that needs that population.
"""
from __future__ import annotations
import numpy as np
from . import lib_eigenstrat as le

DEFAULT_PREFIX = r"C:\Users\benne\aadr_v66\v66.p1_HO"
AUTOSOMES = {str(c) for c in range(1, 23)}


class Panel:
    def __init__(self, prefix: str = DEFAULT_PREFIX, autosomes_only: bool = True):
        self.prefix = prefix
        self.snp = le.read_snp(prefix + ".snp")
        self.ind = le.read_ind(prefix + ".ind")
        self.pg = le.PackedGeno(prefix + ".geno", len(self.ind), len(self.snp))

        # autosomal SNP rows, kept in genomic (.geno) order
        if autosomes_only:
            keep = self.snp["chrom"].isin(AUTOSOMES).to_numpy()
        else:
            keep = np.ones(len(self.snp), dtype=bool)
        self.snp_rows = np.where(keep)[0].astype(np.int64)
        self.n_snp = len(self.snp_rows)

        # fast lookups
        self._id_to_col = {v: i for i, v in enumerate(self.ind["id"].values)}
        self._pop_to_cols = {}
        for i, lab in enumerate(self.ind["pop"].values):
            self._pop_to_cols.setdefault(lab, []).append(i)

    # ------------------------------------------------------------------ select
    def cols_for(self, ids=None, pops=None) -> np.ndarray:
        """Column indices for a set of exact individual-IDs and/or pop labels."""
        cols = []
        for sid in (ids or []):
            c = self._id_to_col.get(sid)
            if c is not None:
                cols.append(c)
        for pop in (pops or []):
            cols.extend(self._pop_to_cols.get(pop, []))
        return np.array(sorted(set(cols)), dtype=np.int64)

    # ----------------------------------------------------------- frequencies
    def frequencies(self, popspec: dict):
        """popspec: name -> dict(ids=[...], pops=[...]).

        Returns (freq, info):
          freq[name] = float64 array (n_snp,) of allele-1 frequency, NaN if the
                       population has no genotype at that SNP.
          info[name] = dict(n_ind, n_snp_covered)
        A single consistent allele coding per SNP is used throughout, which is all
        the f4/D statistics require (they are polarisation-invariant).
        """
        # gather every column we need, read once
        all_cols = []
        cols_by_name = {}
        for name, sel in popspec.items():
            cc = self.cols_for(sel.get("ids"), sel.get("pops"))
            cols_by_name[name] = cc
            all_cols.extend(cc.tolist())
        all_cols = np.array(sorted(set(all_cols)), dtype=np.int64)
        if len(all_cols) == 0:
            raise ValueError("No individuals selected for any population.")
        colpos = {c: i for i, c in enumerate(all_cols)}

        G = self.pg.read(self.snp_rows, all_cols)        # int8 (n_snp, n_ind)
        Gf = G.astype(np.float32)
        Gf[G < 0] = np.nan

        freq, info = {}, {}
        for name, cc in cols_by_name.items():
            if len(cc) == 0:
                freq[name] = np.full(self.n_snp, np.nan)
                info[name] = dict(n_ind=0, n_snp_covered=0)
                continue
            idx = [colpos[c] for c in cc]
            with np.errstate(invalid="ignore"):
                p = np.nanmean(Gf[:, idx], axis=1).astype(np.float64) / 2.0
            freq[name] = p
            info[name] = dict(n_ind=len(cc),
                              n_snp_covered=int(np.isfinite(p).sum()))
        return freq, info
