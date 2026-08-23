"""
bcf_reader.py — a minimal, dependency-free sequential reader for BCF2 files.

The high-coverage archaic genomes that would extend this pipeline's contrast
axis — Chagyrskaya and Mezmaiskaya alongside Altai and Vindija — are published
as BCF. Reading BCF normally means htslib, through `pysam`, `cyvcf2` or the
`bcftools` binary, and none of the three builds on this project's Windows box:
`pysam` has no Windows wheel and fails at the build step, and there is no
WSL or conda to fall back on. That is a packaging accident, not a real obstacle:
BCF2 is a documented, stable binary format, and everything this pipeline needs
from it is four fields.

SCOPE, DELIBERATELY SMALL

Yielded per record: contig name, 1-based position, reference allele, alternate
alleles, and the per-sample genotype calls. INFO fields, FILTER, QUAL and every
FORMAT field other than GT are skipped without being decoded — parsed far enough
to know their length, then stepped over. A general-purpose BCF library is not
wanted here and would be a much larger thing to get right.

Random access is not supported either. The `.csi` index is not read; files are
streamed from the start. Every use in this pipeline is a full pass that
intersects an archaic genome against a fixed SNP panel, so an index would buy
nothing and would be one more binary format to implement correctly.

FORMAT NOTES THAT MATTER

  * BCF2 is BGZF-compressed, and BGZF is a valid multi-member gzip stream, so
    the standard library's `gzip` decompresses it with no special handling.
  * Values are "typed": a descriptor byte packs a count in the high nibble and a
    type code in the low nibble. A count of 15 means the real count follows as
    its own typed integer, which is the case that a naive reader gets wrong on
    long allele lists and multi-sample records.
  * A genotype value encodes the allele as `(allele + 1) << 1 | phased`, so 0 is
    the end-of-vector marker and 1 is a missing allele. Both are returned as -1.
  * Positions on disk are 0-based; VCF and this pipeline are 1-based, so 1 is
    added on the way out.

CORRECTNESS

`tests/test_bcf_reader.py` checks the typed-value decoder against hand-built
byte strings, including the count-15 escape. The stronger check is external and
lives in `scripts/import_archaic_bcf.py`: the AADR 1240K panel already contains
Altai and Vindija, so their genotypes read out of a BCF must agree with the
panel site by site. A parser that is subtly wrong will not survive that.
"""
from __future__ import annotations

import gzip
import struct
from typing import Iterator, NamedTuple

BCF_MAGIC = b"BCF\x02"

# type codes from the BCF2 specification
_T_MISSING, _T_INT8, _T_INT16, _T_INT32, _T_FLOAT, _T_CHAR = 0, 1, 2, 3, 5, 7
_INT_SIZE = {_T_INT8: 1, _T_INT16: 2, _T_INT32: 4}
_INT_FMT = {_T_INT8: "<b", _T_INT16: "<h", _T_INT32: "<i"}
# "missing" sentinels, one per integer width
_INT_MISSING = {_T_INT8: -128, _T_INT16: -32768, _T_INT32: -2147483648}


class Record(NamedTuple):
    chrom: str
    pos: int                 # 1-based
    ref: str
    alts: tuple
    genotypes: tuple         # one tuple of allele indices per sample, -1 missing


def _read_typed_descriptor(buf: memoryview, off: int):
    """Return (count, type_code, new_offset) for a typed value."""
    b = buf[off]
    off += 1
    n, t = b >> 4, b & 0x0F
    if n == 15:                      # escape: the true count is a typed integer
        n, off = _read_typed_int(buf, off)
    return n, t, off


def _read_typed_int(buf: memoryview, off: int):
    """Read a single typed integer (used for the count-15 escape)."""
    b = buf[off]
    off += 1
    t = b & 0x0F
    size = _INT_SIZE[t]
    v = struct.unpack_from(_INT_FMT[t], buf, off)[0]
    return v, off + size


def _skip_typed(buf: memoryview, off: int) -> int:
    """Step over a typed value without decoding it."""
    n, t, off = _read_typed_descriptor(buf, off)
    if t == _T_MISSING:
        return off
    if t == _T_CHAR:
        return off + n
    if t == _T_FLOAT:
        return off + 4 * n
    return off + _INT_SIZE[t] * n


def _read_typed_string(buf: memoryview, off: int):
    n, t, off = _read_typed_descriptor(buf, off)
    if t == _T_MISSING:
        return "", off
    if t != _T_CHAR:
        raise ValueError(f"expected a string, got type code {t}")
    s = bytes(buf[off:off + n]).decode("ascii", "replace")
    return s, off + n


def _parse_header(text: str):
    """Sample names, contig names in IDX order, and the FORMAT index of GT.

    BCF refers to contigs and FORMAT fields by integer index. Those indices come
    from the header's IDX attributes when present; when absent, the spec says
    they are assigned in order of appearance, counting the three dictionaries
    separately. Both cases occur in published archaic files, so both are handled.
    """
    contigs, fmt_idx, samples = {}, {}, []
    auto_contig = auto_fmt = 0
    for line in text.split("\n"):
        line = line.strip()
        if line.startswith("#CHROM"):
            parts = line.split("\t")
            if len(parts) > 9:
                samples = parts[9:]
            continue
        if not line.startswith("##"):
            continue
        if line.startswith("##contig=<") or line.startswith("##FORMAT=<"):
            body = line[line.index("<") + 1:line.rindex(">")]
            fields, depth, cur = {}, 0, ""
            for ch in body:                     # split on commas outside quotes
                if ch == '"':
                    depth ^= 1
                if ch == "," and not depth:
                    if "=" in cur:
                        k, v = cur.split("=", 1)
                        fields[k.strip()] = v.strip().strip('"')
                    cur = ""
                else:
                    cur += ch
            if "=" in cur:
                k, v = cur.split("=", 1)
                fields[k.strip()] = v.strip().strip('"')
            name = fields.get("ID")
            if name is None:
                continue
            if line.startswith("##contig=<"):
                idx = int(fields["IDX"]) if "IDX" in fields else auto_contig
                auto_contig += 1
                contigs[idx] = name
            else:
                idx = int(fields["IDX"]) if "IDX" in fields else auto_fmt
                auto_fmt += 1
                fmt_idx[name] = idx
    return samples, contigs, fmt_idx


class BCFReader:
    """Sequential BCF2 reader. Use as a context manager and iterate."""

    def __init__(self, path):
        self._fh = gzip.open(path, "rb")
        magic = self._fh.read(5)
        if not magic.startswith(BCF_MAGIC):
            self._fh.close()
            raise ValueError(f"{path}: not a BCF2 file (magic {magic!r})")
        l_text = struct.unpack("<I", self._fh.read(4))[0]
        text = self._fh.read(l_text).split(b"\x00", 1)[0].decode("utf-8",
                                                                "replace")
        self.header_text = text
        self.samples, self._contigs, self._fmt_idx = _parse_header(text)
        self.n_samples = len(self.samples)
        if "GT" not in self._fmt_idx:
            raise ValueError(f"{path}: header declares no GT FORMAT field")
        self._gt_key = self._fmt_idx["GT"]

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def close(self):
        self._fh.close()

    def __iter__(self) -> Iterator[Record]:
        read = self._fh.read
        while True:
            head = read(8)
            if len(head) < 8:
                return
            l_shared, l_indiv = struct.unpack("<II", head)
            buf = memoryview(read(l_shared + l_indiv))
            if len(buf) < l_shared + l_indiv:
                return                      # truncated tail; stop cleanly
            yield self._parse_record(buf, l_shared)

    def iter_at(self, wanted) -> Iterator[Record]:
        """Yield only records whose position is in `wanted`.

        `wanted` maps contig name -> a set of 1-based positions.

        The saving is real but bounded, and it is worth being clear about
        which half it is. Decompression still happens for every record — the
        stream is sequential and there is no index — but a record that is not
        wanted is rejected after reading eight bytes of its shared block, so
        the allele strings, the INFO block and the per-sample genotypes are
        never decoded. Against a SNP panel that selects roughly one position in
        thirty, that skips the large majority of the per-record Python work.
        """
        read = self._fh.read
        by_idx = {i: wanted.get(name, ()) for i, name in self._contigs.items()}
        while True:
            head = read(8)
            if len(head) < 8:
                return
            l_shared, l_indiv = struct.unpack("<II", head)
            buf = memoryview(read(l_shared + l_indiv))
            if len(buf) < l_shared + l_indiv:
                return
            chrom_i, pos0 = struct.unpack_from("<ii", buf, 0)
            hits = by_idx.get(chrom_i)
            if not hits or (pos0 + 1) not in hits:
                continue
            yield self._parse_record(buf, l_shared)

    def _parse_record(self, buf: memoryview, l_shared: int) -> Record:
        chrom_i, pos, _rlen, _qual, n_ai, n_fs = struct.unpack_from(
            "<iiifII", buf, 0)
        n_allele = n_ai >> 16
        n_info = n_ai & 0xFFFF
        n_fmt = n_fs >> 24
        n_sample = n_fs & 0xFFFFFF
        off = 24
        off = _skip_typed(buf, off)                     # ID
        alleles = []
        for _ in range(n_allele):
            s, off = _read_typed_string(buf, off)
            alleles.append(s)
        off = _skip_typed(buf, off)                     # FILTER
        for _ in range(n_info):                         # INFO key/value pairs
            off = _skip_typed(buf, off)
            off = _skip_typed(buf, off)

        gts = tuple(tuple() for _ in range(n_sample))
        off = l_shared
        for _ in range(n_fmt):
            key, off = _read_typed_int(buf, off)
            n, t, off = _read_typed_descriptor(buf, off)
            if key == self._gt_key and t in _INT_SIZE:
                gts = self._decode_gt(buf, off, n, t, n_sample)
            # step over this field's values for every sample
            if t == _T_MISSING:
                pass
            elif t == _T_CHAR:
                off += n * n_sample
            elif t == _T_FLOAT:
                off += 4 * n * n_sample
            else:
                off += _INT_SIZE[t] * n * n_sample
        return Record(self._contigs.get(chrom_i, str(chrom_i)), pos + 1,
                      alleles[0] if alleles else "",
                      tuple(alleles[1:]), gts)

    @staticmethod
    def _decode_gt(buf, off, n, t, n_sample):
        size, fmt = _INT_SIZE[t], _INT_FMT[t]
        miss = _INT_MISSING[t]
        out = []
        for s in range(n_sample):
            calls = []
            for a in range(n):
                v = struct.unpack_from(fmt, buf, off + (s * n + a) * size)[0]
                # 0 ends the vector, 1 is a missing allele, miss is "no value"
                calls.append(-1 if v in (0, 1, miss) else (v >> 1) - 1)
            out.append(tuple(calls))
        return tuple(out)


def read_bcf(path) -> Iterator[Record]:
    """Convenience wrapper: iterate a BCF file's records."""
    with BCFReader(path) as r:
        yield from r
