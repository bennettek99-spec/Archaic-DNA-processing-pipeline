"""Tests for the dependency-free BCF2 reader.

The typed-value decoder is the part most likely to be subtly wrong, and wrong in
a way that still produces plausible-looking records, so it is exercised against
hand-built byte strings rather than only against a real file. The count-15
escape in particular is the case a naive reader gets wrong: it appears only on
long allele lists and wide multi-sample records, so a file that happens not to
contain one will pass a smoke test while the parser is broken.

The end-to-end check against real data lives in scripts/import_archaic_bcf.py,
where Altai and Vindija read out of a published BCF are compared site by site
against the same two genomes in the AADR panel.
"""
import gzip
import os
import struct
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from archaic import bcf_reader as br


# ------------------------------------------------------------ typed values ---
def test_descriptor_packs_count_and_type():
    buf = memoryview(bytes([(3 << 4) | br._T_INT8]))
    n, t, off = br._read_typed_descriptor(buf, 0)
    assert (n, t, off) == (3, br._T_INT8, 1)


def test_descriptor_count_15_escape_reads_a_typed_integer():
    """A count of 15 means the real count follows as its own typed integer."""
    buf = memoryview(bytes([(15 << 4) | br._T_INT8,
                            (1 << 4) | br._T_INT16]) + struct.pack("<h", 300))
    n, t, off = br._read_typed_descriptor(buf, 0)
    assert n == 300 and t == br._T_INT8
    assert off == 4


def test_skip_typed_steps_over_each_width_correctly():
    for t, width in ((br._T_INT8, 1), (br._T_INT16, 2), (br._T_INT32, 4),
                     (br._T_FLOAT, 4), (br._T_CHAR, 1)):
        buf = memoryview(bytes([(4 << 4) | t]) + b"\x00" * (4 * width))
        assert br._skip_typed(buf, 0) == 1 + 4 * width
    # a missing value occupies only its descriptor byte
    assert br._skip_typed(memoryview(bytes([br._T_MISSING])), 0) == 1


def test_read_typed_string_and_missing():
    buf = memoryview(bytes([(3 << 4) | br._T_CHAR]) + b"chr")
    s, off = br._read_typed_string(buf, 0)
    assert s == "chr" and off == 4
    s, off = br._read_typed_string(memoryview(bytes([br._T_MISSING])), 0)
    assert s == "" and off == 1


# ---------------------------------------------------------------- genotypes --
def test_decode_gt_maps_encoding_to_allele_indices():
    """A genotype value is (allele + 1) << 1 | phased."""
    vals = [(0 + 1) << 1, (1 + 1) << 1 | 1]        # 0/1, second allele phased
    buf = memoryview(bytes(bytearray([v & 0xFF for v in vals])))
    gts = br.BCFReader._decode_gt(buf, 0, 2, br._T_INT8, 1)
    assert gts == ((0, 1),)


def test_decode_gt_treats_end_of_vector_and_missing_as_no_call():
    # 0 ends the vector, 1 is a missing allele, -128 is the int8 "no value"
    buf = memoryview(bytes(bytearray([0, 1, 0x80, (1 + 1) << 1])))
    gts = br.BCFReader._decode_gt(buf, 0, 2, br._T_INT8, 2)
    assert gts == ((-1, -1), (-1, 1))


def test_decode_gt_separates_samples():
    vals = [2, 2, 4, 4]                             # sample0 0/0, sample1 1/1
    buf = memoryview(bytes(bytearray(vals)))
    assert br.BCFReader._decode_gt(buf, 0, 2, br._T_INT8, 2) == ((0, 0), (1, 1))


# ------------------------------------------------------------------ header ---
def test_header_uses_explicit_idx_when_present():
    text = ('##fileformat=VCFv4.2\n'
            '##contig=<ID=22,length=51304566,IDX=7>\n'
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="d",IDX=3>\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="g",IDX=1>\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n')
    samples, contigs, fmt = br._parse_header(text)
    assert samples == ["A", "B"]
    assert contigs == {7: "22"}
    assert fmt["GT"] == 1 and fmt["DP"] == 3


def test_header_falls_back_to_order_of_appearance_without_idx():
    """The spec assigns indices per dictionary, counted separately."""
    text = ('##contig=<ID=1,length=10>\n'
            '##contig=<ID=2,length=10>\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n')
    samples, contigs, fmt = br._parse_header(text)
    assert contigs == {0: "1", 1: "2"}
    assert fmt["GT"] == 0 and samples == ["S"]


def test_header_tolerates_commas_inside_quoted_descriptions():
    text = ('##FORMAT=<ID=GT,Number=1,Type=String,'
            'Description="Genotype, comma inside",IDX=0>\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n')
    _, _, fmt = br._parse_header(text)
    assert fmt["GT"] == 0


def test_rejects_a_file_that_is_not_bcf(tmp_path):
    p = tmp_path / "not.bcf"
    with gzip.open(p, "wb") as fh:
        fh.write(b"VCFv4.2 plain text, not a BCF at all")
    with pytest.raises(ValueError, match="not a BCF2 file"):
        br.BCFReader(p)


# --------------------------------------------------- real file, if present ---
_REAL = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                     "data", "archaic_hg19", "highcov_ind_22.bcf")


@pytest.mark.skipif(not os.path.exists(_REAL),
                    reason="archaic BCF not downloaded in this checkout")
def test_reads_the_published_archaic_bcf():
    with br.BCFReader(_REAL) as r:
        assert r.samples == ["AltaiNeandertal", "Vindija33.19", "Denisova",
                             "Chagyrskaya-Phalanx"]
        n = 0
        for rec in r:
            assert rec.chrom == "22"
            assert rec.pos > 0
            assert len(rec.genotypes) == 4
            assert all(len(g) == 2 for g in rec.genotypes)
            n += 1
            if n >= 500:
                break
        assert n == 500


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        try:
            fn()
            print(f"ok  {fn.__name__}")
        except TypeError:
            print(f"skip {fn.__name__} (needs a fixture)")
    print(f"\n{len(fns)} tests")
