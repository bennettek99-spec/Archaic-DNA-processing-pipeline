#!/usr/bin/env bash
# =============================================================================
# Oase1 BAM-level archaic-haplotype pipeline (hmmix).
# Reproduces Fu et al. 2015's recent-Neanderthal-ancestor segment analysis at
# read resolution. Run on Linux/macOS/WSL with the environment.yml conda env
# (or: pip install hmmix + conda install -c bioconda samtools bcftools vcftools).
#
# Usage:  bash run_oase1_hmmix.sh
# All large downloads are cached; re-running resumes. Set THREADS as needed.
# =============================================================================
set -euo pipefail
THREADS=${THREADS:-4}
WORK=work; REF=$WORK/ref; DL=$WORK/dl
mkdir -p "$WORK" "$REF" "$DL"

# ---- tool check -------------------------------------------------------------
for t in hmmix samtools bcftools curl; do
  command -v "$t" >/dev/null 2>&1 || { echo "ERROR: '$t' not found. See README / environment.yml"; exit 1; }
done

# =============================================================================
# 1. Oase1 aligned reads from ENA study PRJEB8987
#    (ENA filereport API returns the FTP paths for the run/analysis objects.)
# =============================================================================
echo "== [1/6] Resolving Oase1 files on ENA (PRJEB8987) =="
FR="$DL/ena_filereport.tsv"
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB8987&result=analysis&fields=analysis_accession,submitted_ftp,generated_ftp&format=tsv" > "$FR" || true
# fall back to read_run if the analysis query is empty
if [ "$(wc -l < "$FR")" -le 1 ]; then
  curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB8987&result=read_run&fields=run_accession,submitted_ftp,bam_ftp,fastq_ftp&format=tsv" > "$FR"
fi
echo "   filereport:"; cat "$FR"
# extract every ftp path ending in .bam, download if absent
grep -oE 'ftp\.sra\.ebi\.ac\.uk[^;[:space:]]*\.bam' "$FR" | sort -u | while read -r path; do
  fn="$DL/$(basename "$path")"
  [ -s "$fn" ] || { echo "   downloading $path"; curl -s -o "$fn" "https://$path"; }
done
BAM=$(ls "$DL"/*.bam 2>/dev/null | head -1 || true)
[ -n "$BAM" ] || { echo "No BAM found. Inspect $FR and set BAM= manually (submitted_ftp may host CRAM/FASTQ)."; exit 1; }
echo "   using BAM: $BAM"
samtools index -@ "$THREADS" "$BAM" 2>/dev/null || true

# =============================================================================
# 2. Reference genome, ancestral alleles, strict mask, pre-computed outgroup +
#    mutationrate, and archaic BCFs (from the hmmix Zenodo release).
#    NOTE: fill ZENODO_BASE with the current hmmix data-release record; the tool
#    README (github.com/LauritsSkov/Introgression-detection) links it.
# =============================================================================
echo "== [2/6] Fetching references (cached in $REF) =="
ZENODO_BASE=${ZENODO_BASE:-"https://zenodo.org/record/REPLACE_WITH_HMMIX_RECORD/files"}
fetch () { # url outpath
  [ -s "$2" ] || { echo "   get $(basename "$2")"; curl -sL -o "$2" "$1"; }
}
# hg19 reference + ancestral (Ensembl); strict mask + outgroup + mutrate + archaic (Zenodo)
fetch "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz" "$REF/hg19.fa.gz"
[ -s "$REF/hg19.fa" ] || gunzip -k "$REF/hg19.fa.gz"
fetch "$ZENODO_BASE/strickmask.bed"       "$REF/strickmask.bed"
fetch "$ZENODO_BASE/outgroup.txt"         "$REF/outgroup.txt"
fetch "$ZENODO_BASE/mutationrate.bed"     "$REF/mutationrate.bed"
fetch "$ZENODO_BASE/hg19_ancestral.fa"    "$REF/hg19_ancestral.fa"
mkdir -p "$REF/archaic"
for a in Altai Vindija Chagyrskaya Denisova; do
  fetch "$ZENODO_BASE/${a}.bcf" "$REF/archaic/${a}.bcf"
done
echo "   (if any Zenodo file 404s, set ZENODO_BASE to the current hmmix record — see README)"

# =============================================================================
# 3. BAM -> VCF/BCF on the autosomes (bcftools mpileup + call)
# =============================================================================
echo "== [3/6] Calling variants from the BAM =="
VCF="$WORK/oase1.bcf"
if [ ! -s "$VCF" ]; then
  echo "Oase1" > "$WORK/sample.txt"   # force sample name to match individuals.json
  bcftools mpileup -f "$REF/hg19.fa" -R "$REF/strickmask.bed" -q 20 -Q 20 \
      --threads "$THREADS" -a FORMAT/DP "$BAM" \
    | bcftools call -m -Ou --threads "$THREADS" \
    | bcftools reheader -s "$WORK/sample.txt" \
    | bcftools view -Ob -o "$VCF"
  bcftools index "$VCF"
fi

# =============================================================================
# 4-5. hmmix: observations -> train (haploid) -> decode with archaic annotation
# =============================================================================
echo "== [4/6] hmmix create_ingroup =="
hmmix create_ingroup -ind=individuals.json -vcf="$VCF" \
  -weights="$REF/strickmask.bed" -out="$WORK/obs" \
  -outgroup="$REF/outgroup.txt" -ancestral="$REF/hg19_ancestral.fa"

echo "== [5/6] hmmix train (haploid) + decode =="
hmmix train -obs="$WORK/obs.Oase1.txt" -weights="$REF/strickmask.bed" \
  -mutrates="$REF/mutationrate.bed" -haploid -out="$WORK/trained.Oase1.json"

hmmix decode -obs="$WORK/obs.Oase1.txt" -weights="$REF/strickmask.bed" \
  -mutrates="$REF/mutationrate.bed" -param="$WORK/trained.Oase1.json" -haploid \
  -admixpop="$REF/archaic/*.bcf" -extrainfo > "$WORK/segments.Oase1.txt"

# =============================================================================
# 6. Summarise: segment lengths (cM), >50 cM count, generations estimate
# =============================================================================
echo "== [6/6] Summarising segments =="
python summarize_segments.py "$WORK/segments.Oase1.txt" > "$WORK/oase1_segment_summary.txt"
cat "$WORK/oase1_segment_summary.txt"
echo "Done. See $WORK/segments.Oase1.txt and $WORK/oase1_segment_summary.txt"
