# Oase1 BAM-level archaic-haplotype pipeline (hmmix)

Full read-level reproduction of the Oase1 recent-Neanderthal-ancestor analysis
(Fu et al. 2015, *Nature* 524:216), using **hmmix** (Skov et al. 2018, *PLoS
Genetics* 14:e1007641) to call archaic-introgression segments directly from the
aligned reads and annotate each as Neanderthal- or Denisovan-derived.

This is the read-level counterpart to the array-resolution analysis in
`../oase1_haplotype.py`. The array analysis already recovers the qualitative
result (Oase1's Neanderthal ancestry sits in a few long blocks → a recent
ancestor); this pipeline gives full genomic resolution, precise segment
boundaries, and calibrated segment lengths for dating the admixture.

> **Why it is not run on the project's Windows box:** hmmix's ingest needs
> `samtools`/`bcftools`, and `pysam` has no Windows build; the AADR box also has
> no WSL/conda. Run this on Linux, macOS, or WSL. Everything below is a single
> `bash run_oase1_hmmix.sh` once the prerequisites are installed.

---

## What it does

1. **Fetch Oase1 reads** — ENA study **PRJEB8987** (the accession given in Fu et
   al. 2015's data-availability statement). The script queries the ENA filereport
   API and downloads the aligned BAM(s). Oase1 is low-coverage shotgun (nuclear
   coverage well under 1×), so expect a sparse call set.
2. **Fetch references** — hg19 reference genome + ancestral-allele FASTA, the
   strict-callability mask, the pre-computed 1000G-African **outgroup** and
   **mutation-rate** files from the hmmix Zenodo release, and the **archaic
   BCFs** (Altai, Vindija, Chagyrskaya Neanderthals; Denisova).
3. **BAM → VCF** — `bcftools mpileup | bcftools call` on the autosomes.
4. **hmmix** — `create_ingroup` → `train` (haploid mode, appropriate for a
   low-coverage pseudo-haploid ancient genome) → `decode -admixpop ... -extrainfo`
   to get segments tagged Neanderthal vs Denisovan.
5. **Summarise** — `summarize_segments.py` tabulates segment lengths in cM,
   counts segments >50 cM (Fu et al.'s recent-ancestry criterion), and estimates
   the number of generations to the Neanderthal ancestor from the segment-length
   distribution (mean introgressed length ≈ 100/g cM).

## Expected result (to be confirmed by the run)

Fu et al. (2015) identified **seven** clearly recent Neanderthal segments, several
**> 50 cM**, and a total of ~6–9% of the genome of Neanderthal origin, implying a
Neanderthal ancestor **4–6 generations** before Oase1 lived. Our own genome-wide
*f₄*-ratio independently puts Oase1 at **9.8% ± 2.2%** Neanderthal (damage-
restricted), and the array-resolution segment scan already shows Oase1's longest
blocks (≈ 36 cM at 1240K resolution, a lower bound) far exceeding those of
penecontemporaneous Eurasians. This pipeline should reproduce the long-segment
signature at full resolution and refine the generation estimate.

## Prerequisites

```bash
conda env create -f environment.yml && conda activate oase1-hmmix
# or:  pip install hmmix   &&   conda install -c bioconda samtools bcftools vcftools
```

## Run

```bash
bash run_oase1_hmmix.sh
```

Outputs land in `work/`:
- `obs.Oase1.txt` — per-site archaic-informative observations
- `trained.Oase1.json` — fitted HMM parameters
- `segments.Oase1.txt` — decoded segments (state, coords, archaic annotation)
- `oase1_segment_summary.txt` — lengths in cM, >50 cM count, generation estimate

## Caveats (ancient, low-coverage input)

hmmix was designed for present-day genomes. For Oase1: (i) run **haploid** mode —
low coverage yields pseudo-haploid calls; (ii) restrict to strict-mask regions;
(iii) consider damage-aware genotyping (e.g. trimming terminal bases or using a
damage-restricted BAM) to avoid C→T/G→A ancient-DNA artifacts inflating apparent
derived alleles; (iv) the low SNP density widens segment-boundary uncertainty
relative to a high-coverage genome. The `-haploid` flag and the strict mask are
already set in `run_oase1_hmmix.sh`.

## References
- Fu Q. et al. (2015) *An early modern human from Romania with a recent Neanderthal ancestor.* Nature 524:216. ENA **PRJEB8987**.
- Skov L. et al. (2018) *Detecting archaic introgression using an unadmixed outgroup.* PLoS Genet 14:e1007641. hmmix: https://github.com/LauritsSkov/Introgression-detection
