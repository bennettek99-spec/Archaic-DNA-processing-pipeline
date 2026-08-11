# Selection-candidate sensitivity input

`gower_2021_melanesian_ai_candidates_grch37.bed` contains the 22 top-ranking
candidate intervals in Table 2 of Gower et al. (2021), *Detecting adaptive
introgression in human evolution using convolutional neural networks*, eLife
10:e64669, <https://doi.org/10.7554/eLife.64669>.

The paper states that its empirical VCFs and gene annotations use GRCh37/hg19.
The table reports 1-based inclusive coordinates; this file converts them to
0-based, half-open BED intervals by subtracting one from each start coordinate.

These regions are used only for an exclusion sensitivity analysis. They are
model-ranked candidates, not a definitive catalogue of selected Denisovan loci;
the paper itself requires further assessment and warns that scan candidates can
include false positives.
