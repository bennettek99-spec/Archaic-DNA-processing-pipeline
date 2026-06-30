"""
archaic — genome-wide detection of unexpected archaic introgression in ancient
Eurasian genomes (AADR).

Package layout:
  lib_eigenstrat.py  — memory-mapped reader for AADR EIGENSTRAT/TGENO files.
  stats.py           — D-statistics, f4-ratios, block-jackknife standard errors.
  panel.py           — load the panel, select samples, compute per-population
                       allele frequencies on the autosomes.

See README.md for the method, citations, and validation results.
"""
__version__ = "0.1.0"
