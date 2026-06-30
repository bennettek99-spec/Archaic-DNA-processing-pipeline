#!/usr/bin/env Rscript
# Concordance check: compute the same archaic statistics with ADMIXTOOLS 2 on the
# same AADR packedancestrymap files, so the pure-Python estimates can be validated
# against the field-standard tool.
#
# Install once:
#   install.packages("admixtools")   # or: remotes::install_github("uqrmaie1/admixtools")
#
# Run:
#   Rscript tools/admixtools_concordance.R  C:/Users/benne/aadr_v66/v66.p1_1240K  results/admixtools_results.csv
#
# Output CSV columns: stat,pops,est,se  — compared by tools/compare_admixtools.py.

args <- commandArgs(trailingOnly = TRUE)
prefix  <- ifelse(length(args) >= 1, args[1], "C:/Users/benne/aadr_v66/v66.p1_1240K")
outfile <- ifelse(length(args) >= 2, args[2], "results/admixtools_results.csv")

suppressMessages(library(admixtools))

ALT <- "AltaiNeanderthal.DG"; VIN <- "VindijaG1_final.SG"; CHI <- "Chimp.REF"
MB  <- "Mbuti"; YO <- "Yoruba"
tests <- c("French", "Han", "Sardinian", "Papuan", "Yoruba")

rows <- list()
add <- function(stat, pops, est, se) rows[[length(rows)+1]] <<-
  data.frame(stat=stat, pops=pops, est=est, se=se, stringsAsFactors=FALSE)

# Neanderthal f4-ratio alpha = f4(Altai,Chimp; X,Mbuti) / f4(Altai,Chimp; Vindija,Mbuti)
den <- f4(prefix, ALT, CHI, VIN, MB)
for (X in tests) {
  num <- f4(prefix, ALT, CHI, X, MB)
  add("alpha_Nea", X, num$est / den$est, NA)
  d <- qpdstat(prefix, X, MB, ALT, CHI)          # D(X,Mbuti; Altai,Chimp)
  add("D_Nea", X, d$est, d$se)
}
# differential East-Asian excess: D(French,Han; Altai, Yoruba)
d <- qpdstat(prefix, "French", "Han", ALT, YO)
add("D_diff_FrenchHan", "French;Han", d$est, d$se)

res <- do.call(rbind, rows)
write.csv(res, outfile, row.names = FALSE)
cat("Wrote", outfile, "\n")
