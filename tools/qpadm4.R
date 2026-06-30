#!/usr/bin/env Rscript
# 4-source qpAdm (Anatolia_N + Yamnaya + WHG + Iran_N) for the Italian targets,
# to obtain a fitting ancestry model (the 3-way is rejected at full SNP density).
# Reuses the cached f2 from tools/admixtools_concordance.R.
#   Rscript tools/qpadm4.R results/f2_admixtools results/admixtools_qpadm4.csv
args   <- commandArgs(trailingOnly = TRUE)
f2dir  <- ifelse(length(args) >= 1, args[1], "results/f2_admixtools")
out    <- ifelse(length(args) >= 2, args[2], "results/admixtools_qpadm4.csv")
suppressMessages(library(admixtools))
f2 <- f2_from_precomp(f2dir, verbose = FALSE)

left  <- c("Anatolia_N", "Yamnaya", "WHG", "Iran_N")
right <- c("Mbuti", "Han", "Papuan", "Karitiana", "Natufian", "Ust_Ishim", "MA1")
targets <- c("Etruscan", "Latin", "ImperialRoman", "ItalyBA")

rows <- list()
for (tg in targets) {
  r <- tryCatch(qpadm(f2, left = left, right = right, target = tg, verbose = FALSE),
                error = function(e) NULL)
  if (is.null(r)) next
  w <- r$weights
  g <- function(p) { v <- w$weight[w$left == p]; if (length(v)) v else NA }
  rows[[length(rows)+1]] <- data.frame(
    target = tg, Anatolia_N = g("Anatolia_N"), Yamnaya = g("Yamnaya"),
    WHG = g("WHG"), Iran_N = g("Iran_N"), p = r$popdrop$p[1], stringsAsFactors = FALSE)
  cat(sprintf("%-14s AnatN=%.0f%% Steppe=%.0f%% WHG=%.0f%% Iran=%.0f%%  p=%.3f\n",
              tg, g("Anatolia_N")*100, g("Yamnaya")*100, g("WHG")*100, g("Iran_N")*100, r$popdrop$p[1]))
}
write.csv(do.call(rbind, rows), out, row.names = FALSE)
cat("Wrote", out, "\n")
