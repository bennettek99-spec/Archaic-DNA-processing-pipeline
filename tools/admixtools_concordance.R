#!/usr/bin/env Rscript
# ADMIXTOOLS 2 concordance: compute the same archaic statistics + qpAdm on the
# PLINK export (export_plink.py) so the pure-Python pipeline can be validated
# against the field-standard tool. (ADMIXTOOLS reads PLINK/packedancestrymap, not
# the AADR TGENO format, hence the PLINK export.)
#
# Run (with Rtools on PATH):
#   Rscript tools/admixtools_concordance.R results/plink/aset results
#
# Writes results/admixtools_fstats.csv and results/admixtools_qpadm.csv.

args   <- commandArgs(trailingOnly = TRUE)
pref   <- ifelse(length(args) >= 1, args[1], "results/plink/aset")
outdir <- ifelse(length(args) >= 2, args[2], "results")
suppressMessages(library(admixtools))

ALT <- "Altai"; VIN <- "Vindija"; CHI <- "Chimp"; MB <- "Mbuti"; YO <- "Yoruba"
tests <- c("French", "Han", "Sardinian", "Papuan", "Yoruba")
src   <- c("Anatolia_N", "Yamnaya", "WHG")
right <- c("Mbuti", "Han", "Papuan", "Karitiana", "Iran_N", "Natufian", "Ust_Ishim", "MA1")
targets <- c("Etruscan", "Latin", "ImperialRoman", "ItalyBA")

allpops <- unique(c(ALT, VIN, CHI, MB, YO, tests, src, right, targets))
f2dir <- file.path(outdir, "f2_admixtools")
if (!dir.exists(f2dir)) {
  cat("extract_f2 (this reads the PLINK once)...\n")
  extract_f2(pref, f2dir, pops = allpops, maxmiss = 1, overwrite = TRUE, verbose = FALSE)
}
f2 <- f2_from_precomp(f2dir, verbose = FALSE)

# ---- f-statistics ----
rows <- list(); add <- function(...) rows[[length(rows)+1]] <<- data.frame(..., stringsAsFactors = FALSE)
den <- f4(f2, ALT, CHI, VIN, MB)$est                  # "100% Neanderthal" scale
for (X in tests) {
  num <- f4(f2, ALT, CHI, X, MB)$est
  add(stat = "alpha_Nea", pops = X, est = num / den, se = NA)
  d <- qpdstat(f2, X, MB, ALT, CHI, f4mode = FALSE)    # normalised D(X,Mbuti; Altai,Chimp)
  add(stat = "D_Nea", pops = X, est = d$est, se = d$se)
}
d <- qpdstat(f2, "French", "Han", ALT, YO, f4mode = FALSE)   # East-Asian-excess contrast
add(stat = "D_diff_FrenchHan", pops = "French;Han", est = d$est, se = d$se)
write.csv(do.call(rbind, rows), file.path(outdir, "admixtools_fstats.csv"), row.names = FALSE)

# ---- qpAdm (Anatolia_N + Yamnaya + WHG) ----
qrows <- list()
for (tg in targets) {
  out <- tryCatch(qpadm(f2, left = src, right = right, target = tg, verbose = FALSE),
                  error = function(e) NULL)
  if (!is.null(out)) {
    w <- out$weights; p <- out$popdrop$p[1]
    qrows[[length(qrows)+1]] <- data.frame(
      target = tg,
      Anatolia_N = w$weight[w$left == "Anatolia_N"],
      Yamnaya    = w$weight[w$left == "Yamnaya"],
      WHG        = w$weight[w$left == "WHG"],
      p = p, stringsAsFactors = FALSE)
  }
}
write.csv(do.call(rbind, qrows), file.path(outdir, "admixtools_qpadm.csv"), row.names = FALSE)
cat("Wrote admixtools_fstats.csv and admixtools_qpadm.csv\n")
