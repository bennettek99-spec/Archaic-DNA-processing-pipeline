#!/usr/bin/env Rscript
# Etruscan qpWave/qpAdm cross-check for the V2 publication gap.
#
# Reuses the cached f2 from tools/admixtools_concordance.R (results/f2_admixtools).
# Produces, for each Italian target cohort:
#   1. qpWave rank tests over the source+target set (how many ancestry streams
#      the left/right sets are connected by), and
#   2. the full qpAdm nested-model "popdrop" ladder (dropping each source in
#      turn), so every reduced model's fit p-value is tabled, not just the
#      top-level 3-way mixture.
#
# The 3-way (Anatolia_N + Yamnaya + WHG) is rejected at full SNP density; the
# 4-way (adding Iran_N) is the accepted model (tools/qpadm4.R). This script
# records the rejection evidence formally for the paper.
#
# Usage (Rtools45 on PATH first):
#   Rscript tools/etruscan_qpwave_qpadm.R results/f2_admixtools results/etruscan_admixtools

args   <- commandArgs(trailingOnly = TRUE)
f2dir  <- ifelse(length(args) >= 1, args[1], "results/f2_admixtools")
outdir <- ifelse(length(args) >= 2, args[2], "results/etruscan_admixtools")
suppressMessages(library(admixtools))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
f2 <- f2_from_precomp(f2dir, verbose = FALSE)

src3 <- c("Anatolia_N", "Yamnaya", "WHG")
src4 <- c("Anatolia_N", "Yamnaya", "WHG", "Iran_N")
right <- c("Mbuti", "Han", "Papuan", "Karitiana", "Natufian", "Ust_Ishim", "MA1")
targets <- c("Etruscan", "Latin", "ImperialRoman", "ItalyBA")

# ---- qpWave rank tests ------------------------------------------------------
wave_rows <- list()
for (tg in targets) {
  for (src in list(src3 = src3, src4 = src4)) {
    r <- tryCatch(
      qpwave(f2, left = c(tg, src), right = right, verbose = FALSE),
      error = function(e) NULL)
    if (is.null(r)) next
    rd <- r$rankdrop
    for (k in seq_len(nrow(rd))) {
      wave_rows[[length(wave_rows) + 1]] <- data.frame(
        target = tg,
        model = paste(src, collapse = "+"),
        f4rank = rd$f4rank[k],
        dof = rd$dof[k],
        chi2 = rd$chisq[k],
        p = rd$p[k],
        p_nested = rd$p_nested[k],
        stringsAsFactors = FALSE)
    }
  }
}
wave_df <- do.call(rbind, wave_rows)
write.csv(wave_df, file.path(outdir, "qpwave.csv"), row.names = FALSE)
cat("Wrote", file.path(outdir, "qpwave.csv"), "\n")

# ---- qpAdm nested-model (popdrop) ladder ------------------------------------
drop_rows <- list()
for (tg in targets) {
  for (src in list(src3 = src3, src4 = src4)) {
    full <- paste(src, collapse = "+")
    r <- tryCatch(
      qpadm(f2, left = src, right = right, target = tg, verbose = FALSE),
      error = function(e) NULL)
    if (is.null(r)) next
    pd <- r$popdrop
    for (k in seq_len(nrow(pd))) {
      drop_rows[[length(drop_rows) + 1]] <- data.frame(
        target = tg,
        full_model = full,
        pattern = pd$pat[k],
        dropped = ifelse(pd$pat[k] == "000" | nchar(pd$pat[k]) != length(src),
                         "(full)", pd$pat[k]),
        dof = pd$dof[k],
        chi2 = pd$chisq[k],
        p = pd$p[k],
        feasible = pd$feasible[k],
        stringsAsFactors = FALSE)
    }
  }
}
drop_df <- do.call(rbind, drop_rows)
write.csv(drop_df, file.path(outdir, "qpadm_popdrop.csv"), row.names = FALSE)
cat("Wrote", file.path(outdir, "qpadm_popdrop.csv"), "\n")

# ---- concise console summary of the acceptance/rejection boundary -------------
cat("\n=== qpWave rank tests (p < 0.05 = that rank still needed) ===\n")
for (tg in targets) {
  w <- wave_df[wave_df$target == tg, ]
  if (nrow(w)) {
    cat(sprintf("\n%s:\n", tg))
    for (k in seq_len(nrow(w))) {
      cat(sprintf("  %-24s f4rank=%d chi2=%.1f p=%.4g p_nested=%.4g\n",
                  w$model[k], w$f4rank[k], w$chi2[k], w$p[k], w$p_nested[k]))
    }
  }
}

cat("\n=== qpAdm popdrop ladder (p < 0.05 = REJECTED) ===\n")
for (tg in targets) {
  d <- drop_df[drop_df$target == tg, ]
  if (nrow(d)) {
    cat(sprintf("\n%s:\n", tg))
    for (k in seq_len(nrow(d))) {
      p <- d$p[k]
      flag <- ifelse(is.na(p), "", ifelse(p < 0.05, "(REJECTED)", ""))
      cat(sprintf("  %-24s drop=%-8s p=%.4g %s\n",
                  d$full_model[k], d$dropped[k], p, flag))
    }
  }
}
cat("\nDone.\n")
