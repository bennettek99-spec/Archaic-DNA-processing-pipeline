suppressMessages(library(admixtools))
f2 <- f2_from_precomp("results/f2_admixtools", verbose = FALSE)
src <- c("Anatolia_N", "Yamnaya", "WHG")
rights <- list(
  r4 = c("Mbuti", "Han", "Papuan", "Karitiana"),
  r5 = c("Mbuti", "Han", "Papuan", "Karitiana", "Natufian"),
  r6 = c("Mbuti", "Han", "Papuan", "Karitiana", "Natufian", "Ust_Ishim"),
  r8 = c("Mbuti", "Han", "Papuan", "Karitiana", "Natufian", "Ust_Ishim", "MA1", "Iran_N"))
for (nm in names(rights)) {
  r <- tryCatch(qpadm(f2, left = src, right = rights[[nm]], target = "Etruscan", verbose = FALSE),
                error = function(e) { cat(nm, "ERROR", conditionMessage(e), "\n"); NULL })
  if (!is.null(r)) {
    w <- r$weights
    cat(sprintf("%s |right|=%d  AnatN=%.0f Steppe=%.0f WHG=%.0f  p=%.4f\n",
                nm, length(rights[[nm]]),
                100 * w$weight[w$left == "Anatolia_N"],
                100 * w$weight[w$left == "Yamnaya"],
                100 * w$weight[w$left == "WHG"], r$popdrop$p[1]))
  }
}
