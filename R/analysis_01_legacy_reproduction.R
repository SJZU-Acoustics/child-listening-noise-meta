source("R/helpers.R")
OUT <- file.path(p26_outdir("analysis_01_legacy_reproduction"), "outcome")

D <- p26_data()
d <- D$effects
p <- D$pairs
p$dual <- p$task_group

sink(file.path(OUT, "console.txt"), split = TRUE)
cat("Primary models\n")
cat("metafor", as.character(packageVersion("metafor")),
    "| clubSandwich", as.character(packageVersion("clubSandwich")),
    "| wildmeta", as.character(packageVersion("wildmeta")), "\n")

res <- list()
add <- function(quantity, est, lo, hi, p_model, p_robust = NA, note = "") {
  res[[length(res) + 1]] <<- data.frame(
    quantity = quantity, estimate = est, ci_lo = lo, ci_hi = hi,
    p_model = p_model, p_robust = p_robust, note = note, stringsAsFactors = FALSE)
}

cat("\n=== 1. empty models (audited form: random = study + effect) ===\n")
df_note <- data.frame()
for (nm in c("ALL", "IS", "RT", "PE")) {
  dd <- if (nm == "ALL") d else d[d$outcome == nm, ]
  f <- p26_fit(~ 1, dd); rb <- p26_robust(f)
  aud <- rma.mv(yi_al, vi, random = list(~ 1 | study, ~ 1 | effect_ID),
                data = dd, method = "REML", test = "t")
  cat(sprintf("%-4s g = %+.3f | audited [%+.3f, %+.3f] df %d | adopted [%+.3f, %+.3f] df %d | CR2 p %.3g  k=%d clusters=%d\n",
              nm, f$b[1], aud$ci.lb[1], aud$ci.ub[1], round(aud$ddf[1]),
              f$ci.lb[1], f$ci.ub[1], round(f$ddf[1]), rb$p_CR2[1],
              f$k, length(unique(dd$study))))
  add(paste("empty", nm), f$b[1], f$ci.lb[1], f$ci.ub[1], f$pval[1], rb$p_CR2[1],
      sprintf("adopted spec (dfs=contain); sigma2 study %.3f / effect %.3f",
              f$sigma2[1], f$sigma2[2]))
  df_note <- rbind(df_note, data.frame(
    family = nm, estimate = round(f$b[1], 4),
    ci_audited_df_k = sprintf("[%.3f, %.3f]", aud$ci.lb[1], aud$ci.ub[1]), df_audited = round(aud$ddf[1]),
    ci_adopted_contain = sprintf("[%.3f, %.3f]", f$ci.lb[1], f$ci.ub[1]), df_adopted = round(f$ddf[1]),
    p_adopted = signif(f$pval[1], 3)))
}
write.csv(df_note, file.path(OUT, "empty_model_df_comparison.csv"), row.names = FALSE)

cat("\n=== 2. SNR association per outcome family, RQ1 (audited form:",
    "adjusted for age and masking) ===\n")
snr_p <- c()
for (nm in c("IS", "RT", "PE")) {
  dd <- d[d$outcome == nm, ]
  f <- p26_fit(~ SNR_c + age_c + masking, dd); rb <- p26_robust(f)
  i <- which(rownames(f$b) == "SNR_c")
  cat(sprintf("%-3s slope %+.4f [%+.4f, %+.4f]  p_model %.3g | CR2 p %.3g (df %.1f)\n",
              nm, f$b[i], f$ci.lb[i], f$ci.ub[i], f$pval[i], rb$p_CR2[i], rb$df_CR2[i]))
  add(paste("slope", nm), f$b[i], f$ci.lb[i], f$ci.ub[i], f$pval[i], rb$p_CR2[i],
      "adjusted for age and masking")
  snr_p <- c(snr_p, rb$p_CR2[i]); names(snr_p)[length(snr_p)] <- nm
  if (!is.null(p26_disagree(rb[i, , drop = FALSE])))
    cat("    *** model-based and robust disagree at .05 -> exploratory\n")
  fm <- p26_fit(~ SNR_c, dd); rbm <- p26_robust(fm)
  cat(sprintf("    (marginal, unadjusted: %+.4f  p_model %.3g | CR2 p %.3g)\n",
              fm$b[2], fm$pval[2], rbm$p_CR2[2]))
  add(paste("slope", nm, "(marginal)"), fm$b[2], fm$ci.lb[2], fm$ci.ub[2],
      fm$pval[2], rbm$p_CR2[2], "unadjusted; not the audited quantity")
}
cat("  BH across the three-family RQ1 screen (CR2 p): ",
    paste(sprintf("%s q=%.4f", names(snr_p), p.adjust(snr_p, "BH")), collapse = " | "), "\n")

cat("\n=== 3. moderators (audited form: SNR_c * age_c + Outcome_Type + masking) ===\n")
f <- p26_fit(~ SNR_c * age_c + Outcome_Type + masking, d); rb <- p26_robust(f)
rownames(rb) <- rownames(f$b)
for (term in c("SNR_c:age_c", "age_c", "maskingIM", "Outcome_TypeRT", "Outcome_TypePE")) {
  i <- which(rownames(f$b) == term)
  cat(sprintf("%-16s %+.5f  p_model %.4g | CR2 p %.4g (df %.1f)\n",
              term, f$b[i], f$pval[i], rb$p_CR2[i], rb$df_CR2[i]))
  add(term, f$b[i], f$ci.lb[i], f$ci.ub[i], f$pval[i], rb$p_CR2[i], "full moderator model")
}
dis <- p26_disagree(rb)
cat(if (is.null(dis)) "  no model/robust disagreement in this model\n" else "")
if (!is.null(dis)) { cat("  *** disagreement -> exploratory:\n"); print(dis) }

cat("\n=== 4. Same-task versus cross-task aligned pairs ===\n")
m2 <- rma.mv(diff ~ dual, v_diff, random = list(~ 1 | study, ~ 1 | pair_ID),
             data = p, method = "REML", test = "t")
bt <- Wald_test_cwb(full_model = m2, constraints = constrain_zero(2),
                    R = BOOT_R, seed = 26)
add("RQ2 dual-task vs same-task", m2$b[2], m2$ci.lb[2], m2$ci.ub[2], m2$pval[2], bt$p_val,
    paste0("robust p = cluster wild bootstrap, ", length(unique(p$study)), " clusters"))
out <- do.call(rbind, res)
write.csv(out, file.path(OUT, "legacy_reproduction.csv"), row.names = FALSE)
sink()
