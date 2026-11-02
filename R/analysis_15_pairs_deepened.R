
source(file.path("R", "helpers.R"))
OUT <- file.path(p26_outdir("analysis_15_pairs_deepened"), "outcome")

D <- p26_data()
p <- D$pairs
p <- p[order(p$study, p$pair_ID), ]

sink(file.path(OUT, "console.txt"), split = TRUE)
cat("A15 - the aligned IS-RT pairs, deepened\n\n")

cat("=== 0. what the pairs can identify ===\n")
sw <- tapply(p$SNR, p$study, function(x) length(unique(x)) > 1)
cat(sprintf("pairs: %d, studies: %d\n", nrow(p), length(unique(p$study))))
cat(sprintf("studies varying SNR within the pairs: %d of %d\n", sum(sw), length(sw)))
cat(sprintf("within-study share of pair SNR variation: %.1f%%\n",
            100 * sum((p$SNR - ave(p$SNR, p$study, FUN = mean))^2) /
              sum((p$SNR - mean(p$SNR))^2)))
bw <- tapply(p$IS_m2, p$study, function(x) length(unique(round(x, 6))) > 1)
cat(sprintf("studies varying IS baseline accuracy within the pairs: %d of %d\n",
            sum(bw), length(bw)))

cat("\n=== 1. ISRT_diff ~ SNR, within/between (campaign rule), CWB inference ===\n")
p$SNR_between <- ave(p$SNR, p$study, FUN = mean)
p$SNR_within  <- p$SNR - p$SNR_between
p$SNR_between_c <- p$SNR_between - mean(p$SNR_between)
m_dec <- rma.mv(diff ~ SNR_within + SNR_between_c, v_diff,
                random = list(~ 1 | study, ~ 1 | pair_ID),
                data = p, method = "REML", test = "t")
print(data.frame(term = rownames(m_dec$b), est = as.numeric(m_dec$b),
                 p_model = m_dec$pval), row.names = FALSE, digits = 4)
cwb_dec <- rep(NA_real_, 3)
for (i in 2:3) {
  bt <- try(Wald_test_cwb(full_model = m_dec, constraints = constrain_zero(i),
                          R = BOOT_R, seed = 26), silent = TRUE)
  if (inherits(bt, "try-error")) {
    cat(sprintf("  CWB failed for %s: %s\n", rownames(m_dec$b)[i],
                conditionMessage(attr(bt, "condition"))))
  } else {
    cwb_dec[i] <- bt$p_val
    cat(sprintf("  %-16s est %+.4f  CWB p = %.4f\n",
                rownames(m_dec$b)[i], m_dec$b[i], bt$p_val))
  }
}

cat("\n=== 2. the quadrant structure of (g_IS, g_RT) ===\n")
cat("positive = more adverse in noise on both axes\n")
p$quad <- ifelse(p$g_IS > 0 & p$g_RT > 0, "both adverse",
          ifelse(p$g_IS > 0 & p$g_RT <= 0, "accuracy only",
          ifelse(p$g_IS <= 0 & p$g_RT > 0, "speed only", "neither")))
qq <- table(p$quad)
print(qq)
cat(sprintf("\nby task_group:\n"))
print(table(p$quad, p$task_group))
cat(sprintf("\nmedian g_IS %+.3f, median g_RT %+.3f; concordant sign pairs: %d of %d\n",
            median(p$g_IS), median(p$g_RT),
            sum(p$quad %in% c("both adverse", "neither")), nrow(p)))

cat("\n=== 3. cor(g_IS, g_RT): common cost (positive) or trade-off (negative)? ===\n")
cat(sprintf("Pearson %+.3f, Spearman %+.3f (naive, k=%d pairs)\n",
            cor(p$g_IS, p$g_RT), cor(p$g_IS, p$g_RT, method = "spearman"), nrow(p)))
set.seed(26)
studies <- unique(p$study)
boot <- replicate(BOOT_R, {
  ss <- sample(studies, replace = TRUE)
  idx <- unlist(lapply(seq_along(ss), function(i) which(p$study == ss[i])))
  if (length(unique(p$study[idx])) < 3 || sd(p$g_IS[idx]) == 0 || sd(p$g_RT[idx]) == 0)
    return(NA_real_)
  cor(p$g_IS[idx], p$g_RT[idx])
})
boot <- boot[!is.na(boot)]
cat(sprintf("cluster bootstrap (resampling %d studies, %d draws): median %+.3f, 95%% CI [%+.3f, %+.3f]\n",
            length(studies), length(boot), median(boot),
            quantile(boot, .025), quantile(boot, .975)))

cat("\n=== 4. does the dissociation depend on task difficulty (IS baseline)? ===\n")
cat("IS baseline = quiet-condition accuracy of the pair's IS side (percent)\n")
p$base_between <- ave(p$IS_m2, p$study, FUN = mean)
p$base_within  <- p$IS_m2 - p$base_between
p$base_between_c <- p$base_between - mean(p$base_between)
cat(sprintf("IS baseline: median %.1f%%, range %.1f-%.1f\n",
            median(p$IS_m2), min(p$IS_m2), max(p$IS_m2)))
m_base <- rma.mv(diff ~ base_within + base_between_c, v_diff,
                 random = list(~ 1 | study, ~ 1 | pair_ID),
                 data = p, method = "REML", test = "t")
print(data.frame(term = rownames(m_base$b), est = as.numeric(m_base$b),
                 p_model = m_base$pval), row.names = FALSE, digits = 4)
cwb_base <- rep(NA_real_, 3)
for (i in 2:3) {
  bt <- try(Wald_test_cwb(full_model = m_base, constraints = constrain_zero(i),
                          R = BOOT_R, seed = 26), silent = TRUE)
  if (inherits(bt, "try-error")) {
    cat(sprintf("  CWB failed for %s: %s\n", rownames(m_base$b)[i],
                conditionMessage(attr(bt, "condition"))))
  } else {
    cwb_base[i] <- bt$p_val
    cat(sprintf("  %-16s est %+.4f  CWB p = %.4f\n",
                rownames(m_base$b)[i], m_base$b[i], bt$p_val))
  }
}
cat("  prediction of the ceiling account: NEGATIVE base slope (higher quiet\n")
cat("  accuracy -> cost appears in RT -> ISRT_diff falls)\n")

quad_out <- as.data.frame(table(p$quad, p$task_group))
names(quad_out) <- c("quadrant", "task_group", "pairs")
write.csv(quad_out, file.path(OUT, "quadrants_by_task.csv"), row.names = FALSE)
write.csv(data.frame(
  term = rownames(m_dec$b), est = as.numeric(m_dec$b), p_model = m_dec$pval,
  p_CWB = cwb_dec), file.path(OUT, "snr_decomposition_cwb.csv"), row.names = FALSE)
write.csv(data.frame(
  term = rownames(m_base$b), est = as.numeric(m_base$b), p_model = m_base$pval,
  p_CWB = cwb_base), file.path(OUT, "difficulty_moderation_cwb.csv"), row.names = FALSE)
write.csv(data.frame(
  pearson = cor(p$g_IS, p$g_RT), spearman = cor(p$g_IS, p$g_RT, method = "spearman"),
  boot_median = median(boot), boot_lo = quantile(boot, .025),
  boot_hi = quantile(boot, .975)),
  file.path(OUT, "pair_correlation.csv"), row.names = FALSE)
write.csv(p[, c("pair_ID", "study", "task_group", "SNR", "IS_m2", "g_IS",
                "g_RT", "diff", "quad")],
          file.path(OUT, "pair_level_data.csv"), row.names = FALSE)
cat("\nwritten to", OUT, "\n")
sink()
