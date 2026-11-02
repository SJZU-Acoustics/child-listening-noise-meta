
source(file.path("R", "helpers.R"))
OUT <- file.path(p26_outdir("analysis_05_heterogeneity"), "outcome")

D <- p26_data()
d <- D$effects
d <- d[order(d$study, d$effect_ID), ]

sink(file.path(OUT, "console.txt"), split = TRUE)
cat("A05 - heterogeneity structure\n\n")

mkV <- function(sub) vcalc(vi = vi, cluster = study, grp1 = experiment_group_id,
                           grp2 = control_group_id, w1 = n1, w2 = n2,
                           data = sub, sparse = FALSE)

cat("=== 1. tau^2 decomposition and multilevel I^2 (shared-arm V) ===\n")
het <- list()
for (nm in c("ALL", "IS", "RT", "PE")) {
  dd <- if (nm == "ALL") d else d[d$outcome == nm, ]
  V <- mkV(dd)
  m <- rma.mv(yi_al, V, data = dd, random = list(~ 1 | study, ~ 1 | effect_ID),
              method = "REML", test = "t", dfs = "contain")
  W <- solve(V); X <- model.matrix(m)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  v_typ <- (m$k - m$p) / sum(diag(P))
  tot <- sum(m$sigma2) + v_typ
  het[[length(het) + 1]] <- data.frame(
    family = nm, k = m$k, clusters = length(unique(dd$study)),
    tau2_study = m$sigma2[1], tau2_within = m$sigma2[2], v_typical = v_typ,
    I2_total = 100 * sum(m$sigma2) / tot,
    I2_study = 100 * m$sigma2[1] / tot,
    I2_within = 100 * m$sigma2[2] / tot,
    QE = m$QE, QE_p = m$QEp)
  cat(sprintf("%-4s k=%3d clusters=%2d | tau2 study %.4f / within %.4f | I2 %.1f%% (study %.1f%%, within %.1f%%) | Q=%.0f p=%.3g\n",
              nm, m$k, length(unique(dd$study)), m$sigma2[1], m$sigma2[2],
              100 * sum(m$sigma2) / tot, 100 * m$sigma2[1] / tot,
              100 * m$sigma2[2] / tot, m$QE, m$QEp))
}
het <- do.call(rbind, het)
write.csv(het, file.path(OUT, "heterogeneity_decomposition.csv"), row.names = FALSE)

cat("\n  Which level dominates decides what a reader should take away:\n")
for (i in seq_len(nrow(het))) {
  dom <- if (het$I2_study[i] > het$I2_within[i]) "BETWEEN studies" else "WITHIN studies"
  cat(sprintf("  %-4s %s dominates (%.0f%% vs %.0f%%)\n", het$family[i], dom,
              max(het$I2_study[i], het$I2_within[i]),
              min(het$I2_study[i], het$I2_within[i])))
}

cat("\n=== 2. prediction intervals - what the pooled effect licenses ===\n")
pi <- list()
for (nm in c("ALL", "IS", "RT", "PE")) {
  dd <- if (nm == "ALL") d else d[d$outcome == nm, ]
  m <- rma.mv(yi_al, mkV(dd), data = dd, random = list(~ 1 | study, ~ 1 | effect_ID),
              method = "REML", test = "t", dfs = "contain")
  pr <- predict(m)
  pi[[length(pi) + 1]] <- data.frame(family = nm, estimate = m$b[1],
                                     ci_lo = m$ci.lb[1], ci_hi = m$ci.ub[1],
                                     pi_lo = pr$pi.lb, pi_hi = pr$pi.ub,
                                     pi_covers_zero = pr$pi.lb < 0)
  cat(sprintf("%-4s g = %+.3f  CI [%+.3f, %+.3f]  PI [%+.3f, %+.3f]%s\n",
              nm, m$b[1], m$ci.lb[1], m$ci.ub[1], pr$pi.lb, pr$pi.ub,
              if (pr$pi.lb < 0) "   <- PI includes zero" else ""))
}
pi <- do.call(rbind, pi)
write.csv(pi, file.path(OUT, "prediction_intervals.csv"), row.names = FALSE)

cat("\n=== 3. the legacy Table 2 transposition, settled on current data ===\n")
cat("The 0605 draft reported its two tau^2 components the wrong way round.\n")
cat("Correct values under the adopted specification, for the paper:\n")
for (i in seq_len(nrow(het))) {
  cat(sprintf("  %-4s tau2(study, level 3) = %.4f | tau2(effect, level 2) = %.4f | ratio %.2f\n",
              het$family[i], het$tau2_study[i], het$tau2_within[i],
              het$tau2_study[i] / het$tau2_within[i]))
}

cat("\n=== 4. does the SNR slope itself vary across studies? (A04's open question) ===\n")
cat("Testing a random slope only where enough studies manipulate SNR.\n")
for (nm in c("ALL", "IS")) {
  dd <- if (nm == "ALL") d else d[d$outcome == nm, ]
  dd <- dd[dd$study %in% names(which(tapply(dd$SNR, dd$study, function(x)
    length(unique(x))) > 1)), ]
  if (length(unique(dd$study)) < 5) { cat(sprintf("  %-4s too few varying studies\n", nm)); next }
  V <- mkV(dd)
  m0 <- try(rma.mv(yi_al, V, mods = ~ SNR_c, data = dd,
                   random = list(~ 1 | study, ~ 1 | effect_ID),
                   method = "ML", test = "t", dfs = "contain"), silent = TRUE)
  m1 <- try(rma.mv(yi_al, V, mods = ~ SNR_c, data = dd,
                   random = list(~ SNR_c | study, ~ 1 | effect_ID), struct = "GEN",
                   method = "ML", test = "t", dfs = "contain"), silent = TRUE)
  if (inherits(m0, "try-error") || inherits(m1, "try-error")) {
    cat(sprintf("  %-4s random-slope model did not converge\n", nm)); next
  }
  lr <- 2 * (logLik(m1) - logLik(m0)); dfd <- attr(logLik(m1), "df") - attr(logLik(m0), "df")
  cat(sprintf("  %-4s %d studies, k=%d | LRT chi2=%.2f df=%d p=%.4f | SD of slope %.4f\n",
              nm, length(unique(dd$study)), nrow(dd), lr, dfd,
              pchisq(lr, dfd, lower.tail = FALSE), sqrt(m1$tau2[2])))
  cat("       (LRT on a variance component is anticonservative here; read with the\n")
  cat("        per-study slope spread in A04, not on its own)\n")
}

cat("\nwritten to", OUT, "\n")
sink()
