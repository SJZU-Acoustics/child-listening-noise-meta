
source(file.path("R", "helpers.R"))
OUT <- file.path(p26_outdir("analysis_06_publication_bias"), "outcome")

D <- p26_data()
d <- D$effects
d <- d[order(d$study, d$effect_ID), ]

sink(file.path(OUT, "console.txt"), split = TRUE)
cat("A06 - publication bias / small-study effects\n\n")

mkV <- function(sub) vcalc(vi = vi, cluster = study, grp1 = experiment_group_id,
                           grp2 = control_group_id, w1 = n1, w2 = n2,
                           data = sub, sparse = FALSE)

d$sqrt_inv_n <- sqrt(1 / pmin(d$n1, d$n2))
d$se <- sqrt(d$vi)

cat("=== 1. two Egger variants, each with CR2 clustering ===\n")
rows <- list()
for (nm in c("ALL", "IS", "RT", "PE")) {
  dd <- if (nm == "ALL") d else d[d$outcome == nm, ]
  V <- mkV(dd)
  for (variant in c("classic (SE)", "sample-size (sqrt(1/n))")) {
    f <- if (variant == "classic (SE)") ~ se else ~ sqrt_inv_n
    m <- try(rma.mv(yi_al, V, mods = f, data = dd,
                    random = list(~ 1 | study, ~ 1 | effect_ID),
                    method = "REML", test = "t", dfs = "contain"), silent = TRUE)
    if (inherits(m, "try-error")) { cat(nm, variant, ": failed\n"); next }
    ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = dd$study)
    rows[[length(rows) + 1]] <- data.frame(
      family = nm, variant = variant, slope = m$b[2],
      ci_lo = m$ci.lb[2], ci_hi = m$ci.ub[2],
      p_model = m$pval[2], p_CR2 = ct$p_Satt[2], df_CR2 = ct$df_Satt[2],
      intercept = m$b[1])
    cat(sprintf("%-4s %-24s slope %+.4f  p_model %-9.3g CR2 p %-8.4f df %4.1f | adjusted intercept %+.3f\n",
                nm, variant, m$b[2], m$pval[2], ct$p_Satt[2], ct$df_Satt[2], m$b[1]))
  }
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUT, "egger_tests.csv"), row.names = FALSE)

cat("\n  The adjusted intercept is the estimate the test implies at infinite\n")
cat("  precision. Compare it with the pooled estimate to judge magnitude, not\n")
cat("  just significance.\n")

cat("\n=== 2. how much would the pooled estimate move? ===\n")
for (nm in c("ALL", "IS", "RT", "PE")) {
  dd <- if (nm == "ALL") d else d[d$outcome == nm, ]
  m <- rma.mv(yi_al, mkV(dd), data = dd, random = list(~ 1 | study, ~ 1 | effect_ID),
              method = "REML", test = "t", dfs = "contain")
  s <- res[res$family == nm & res$variant == "sample-size (sqrt(1/n))", ]
  cat(sprintf("%-4s pooled %+.3f -> adjusted intercept %+.3f  (shift %+.3f, %.0f%%)\n",
              nm, m$b[1], s$intercept, s$intercept - m$b[1],
              100 * (s$intercept - m$b[1]) / abs(m$b[1])))
}

cat("\n=== 3. residual funnel data (exported for inspection) ===\n")
fun <- list()
for (nm in c("IS", "RT", "PE")) {
  dd <- d[d$outcome == nm, ]
  m <- rma.mv(yi_al, mkV(dd), data = dd, random = list(~ 1 | study, ~ 1 | effect_ID),
              method = "REML", test = "t", dfs = "contain")
  fun[[length(fun) + 1]] <- data.frame(
    outcome = nm, study = dd$study, effect_id = dd$effect_ID,
    yi = dd$yi_al, se = sqrt(dd$vi), n = pmin(dd$n1, dd$n2),
    residual = as.numeric(residuals(m)))
}
fun <- do.call(rbind, fun)
write.csv(fun, file.path(OUT, "funnel_data.csv"), row.names = FALSE)
cat(sprintf("  %d rows written; residual mean %.4f, SD %.3f\n",
            nrow(fun), mean(fun$residual), sd(fun$residual)))

cat("\n  asymmetry of residuals about zero, by family:\n")
for (nm in unique(fun$outcome)) {
  s <- fun[fun$outcome == nm, ]
  small <- s$n <= median(s$n)
  cat(sprintf("  %-3s small-sample half mean residual %+.3f | large-sample half %+.3f | gap %+.3f\n",
              nm, mean(s$residual[small]), mean(s$residual[!small]),
              mean(s$residual[small]) - mean(s$residual[!small])))
}

cat("\n=== 4. what these tests can actually detect ===\n")
for (nm in c("ALL", "IS", "RT", "PE")) {
  dd <- if (nm == "ALL") d else d[d$outcome == nm, ]
  cat(sprintf("  %-4s %2d clusters, n range %d-%d (median %d)\n", nm,
              length(unique(dd$study)), min(pmin(dd$n1, dd$n2)),
              max(pmin(dd$n1, dd$n2)), median(pmin(dd$n1, dd$n2))))
}
cat("  Egger-type tests are widely held to need ~10+ independent studies and a\n")
cat(sprintf("  wide spread of precisions. PE has %d clusters; RT has %d.\n",
            length(unique(d$study[d$outcome == "PE"])), length(unique(d$study[d$outcome == "RT"]))))

cat("\nwritten to", OUT, "\n")
sink()
