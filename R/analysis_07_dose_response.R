
source(file.path("R", "helpers.R"))
OUT <- file.path(p26_outdir("analysis_07_dose_response"), "outcome")

D <- p26_data()
d <- D$effects
d <- d[order(d$study, d$effect_ID), ]
d$SNR_between <- ave(d$SNR, d$study, FUN = mean)
d$SNR_within  <- d$SNR - d$SNR_between
d$SNR_between_c <- d$SNR_between - mean(d$SNR_between)
d$SNR_within2 <- d$SNR_within^2

mkV <- function(sub) vcalc(vi = vi, cluster = study, grp1 = experiment_group_id,
                           grp2 = control_group_id, w1 = n1, w2 = n2,
                           data = sub, sparse = FALSE)

sink(file.path(OUT, "console.txt"), split = TRUE)
cat("A07 - dose-response shape\n\n")

cat("=== 1. what can carry a curvature test ===\n")
lv <- tapply(d$SNR, d$study, function(x) length(unique(x)))
cat(sprintf("studies with >=3 SNR levels: %d (%s)\n", sum(lv >= 3),
            paste(names(lv)[lv >= 3], collapse = "; ")))
for (nm in c("IS", "RT", "PE")) {
  s <- d[d$outcome == nm, ]
  l <- tapply(s$SNR, s$study, function(x) length(unique(x)))
  cat(sprintf("  %-3s %d studies with >=3 levels, %d effects in them\n",
              nm, sum(l >= 3), sum(s$study %in% names(l)[l >= 3])))
}

cat("\n=== 2. quadratic term on the within-study SNR component ===\n")
rows <- list()
for (nm in c("ALL", "IS", "RT")) {
  dd <- if (nm == "ALL") d else d[d$outcome == nm, ]
  l <- tapply(dd$SNR, dd$study, function(x) length(unique(x)))
  dd <- dd[dd$study %in% names(l)[l >= 3], ]
  if (length(unique(dd$study)) < 3) { cat(sprintf("  %-4s too few studies\n", nm)); next }
  V <- mkV(dd)
  m <- try(rma.mv(yi_al, V, mods = ~ SNR_within + SNR_within2 + SNR_between_c,
                  data = dd, random = list(~ 1 | study, ~ 1 | effect_ID),
                  method = "REML", test = "t", dfs = "contain"), silent = TRUE)
  if (inherits(m, "try-error")) { cat(sprintf("  %-4s failed\n", nm)); next }
  ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = dd$study)
  i <- which(rownames(m$b) == "SNR_within2")
  cat(sprintf("%-4s %2d studies k=%3d | linear %+.4f | quadratic %+.5f  p_model %.4g  CR2 p %.4f (df %.1f)\n",
              nm, length(unique(dd$study)), nrow(dd),
              m$b[which(rownames(m$b) == "SNR_within")], m$b[i],
              m$pval[i], ct$p_Satt[i], ct$df_Satt[i]))
  rows[[length(rows) + 1]] <- data.frame(
    family = nm, n_studies = length(unique(dd$study)), k = nrow(dd),
    linear = m$b[which(rownames(m$b) == "SNR_within")], quadratic = m$b[i],
    p_model = m$pval[i], p_CR2 = ct$p_Satt[i], df_CR2 = ct$df_Satt[i])
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUT, "curvature_tests.csv"), row.names = FALSE)

cat("\n=== 3. per-study shape, studies with >=3 levels (descriptive) ===\n")
shape <- list()
for (st in names(lv)[lv >= 3]) {
  for (nm in c("IS", "RT", "PE")) {
    ss <- d[d$study == st & d$outcome == nm, ]
    if (length(unique(ss$SNR)) < 3) next
    agg <- aggregate(yi_al ~ SNR, ss, mean)
    agg <- agg[order(agg$SNR), ]
    cat(sprintf("  %-22s %-3s  ", st, nm))
    cat(paste(sprintf("%+.0fdB:%+.2f", agg$SNR, agg$yi_al), collapse = "  "), "\n")
    shape[[length(shape) + 1]] <- data.frame(study = st, outcome = nm,
                                             SNR = agg$SNR, mean_effect = agg$yi_al)
  }
}
shape <- do.call(rbind, shape)
write.csv(shape, file.path(OUT, "per_study_shape.csv"), row.names = FALSE)

cat("\n  monotone in the expected direction (worse effect at worse SNR)?\n")
for (i in unique(paste(shape$study, shape$outcome))) {
  s <- shape[paste(shape$study, shape$outcome) == i, ]
  s <- s[order(s$SNR), ]
  mono <- all(diff(s$mean_effect) <= 0)
  cat(sprintf("  %-28s %s\n", i, if (mono) "monotone decreasing" else "NOT monotone"))
}

cat("\n=== 4. is the IS scale bounded where these studies sit? ===\n")
cat("  (a curvature finding in IS could be a ceiling artefact rather than a\n")
cat("   property of listening -- accuracy cannot exceed 100%)\n")
src <- read.csv(file.path(P26_ROOT, ".cache", "data", "source_rows.csv"),
                stringsAsFactors = FALSE)
unit <- setNames(src$unit, src$source_id)
dis <- d[d$outcome == "IS", ]
dis$ctrl_unit <- unit[dis$control_group_id]
pc <- dis[dis$ctrl_unit != "RAU", ]
cat(sprintf("  IS control-arm accuracy (n=%d, RAU excluded): median %.1f%%, %.0f%% above 90%%\n",
            nrow(pc), median(pc$m2), 100 * mean(pc$m2 > 90)))
cat(sprintf("  noise-arm accuracy: median %.1f%%, %.0f%% above 90%%\n",
            median(pc$m1), 100 * mean(pc$m1 > 90)))
cat("  So the control condition sits near ceiling and the noise condition does not:\n")
cat("  the measurable range is asymmetric, which by itself bends the dose-response.\n")

cat("\nwritten to", OUT, "\n")
sink()
