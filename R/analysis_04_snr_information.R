
source(file.path("R", "helpers.R"))
OUT <- file.path(p26_outdir("analysis_04_snr_information"), "outcome")

D <- p26_data()
d <- D$effects
d <- d[order(d$study, d$effect_ID), ]

sink(file.path(OUT, "console.txt"), split = TRUE)
cat("A04 - where the SNR information lives\n\n")

cat("=== 1. does each study manipulate SNR at all? ===\n")
per <- do.call(rbind, lapply(split(d, d$study), function(s) data.frame(
  study = s$study[1], n_effects = nrow(s),
  n_snr_levels = length(unique(s$SNR)),
  snr_min = min(s$SNR), snr_max = max(s$SNR),
  snr_range = diff(range(s$SNR)),
  outcomes = paste(sort(unique(s$outcome)), collapse = "/"))))
per <- per[order(-per$snr_range, -per$n_snr_levels), ]
print(per, row.names = FALSE)
cat(sprintf("\nstudies with >1 SNR level: %d of %d (%.0f%%)\n",
            sum(per$n_snr_levels > 1), nrow(per), 100 * mean(per$n_snr_levels > 1)))
cat(sprintf("effects in single-SNR studies: %d of %d -- these carry NO within-study\n",
            sum(per$n_effects[per$n_snr_levels == 1]), nrow(d)))
cat("  information about the SNR slope, however many of them there are.\n")
write.csv(per, file.path(OUT, "snr_by_study.csv"), row.names = FALSE)

cat("\nby outcome family:\n")
for (nm in c("IS", "RT", "PE")) {
  s <- d[d$outcome == nm, ]
  pl <- tapply(s$SNR, s$study, function(x) length(unique(x)))
  cat(sprintf("  %-3s %3d effects, %2d studies, %2d of them vary SNR (%d effects in those)\n",
              nm, nrow(s), length(pl), sum(pl > 1),
              sum(s$study %in% names(pl)[pl > 1])))
}

cat("\n=== 2. how much SNR variation is within studies? ===\n")
gm <- ave(d$SNR, d$study, FUN = mean)
ss_w <- sum((d$SNR - gm)^2); ss_b <- sum((gm - mean(d$SNR))^2)
cat(sprintf("total SS %.1f = within-study %.1f (%.1f%%) + between-study %.1f (%.1f%%)\n",
            ss_w + ss_b, ss_w, 100 * ss_w / (ss_w + ss_b), ss_b, 100 * ss_b / (ss_w + ss_b)))
for (nm in c("IS", "RT", "PE")) {
  s <- d[d$outcome == nm, ]
  g <- ave(s$SNR, s$study, FUN = mean)
  w <- sum((s$SNR - g)^2); b <- sum((g - mean(s$SNR))^2)
  cat(sprintf("  %-3s within %.1f%% / between %.1f%%\n", nm, 100 * w / (w + b), 100 * b / (w + b)))
}

cat("\n=== 3. per-study SNR slopes, where estimable ===\n")
sl <- list()
for (nm in c("IS", "RT", "PE")) {
  s <- d[d$outcome == nm, ]
  for (st in unique(s$study)) {
    ss <- s[s$study == st, ]
    if (length(unique(ss$SNR)) < 2 || nrow(ss) < 3) next
    m <- try(rma(yi_al, vi, mods = ~ SNR, data = ss, method = "FE"), silent = TRUE)
    if (inherits(m, "try-error")) next
    sl[[length(sl) + 1]] <- data.frame(outcome = nm, study = st, k = nrow(ss),
                                       levels = length(unique(ss$SNR)),
                                       slope = m$b[2], se = m$se[2])
  }
}
sl <- do.call(rbind, sl)
print(sl[order(sl$outcome, sl$slope), ], row.names = FALSE, digits = 3)
cat("\nper-family summary of the study-specific slopes:\n")
for (nm in unique(sl$outcome)) {
  s <- sl[sl$outcome == nm, ]
  cat(sprintf("  %-3s %2d estimable studies | median %+.3f | range %+.3f to %+.3f | %d negative\n",
              nm, nrow(s), median(s$slope), min(s$slope), max(s$slope), sum(s$slope < 0)))
}
write.csv(sl, file.path(OUT, "per_study_slopes.csv"), row.names = FALSE)

cat("\n=== 4. THE decomposition: within-study vs between-study SNR effect ===\n")
d$SNR_between <- ave(d$SNR, d$study, FUN = mean)
d$SNR_within  <- d$SNR - d$SNR_between
d$SNR_between_c <- d$SNR_between - mean(d$SNR_between)

dec <- list()
for (nm in c("ALL", "IS", "RT", "PE")) {
  dd <- if (nm == "ALL") d else d[d$outcome == nm, ]
  V <- vcalc(vi = vi, cluster = study, grp1 = experiment_group_id,
             grp2 = control_group_id, w1 = n1, w2 = n2, data = dd, sparse = FALSE)
  m <- try(rma.mv(yi_al, V, mods = ~ SNR_within + SNR_between_c, data = dd,
                  random = list(~ 1 | study, ~ 1 | effect_ID),
                  method = "REML", test = "t", dfs = "contain"), silent = TRUE)
  if (inherits(m, "try-error")) { cat(nm, ": model failed\n"); next }
  ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = dd$study)
  for (i in 2:3) {
    dec[[length(dec) + 1]] <- data.frame(
      family = nm, term = rownames(m$b)[i], estimate = m$b[i],
      ci_lo = m$ci.lb[i], ci_hi = m$ci.ub[i],
      p_model = m$pval[i], p_CR2 = ct$p_Satt[i], df_CR2 = ct$df_Satt[i])
  }
  cat(sprintf("%-4s within %+.4f (CR2 p %.4f, df %.1f) | between %+.4f (CR2 p %.4f, df %.1f)\n",
              nm, m$b[2], ct$p_Satt[2], ct$df_Satt[2], m$b[3], ct$p_Satt[3], ct$df_Satt[3]))
}
dec <- do.call(rbind, dec)
write.csv(dec, file.path(OUT, "within_between_decomposition.csv"), row.names = FALSE)
cat("\n  A within-study slope is the interpretable quantity: same children, worse SNR.\n")
cat("  A between-study slope compares different literatures and is confounded.\n")

cat("\n=== 5. leave-one-study-out on the IS slope (is it one study?) ===\n")
dis <- d[d$outcome == "IS", ]
loo <- list()
for (st in c("(none)", unique(dis$study))) {
  ss <- if (st == "(none)") dis else dis[dis$study != st, ]
  if (length(unique(ss$SNR)) < 2) next
  V <- vcalc(vi = vi, cluster = study, grp1 = experiment_group_id,
             grp2 = control_group_id, w1 = n1, w2 = n2, data = ss, sparse = FALSE)
  m <- try(rma.mv(yi_al, V, mods = ~ SNR_c + age_c + masking, data = ss,
                  random = list(~ 1 | study, ~ 1 | effect_ID),
                  method = "REML", test = "t", dfs = "contain"), silent = TRUE)
  if (inherits(m, "try-error")) next
  ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = ss$study)
  loo[[length(loo) + 1]] <- data.frame(dropped = st, k = nrow(ss),
                                       slope = m$b[2], p_CR2 = ct$p_Satt[2],
                                       df = ct$df_Satt[2])
}
loo <- do.call(rbind, loo)
loo <- loo[order(loo$slope), ]
print(loo, row.names = FALSE, digits = 4)
full <- loo$slope[loo$dropped == "(none)"]
cat(sprintf("\nfull-sample slope %+.4f; leave-one-out range %+.4f to %+.4f (%.1f%% swing)\n",
            full, min(loo$slope), max(loo$slope),
            100 * (max(loo$slope) - min(loo$slope)) / abs(full)))
cat(sprintf("still significant at .05 after dropping any single study: %d of %d\n",
            sum(loo$p_CR2[loo$dropped != "(none)"] < .05), sum(loo$dropped != "(none)")))
write.csv(loo, file.path(OUT, "leave_one_out_IS.csv"), row.names = FALSE)

cat("\nwritten to", OUT, "\n")
sink()
