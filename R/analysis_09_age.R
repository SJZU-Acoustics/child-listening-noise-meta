
source(file.path("R", "helpers.R"))
OUT <- file.path(p26_outdir("analysis_09_age"), "outcome")

D <- p26_data()
d <- D$effects
d <- d[order(d$study, d$effect_ID), ]

mkV <- function(sub) vcalc(vi = vi, cluster = study, grp1 = experiment_group_id,
                           grp2 = control_group_id, w1 = n1, w2 = n2,
                           data = sub, sparse = FALSE)
addcols <- function(x) {
  x$age_between <- ave(x$age, x$study, FUN = mean)
  x$age_within  <- x$age - x$age_between
  x$age_between_c <- x$age_between - mean(x$age_between)
  x$age_within2 <- x$age_within^2
  x$age_c2 <- (x$age - mean(x$age))^2
  x
}

sink(file.path(OUT, "console.txt"), split = TRUE)
cat("A09 - age\n\n")

cat("=== 1. where the age information is ===\n")
na <- tapply(d$age, d$study, function(x) length(unique(x)))
cat(sprintf("studies with >1 age value: %d of %d; effects in them: %d of %d\n",
            sum(na > 1), length(na), sum(d$study %in% names(na)[na > 1]), nrow(d)))
g <- ave(d$age, d$study, FUN = mean)
w <- sum((d$age - g)^2); b <- sum((g - mean(d$age))^2)
cat(sprintf("age variation: within-study %.1f%% / between-study %.1f%%\n",
            100 * w / (w + b), 100 * b / (w + b)))
cat("age_basis:\n"); print(table(d$age_basis))
cat(sprintf("age range %.1f-%.1f years\n", min(d$age), max(d$age)))

cat("\n=== 2. decomposed linear age (the A04 correction applied to age) ===\n")
rows <- list()
for (nm in c("ALL", "IS", "RT", "PE")) {
  dd <- addcols(if (nm == "ALL") d else d[d$outcome == nm, ])
  if (length(unique(dd$age_within)) < 2) {
    cat(sprintf("%-4s no within-study age variation at all\n", nm)); next }
  m <- try(rma.mv(yi_al, mkV(dd), mods = ~ age_within + age_between_c, data = dd,
                  random = list(~ 1 | study, ~ 1 | effect_ID),
                  method = "REML", test = "t", dfs = "contain"), silent = TRUE)
  if (inherits(m, "try-error")) { cat(sprintf("%-4s failed\n", nm)); next }
  ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = dd$study)
  cat(sprintf("%-4s within %+.4f (CR2 p %.4f, df %.1f) | between %+.4f (CR2 p %.4f, df %.1f)\n",
              nm, m$b[2], ct$p_Satt[2], ct$df_Satt[2], m$b[3], ct$p_Satt[3], ct$df_Satt[3]))
  for (i in 2:3) rows[[length(rows) + 1]] <- data.frame(
    family = nm, term = rownames(m$b)[i], estimate = m$b[i],
    ci_lo = m$ci.lb[i], ci_hi = m$ci.ub[i],
    p_model = m$pval[i], p_CR2 = ct$p_Satt[i], df_CR2 = ct$df_Satt[i])
}
dec <- do.call(rbind, rows)
write.csv(dec, file.path(OUT, "age_decomposition.csv"), row.names = FALSE)

extra <- list()
cat("\n=== 3. is the age effect non-linear? ===\n")
for (nm in c("ALL", "IS")) {
  dd <- addcols(if (nm == "ALL") d else d[d$outcome == nm, ])
  m <- try(rma.mv(yi_al, mkV(dd), mods = ~ age_within + age_within2 + age_between_c,
                  data = dd, random = list(~ 1 | study, ~ 1 | effect_ID),
                  method = "REML", test = "t", dfs = "contain"), silent = TRUE)
  if (inherits(m, "try-error")) { cat(sprintf("  %-4s failed\n", nm)); next }
  ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = dd$study)
  i <- which(rownames(m$b) == "age_within2")
  cat(sprintf("  %-4s quadratic (within) %+.5f  p_model %.4g  CR2 p %.4f (df %.1f)\n",
              nm, m$b[i], m$pval[i], ct$p_Satt[i], ct$df_Satt[i]))
  extra[[length(extra) + 1]] <- data.frame(test = paste("age^2 within", nm),
                                           estimate = m$b[i], p_CR2 = ct$p_Satt[i])
}
cat("  and as a coarse age band, which needs no functional form:\n")
d$band <- cut(d$age, c(0, 8, 10, 99), labels = c("6-8", "8-10", "10-13"), right = FALSE)
print(table(d$band, d$outcome))
for (nm in c("ALL", "IS")) {
  dd <- if (nm == "ALL") d else d[d$outcome == nm, ]
  m <- try(rma.mv(yi_al, mkV(dd), mods = ~ band, data = dd,
                  random = list(~ 1 | study, ~ 1 | effect_ID),
                  method = "REML", test = "t", dfs = "contain"), silent = TRUE)
  if (inherits(m, "try-error")) next
  ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = dd$study)
  cat(sprintf("  %-4s band 8-10 vs 6-8 %+.3f (CR2 p %.3f) | 10-13 vs 6-8 %+.3f (CR2 p %.3f)\n",
              nm, m$b[2], ct$p_Satt[2], m$b[3], ct$p_Satt[3]))
  extra[[length(extra) + 1]] <- data.frame(test = paste(c("band 8-10", "band 10-13"), nm),
                                           estimate = m$b[2:3], p_CR2 = ct$p_Satt[2:3])
}

cat("\n=== 4. does the null survive restricting to reported ages? ===\n")
for (lab in list(c("all", ""), c("reported only", "reported"))) {
  dd <- if (lab[2] == "") d else d[d$age_basis == "reported", ]
  dd <- addcols(dd)
  m <- try(rma.mv(yi_al, mkV(dd), mods = ~ age_within + age_between_c, data = dd,
                  random = list(~ 1 | study, ~ 1 | effect_ID),
                  method = "REML", test = "t", dfs = "contain"), silent = TRUE)
  if (inherits(m, "try-error")) { cat(sprintf("  %-14s failed\n", lab[1])); next }
  ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = dd$study)
  cat(sprintf("  %-14s k=%3d studies=%2d | within %+.4f (CR2 p %.3f) | between %+.4f (CR2 p %.3f)\n",
              lab[1], nrow(dd), length(unique(dd$study)),
              m$b[2], ct$p_Satt[2], m$b[3], ct$p_Satt[3]))
}

cat("\n=== 5. does age moderate the SNR gradient? (decomposed) ===\n")
dd <- addcols(d)
dd$SNR_between <- ave(dd$SNR, dd$study, FUN = mean)
dd$SNR_within  <- dd$SNR - dd$SNR_between
m <- try(rma.mv(yi_al, mkV(dd), mods = ~ SNR_within * age_within + SNR_between + age_between_c,
                data = dd, random = list(~ 1 | study, ~ 1 | effect_ID),
                method = "REML", test = "t", dfs = "contain"), silent = TRUE)
if (!inherits(m, "try-error")) {
  ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = dd$study)
  i <- grep(":", rownames(m$b))
  cat(sprintf("  SNR_within x age_within %+.5f  p_model %.4g  CR2 p %.4f (df %.1f)\n",
              m$b[i], m$pval[i], ct$p_Satt[i], ct$df_Satt[i]))
  extra[[length(extra) + 1]] <- data.frame(test = "SNR_within x age_within",
                                           estimate = m$b[i], p_CR2 = ct$p_Satt[i])
} else cat("  interaction model failed\n")

write.csv(do.call(rbind, extra), file.path(OUT, "extra_tests.csv"), row.names = FALSE)
cat("\nwritten to", OUT, "\n")
sink()
