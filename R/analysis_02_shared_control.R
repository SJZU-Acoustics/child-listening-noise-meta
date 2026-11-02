
source(file.path("R", "helpers.R"))
OUT <- file.path(p26_outdir("analysis_02_shared_control"), "outcome")

D <- p26_data()
d <- D$effects
d <- d[order(d$study, d$effect_ID), ]

sink(file.path(OUT, "console.txt"), split = TRUE)
cat("A02 - shared-control structure\n\n")

cat("=== 1. how much sharing is there ===\n")
d$arm_key <- paste(d$study, d$control_group_id)
tab <- table(d$arm_key)
d$arm_reuse <- as.integer(tab[d$arm_key])
cat(sprintf("distinct control arms: %d over %d effects\n", length(tab), nrow(d)))
cat(sprintf("effects on a shared control arm: %d (%.1f%%)\n",
            sum(d$arm_reuse > 1), 100 * mean(d$arm_reuse > 1)))
cat(sprintf("studies with any sharing: %d of %d\n",
            length(unique(d$study[d$arm_reuse > 1])), length(unique(d$study))))
cat("reuse distribution (effects per control arm):\n")
print(table(as.integer(tab)))
shr <- aggregate(arm_reuse ~ study, d, function(x) c(n = length(x), max = max(x)))
shr <- data.frame(study = shr$study, n_effects = shr$arm_reuse[, "n"],
                  max_reuse = shr$arm_reuse[, "max"])
write.csv(shr[order(-shr$max_reuse), ], file.path(OUT, "sharing_by_study.csv"), row.names = FALSE)

cat("\n=== 2. building the block-diagonal V ===\n")
V_shared <- vcalc(vi = vi, cluster = study, grp1 = experiment_group_id,
                  grp2 = control_group_id, w1 = n1, w2 = n2, data = d,
                  checkpd = TRUE, sparse = FALSE)
V_diag <- diag(d$vi)
V_rho5 <- vcalc(vi = vi, cluster = study, obs = effect_ID, rho = 0.5, data = d)

od   <- row(V_shared) != col(V_shared)
nz   <- od & (V_shared != 0)
Rmat <- cov2cor(V_shared)
cat(sprintf("V_shared: %d x %d | diagonal equals vi: %s | non-zero off-diagonals: %d\n",
            nrow(V_shared), ncol(V_shared),
            isTRUE(all.equal(diag(V_shared), d$vi, check.attributes = FALSE)),
            sum(nz)))
cat(sprintf("  implied sampling correlations: median %.3f, range %.3f to %.3f\n",
            median(Rmat[nz]), min(Rmat[nz]), max(Rmat[nz])))
stopifnot(all(Rmat[nz] > 0), all(Rmat[nz] < 1))
ev <- eigen(V_shared, symmetric = TRUE, only.values = TRUE)$values
cat(sprintf("  positive definite: %s (min eigenvalue %.3g)\n", all(ev > 0), min(ev)))

lone <- which(d$arm_reuse == 1)
cat(sprintf("  effects on a non-shared arm with zero off-diagonal: %d of %d\n",
            sum(sapply(lone, function(i) all(V_shared[i, -i] == 0))), length(lone)))

cat("\n=== 3. headline models under each working model ===\n")
fit_one <- function(V, sub, mods = ~ 1, idx = 1) {
  m <- rma.mv(yi_al, V, mods = mods, data = sub,
              random = list(~ 1 | study, ~ 1 | effect_ID),
              method = "REML", test = "t", dfs = "contain")
  rb <- clubSandwich::coef_test(m, vcov = "CR2", cluster = sub$study)
  list(est = m$b[idx], lo = m$ci.lb[idx], hi = m$ci.ub[idx],
       p_model = m$pval[idx], se_CR2 = rb$SE[idx], df = rb$df_Satt[idx],
       p_CR2 = rb$p_Satt[idx], s2 = m$sigma2)
}
subV <- function(V, keep) V[keep, keep, drop = FALSE]

rows <- list()
for (nm in c("ALL", "IS", "RT", "PE")) {
  keep <- if (nm == "ALL") rep(TRUE, nrow(d)) else d$outcome == nm
  dd <- d[keep, ]
  for (wm in c("diagonal", "shared-arm", "rho=0.5")) {
    V <- switch(wm, "diagonal" = subV(V_diag, keep),
                "shared-arm" = subV(V_shared, keep), "rho=0.5" = subV(V_rho5, keep))
    r <- fit_one(V, dd)
    rows[[length(rows) + 1]] <- data.frame(
      quantity = paste("empty", nm), working_model = wm, estimate = r$est,
      ci_lo = r$lo, ci_hi = r$hi, p_model = r$p_model,
      se_CR2 = r$se_CR2, df_CR2 = r$df, p_CR2 = r$p_CR2)
  }
}
for (nm in c("IS", "RT", "PE")) {
  keep <- d$outcome == nm
  dd <- d[keep, ]
  for (wm in c("diagonal", "shared-arm", "rho=0.5")) {
    V <- switch(wm, "diagonal" = subV(V_diag, keep),
                "shared-arm" = subV(V_shared, keep), "rho=0.5" = subV(V_rho5, keep))
    r <- fit_one(V, dd, ~ SNR_c + age_c + masking, idx = 2)
    rows[[length(rows) + 1]] <- data.frame(
      quantity = paste("SNR slope", nm), working_model = wm, estimate = r$est,
      ci_lo = r$lo, ci_hi = r$hi, p_model = r$p_model,
      se_CR2 = r$se_CR2, df_CR2 = r$df, p_CR2 = r$p_CR2)
  }
}
res <- do.call(rbind, rows)
for (q in unique(res$quantity)) {
  cat("\n", q, "\n", sep = "")
  s <- res[res$quantity == q, ]
  for (i in seq_len(nrow(s))) {
    cat(sprintf("  %-11s est %+.4f [%+.4f, %+.4f]  p_model %-9.3g CR2 p %-8.4f df %.1f\n",
                s$working_model[i], s$estimate[i], s$ci_lo[i], s$ci_hi[i],
                s$p_model[i], s$p_CR2[i], s$df_CR2[i]))
  }
}

cat("\n=== 4. did modelling V buy df or change a verdict? ===\n")
wide <- reshape(res[, c("quantity", "working_model", "estimate", "df_CR2", "p_CR2")],
                idvar = "quantity", timevar = "working_model", direction = "wide")
names(wide) <- gsub("\\.", "_", names(wide))
wide$df_gain <- round(wide$`df_CR2_shared-arm` - wide$df_CR2_diagonal, 2)
wide$verdict_change <- ifelse((wide$p_CR2_diagonal < .05) !=
                                (wide$`p_CR2_shared-arm` < .05), "CHANGED", "")
print(wide[, c("quantity", "estimate_diagonal", "estimate_shared-arm",
               "df_CR2_diagonal", "df_CR2_shared-arm", "df_gain",
               "p_CR2_diagonal", "p_CR2_shared-arm", "verdict_change")],
      row.names = FALSE, digits = 4)

write.csv(res, file.path(OUT, "working_model_comparison.csv"), row.names = FALSE)
write.csv(wide, file.path(OUT, "df_and_verdict_change.csv"), row.names = FALSE)
cat("\nwritten to", OUT, "\n")
sink()
