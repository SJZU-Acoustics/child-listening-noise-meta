
source(file.path("R", "helpers.R"))
OUT <- file.path(p26_outdir("analysis_03_r_table"), "outcome")

D <- p26_data()
d <- D$effects
p <- D$pairs
d <- d[order(d$study, d$effect_ID), ]

sink(file.path(OUT, "console.txt"), split = TRUE)
cat("A03 - paired-variance r-table across every headline\n\n")

RS <- c(0, 0.3, 0.5, 0.7, 0.9)

cat("=== 0. how far the delivered variances already are from independence ===\n")
d$v_indep <- d$vi
for (r in RS) {
  vp <- p26_within_variance(d, r)
  cat(sprintf("  r = %.1f : mean vi ratio (paired / delivered) = %.3f\n",
              r, mean(vp / d$v_indep)))
}
cat("  (a ratio of 1 marks the r the delivered variances behave like)\n")
rr <- sapply(seq(0, 0.95, 0.01),
             function(r) abs(mean(p26_within_variance(d, r) / d$v_indep) - 1))
cat(sprintf("  implied r = %.2f\n", seq(0, 0.95, 0.01)[which.min(rr)]))

cat("\n=== 1. every headline across the r-table (shared-arm V + CR2) ===\n")

build_V <- function(sub, vvec) {
  s <- sub; s$vi <- vvec
  vcalc(vi = vi, cluster = study, grp1 = experiment_group_id,
        grp2 = control_group_id, w1 = n1, w2 = n2, data = s, sparse = FALSE)
}
fit <- function(sub, vvec, mods = ~ 1, idx = 1) {
  V <- build_V(sub, vvec)
  m <- rma.mv(yi_al, V, mods = mods, data = sub,
              random = list(~ 1 | study, ~ 1 | effect_ID),
              method = "REML", test = "t", dfs = "contain")
  ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = sub$study)
  c(est = m$b[idx], p_model = m$pval[idx], p_CR2 = ct$p_Satt[idx], df = ct$df_Satt[idx])
}

rows <- list()
for (nm in c("ALL", "IS", "RT", "PE")) {
  dd <- if (nm == "ALL") d else d[d$outcome == nm, ]
  r0 <- fit(dd, dd$v_indep)
  rows[[length(rows) + 1]] <- data.frame(quantity = paste("empty", nm), r = "delivered",
                                         t(r0))
  for (r in RS) {
    rr2 <- fit(dd, p26_within_variance(dd, r))
    rows[[length(rows) + 1]] <- data.frame(quantity = paste("empty", nm),
                                           r = sprintf("%.1f", r), t(rr2))
  }
}
for (nm in c("IS", "RT", "PE")) {
  dd <- d[d$outcome == nm, ]
  r0 <- fit(dd, dd$v_indep, ~ SNR_c + age_c + masking, idx = 2)
  rows[[length(rows) + 1]] <- data.frame(quantity = paste("SNR slope", nm),
                                         r = "delivered", t(r0))
  for (r in RS) {
    rr2 <- fit(dd, p26_within_variance(dd, r), ~ SNR_c + age_c + masking, idx = 2)
    rows[[length(rows) + 1]] <- data.frame(quantity = paste("SNR slope", nm),
                                           r = sprintf("%.1f", r), t(rr2))
  }
}
res <- do.call(rbind, rows)
names(res) <- c("quantity", "r", "estimate", "p_model", "p_CR2", "df_CR2")

for (q in unique(res$quantity)) {
  s <- res[res$quantity == q, ]
  cat("\n", q, "\n", sep = "")
  for (i in seq_len(nrow(s)))
    cat(sprintf("  r=%-9s est %+.4f  p_model %-10.3g CR2 p %-8.4f df %.1f%s\n",
                s$r[i], s$estimate[i], s$p_model[i], s$p_CR2[i], s$df_CR2[i],
                if (s$p_CR2[i] < .05) "" else "   <- not significant"))
}

write.csv(res, file.path(OUT, "r_table_headlines.csv"), row.names = FALSE)
sink()
