
source(file.path("R", "helpers.R"))
OUT <- file.path(p26_outdir("analysis_22_scope_sensitivity"), "outcome")

D <- p26_data()
SCOPES <- c("base", "condition_subset", "no_composite_PE", "both")

fit_one <- function(d, family, mods, idx) {
  if (family != "ALL") d <- d[d$outcome == family, ]
  d <- d[order(d$study, d$effect_ID), ]
  d$SNR_between <- ave(d$SNR, d$study, FUN = mean)
  d$SNR_within  <- d$SNR - d$SNR_between
  d$SNR_between_c <- d$SNR_between - mean(d$SNR_between)
  V <- vcalc(vi = vi, cluster = study, grp1 = experiment_group_id,
             grp2 = control_group_id, w1 = n1, w2 = n2, data = d, sparse = FALSE)
  m <- rma.mv(yi_al, V, mods = mods, data = d,
              random = list(~ 1 | study, ~ 1 | effect_ID),
              method = "REML", test = "t", dfs = "contain")
  ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = d$study)
  ci <- clubSandwich::conf_int(m, vcov = "CR2", cluster = d$study)
  data.frame(estimate = as.numeric(m$b)[idx], ci_lo = ci$CI_L[idx], ci_hi = ci$CI_U[idx],
             p_CR2 = ct$p_Satt[idx], df_CR2 = ct$df_Satt[idx],
             k = m$k, studies = length(unique(d$study)))
}

QUANT <- list(
  list(name = "pooled ALL",     family = "ALL", mods = ~ 1, idx = 1),
  list(name = "pooled IS",      family = "IS",  mods = ~ 1, idx = 1),
  list(name = "pooled RT",      family = "RT",  mods = ~ 1, idx = 1),
  list(name = "pooled PE",      family = "PE",  mods = ~ 1, idx = 1),
  list(name = "SNR within ALL", family = "ALL", mods = ~ SNR_within + SNR_between_c, idx = 2),
  list(name = "SNR within IS",  family = "IS",  mods = ~ SNR_within + SNR_between_c, idx = 2))

sink(file.path(OUT, "console.txt"), split = TRUE)
cat("A22 - scope sensitivity: base vs narrower condition set vs composite PE removed\n\n")

rows <- list()
for (s in SCOPES) {
  x <- p26_scope(D, s)
  cat(sprintf("scope %-17s effects %3d (IS %d / RT %d / PE %d) | studies %d | pairs %d\n", s,
              nrow(x$effects), sum(x$effects$outcome == "IS"), sum(x$effects$outcome == "RT"),
              sum(x$effects$outcome == "PE"), length(unique(x$effects$study)), nrow(x$pairs)))
  for (q in QUANT) {
    r <- fit_one(x$effects, q$family, q$mods, q$idx)
    rows[[length(rows) + 1]] <- cbind(scope = s, quantity = q$name, r)
  }
  p <- x$pairs
  rr <- cor(p$g_IS, p$g_RT)
  rows[[length(rows) + 1]] <- data.frame(
    scope = s, quantity = "IS-RT pair r", estimate = rr, ci_lo = NA, ci_hi = NA,
    p_CR2 = NA, df_CR2 = NA, k = nrow(p), studies = length(unique(p$study)))
  rows[[length(rows) + 1]] <- data.frame(
    scope = s, quantity = "pairs adverse on both", estimate = sum(p$g_IS > 0 & p$g_RT > 0),
    ci_lo = NA, ci_hi = NA, p_CR2 = NA, df_CR2 = NA, k = nrow(p),
    studies = length(unique(p$study)))
}
res <- do.call(rbind, rows)
rownames(res) <- NULL

cat("\n=== estimates by scope (shared-arm V, CR2) ===\n")
print(res, row.names = FALSE, digits = 4)

base <- res[res$scope == "base", c("quantity", "estimate")]
cat("\n=== largest change from the base, by quantity ===\n")
for (q in unique(res$quantity)) {
  s <- res[res$quantity == q & res$scope != "base", ]
  b <- base$estimate[base$quantity == q]
  sig <- if (all(is.na(s$p_CR2))) "" else
    sprintf(" | significant in %d/3 alternatives", sum(s$p_CR2 < .05))
  cat(sprintf("%-22s base %+.3f | change %+.3f to %+.3f%s\n", q, b,
              min(s$estimate - b), max(s$estimate - b), sig))
}

write.csv(res, file.path(OUT, "scope_sensitivity.csv"), row.names = FALSE)
cat("\nwritten to", OUT, "\n")
sink()
