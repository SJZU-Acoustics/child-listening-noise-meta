
source(file.path("R", "helpers.R"))
OUT <- file.path(p26_outdir("analysis_13_specification_curve"), "outcome")

D <- p26_data()
base <- D$effects
base <- base[order(base$study, base$effect_ID), ]

sink(file.path(OUT, "console.txt"), split = TRUE)
cat("A13 - exit specification curve\n\n")
cat("48 specifications = 3 working models x 4 variance assumptions x\n")
cat("                    Holt in/out x study 13 merged/split\n\n")

specs <- expand.grid(
  working = c("shared-arm", "diagonal", "rho=0.5"),
  variance = c("delivered", "paired r=0", "paired r=0.5", "paired r=0.9"),
  holt = c("kept", "dropped"),
  study13 = c("merged", "split"),
  stringsAsFactors = FALSE)
cat(sprintf("%d specifications\n\n", nrow(specs)))

build <- function(sp, family) {
  d <- base
  if (sp$study13 == "split") {
    i <- d$report_source == "25.pdf"
    d$study[i] <- paste0(d$study[i], " (25.pdf)")
  }
  if (sp$holt == "dropped") d <- d[d$design == "within-subject", ]
  if (family != "ALL") d <- d[d$outcome == family, ]
  if (nrow(d) < 10 || length(unique(d$study)) < 4) return(NULL)
  if (sp$variance != "delivered") {
    r <- as.numeric(sub("paired r=", "", sp$variance))
    d$vi <- p26_within_variance(d, r)
  }
  d$SNR_between <- ave(d$SNR, d$study, FUN = mean)
  d$SNR_within  <- d$SNR - d$SNR_between
  d$SNR_between_c <- d$SNR_between - mean(d$SNR_between)
  V <- switch(sp$working,
    "diagonal"   = diag(d$vi),
    "shared-arm" = vcalc(vi = vi, cluster = study, grp1 = experiment_group_id,
                         grp2 = control_group_id, w1 = n1, w2 = n2, data = d, sparse = FALSE),
    "rho=0.5"    = vcalc(vi = vi, cluster = study, obs = effect_ID, rho = 0.5, data = d))
  list(d = d, V = V)
}

run_one <- function(sp, family, mods, idx) {
  b <- build(sp, family)
  if (is.null(b)) return(NULL)
  m <- try(rma.mv(yi_al, b$V, mods = mods, data = b$d,
                  random = list(~ 1 | study, ~ 1 | effect_ID),
                  method = "REML", test = "t", dfs = "contain"), silent = TRUE)
  if (inherits(m, "try-error")) return(NULL)
  ct <- try(clubSandwich::coef_test(m, vcov = "CR2", cluster = b$d$study), silent = TRUE)
  if (inherits(ct, "try-error")) return(NULL)
  if (idx > length(m$b)) return(NULL)
  data.frame(estimate = as.numeric(m$b)[idx], p_model = m$pval[idx],
             p_CR2 = ct$p_Satt[idx], df_CR2 = ct$df_Satt[idx],
             k = m$k, clusters = length(unique(b$d$study)))
}

QUANT <- list(
  list(name = "pooled ALL",       family = "ALL", mods = ~ 1,                                idx = 1),
  list(name = "pooled IS",        family = "IS",  mods = ~ 1,                                idx = 1),
  list(name = "pooled RT",        family = "RT",  mods = ~ 1,                                idx = 1),
  list(name = "pooled PE",        family = "PE",  mods = ~ 1,                                idx = 1),
  list(name = "SNR within ALL",   family = "ALL", mods = ~ SNR_within + SNR_between_c,       idx = 2),
  list(name = "SNR within IS",    family = "IS",  mods = ~ SNR_within + SNR_between_c,       idx = 2))

all_rows <- list()
for (q in QUANT) {
  cat(sprintf("running %-16s ", q$name)); flush.console()
  for (i in seq_len(nrow(specs))) {
    r <- run_one(specs[i, ], q$family, q$mods, q$idx)
    if (is.null(r)) next
    all_rows[[length(all_rows) + 1]] <- cbind(quantity = q$name, specs[i, ], r)
  }
  cat("done\n")
}
res <- do.call(rbind, all_rows)
rownames(res) <- NULL
write.csv(res, file.path(OUT, "specification_curve.csv"), row.names = FALSE)

cat("\n=== survival across the 48 specifications ===\n")
summ <- list()
for (q in unique(res$quantity)) {
  s <- res[res$quantity == q, ]
  sig <- sum(s$p_CR2 < .05)
  same_sign <- all(sign(s$estimate) == sign(s$estimate[1]))
  cat(sprintf("%-16s %2d/%2d fitted | significant %2d/%2d | est %+.4f to %+.4f | sign stable %s | df %.1f-%.1f\n",
              q, nrow(s), nrow(specs), sig, nrow(s),
              min(s$estimate), max(s$estimate), same_sign,
              min(s$df_CR2), max(s$df_CR2)))
  summ[[length(summ) + 1]] <- data.frame(
    quantity = q, fitted = nrow(s), significant = sig,
    est_min = min(s$estimate), est_max = max(s$estimate),
    median_est = median(s$estimate), sign_stable = same_sign,
    p_min = min(s$p_CR2), p_max = max(s$p_CR2))
}
summ <- do.call(rbind, summ)
write.csv(summ, file.path(OUT, "specification_summary.csv"), row.names = FALSE)

cat("\n=== which discretionary choice moves the estimate most? ===\n")
for (q in unique(res$quantity)) {
  s <- res[res$quantity == q, ]
  cat(sprintf("\n%s (median %+.4f)\n", q, median(s$estimate)))
  for (fac in c("working", "variance", "holt", "study13")) {
    a <- aggregate(s$estimate, list(s[[fac]]), median)
    rng <- diff(range(a$x))
    cat(sprintf("  %-9s range of medians %.4f   (%s)\n", fac, rng,
                paste(sprintf("%s %+.3f", a$Group.1, a$x), collapse = ", ")))
  }
}

cat("\n=== worst case per quantity (largest CR2 p) ===\n")
for (q in unique(res$quantity)) {
  s <- res[res$quantity == q, ]
  w <- s[which.max(s$p_CR2), ]
  cat(sprintf("%-16s p=%.4f under %s / %s / Holt %s / study13 %s\n",
              q, w$p_CR2, w$working, w$variance, w$holt, w$study13))
}

cat("\nwritten to", OUT, "\n")
sink()
