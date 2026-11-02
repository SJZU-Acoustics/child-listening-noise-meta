
suppressMessages({
  library(metafor)
  library(clubSandwich)
  library(wildmeta)
})

P26_ROOT <- getwd()

BOOT_R <- 1000  # house default

p26_data <- function(check_freeze = TRUE, base = "freeze") {
  if (check_freeze) p26_assert_freeze()
  src <- file.path(P26_ROOT, ".cache", "data")

  d <- read.csv(file.path(src, "analysis_table.csv"), stringsAsFactors = FALSE)
  p <- read.csv(file.path(src, "aligned_pairs.csv"), stringsAsFactors = FALSE)

  d <- escalc(measure = "SMD", m1i = m1, sd1i = sd1, n1i = n1,
              m2i = m2, sd2i = sd2, n2i = n2, data = d)

  stopifnot(all((d$direction_raw == "higher_better") == (d$outcome == "IS")))
  d$yi_al <- d$yi * ifelse(d$outcome == "IS", -1, 1)

  d$SNR_c <- d$SNR - mean(d$SNR)
  d$age_c <- d$age - mean(d$age)
  d$masking <- factor(d$noise_masking_type,
                      levels = c("Energetic Masking", "Informational Masking"),
                      labels = c("EM", "IM"))
  d$Outcome_Type <- factor(d$outcome, levels = c("IS", "RT", "PE"))

  gI <- escalc(measure = "SMD", m1i = IS_m1, sd1i = IS_sd1, n1i = IS_n1,
               m2i = IS_m2, sd2i = IS_sd2, n2i = IS_n2, data = p)
  gR <- escalc(measure = "SMD", m1i = RT_m1, sd1i = RT_sd1, n1i = RT_n1,
               m2i = RT_m2, sd2i = RT_sd2, n2i = RT_n2, data = p)
  p$g_IS   <- -gI$yi                      # aligned: positive = more adverse
  p$g_RT   <-  gR$yi
  p$diff   <- p$g_IS - p$g_RT             # ISRT_diff (RQ2 outcome)
  p$v_diff <- gI$vi + gR$vi
  p$SNR_c  <- p$SNR - mean(p$SNR)
  p$task_group <- factor(p$task_group, levels = c("same-task", "dual-task"))

  list(effects = d, pairs = p)
}

p26_scope <- function(dat, scope = c("base", "condition_subset", "no_composite_PE", "both")) {
  scope <- match.arg(scope)
  e <- dat$effects
  keep <- rep(TRUE, nrow(e))
  if (scope %in% c("condition_subset", "both")) keep <- keep & e$in_condition_subset
  if (scope %in% c("no_composite_PE", "both")) keep <- keep & !e$composite_PE
  e <- e[keep, ]
  p <- dat$pairs[dat$pairs$effect_ID_IS %in% e$effect_ID &
                 dat$pairs$effect_ID_RT %in% e$effect_ID, ]
  list(effects = e, pairs = p)
}

p26_assert_freeze <- function() {
  d <- read.csv(file.path(P26_ROOT, ".cache", "data", "analysis_table.csv"))
  p <- read.csv(file.path(P26_ROOT, ".cache", "data", "aligned_pairs.csv"))
  stopifnot(nrow(d) == 269L, length(unique(d$study)) == 24L,
            nrow(p) == 47L, sum(d$outcome == "IS") == 166L,
            sum(d$outcome == "RT") == 67L, sum(d$outcome == "PE") == 36L)
  invisible(TRUE)
}

p26_fit <- function(formula = ~ 1, data, ...) {
  rma.mv(yi_al, vi, mods = formula, data = data,
         random = list(~ 1 | study, ~ 1 | effect_ID),
         method = "REML", test = "t", dfs = "contain", ...)
}

p26_robust <- function(fit, cluster = NULL) {
  if (is.null(cluster)) cluster <- fit$data$study
  ct <- clubSandwich::coef_test(fit, vcov = "CR2", cluster = cluster)
  mb <- data.frame(coef = rownames(fit$b), est = as.numeric(fit$b),
                   se_model = fit$se, p_model = fit$pval)
  cbind(mb, se_CR2 = ct$SE, df_CR2 = ct$df_Satt, p_CR2 = ct$p_Satt)
}

p26_report <- function(fit, label = "") {
  rb <- p26_robust(fit)
  cat(if (nzchar(label)) paste0("\n", label, "\n") else "\n")
  print(round(rb[, c("est", "se_model", "p_model", "se_CR2", "df_CR2", "p_CR2")], 4))
  cat(sprintf("  sigma^2: study %.4f | effect %.4f   k=%d, clusters=%d\n",
              fit$sigma2[1], fit$sigma2[2], fit$k, length(unique(fit$data$study))))
  invisible(rb)
}

p26_disagree <- function(rb) {
  d <- rb[(rb$p_model < .05) != (rb$p_CR2 < .05), , drop = FALSE]
  if (nrow(d)) d[, c("coef", "p_model", "p_CR2")] else NULL
}

p26_outdir <- function(analysis) {
  for (s in "outcome") {
    dir.create(file.path(P26_ROOT, ".cache", "analysis", analysis, s),
               recursive = TRUE, showWarnings = FALSE)
  }
  file.path(P26_ROOT, ".cache", "analysis", analysis)
}

p26_within_variance <- function(data, r) {
  stopifnot(length(r) == 1L, is.finite(r), r >= 0, r < 1,
            all(data$design %in% c("within-subject", "between-subjects")))
  v <- data$vi
  within <- data$design == "within-subject"
  stopifnot(all(data$n1[within] == data$n2[within]))
  n <- data$n1[within]
  v[within] <- 2 * (1 - r) / n + data$yi[within]^2 / (2 * n)
  stopifnot(identical(v[!within], data$vi[!within]))
  v
}
