
source(file.path("R", "helpers.R"))
OUT <- file.path(p26_outdir("analysis_16_log_ratio"), "outcome")

D <- p26_data()
d <- D$effects
d <- d[order(d$study, d$effect_ID), ]
src <- read.csv(file.path(P26_ROOT, ".cache", "data", "source_rows.csv"),
                stringsAsFactors = FALSE)
unit <- setNames(src$unit, src$source_id)
d$ctrl_unit <- unit[d$control_group_id]

sink(file.path(OUT, "console.txt"), split = TRUE)
cat("A16 - log ratio of means: a third scale\n\n")

rom <- escalc(measure = "ROM", m1i = m1, sd1i = sd1, n1i = n1,
              m2i = m2, sd2i = sd2, n2i = n2, data = d)
d$rom    <- ifelse(d$outcome == "IS", -rom$yi, rom$yi)
d$rom_vi <- rom$vi
bad <- is.na(d$rom) | d$m1 <= 0 | d$m2 <= 0
cat(sprintf("ROM undefined or non-positive means: %d rows (excluded)\n", sum(bad)))
d <- d[!bad, ]

mkV_rom <- function(sub) {
  V <- diag(sub$rom_vi)
  key <- paste(sub$study, sub$control_group_id)
  for (arm in unique(key)) {
    idx <- which(key == arm)
    if (length(idx) > 1) {
      cv <- sub$sd2[idx[1]]^2 / (sub$n2[idx[1]] * sub$m2[idx[1]]^2)
      for (a in idx) for (b in idx) if (a < b) V[a, b] <- V[b, a] <- cv
    }
  }
  V
}
V_all <- mkV_rom(d)
cat(sprintf("V validation: max |diag(V) - escalc vi| = %.3g ; PD: %s\n",
            max(abs(diag(V_all) - d$rom_vi)),
            min(eigen(V_all, symmetric = TRUE, only.values = TRUE)$values) > 0))

fit_all <- function(sub, label) {
  V <- mkV_rom(sub)
  sub$SNR_between <- ave(sub$SNR, sub$study, FUN = mean)
  sub$SNR_within  <- sub$SNR - sub$SNR_between
  sub$SNR_between_c <- sub$SNR_between - mean(sub$SNR_between)
  m0 <- rma.mv(rom, V, data = sub, random = list(~ 1 | study, ~ 1 | effect_ID),
               method = "REML", test = "t", dfs = "contain")
  ct0 <- clubSandwich::coef_test(m0, vcov = "CR2", cluster = sub$study)
  m1 <- rma.mv(rom, V, mods = ~ SNR_within + SNR_between_c, data = sub,
               random = list(~ 1 | study, ~ 1 | effect_ID),
               method = "REML", test = "t", dfs = "contain")
  ct1 <- clubSandwich::coef_test(m1, vcov = "CR2", cluster = sub$study)
  pct <- if (label %in% c("IS")) 100 * (1 - exp(-m0$b[1])) else NA
  cat(sprintf("\n%-3s pooled ROM %+.4f [%.4f, %.4f]  CR2 p %.4g (df %.1f)%s\n",
              label, m0$b[1], m0$ci.lb[1], m0$ci.ub[1], ct0$p_Satt[1],
              ct0$df_Satt[1],
              ifelse(is.na(pct), "", sprintf("  = %.1f%% of quiet performance lost", pct))))
  for (i in 2:3)
    cat(sprintf("    %-16s %+.5f  CR2 p %.4g (df %.1f)\n",
                rownames(m1$b)[i], m1$b[i], ct1$p_Satt[i], ct1$df_Satt[i]))
  data.frame(family = label,
             pooled_est = m0$b[1], pooled_lo = m0$ci.lb[1],
             pooled_hi = m0$ci.ub[1], pooled_p = ct0$p_Satt[1],
             pooled_df = ct0$df_Satt[1],
             snr_within_est = m1$b[2], snr_within_p = ct1$p_Satt[2],
             snr_within_df = ct1$df_Satt[2],
             snr_between_est = m1$b[3], snr_between_p = ct1$p_Satt[3],
             snr_between_df = ct1$df_Satt[3])
}

cat("=== 1. ROM per family (RT all rows; IS percent rows only) ===\n")
res <- rbind(
  fit_all(d[d$outcome == "RT", ], "RT"),
  fit_all(d[d$outcome == "IS" & d$ctrl_unit != "RAU", ], "IS"),
  fit_all(d, "ALL"))
cat(sprintf("\nIS RAU rows excluded: %d\n",
            sum(d$outcome == "IS" & d$ctrl_unit == "RAU")))

cat("\n=== 2. three-scale concordance (CR2 p values) ===\n")
cat("quantity                         SMD(A02/A04)  raw(A11)     ROM(A16)\n")
OC <- function(a, f) read.csv(file.path(P26_ROOT, ".cache", "analysis", a, "outcome", f))
wm <- OC("analysis_02_shared_control", "working_model_comparison.csv")
wm <- wm[wm$working_model == "shared-arm", ]
wb <- OC("analysis_04_snr_information", "within_between_decomposition.csv")
wb <- wb[wb$term == "SNR_within", ]
rh <- OC("analysis_11_raw_scale", "raw_scale_headlines.csv")
conc <- data.frame(
  quantity = c("pooled IS", "pooled RT", "IS SNR within", "RT SNR within"),
  SMD = sprintf("%.3g", c(wm$p_CR2[wm$quantity == "empty IS"], wm$p_CR2[wm$quantity == "empty RT"],
                          wb$p_CR2[wb$family == "IS"], wb$p_CR2[wb$family == "RT"])),
  raw = sprintf("%.3g", c(rh$p_CR2[1], rh$p_CR2[2], rh$snr_within_p[1], rh$snr_within_p[2])),
  ROM = c(sprintf("%.3g", res$pooled_p[res$family == "IS"]),
          sprintf("%.3g", res$pooled_p[res$family == "RT"]),
          sprintf("%.3g", res$snr_within_p[res$family == "IS"]),
          sprintf("%.3g", res$snr_within_p[res$family == "RT"])))
print(conc, row.names = FALSE)

write.csv(res, file.path(OUT, "rom_results.csv"), row.names = FALSE)
write.csv(conc, file.path(OUT, "three_scale_concordance.csv"), row.names = FALSE)
write.csv(d[, c("effect_ID", "study", "outcome", "m1", "m2", "rom", "rom_vi")],
          file.path(OUT, "rom_effect_data.csv"), row.names = FALSE)
cat("\nwritten to", OUT, "\n")
sink()
