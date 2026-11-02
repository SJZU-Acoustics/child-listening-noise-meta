
source(file.path("R", "helpers.R"))
OUT <- file.path(p26_outdir("analysis_14_variability"), "outcome")

D <- p26_data()
d <- D$effects
d <- d[order(d$study, d$effect_ID), ]

sink(file.path(OUT, "console.txt"), split = TRUE)
cat("A14 - dispersion: does noise make children more variable?\n\n")

vr <- escalc(measure = "VR", m1i = m1, sd1i = sd1, n1i = n1,
             m2i = m2, sd2i = sd2, n2i = n2, data = d)
d$lnVR <- vr$yi; d$lnVR_vi <- vr$vi

mkV_vr <- function(sub) {
  k <- nrow(sub)
  V <- diag(sub$lnVR_vi)
  key <- paste(sub$study, sub$control_group_id)
  for (arm in unique(key)) {
    idx <- which(key == arm)
    if (length(idx) > 1) {
      cv <- 1 / (2 * (sub$n2[idx[1]] - 1))   # shared-control covariance
      for (a in idx) for (b in idx) if (a < b) V[a, b] <- V[b, a] <- cv
    }
  }
  V
}
V_all <- mkV_vr(d)
cat(sprintf("V validation: max |diag(V) - escalc vi| = %.3g\n",
            max(abs(diag(V_all) - d$lnVR_vi))))
cat(sprintf("V positive definite: %s (min eigenvalue %.4g)\n",
            min(eigen(V_all, symmetric = TRUE, only.values = TRUE)$values) > 0,
            min(eigen(V_all, symmetric = TRUE, only.values = TRUE)$values)))
shared <- paste(d$study, d$control_group_id) %in%
  names(which(table(paste(d$study, d$control_group_id)) > 1))
cat(sprintf("effects on a shared control arm: %d of %d\n", sum(shared), nrow(d)))

fit_pool <- function(sub, V, label) {
  m <- rma.mv(lnVR, V, data = sub, random = list(~ 1 | study, ~ 1 | effect_ID),
              method = "REML", test = "t", dfs = "contain")
  ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = sub$study)
  cat(sprintf("%-5s lnVR = %+.4f [%.4f, %.4f]  p_model %.4g | CR2 p %.4g (df %.1f) | k=%d, %d studies\n",
              label, m$b[1], m$ci.lb[1], m$ci.ub[1], m$pval[1], ct$p_Satt[1],
              ct$df_Satt[1], m$k, length(unique(sub$study))))
  data.frame(family = label, k = m$k, studies = length(unique(sub$study)),
             est = m$b[1], ci_lo = m$ci.lb[1], ci_hi = m$ci.ub[1],
             p_model = m$pval[1], p_CR2 = ct$p_Satt[1], df_CR2 = ct$df_Satt[1])
}

cat("\n=== 2. pooled lnVR per outcome family (shared-arm V + CR2) ===\n")
cat("positive = more variable in noise\n")
pool <- rbind(
  fit_pool(d, V_all, "ALL"),
  fit_pool(d[d$outcome == "IS", ], mkV_vr(d[d$outcome == "IS", ]), "IS"),
  fit_pool(d[d$outcome == "RT", ], mkV_vr(d[d$outcome == "RT", ]), "RT"),
  fit_pool(d[d$outcome == "PE", ], mkV_vr(d[d$outcome == "PE", ]), "PE"))
print(pool, row.names = FALSE, digits = 4)

cat("\n=== 3. SNR decomposition on lnVR (campaign rule: within/between) ===\n")
d$SNR_between <- ave(d$SNR, d$study, FUN = mean)
d$SNR_within  <- d$SNR - d$SNR_between
d$SNR_between_c <- d$SNR_between - mean(d$SNR_between)
cat(sprintf("within-study share of SNR variation: %.1f%%\n",
            100 * sum(d$SNR_within^2) / (sum(d$SNR_within^2) + sum((d$SNR_between - mean(d$SNR_between))^2))))

fit_snr <- function(sub, label) {
  sub$SNR_between <- ave(sub$SNR, sub$study, FUN = mean)
  sub$SNR_within  <- sub$SNR - sub$SNR_between
  sub$SNR_between_c <- sub$SNR_between - mean(sub$SNR_between)
  m <- rma.mv(lnVR, mkV_vr(sub), mods = ~ SNR_within + SNR_between_c,
              data = sub, random = list(~ 1 | study, ~ 1 | effect_ID),
              method = "REML", test = "t", dfs = "contain")
  ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = sub$study)
  out <- data.frame(family = label, term = rownames(m$b),
                    est = as.numeric(m$b), p_model = m$pval,
                    p_CR2 = ct$p_Satt, df_CR2 = ct$df_Satt)
  print(out[out$term != "intrcpt", ], row.names = FALSE, digits = 4)
  out
}
snr_res <- rbind(fit_snr(d, "ALL"), fit_snr(d[d$outcome == "IS", ], "IS"),
                 fit_snr(d[d$outcome == "RT", ], "RT"))
n_snr_manip_PE <- length(unique(d$study[d$outcome == "PE" &
                       ave(d$SNR, d$study, d$outcome == "PE", FUN = function(x) 0) == 0]))
cat("PE SNR decomposition not run: only one PE study manipulates SNR (A04).\n")

cat("\n=== 4. IS: is lnVR just the mechanical shadow of the accuracy loss? ===\n")
src <- read.csv(file.path(P26_ROOT, ".cache", "data", "source_rows.csv"),
                stringsAsFactors = FALSE)
unit <- setNames(src$unit, src$source_id)
d$ctrl_unit <- unit[d$control_group_id]
dis <- d[d$outcome == "IS" & d$ctrl_unit != "RAU", ]
dis$raw_drop <- dis$m2 - dis$m1
dis$base <- dis$m2
dis$drop_c <- dis$raw_drop - mean(dis$raw_drop)
dis$base_c <- dis$base - mean(dis$base)
m_art <- rma.mv(lnVR, mkV_vr(dis), mods = ~ drop_c + base_c, data = dis,
                random = list(~ 1 | study, ~ 1 | effect_ID),
                method = "REML", test = "t", dfs = "contain")
ct_art <- clubSandwich::coef_test(m_art, vcov = "CR2", cluster = dis$study)
art <- data.frame(term = rownames(m_art$b), est = as.numeric(m_art$b),
                  p_model = m_art$pval, p_CR2 = ct_art$p_Satt, df_CR2 = ct_art$df_Satt)
print(art, row.names = FALSE, digits = 4)
cat("  if lnVR were mechanical, raw drop would absorb it: residual intercept\n")
cat("  term and a strong drop coefficient = scale shadow, not dispersion.\n")
cat(sprintf("cor(lnVR, raw drop) = %+.3f ; cor(lnVR, baseline) = %+.3f\n",
            cor(dis$lnVR, dis$raw_drop), cor(dis$lnVR, dis$base)))

cat("\n=== 5. IS spread index: SD relative to the binomial bound ===\n")
bnd <- function(p) sqrt(p * (100 - p))
bad <- dis$m1 <= 0 | dis$m1 >= 100 | dis$m2 <= 0 | dis$m2 >= 100
cat(sprintf("rows outside (0,100) on either arm, excluded from the index: %d\n",
            sum(bad)))
dis <- dis[!bad, ]
dis$spread1 <- dis$sd1 / bnd(dis$m1)
dis$spread2 <- dis$sd2 / bnd(dis$m2)
dis$lnSR <- log(dis$spread1 / dis$spread2)
cat(sprintf("control-arm spread index: median %.3f (IQR %.3f-%.3f)\n",
            median(dis$spread2), quantile(dis$spread2, .25), quantile(dis$spread2, .75)))
cat(sprintf("noise-arm spread index  : median %.3f (IQR %.3f-%.3f)\n",
            median(dis$spread1), quantile(dis$spread1, .25), quantile(dis$spread1, .75)))
cat(sprintf("ln spread ratio: median %+.3f (IQR %+.3f to %+.3f); fraction > 0: %.0f%%\n",
            median(dis$lnSR), quantile(dis$lnSR, .25, na.rm = TRUE),
            quantile(dis$lnSR, .75, na.rm = TRUE), 100 * mean(dis$lnSR > 0)))
w <- 1 / dis$lnVR_vi
m_sr <- rma.mv(lnSR, V = diag(dis$lnVR_vi), data = dis,
               random = list(~ 1 | study, ~ 1 | effect_ID),
               method = "REML", test = "t", dfs = "contain")
ct_sr <- clubSandwich::coef_test(m_sr, vcov = "CR2", cluster = dis$study)
cat(sprintf("pooled ln spread ratio = %+.4f [%.4f, %.4f]  CR2 p %.4g (df %.1f)\n",
            m_sr$b[1], m_sr$ci.lb[1], m_sr$ci.ub[1], ct_sr$p_Satt[1], ct_sr$df_Satt[1]))
cat("  > 0 means: noise spreads children MORE than the bounded scale allows\n")
cat("  mechanically -- a dispersion effect beyond the mean-shift coupling.\n")

cat("\n=== 6. RT: does variability rise BEYOND the proportional coupling? ===\n")
drt <- d[d$outcome == "RT", ]
cvr <- escalc(measure = "CVR", m1i = m1, sd1i = sd1, n1i = n1,
              m2i = m2, sd2i = sd2, n2i = n2, data = drt)
drt$lnCVR <- cvr$yi; drt$lnCVR_vi <- cvr$vi
drt$SNR_between <- ave(drt$SNR, drt$study, FUN = mean)
drt$SNR_within  <- drt$SNR - drt$SNR_between
drt$SNR_between_c <- drt$SNR_between - mean(drt$SNR_between)
m_cvr <- rma.mv(lnCVR, lnCVR_vi, data = drt,
                random = list(~ 1 | study, ~ 1 | effect_ID),
                method = "REML", test = "t", dfs = "contain")
ct_cvr <- clubSandwich::coef_test(m_cvr, vcov = "CR2", cluster = drt$study)
cat(sprintf("pooled lnCVR = %+.4f [%.4f, %.4f]  p_model %.4g | CR2 p %.4g (df %.1f)\n",
            m_cvr$b[1], m_cvr$ci.lb[1], m_cvr$ci.ub[1], m_cvr$pval[1],
            ct_cvr$p_Satt[1], ct_cvr$df_Satt[1]))
m_cvr2 <- rma.mv(lnCVR, lnCVR_vi, mods = ~ SNR_within + SNR_between_c, data = drt,
                 random = list(~ 1 | study, ~ 1 | effect_ID),
                 method = "REML", test = "t", dfs = "contain")
ct_cvr2 <- clubSandwich::coef_test(m_cvr2, vcov = "CR2", cluster = drt$study)
cvr_snr <- data.frame(term = rownames(m_cvr2$b), est = as.numeric(m_cvr2$b),
                      p_model = m_cvr2$pval, p_CR2 = ct_cvr2$p_Satt,
                      df_CR2 = ct_cvr2$df_Satt)
print(cvr_snr[cvr_snr$term != "intrcpt", ], row.names = FALSE, digits = 4)

cat("\n=== 7. FDR across family R12 ===\n")
snr_nb <- snr_res[snr_res$term != "intrcpt", ]
fam <- rbind(
  data.frame(q = paste0("pooled lnVR ", pool$family), p = pool$p_CR2),
  data.frame(q = paste0("lnVR SNR ", snr_nb$family, " ", snr_nb$term),
             p = snr_nb$p_CR2),
  data.frame(q = paste0("IS artefact model: ", art$term[art$term != "intrcpt"]),
             p = art$p_CR2[art$term != "intrcpt"]))
fam <- fam[!is.na(fam$p), ]
fam$q_BH <- p.adjust(fam$p, method = "BH")
fam <- fam[order(fam$p), ]
print(fam, row.names = FALSE, digits = 4)

write.csv(pool, file.path(OUT, "pooled_lnVR.csv"), row.names = FALSE)
write.csv(snr_res, file.path(OUT, "lnVR_SNR_decomposition.csv"), row.names = FALSE)
write.csv(art, file.path(OUT, "IS_artefact_model.csv"), row.names = FALSE)
write.csv(dis[, c("effect_ID", "study", "m1", "m2", "sd1", "sd2", "lnVR",
                  "raw_drop", "base", "spread1", "spread2", "lnSR")],
          file.path(OUT, "IS_dispersion_data.csv"), row.names = FALSE)
write.csv(data.frame(term = c("pooled", cvr_snr$term),
                     est = c(m_cvr$b[1], cvr_snr$est),
                     p_CR2 = c(ct_cvr$p_Satt[1], cvr_snr$p_CR2),
                     df_CR2 = c(ct_cvr$df_Satt[1], cvr_snr$df_CR2)),
          file.path(OUT, "RT_lnCVR.csv"), row.names = FALSE)
write.csv(fam, file.path(OUT, "family_R12_fdr.csv"), row.names = FALSE)
cat("\nwritten to", OUT, "\n")
sink()
