
source(file.path("R", "helpers.R"))
OUT <- file.path(p26_outdir("analysis_18_inference_layers"), "outcome")

D <- p26_data()
d <- D$effects
d <- d[order(d$study, d$effect_ID), ]

sink(file.path(OUT, "console.txt"), split = TRUE)
cat("A18 - CWB on the within-slopes, and the cell-level check\n\n")

mkV <- function(sub) vcalc(vi = vi, cluster = study, grp1 = experiment_group_id,
                           grp2 = control_group_id, w1 = n1, w2 = n2,
                           data = sub, sparse = FALSE)

cat("=== 1. cluster wild bootstrap on the within-study SNR slopes ===\n")
a15 <- read.csv(file.path(P26_ROOT, ".cache", "analysis", "analysis_15_pairs_deepened",
                          "outcome", "snr_decomposition_cwb.csv"))
a15 <- a15[a15$term == "SNR_within", ]
cat(sprintf("(%d resamples; the RQ1 slopes' CR2 df are small -- where Satterthwaite\n", BOOT_R))
cat(sprintf("is most fragile, and where A15's pair slope went from model p %.2g\n", a15$p_model))
cat(sprintf("to CWB p %.3f)\n\n", a15$p_CWB))
cwb_rows <- list()
for (fam in c("ALL", "IS", "RT")) {
  sub <- if (fam == "ALL") d else d[d$outcome == fam, ]
  sub$SNR_between <- ave(sub$SNR, sub$study, FUN = mean)
  sub$SNR_within  <- sub$SNR - sub$SNR_between
  sub$SNR_between_c <- sub$SNR_between - mean(sub$SNR_between)
  V <- mkV(sub)
  m <- rma.mv(yi_al, V, mods = ~ SNR_within + SNR_between_c, data = sub,
              random = list(~ 1 | study, ~ 1 | effect_ID),
              method = "REML", test = "t", dfs = "contain")
  ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = sub$study)
  bt <- try(Wald_test_cwb(full_model = m, constraints = constrain_zero(2),
                          R = BOOT_R, seed = 26), silent = TRUE)
  p_cwb <- if (inherits(bt, "try-error")) NA_real_ else bt$p_val
  cat(sprintf("%-4s within-slope %+.4f: model p %.4g | CR2 p %.4g (df %.1f) | CWB p %s\n",
              fam, m$b[2], m$pval[2], ct$p_Satt[2], ct$df_Satt[2],
              ifelse(is.na(p_cwb), paste("FAILED -", conditionMessage(attr(bt, "condition"))),
                     sprintf("%.4f", p_cwb))))
  cwb_rows[[fam]] <- data.frame(family = fam, est = m$b[2],
                                p_model = m$pval[2], p_CR2 = ct$p_Satt[2],
                                df_CR2 = ct$df_Satt[2], p_CWB = p_cwb)
}
cwb_res <- do.call(rbind, cwb_rows)
cat("PE not run: only one PE study manipulates SNR (A04).\n")

cat("\n=== 2. study x outcome cell aggregates: does the pooled effect move? ===\n")
cells <- list()
for (cx in unique(paste(d$study, d$outcome))) {
  idx <- which(paste(d$study, d$outcome) == cx)
  sub <- d[idx, ]
  Vc <- if (nrow(sub) == 1) matrix(sub$vi, 1, 1) else mkV(sub)
  prec <- solve(Vc); one <- rep(1, nrow(sub))
  cells[[cx]] <- data.frame(study = sub$study[1], outcome = sub$outcome[1],
                            yi = as.numeric((t(one) %*% prec %*% sub$yi_al) / (t(one) %*% prec %*% one)),
                            vi = as.numeric(1 / (t(one) %*% prec %*% one)),
                            k_cell = nrow(sub))
}
cell <- do.call(rbind, cells)
cat(sprintf("%d cells from %d effects, %d studies\n", nrow(cell), nrow(d),
            length(unique(cell$study))))
cell_rows <- list()
for (fam in c("ALL", "IS", "RT", "PE")) {
  sub <- if (fam == "ALL") cell else cell[cell$outcome == fam, ]
  m <- rma.mv(yi, vi, random = ~ 1 | study, data = sub,
              method = "REML", test = "t", dfs = "contain")
  ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = sub$study)
  cat(sprintf("%-4s cell-level pooled %+.3f [%.3f, %.3f]  CR2 p %.4g (df %.1f, %d cells, %d studies)\n",
              fam, m$b[1], m$ci.lb[1], m$ci.ub[1], ct$p_Satt[1], ct$df_Satt[1],
              m$k, length(unique(sub$study))))
  cell_rows[[fam]] <- data.frame(family = fam, est = m$b[1], lo = m$ci.lb[1],
                                 hi = m$ci.ub[1], p_CR2 = ct$p_Satt[1],
                                 df_CR2 = ct$df_Satt[1], cells = m$k,
                                 studies = length(unique(sub$study)))
}
cell_res <- do.call(rbind, cell_rows)
A02REF <- local({
  w <- read.csv(file.path(P26_ROOT, ".cache", "analysis", "analysis_02_shared_control",
                          "outcome", "working_model_comparison.csv"))
  w <- w[w$working_model == "shared-arm" & grepl("^empty", w$quantity), ]
  setNames(w$estimate, sub("empty ", "", w$quantity))
})
cat(sprintf("\nreference, effect-level adopted estimates: ALL %+.3f, IS %+.3f,\nRT %+.3f, PE %+.3f (A02)\n",
            A02REF["ALL"], A02REF["IS"], A02REF["RT"], A02REF["PE"]))

write.csv(cwb_res, file.path(OUT, "cwb_within_slopes.csv"), row.names = FALSE)
write.csv(cell_res, file.path(OUT, "cell_level_pooled.csv"), row.names = FALSE)
write.csv(cell, file.path(OUT, "cell_aggregates.csv"), row.names = FALSE)
cat("\nwritten to", OUT, "\n")
sink()
