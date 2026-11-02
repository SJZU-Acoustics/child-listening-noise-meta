
source(file.path("R", "helpers.R"))
OUT <- file.path(p26_outdir("analysis_17_small_study"), "outcome")

D <- p26_data()
d <- D$effects
d <- d[order(d$study, d$effect_ID), ]

sink(file.path(OUT, "console.txt"), split = TRUE)
cat("A17 - small-study effects: correction, not just detection\n\n")

mkV <- function(sub) vcalc(vi = vi, cluster = study, grp1 = experiment_group_id,
                           grp2 = control_group_id, w1 = n1, w2 = n2,
                           data = sub, sparse = FALSE)

d$sqrt_inv_n <- sqrt(1 / pmin(d$n1, d$n2))
d$inv_n      <- 1 / pmin(d$n1, d$n2)

fit_pp <- function(sub, label) {
  V <- mkV(sub)
  m_pet   <- rma.mv(yi_al, V, mods = ~ sqrt_inv_n, data = sub,
                    random = list(~ 1 | study, ~ 1 | effect_ID),
                    method = "REML", test = "t", dfs = "contain")
  ct_pet  <- clubSandwich::coef_test(m_pet, vcov = "CR2", cluster = sub$study)
  m_peese <- rma.mv(yi_al, V, mods = ~ inv_n, data = sub,
                    random = list(~ 1 | study, ~ 1 | effect_ID),
                    method = "REML", test = "t", dfs = "contain")
  ct_peese <- clubSandwich::coef_test(m_peese, vcov = "CR2", cluster = sub$study)
  m_pool  <- rma.mv(yi_al, V, data = sub,
                    random = list(~ 1 | study, ~ 1 | effect_ID),
                    method = "REML", test = "t", dfs = "contain")
  cat(sprintf("\n%s (k=%d, %d studies; pooled %+.3f)\n", label, m_pet$k,
              length(unique(sub$study)), m_pool$b[1]))
  cat(sprintf("  PET   : slope %+.3f (CR2 p %.4g, df %.1f) | adjusted (intercept) %+.3f [%.3f, %.3f]\n",
              m_pet$b[2], ct_pet$p_Satt[2], ct_pet$df_Satt[2],
              m_pet$b[1], m_pet$ci.lb[1], m_pet$ci.ub[1]))
  cat(sprintf("  PEESE : slope %+.3f (CR2 p %.4g, df %.1f) | adjusted (intercept) %+.3f [%.3f, %.3f]\n",
              m_peese$b[2], ct_peese$p_Satt[2], ct_peese$df_Satt[2],
              m_peese$b[1], m_peese$ci.lb[1], m_peese$ci.ub[1]))
  data.frame(family = label, pooled = m_pool$b[1],
             pet_slope = m_pet$b[2], pet_slope_p = ct_pet$p_Satt[2],
             pet_adj = m_pet$b[1], pet_adj_lo = m_pet$ci.lb[1],
             pet_adj_hi = m_pet$ci.ub[1],
             peese_slope = m_peese$b[2], peese_slope_p = ct_peese$p_Satt[2],
             peese_adj = m_peese$b[1], peese_adj_lo = m_peese$ci.lb[1],
             peese_adj_hi = m_peese$ci.ub[1])
}

cat("=== 1. PET-PEESE on the artefact-free precision scale ===\n")
pp <- rbind(fit_pp(d, "ALL"),
            fit_pp(d[d$outcome == "IS", ], "IS"),
            fit_pp(d[d$outcome == "RT", ], "RT"),
            fit_pp(d[d$outcome == "PE", ], "PE"))
cat(sprintf("\nRT has %d clusters, PE %d -- underpowered legs, read with A06's caution.\n",
            length(unique(d$study[d$outcome == "RT"])), length(unique(d$study[d$outcome == "PE"]))))
cat(sprintf("\nimplied shrinkage of the pooled estimate (conditional rule):\n"))
for (i in seq_len(nrow(pp))) {
  adj <- if (pp$pet_slope_p[i] < 0.05) pp$peese_adj[i] else pp$pet_adj[i]
  cat(sprintf("  %-4s pooled %+.3f -> adjusted %+.3f (%.0f%% shrinkage)\n",
              pp$family[i], pp$pooled[i], adj,
              100 * (pp$pooled[i] - adj) / pp$pooled[i]))
}

cat("\n=== 2. selection model attempt ===\n")
cat("selmodel has no rma.mv method (metafor 5.0.1); run on study x outcome\n")
cat("cell aggregates under rma.uni instead -- approximate: it ignores the\n")
cat("clustering of cells within studies, and is reported as such.\n")
cells <- list()
for (cx in unique(paste(d$study, d$outcome))) {
  idx <- which(paste(d$study, d$outcome) == cx)
  sub <- d[idx, ]
  Vc <- if (nrow(sub) == 1) matrix(sub$vi, 1, 1) else mkV(sub)
  prec <- solve(Vc)
  one <- rep(1, nrow(sub))
  beta <- as.numeric((t(one) %*% prec %*% sub$yi_al) / (t(one) %*% prec %*% one))
  varb <- as.numeric(1 / (t(one) %*% prec %*% one))
  cells[[cx]] <- data.frame(study = sub$study[1], outcome = sub$outcome[1],
                            yi = beta, vi = varb, k_cell = nrow(sub))
}
cell <- do.call(rbind, cells)
cat(sprintf("%d cells from %d effects (%d studies)\n", nrow(cell), nrow(d),
            length(unique(cell$study))))
sel_out <- list()
for (fam in c("ALL", "IS")) {
  sub <- if (fam == "ALL") cell else cell[cell$outcome == "IS", ]
  base_fit <- rma.uni(yi, vi, data = sub, method = "REML", test = "t")
  sm <- try(selmodel(base_fit, type = "stepfun", steps = c(.025, .10, .50, 1),
                     verbose = FALSE), silent = TRUE)
  if (inherits(sm, "try-error")) {
    cat(sprintf("%s: selection model FAILED - %s\n", fam,
                conditionMessage(attr(sm, "condition"))))
    sel_out[[fam]] <- data.frame(family = fam, converged = FALSE)
  } else {
    cat(sprintf("%s: adjusted %+.3f [%.3f, %.3f]; LRT p %.4g; deltas %s\n",
                fam, sm$b[1], sm$ci.lb[1], sm$ci.ub[1], sm$LRTp,
                paste(round(sm$delta, 3), collapse = ", ")))
    sel_out[[fam]] <- data.frame(family = fam, converged = TRUE,
                                 adj = sm$b[1], lo = sm$ci.lb[1], hi = sm$ci.ub[1],
                                 LRT_p = sm$LRTp,
                                 deltas = paste(round(sm$delta, 3), collapse = ";"))
  }
}
sel <- do.call(rbind, sel_out)
i2 <- read.csv(file.path(P26_ROOT, ".cache", "analysis", "analysis_05_heterogeneity", "outcome",
                         "heterogeneity_decomposition.csv"))$I2_total
cat(sprintf("\ncaution: the step function is fit on %d/%d cells whose effects are\n",
            nrow(cell), sum(cell$outcome == "IS")))
cat(sprintf("highly heterogeneous (I2 %.0f-%.0f%%); a one-parameter selection story that\n",
            min(i2), max(i2)))
cat("cannot absorb that heterogeneity overcorrects -- read the adjusted\n")
cat("intercept as a bound on how bad it could be, not as an estimate.\n")
write.csv(cell, file.path(OUT, "cell_aggregates.csv"), row.names = FALSE)

cat("\n=== 3. reconciliation with A06 and A12 ===\n")
eg6 <- read.csv(file.path(P26_ROOT, ".cache", "analysis", "analysis_06_publication_bias", "outcome", "egger_tests.csv"))
dd12 <- read.csv(file.path(P26_ROOT, ".cache", "analysis", "analysis_12_design_descriptors", "outcome", "design_descriptors.csv"))
n12 <- dd12[dd12$term == "n_between_c", ]
eg6 <- eg6[grepl("sample-size", eg6$variant), ]
npart <- pmin(d$n1, d$n2)
cat(sprintf("A06 sample-size Egger: CR2 p %.3f to %.3f; A12 moderator: %+.4f/participant (CR2 p %.3f).\n",
            min(eg6$p_CR2), max(eg6$p_CR2), n12$estimate, n12$p_CR2))
cat(sprintf("A12 implies %.2f SMD across the observed %d-to-%d participant range.\n",
            abs(n12$estimate) * (max(npart) - min(npart)), min(npart), max(npart)))

write.csv(pp, file.path(OUT, "pet_peese.csv"), row.names = FALSE)
write.csv(sel, file.path(OUT, "selection_model.csv"), row.names = FALSE)
cat("\nwritten to", OUT, "\n")
sink()
