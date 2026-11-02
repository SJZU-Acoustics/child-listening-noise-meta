
source(file.path("R", "helpers.R"))
OUT <- file.path(p26_outdir("analysis_10_moderator_estimability"), "outcome")

D <- p26_data()
d <- D$effects
d <- d[order(d$study, d$effect_ID), ]

mkV <- function(sub) vcalc(vi = vi, cluster = study, grp1 = experiment_group_id,
                           grp2 = control_group_id, w1 = n1, w2 = n2,
                           data = sub, sparse = FALSE)

sink(file.path(OUT, "console.txt"), split = TRUE)
cat("A10 - moderator estimability\n\n")

cat("=== 1. can each moderator be identified within studies? ===\n")
mods <- c("SNR", "age", "masking", "noise_type", "task", "outcome")
aud <- list()
for (v in mods) {
  n <- tapply(d[[v]], d$study, function(x) length(unique(x)))
  eff_in <- sum(d$study %in% names(n)[n > 1])
  aud[[length(aud) + 1]] <- data.frame(
    moderator = v, studies_varying = sum(n > 1), studies_total = length(n),
    effects_in_varying = eff_in,
    verdict = if (sum(n > 1) >= 8) "within-study estimable"
              else if (sum(n > 1) >= 3) "weakly identified"
              else "BETWEEN-STUDY ONLY")
  cat(sprintf("  %-11s varies within %2d/%d studies (%3d effects)  ->  %s\n",
              v, sum(n > 1), length(n), eff_in, aud[[length(aud)]]$verdict))
}
aud <- do.call(rbind, aud)
write.csv(aud, file.path(OUT, "moderator_estimability.csv"), row.names = FALSE)

cat("\n  Reading: a moderator that never varies within a study cannot be\n")
cat("  distinguished from everything else that differs between those studies.\n")
cat("  Its pooled coefficient is a between-literature comparison, whatever the\n")
cat("  model prints.\n")

cat("\n=== 2. the one well-identified contrast: outcome type, within study ===\n")
n_out <- tapply(d$outcome, d$study, function(x) length(unique(x)))
dw <- d[d$study %in% names(n_out)[n_out > 1], ]
cat(sprintf("%d studies measuring >1 outcome family, %d effects\n",
            length(unique(dw$study)), nrow(dw)))
cat("  within-study composition:\n")
print(table(dw$study, dw$outcome))

V <- mkV(dw)
m <- rma.mv(yi_al, V, mods = ~ Outcome_Type, data = dw,
            random = list(~ 1 | study, ~ 1 | effect_ID),
            method = "REML", test = "t", dfs = "contain")
ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = dw$study)
cat("\n  random-effects (between+within information):\n")
oc <- data.frame(term = rownames(m$b), estimate = as.numeric(m$b),
                 ci_lo = m$ci.lb, ci_hi = m$ci.ub,
                 p_model = m$pval, p_CR2 = ct$p_Satt, df_CR2 = ct$df_Satt)
print(oc, row.names = FALSE, digits = 4)

cat("\n  with study FIXED effects (purely within-study identification):\n")
mfe <- rma.mv(yi_al, V, mods = ~ Outcome_Type + factor(study), data = dw,
              random = ~ 1 | effect_ID, method = "REML", test = "t", dfs = "contain")
ctfe <- clubSandwich::coef_test(mfe, vcov = "CR2", cluster = dw$study)
keep <- grep("Outcome_Type", rownames(mfe$b))
ocfe <- data.frame(term = rownames(mfe$b)[keep], estimate = as.numeric(mfe$b)[keep],
                   p_model = mfe$pval[keep], p_CR2 = ctfe$p_Satt[keep],
                   df_CR2 = ctfe$df_Satt[keep])
print(ocfe, row.names = FALSE, digits = 4)
write.csv(rbind(cbind(model = "random effects", oc[, c("term","estimate","p_model","p_CR2","df_CR2")]),
                cbind(model = "study fixed effects", ocfe)),
          file.path(OUT, "outcome_type_within_study.csv"), row.names = FALSE)

cat("\n=== 3. the audit's claim, tested: masking under study fixed effects ===\n")
a01 <- read.csv(file.path(P26_ROOT, ".cache", "analysis", "analysis_01_legacy_reproduction",
                          "outcome", "legacy_reproduction.csv"))
a01 <- a01[a01$quantity == "maskingIM", ]
cat(sprintf("A01 reported masking IM vs EM at %+.3f, p_model %.3f / CR2 %.3f on pooled data.\n",
            a01$estimate, a01$p_model, a01$p_robust))
cat("If that contrast is between-study, adding study fixed effects should leave\n")
cat("almost nothing to estimate it from.\n\n")
nm_ <- tapply(d$masking, d$study, function(x) length(unique(x)))
cat(sprintf("studies where masking varies: %d (%s), %d effects\n",
            sum(nm_ > 1), paste(names(nm_)[nm_ > 1], collapse = "; "),
            sum(d$study %in% names(nm_)[nm_ > 1])))
dm <- d[d$study %in% names(nm_)[nm_ > 1], ]
mm <- try(rma.mv(yi_al, mkV(dm), mods = ~ masking, data = dm,
                 random = list(~ 1 | study, ~ 1 | effect_ID),
                 method = "REML", test = "t", dfs = "contain"), silent = TRUE)
if (!inherits(mm, "try-error")) {
  cat(sprintf("\nmasking IM vs EM, within the single study that varies it: %+.4f  p_model %.4g\n",
              mm$b[2], mm$pval[2]))
  cat("  (one cluster: no cluster-robust inference is possible, and no\n")
  cat("   generalisation beyond that study is available)\n")
}
cat("\nVerdict on masking: it is a between-study contrast in this evidence base.\n")
cat("It should be reported descriptively, or not at all, but not as a moderator\n")
cat("of the noise effect.\n")

sink()
