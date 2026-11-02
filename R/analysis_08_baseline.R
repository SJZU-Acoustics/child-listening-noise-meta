
source(file.path("R", "helpers.R"))
OUT <- file.path(p26_outdir("analysis_08_baseline"), "outcome")

D <- p26_data()
d <- D$effects
d <- d[order(d$study, d$effect_ID), ]
src <- read.csv(file.path(P26_ROOT, ".cache", "data", "source_rows.csv"),
                stringsAsFactors = FALSE)
unit <- setNames(src$unit, src$source_id)
d$ctrl_unit <- unit[d$control_group_id]

mkV <- function(sub) vcalc(vi = vi, cluster = study, grp1 = experiment_group_id,
                           grp2 = control_group_id, w1 = n1, w2 = n2,
                           data = sub, sparse = FALSE)

sink(file.path(OUT, "console.txt"), split = TRUE)
cat("A08 - the control arm as information\n\n")

cat("=== 1. what is on a common scale ===\n")
cat(sprintf("IS  : %d effects on percent (%d RAU excluded)\n",
            sum(d$outcome == "IS" & d$ctrl_unit != "RAU"),
            sum(d$outcome == "IS" & d$ctrl_unit == "RAU")))
cat(sprintf("RT  : %d effects in ms (T08 unified)\n", sum(d$outcome == "RT")))
cat("PE  : mixed rating scales, no common metric -- excluded\n")

cat("\n=== 2. IS: does the penalty depend on baseline accuracy? ===\n")
dis <- d[d$outcome == "IS" & d$ctrl_unit != "RAU", ]
dis$base <- dis$m2
dis$base_between <- ave(dis$base, dis$study, FUN = mean)
dis$base_within  <- dis$base - dis$base_between
dis$base_between_c <- dis$base_between - mean(dis$base_between)

cat(sprintf("baseline accuracy: median %.1f%%  IQR %.1f-%.1f  range %.1f-%.1f  (%d studies)\n",
            median(dis$base), quantile(dis$base, .25), quantile(dis$base, .75),
            min(dis$base), max(dis$base), length(unique(dis$study))))
cat(sprintf("within-study variation in baseline: %.1f%% of total SS\n",
            100 * sum(dis$base_within^2) / (sum(dis$base_within^2) + sum((dis$base_between - mean(dis$base))^2))))

V <- mkV(dis)
m <- rma.mv(yi_al, V, mods = ~ base_within + base_between_c, data = dis,
            random = list(~ 1 | study, ~ 1 | effect_ID),
            method = "REML", test = "t", dfs = "contain")
ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = dis$study)
base_res <- data.frame(term = rownames(m$b), estimate = as.numeric(m$b),
                       ci_lo = m$ci.lb, ci_hi = m$ci.ub,
                       p_model = m$pval, p_CR2 = ct$p_Satt, df_CR2 = ct$df_Satt)
print(base_res, row.names = FALSE, digits = 4)
cat("  a positive slope means: the higher the quiet-condition accuracy, the\n")
cat("  larger the standardised penalty from noise.\n")

cat("\n=== 3. is the IS effect tracking the control SD rather than the raw change? ===\n")
dis$raw_drop <- dis$m2 - dis$m1                # percentage points lost in noise
dis$ctrl_sd  <- dis$sd2
cat(sprintf("raw drop (percentage points): median %.1f, IQR %.1f-%.1f\n",
            median(dis$raw_drop), quantile(dis$raw_drop, .25), quantile(dis$raw_drop, .75)))
cat(sprintf("control SD: median %.2f, IQR %.2f-%.2f\n",
            median(dis$ctrl_sd), quantile(dis$ctrl_sd, .25), quantile(dis$ctrl_sd, .75)))
cat(sprintf("\ncorrelation of control SD with baseline accuracy: r = %+.3f  (Spearman %+.3f)\n",
            cor(dis$ctrl_sd, dis$base), cor(dis$ctrl_sd, dis$base, method = "spearman")))
cat(sprintf("correlation of |g| with control SD              : r = %+.3f  (Spearman %+.3f)\n",
            cor(abs(dis$yi_al), dis$ctrl_sd), cor(abs(dis$yi_al), dis$ctrl_sd, method = "spearman")))
cat(sprintf("correlation of |g| with raw drop                : r = %+.3f  (Spearman %+.3f)\n",
            cor(abs(dis$yi_al), dis$raw_drop), cor(abs(dis$yi_al), dis$raw_drop, method = "spearman")))

cat("\n  Formal test: does g depend on control SD once the raw drop is held fixed?\n")
dis$sd_c   <- dis$ctrl_sd - mean(dis$ctrl_sd)
dis$drop_c <- dis$raw_drop - mean(dis$raw_drop)
m2f <- rma.mv(yi_al, V, mods = ~ drop_c + sd_c, data = dis,
              random = list(~ 1 | study, ~ 1 | effect_ID),
              method = "REML", test = "t", dfs = "contain")
ct2 <- clubSandwich::coef_test(m2f, vcov = "CR2", cluster = dis$study)
art <- data.frame(term = rownames(m2f$b), estimate = as.numeric(m2f$b),
                  p_model = m2f$pval, p_CR2 = ct2$p_Satt, df_CR2 = ct2$df_Satt)
print(art, row.names = FALSE, digits = 4)
cat("  A strong negative sd_c term is the artefact signature: the same raw drop\n")
cat("  yields a bigger g when the control SD happens to be small.\n")

m_null <- rma.mv(yi_al, V, data = dis, random = list(~ 1 | study, ~ 1 | effect_ID),
                 method = "REML", test = "t", dfs = "contain")
cat(sprintf("\n  tau2 total: null %.3f -> with drop+SD %.3f  (%.0f%% of between-effect variance explained)\n",
            sum(m_null$sigma2), sum(m2f$sigma2),
            100 * (sum(m_null$sigma2) - sum(m2f$sigma2)) / sum(m_null$sigma2)))

cat("\n=== 4. RT: does the penalty depend on baseline latency? ===\n")
drt <- d[d$outcome == "RT", ]
drt$base <- drt$m2
drt$base_between <- ave(drt$base, drt$study, FUN = mean)
drt$base_within  <- drt$base - drt$base_between
drt$base_between_c <- drt$base_between - mean(drt$base_between)
cat(sprintf("baseline RT (ms): median %.0f  range %.0f-%.0f  (%d studies)\n",
            median(drt$base), min(drt$base), max(drt$base), length(unique(drt$study))))
Vr <- mkV(drt)
mr <- rma.mv(yi_al, Vr, mods = ~ base_within + base_between_c, data = drt,
             random = list(~ 1 | study, ~ 1 | effect_ID),
             method = "REML", test = "t", dfs = "contain")
ctr <- clubSandwich::coef_test(mr, vcov = "CR2", cluster = drt$study)
rt_res <- data.frame(term = rownames(mr$b), estimate = as.numeric(mr$b),
                     p_model = mr$pval, p_CR2 = ctr$p_Satt, df_CR2 = ctr$df_Satt)
print(rt_res, row.names = FALSE, digits = 4)

write.csv(base_res, file.path(OUT, "IS_baseline_moderation.csv"), row.names = FALSE)
write.csv(art, file.path(OUT, "IS_scale_artefact_test.csv"), row.names = FALSE)
write.csv(rt_res, file.path(OUT, "RT_baseline_moderation.csv"), row.names = FALSE)
write.csv(dis[, c("effect_ID", "study", "base", "raw_drop", "ctrl_sd", "yi_al", "vi")],
          file.path(OUT, "IS_baseline_data.csv"), row.names = FALSE)
cat("\nwritten to", OUT, "\n")
sink()
