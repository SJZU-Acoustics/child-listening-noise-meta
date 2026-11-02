
source(file.path("R", "helpers.R"))
OUT <- file.path(p26_outdir("analysis_11_raw_scale"), "outcome")

D <- p26_data()
d <- D$effects
d <- d[order(d$study, d$effect_ID), ]
src <- read.csv(file.path(P26_ROOT, ".cache", "data", "source_rows.csv"),
                stringsAsFactors = FALSE)
d$ctrl_unit <- setNames(src$unit, src$source_id)[d$control_group_id]

sink(file.path(OUT, "console.txt"), split = TRUE)
cat("A11 - the raw scale as a parallel outcome\n\n")

prep <- function(sub, adverse_is_higher) {
  if (adverse_is_higher) {   # RT: noise - control
    e <- escalc(measure = "MD", m1i = m1, sd1i = sd1, n1i = n1,
                m2i = m2, sd2i = sd2, n2i = n2, data = sub)
  } else {                   # IS: control - noise
    e <- escalc(measure = "MD", m1i = m2, sd1i = sd2, n1i = n2,
                m2i = m1, sd2i = sd1, n2i = n1, data = sub)
  }
  e
}
mkV <- function(sub) vcalc(vi = vi, cluster = study, grp1 = experiment_group_id,
                           grp2 = control_group_id, w1 = n1, w2 = n2,
                           data = sub, sparse = FALSE)
fit <- function(sub, mods = ~ 1, idx = 1) {
  V <- mkV(sub)
  m <- rma.mv(yi, V, mods = mods, data = sub,
              random = list(~ 1 | study, ~ 1 | effect_ID),
              method = "REML", test = "t", dfs = "contain")
  ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = sub$study)
  pr <- predict(m)
  list(m = m, est = m$b[idx], lo = m$ci.lb[idx], hi = m$ci.ub[idx],
       p_model = m$pval[idx], p_CR2 = ct$p_Satt[idx], df = ct$df_Satt[idx],
       pi_lo = pr$pi.lb[1], pi_hi = pr$pi.ub[1])
}

cat("=== 1. IS: accuracy lost, in percentage points ===\n")
dis <- prep(d[d$outcome == "IS" & d$ctrl_unit != "RAU", ], FALSE)
cat(sprintf("%d effects, %d studies (%d RAU rows excluded)\n",
            nrow(dis), length(unique(dis$study)),
            sum(d$outcome == "IS" & d$ctrl_unit == "RAU")))
r <- fit(dis)
cat(sprintf("pooled loss  %+.2f points [%+.2f, %+.2f]  p_model %.3g | CR2 p %.4g (df %.1f)\n",
            r$est, r$lo, r$hi, r$p_model, r$p_CR2, r$df))
cat(sprintf("prediction interval [%+.2f, %+.2f] points\n", r$pi_lo, r$pi_hi))

dis$SNR_between <- ave(dis$SNR, dis$study, FUN = mean)
dis$SNR_within  <- dis$SNR - dis$SNR_between
dis$SNR_between_c <- dis$SNR_between - mean(dis$SNR_between)
rs <- fit(dis, ~ SNR_within + SNR_between_c, idx = 2)
rb <- fit(dis, ~ SNR_within + SNR_between_c, idx = 3)
cat(sprintf("SNR within-study  %+.3f points per dB  CR2 p %.4g (df %.1f)\n",
            rs$est, rs$p_CR2, rs$df))
cat(sprintf("SNR between-study %+.3f points per dB  CR2 p %.4g (df %.1f)\n",
            rb$est, rb$p_CR2, rb$df))

cat("\n=== 2. RT: latency added, in milliseconds ===\n")
drt <- prep(d[d$outcome == "RT", ], TRUE)
cat(sprintf("%d effects, %d studies\n", nrow(drt), length(unique(drt$study))))
rr <- fit(drt)
cat(sprintf("pooled slowing %+.0f ms [%+.0f, %+.0f]  p_model %.3g | CR2 p %.4g (df %.1f)\n",
            rr$est, rr$lo, rr$hi, rr$p_model, rr$p_CR2, rr$df))
cat(sprintf("prediction interval [%+.0f, %+.0f] ms\n", rr$pi_lo, rr$pi_hi))
cat(sprintf("relative to a median baseline of %.0f ms, that is %.1f%%\n",
            median(drt$m2), 100 * rr$est / median(drt$m2)))
drt$SNR_between <- ave(drt$SNR, drt$study, FUN = mean)
drt$SNR_within  <- drt$SNR - drt$SNR_between
drt$SNR_between_c <- drt$SNR_between - mean(drt$SNR_between)
rrs <- fit(drt, ~ SNR_within + SNR_between_c, idx = 2)
cat(sprintf("SNR within-study %+.1f ms per dB  CR2 p %.4g (df %.1f)\n",
            rrs$est, rrs$p_CR2, rrs$df))

cat("\n=== 3. do the two scales tell the same story? ===\n")
smd <- function(sub, mods = ~ 1, idx = 1) {
  V <- mkV(sub)
  m <- rma.mv(yi_al, V, mods = mods, data = sub,
              random = list(~ 1 | study, ~ 1 | effect_ID),
              method = "REML", test = "t", dfs = "contain")
  ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = sub$study)
  c(est = m$b[idx], p = ct$p_Satt[idx])
}
cmp <- data.frame()
for (nm in c("IS", "RT")) {
  sub_raw <- if (nm == "IS") dis else drt
  sub_smd <- d[d$outcome == nm & (nm == "RT" | d$ctrl_unit != "RAU"), ]
  sub_smd$SNR_between <- ave(sub_smd$SNR, sub_smd$study, FUN = mean)
  sub_smd$SNR_within  <- sub_smd$SNR - sub_smd$SNR_between
  sub_smd$SNR_between_c <- sub_smd$SNR_between - mean(sub_smd$SNR_between)
  a <- smd(sub_smd); b <- smd(sub_smd, ~ SNR_within + SNR_between_c, 2)
  ra <- fit(sub_raw); rbb <- fit(sub_raw, ~ SNR_within + SNR_between_c, 2)
  cmp <- rbind(cmp,
    data.frame(family = nm, quantity = "pooled effect",
               smd_est = a["est"], smd_p = a["p"], raw_est = ra$est, raw_p = ra$p_CR2),
    data.frame(family = nm, quantity = "SNR within-study",
               smd_est = b["est"], smd_p = b["p"], raw_est = rbb$est, raw_p = rbb$p_CR2))
}
cmp$same_verdict <- ifelse((cmp$smd_p < .05) == (cmp$raw_p < .05), "agree", "DISAGREE")
print(cmp, row.names = FALSE, digits = 4)

cat(sprintf("\nrank correlation of the two scales across IS effects: rho = %+.3f\n",
            cor(dis$yi, dis$yi_al, method = "spearman")))
cat(sprintf("rank correlation across RT effects:                    rho = %+.3f\n",
            cor(drt$yi, drt$yi_al, method = "spearman")))

cat("\n=== 4. per-study raw losses (for a forest-style display) ===\n")
per <- do.call(rbind, lapply(split(dis, dis$study), function(s) {
  if (nrow(s) == 1) {
    return(data.frame(study = s$study[1], k = 1L, points_lost = s$yi[1],
                      lo = s$yi[1] - 1.96 * sqrt(s$vi[1]),
                      hi = s$yi[1] + 1.96 * sqrt(s$vi[1])))
  }
  Vs <- mkV(s)   # assign first: rma.mv cannot resolve a call in the V argument
  m <- rma.mv(yi, Vs, data = s, random = ~ 1 | effect_ID, method = "REML")
  data.frame(study = s$study[1], k = nrow(s), points_lost = m$b[1],
             lo = m$ci.lb[1], hi = m$ci.ub[1])
}))
per <- per[order(-per$points_lost), ]
print(per, row.names = FALSE, digits = 3)
write.csv(per, file.path(OUT, "IS_per_study_points_lost.csv"), row.names = FALSE)
write.csv(cmp, file.path(OUT, "scale_agreement.csv"), row.names = FALSE)
write.csv(data.frame(
  scale = c("IS percentage points", "RT milliseconds"),
  k = c(nrow(dis), nrow(drt)),
  pooled = c(r$est, rr$est), ci_lo = c(r$lo, rr$lo), ci_hi = c(r$hi, rr$hi),
  pi_lo = c(r$pi_lo, rr$pi_lo), pi_hi = c(r$pi_hi, rr$pi_hi),
  p_CR2 = c(r$p_CR2, rr$p_CR2),
  snr_within = c(rs$est, rrs$est), snr_within_p = c(rs$p_CR2, rrs$p_CR2),
  snr_between = c(rb$est, NA), snr_between_p = c(rb$p_CR2, NA)),
  file.path(OUT, "raw_scale_headlines.csv"), row.names = FALSE)

cat("\nwritten to", OUT, "\n")
sink()
