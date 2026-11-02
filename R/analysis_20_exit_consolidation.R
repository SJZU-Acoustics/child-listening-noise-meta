
source(file.path("R", "helpers.R"))
OUT <- file.path(p26_outdir("analysis_20_exit_consolidation"), "outcome")

D <- p26_data()
d <- D$effects
d <- d[order(d$study, d$effect_ID), ]

sink(file.path(OUT, "console.txt"), split = TRUE)
cat("A20 - exit consolidation: R2/R4 family decision + families at exit\n\n")

mkV <- function(sub) vcalc(vi = vi, cluster = study, grp1 = experiment_group_id,
                           grp2 = control_group_id, w1 = n1, w2 = n2,
                           data = sub, sparse = FALSE)

clamp_cwb <- function(p) pmax(p, 1 / BOOT_R)

cat("=== 1. cluster wild bootstrap on the BETWEEN-study SNR slopes ===\n")
cat("(completing A18's layer: same decomposition model, shared-arm V,\n")
cat(sprintf(" %d resamples, seed 26; constraint = SNR_between_c)\n\n", BOOT_R))

btw_rows <- list()
for (fam in c("ALL", "IS", "RT", "PE")) {
  sub <- if (fam == "ALL") d else d[d$outcome == fam, ]
  sub$SNR_between <- ave(sub$SNR, sub$study, FUN = mean)
  sub$SNR_within  <- sub$SNR - sub$SNR_between
  sub$SNR_between_c <- sub$SNR_between - mean(sub$SNR_between)
  V <- mkV(sub)
  m <- rma.mv(yi_al, V, mods = ~ SNR_within + SNR_between_c, data = sub,
              random = list(~ 1 | study, ~ 1 | effect_ID),
              method = "REML", test = "t", dfs = "contain")
  ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = sub$study)
  bt <- try(Wald_test_cwb(full_model = m, constraints = constrain_zero(3),
                          R = BOOT_R, seed = 26), silent = TRUE)
  p_cwb <- if (inherits(bt, "try-error")) NA_real_ else bt$p_val
  cat(sprintf("%-4s between-slope %+.4f: model p %.4g | CR2 p %.4g (df %.1f) | CWB p %s\n",
              fam, m$b[3], m$pval[3], ct$p_Satt[3], ct$df_Satt[3],
              ifelse(is.na(p_cwb), "FAILED", sprintf("%.4f", p_cwb))))
  btw_rows[[fam]] <- data.frame(family = fam, est = m$b[3],
                                p_model = m$pval[3], p_CR2 = ct$p_Satt[3],
                                df_CR2 = ct$df_Satt[3], p_CWB = p_cwb)
}
btw_res <- do.call(rbind, btw_rows)
cat("PE within-slope: not CWB-testable (1 manipulating study, A04) --\n")
cat("excluded from the family on estimability, reported descriptively.\n")

OC <- function(a, f) read.csv(file.path(P26_ROOT, ".cache", "analysis", a, "outcome", f),
                              check.names = FALSE)
pick <- function(tab, col, ...) {
  key <- list(...)
  i <- Reduce(`&`, Map(function(k, v) tab[[k]] == v, names(key), key))
  stopifnot(sum(i) == 1)
  tab[[col]][i]
}
cwb18 <- OC("analysis_18_inference_layers", "cwb_within_slopes.csv")
wm02  <- OC("analysis_02_shared_control", "working_model_comparison.csv")
leg01 <- OC("analysis_01_legacy_reproduction", "legacy_reproduction.csv")
dec04 <- OC("analysis_04_snr_information", "within_between_decomposition.csv")
eg06  <- OC("analysis_06_publication_bias", "egger_tests.csv")
cv07  <- OC("analysis_07_dose_response", "curvature_tests.csv")
bs08  <- OC("analysis_08_baseline", "IS_baseline_moderation.csv")
ar08  <- OC("analysis_08_baseline", "IS_scale_artefact_test.csv")
rt08  <- OC("analysis_08_baseline", "RT_baseline_moderation.csv")
ag09  <- OC("analysis_09_age", "age_decomposition.csv")
ex09  <- OC("analysis_09_age", "extra_tests.csv")
ot10  <- OC("analysis_10_moderator_estimability", "outcome_type_within_study.csv")
rw11  <- OC("analysis_11_raw_scale", "raw_scale_headlines.csv")
dd12  <- OC("analysis_12_design_descriptors", "design_descriptors.csv")
vr14  <- OC("analysis_14_variability", "family_R12_fdr.csv")
sn15  <- OC("analysis_15_pairs_deepened", "snr_decomposition_cwb.csv")
df15  <- OC("analysis_15_pairs_deepened", "difficulty_moderation_cwb.csv")
pp17  <- OC("analysis_17_small_study", "pet_peese.csv")
wm <- function(q) pick(wm02, "p_CR2", quantity = q, working_model = "shared-arm")
lg <- function(q) pick(leg01, "p_robust", quantity = q)
eg <- function(f) pick(eg06, "p_CR2", family = f, variant = "sample-size (sqrt(1/n))")
ag <- function(f, t) pick(ag09, "p_CR2", family = f, term = t)
vr <- function(q) pick(vr14, "p", q = q)

cat("\n=== 2. family R2 under each defensible reading (BH q) ===\n\n")

r2 <- list()

pa <- data.frame(reading = "A: pooled-3, CR2",
                 member = c("SNR slope IS", "SNR slope RT", "SNR slope PE"),
                 p = sapply(c("SNR slope IS", "SNR slope RT", "SNR slope PE"), wm))
pa$q <- p.adjust(pa$p, "BH")
r2[["A"]] <- pa

pb <- data.frame(reading = "B: decomposition-7, CR2",
                 member = c("within ALL", "within IS", "within RT",
                            "between ALL", "between IS", "between RT", "between PE"),
                 p = c(cwb18$p_CR2[match(c("ALL", "IS", "RT"), cwb18$family)],
                       btw_res$p_CR2))
pb$q <- p.adjust(pb$p, "BH")
r2[["B"]] <- pb

pc <- data.frame(reading = "C: decomposition-7, CWB",
                 member = pb$member,
                 p = clamp_cwb(c(cwb18$p_CWB[match(c("ALL", "IS", "RT"), cwb18$family)],
                                 btw_res$p_CWB)))
pc$q <- p.adjust(pc$p, "BH")
r2[["C"]] <- pc

pb8 <- data.frame(reading = "B8: decomposition-8 incl PE-within, CR2",
                  member = c(pb$member, "within PE (1 study)"),
                  p = c(pb$p, pick(dec04, "p_CR2", family = "PE", term = "SNR_within")))
pb8$q <- p.adjust(pb8$p, "BH")
r2[["B8"]] <- pb8

r2_tab <- do.call(rbind, r2)
print(r2_tab, row.names = FALSE, digits = 3)

cat("\n=== 3. families at exit (BH within each; ledger-recorded p) ===\n\n")

fam_list <- list(
  R1 = data.frame(member = c("pooled g ALL", "pooled g IS", "pooled g RT", "pooled g PE"),
                  source = "A01/A02",
                  p = sapply(c("empty ALL", "empty IS", "empty RT", "empty PE"), wm)),
  R3 = data.frame(member = c("SNR x age", "age main", "outcome RT vs IS", "outcome PE vs IS"),
                  source = "A01",
                  p = sapply(c("SNR_c:age_c", "age_c", "Outcome_TypeRT", "Outcome_TypePE"), lg)),
  R4 = data.frame(member = c("pairs SNR within", "pairs SNR between", "dual vs same task",
                             "pairs IS-baseline within", "pairs IS-baseline between"),
                  source = "A15/A01",
                  p = c(pick(sn15, "p_CWB", term = "SNR_within"),
                        pick(sn15, "p_CWB", term = "SNR_between_c"),
                        lg("RQ2 dual-task vs same-task"),
                        pick(df15, "p_CWB", term = "base_within"),
                        pick(df15, "p_CWB", term = "base_between_c"))),
  R5 = data.frame(member = c("Egger-n/PET ALL", "Egger-n/PET IS", "Egger-n/PET RT", "Egger-n/PET PE",
                             "PEESE ALL", "PEESE IS", "PEESE RT", "PEESE PE",
                             "sample size between-study"),
                  source = "A06/A17/A12",
                  p = c(sapply(c("ALL", "IS", "RT", "PE"), eg),
                        pp17$peese_slope_p[match(c("ALL", "IS", "RT", "PE"), pp17$family)],
                        pick(dd12, "p_CR2", term = "n_between_c"))),
  R6 = data.frame(member = c("SNR^2 within ALL", "SNR^2 within IS"),
                  source = "A07", p = cv07$p_CR2[match(c("ALL", "IS"), cv07$family)]),
  R7 = data.frame(member = c("IS baseline within", "IS baseline between", "IS raw drop",
                             "IS control SD (artefact term)", "RT baseline within", "RT baseline between"),
                  source = "A08",
                  p = c(pick(bs08, "p_CR2", term = "base_within"),
                        pick(bs08, "p_CR2", term = "base_between_c"),
                        pick(ar08, "p_CR2", term = "drop_c"),
                        pick(ar08, "p_CR2", term = "sd_c"),
                        pick(rt08, "p_CR2", term = "base_within"),
                        pick(rt08, "p_CR2", term = "base_between_c"))),
  R8 = data.frame(member = c("age within ALL", "age within IS", "age within RT", "age within PE",
                             "age between ALL", "age between IS", "age between RT", "age between PE",
                             "age^2 within ALL", "age^2 within IS",
                             "band 8-10 ALL", "band 8-10 IS", "band 10-13 ALL", "band 10-13 IS",
                             "SNR_within x age_within"),
                  source = "A09",
                  p = c(sapply(c("ALL", "IS", "RT", "PE"), ag, t = "age_within"),
                        sapply(c("ALL", "IS", "RT", "PE"), ag, t = "age_between_c"),
                        ex09$p_CR2[match(c("age^2 within ALL", "age^2 within IS",
                                           "band 8-10 ALL", "band 8-10 IS",
                                           "band 10-13 ALL", "band 10-13 IS",
                                           "SNR_within x age_within"), ex09$test)])),
  R9 = data.frame(member = c("RT vs IS (RE)", "RT vs IS (study FE)", "PE vs IS (RE)", "PE vs IS (study FE)"),
                  source = "A10",
                  p = c(pick(ot10, "p_CR2", model = "random effects", term = "Outcome_TypeRT"),
                        pick(ot10, "p_CR2", model = "study fixed effects", term = "Outcome_TypeRT"),
                        pick(ot10, "p_CR2", model = "random effects", term = "Outcome_TypePE"),
                        pick(ot10, "p_CR2", model = "study fixed effects", term = "Outcome_TypePE"))),
  R10 = data.frame(member = c("IS pooled points lost", "IS SNR within (pts/dB)", "IS SNR between (pts/dB)",
                              "RT pooled ms", "RT SNR within (ms/dB)"),
                   source = "A11",
                   p = c(rw11$p_CR2[1], rw11$snr_within_p[1], rw11$snr_between_p[1],
                         rw11$p_CR2[2], rw11$snr_within_p[2])),
  R12 = data.frame(member = c("lnVR ALL", "lnVR IS", "lnVR RT", "lnVR PE",
                              "lnVR SNR within ALL", "lnVR SNR within IS", "lnVR SNR within RT",
                              "lnVR SNR between ALL", "lnVR SNR between IS", "lnVR SNR between RT",
                              "lnVR artefact: raw drop", "lnVR artefact: baseline"),
                   source = "A14",
                   p = sapply(c("pooled lnVR ALL", "pooled lnVR IS", "pooled lnVR RT", "pooled lnVR PE",
                                "lnVR SNR ALL SNR_within", "lnVR SNR IS SNR_within", "lnVR SNR RT SNR_within",
                                "lnVR SNR ALL SNR_between_c", "lnVR SNR IS SNR_between_c",
                                "lnVR SNR RT SNR_between_c",
                                "IS artefact model: drop_c", "IS artefact model: base_c"), vr))
)
stopifnot(!anyNA(unlist(lapply(fam_list, `[[`, "p"))))

r2c <- r2[["C"]]
fam_list$R2 <- data.frame(member = r2c$member, source = "A04/A18/A20", p = r2c$p)

fam_tab <- do.call(rbind, lapply(names(fam_list), function(f) {
  x <- fam_list[[f]]
  x$family <- f
  x$n_family <- nrow(x)
  x$q <- p.adjust(x$p, "BH")
  x$sig_q05 <- x$q < 0.05
  x[, c("family", "n_family", "member", "source", "p", "q", "sig_q05")]
}))
fam_tab <- fam_tab[order(fam_tab$family, fam_tab$q), ]
print(fam_tab, row.names = FALSE, digits = 3)

cat("\nsurvivors at q < 0.05, by family:\n")
sv <- fam_tab[fam_tab$sig_q05, c("family", "member", "q")]
if (nrow(sv)) print(sv, row.names = FALSE, digits = 3) else cat("(none)\n")

write.csv(btw_res, file.path(OUT, "cwb_between_slopes.csv"), row.names = FALSE)
write.csv(r2_tab, file.path(OUT, "r2_decision.csv"), row.names = FALSE)
write.csv(fam_tab, file.path(OUT, "families_at_exit.csv"), row.names = FALSE)
cat("\nwritten to", OUT, "\n")
sink()
