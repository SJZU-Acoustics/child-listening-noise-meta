
source(file.path("R", "helpers.R"))
OUT <- file.path(p26_outdir("analysis_19_reuse_effort"), "outcome")

D <- p26_data()
d <- D$effects
d <- d[order(d$study, d$effect_ID), ]

sink(file.path(OUT, "console.txt"), split = TRUE)
cat("A19 - reuse sensitivity and the effort side\n\n")

cat("=== 1. one-effect-per-control-arm resampling (B = 500) ===\n")
key <- paste(d$study, d$control_group_id)
arm_counts <- table(key)
cat(sprintf("distinct control arms: %d, reused: %d (up to %dx); effects on reused arms: %d\n",
            length(arm_counts), sum(arm_counts > 1), max(arm_counts),
            sum(arm_counts[key] > 1)))
set.seed(26)
B <- 500
draw <- replicate(B, {
  keep <- unlist(lapply(split(seq_len(nrow(d)), key), function(ix)
    if (length(ix) == 1) ix else sample(ix, 1)))
  sub <- d[sort(keep), ]
  out <- c(ALL = NA_real_, IS = NA_real_)
  m <- try(rma.mv(yi_al, vi, random = list(~ 1 | study, ~ 1 | effect_ID),
                  data = sub, method = "REML", test = "t", dfs = "contain"),
           silent = TRUE)
  if (!inherits(m, "try-error")) {
    ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = sub$study)
    out["ALL"] <- m$b[1]
  }
  subI <- sub[sub$outcome == "IS", ]
  mI <- try(rma.mv(yi_al, vi, random = list(~ 1 | study, ~ 1 | effect_ID),
                   data = subI, method = "REML", test = "t", dfs = "contain"),
            silent = TRUE)
  if (!inherits(mI, "try-error")) out["IS"] <- mI$b[1]
  out
})
draw <- t(draw)
A02REF <- local({
  w <- read.csv(file.path(P26_ROOT, ".cache", "analysis", "analysis_02_shared_control",
                          "outcome", "working_model_comparison.csv"))
  w <- w[w$working_model == "shared-arm" & grepl("^empty", w$quantity), ]
  setNames(w$estimate, sub("empty ", "", w$quantity))
})
cat(sprintf("ALL: median %+.3f, range %+.3f to %+.3f (adopted: %+.3f)\n",
            median(draw[, "ALL"], na.rm = TRUE), min(draw[, "ALL"], na.rm = TRUE),
            max(draw[, "ALL"], na.rm = TRUE), A02REF["ALL"]))
cat(sprintf("IS : median %+.3f, range %+.3f to %+.3f (adopted: %+.3f)\n",
            median(draw[, "IS"], na.rm = TRUE), min(draw[, "IS"], na.rm = TRUE),
            max(draw[, "IS"], na.rm = TRUE), A02REF["IS"]))
cat(sprintf("max deviation from adopted: ALL %.3f, IS %.3f\n",
            max(abs(draw[, "ALL"] - A02REF["ALL"]), na.rm = TRUE),
            max(abs(draw[, "IS"] - A02REF["IS"]), na.rm = TRUE)))

cat("\n=== 2. study-level PE profile and IS-PE concordance ===\n")
mkV <- function(sub) vcalc(vi = vi, cluster = study, grp1 = experiment_group_id,
                           grp2 = control_group_id, w1 = n1, w2 = n2,
                           data = sub, sparse = FALSE)
cell_fe <- function(sub) {
  Vc <- if (nrow(sub) == 1) matrix(sub$vi, 1, 1) else mkV(sub)
  prec <- solve(Vc); one <- rep(1, nrow(sub))
  c(beta = as.numeric((t(one) %*% prec %*% sub$yi_al) / (t(one) %*% prec %*% one)),
    se = sqrt(as.numeric(1 / (t(one) %*% prec %*% one))), k = nrow(sub))
}
pe <- d[d$outcome == "PE", ]
pe_prof <- do.call(rbind, lapply(split(pe, pe$study), function(s) {
  r <- cell_fe(s)
  data.frame(study = s$study[1], k = r["k"], g_PE = r["beta"], se = r["se"],
             n_adverse = sum(s$yi_al > 0))
}))
pe_prof <- pe_prof[order(-pe_prof$g_PE), ]
print(pe_prof, row.names = FALSE, digits = 3)
cat(sprintf("\nPE effects adverse in sign: %d of %d\n", sum(pe$yi_al > 0), nrow(pe)))

both <- intersect(unique(pe$study), unique(d$study[d$outcome == "IS"]))
cat(sprintf("\nstudies with both IS and PE: %d\n", length(both)))
conc <- do.call(rbind, lapply(both, function(s) {
  data.frame(study = s,
             g_IS = cell_fe(d[d$study == s & d$outcome == "IS", ])["beta"],
             g_PE = cell_fe(d[d$study == s & d$outcome == "PE", ])["beta"])
}))
conc <- conc[order(-conc$g_IS), ]
print(conc, row.names = FALSE, digits = 3)
cat(sprintf("\nsign concordance (both adverse): %d of %d\n",
            sum(conc$g_IS > 0 & conc$g_PE > 0), nrow(conc)))
cat(sprintf("rank correlation (Spearman, %d studies -- descriptive only): %+.2f\n",
            nrow(conc), cor(conc$g_IS, conc$g_PE, method = "spearman")))

write.csv(data.frame(ALL = draw[, "ALL"], IS = draw[, "IS"]),
          file.path(OUT, "reuse_resamples.csv"), row.names = FALSE)
write.csv(pe_prof, file.path(OUT, "PE_study_profile.csv"), row.names = FALSE)
write.csv(conc, file.path(OUT, "IS_PE_concordance.csv"), row.names = FALSE)
cat("\nwritten to", OUT, "\n")
sink()
