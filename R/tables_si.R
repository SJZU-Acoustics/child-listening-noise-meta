
source(file.path("R", "helpers.R"))
OUTC <- function(a, f) file.path(".cache", "analysis", a, "outcome", f)
TDIR <- file.path("output", "tables", "si")
dir.create(TDIR, recursive = TRUE, showWarnings = FALSE)

D <- p26_data()
d <- D$effects
fmt <- function(x, k = 2) formatC(x, format = "f", digits = k)
num <- function(x) gsub("-", "$-$", x)
fp <- function(p) ifelse(p < 0.001, "$<$0.001", fmt(p, 3))
lab_study <- function(s) {
  s <- gsub("&", "\\\\&", sub(",(\\d)", ", \\1", s))
  s <- gsub("\u00e4", "{\\\"a}", s, fixed = TRUE)
  gsub("\u00f6", "{\\\"o}", s, fixed = TRUE)
}

tex_table <- function(file, header, body_lines, colspec) {
  con <- file(file.path(TDIR, file), "w")
  writeLines(c(
    sprintf("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}%s@{}}",
            colspec),
    "\\toprule", header, "\\midrule", body_lines,
    "\\bottomrule", "\\end{tabular*}"), con)
  close(con)
}

cells <- read.csv(OUTC("analysis_18_inference_layers", "cell_aggregates.csv"))
cells <- cells[order(cells$study, cells$outcome), ]
age_codes <- c(reported = "R", computed = "C", estimated = "E",
               "group label" = "G", "recruited-sample proxy" = "P")
coding <- lapply(seq_len(nrow(cells)), function(i) {
  z <- cells[i, ]
  x <- d[d$study == z$study & d$outcome == z$outcome, ]
  stopifnot(nrow(x) == z$k_cell, !anyNA(age_codes[x$age_basis]))
  ages <- range(x$age)
  data.frame(study = z$study, outcome = z$outcome, k = nrow(x),
    age = if (diff(ages) == 0) fmt(ages[1], 2) else paste(fmt(ages, 2), collapse = "--"),
    age_basis = paste(sort(unique(age_codes[x$age_basis])), collapse = "/"),
    task = paste(sort(unique(x$task)), collapse = "; "),
    masker = paste(sort(unique(x$noise_type)), collapse = ", "),
    masking = paste(sort(unique(as.character(x$masking))), collapse = "/"),
    estimate = z$yi, lower = z$yi - 1.96 * sqrt(z$vi),
    upper = z$yi + 1.96 * sqrt(z$vi))
})
coding <- do.call(rbind, coding)
stopifnot(nrow(coding) == 39, sum(coding$k) == 269)
write.csv(coding, file.path("output", "numbers", "study_outcome_coding.csv"),
          row.names = FALSE)
body <- sprintf("%s & %s & %d & %s (%s) & %s & %s & %s & %s [%s, %s] \\\\",
  lab_study(coding$study), coding$outcome, coding$k, coding$age, coding$age_basis,
  coding$task, coding$masker, coding$masking, num(fmt(coding$estimate)),
  num(fmt(coding$lower)), num(fmt(coding$upper)))
writeLines(body, file.path(TDIR, "si_table_cells.tex"))

oh <- read.csv(file.path("output", "numbers", "ohat_domain_matrix.csv"),
               check.names = FALSE)
oh[-(1:2)] <- lapply(oh[-(1:2)], function(x) ifelse(x == "--", "$-$$-$", sub("^-$", "$-$", x)))
body <- sprintf("%s & %s & %s & %s & %s & %s & %s & %s \\\\",
                lab_study(gsub(", ", ",", oh$Study)),
                oh$`Selection bias`, oh$`Confounding bias`,
                oh$`Performance bias`, oh$`Attrition/exclusion bias`,
                oh$`Detection bias`, oh$`Selective reporting bias`,
                oh$`Other sources of bias`)
tex_table("si_table_ohat.tex",
  "Study & Selection & Confounding & Performance & Attrition & Detection & Reporting & Other \\\\",
  body, "llllllll")

loo <- read.csv(OUTC("analysis_04_snr_information", "leave_one_out_IS.csv"))
loo <- loo[order(loo$slope), ]
body <- sprintf("%s & %d & %s & %s \\\\",
                lab_study(sub("^\\(none\\)$", "None (all studies)", loo$dropped)), loo$k,
                num(fmt(loo$slope, 3)), fp(loo$p_CR2))
tex_table("si_table_loso.tex",
  "Study dropped & $k$ & SNR slope ($g$/dB) & $p$ \\\\",
  body, "lrll")

sc <- read.csv(OUTC("analysis_13_specification_curve", "specification_curve.csv"))
sc$cellv <- sprintf("%s (%s)", num(fmt(sc$estimate)), fp(sc$p_CR2))
wide <- reshape(sc[, c("quantity", "working", "variance", "holt", "study13",
                       "cellv")],
                idvar = c("working", "variance", "holt", "study13"),
                timevar = "quantity", direction = "wide")
names(wide) <- sub("cellv.", "", names(wide))
wide$working <- c("shared-arm" = "Shared arm", "diagonal" = "Diagonal",
                  "rho=0.5" = "$\\rho = 0.5$")[wide$working]
wide$variance <- sub("^paired r=(.*)$", "Paired, $r = \\1$",
                     sub("^delivered$", "Independent", wide$variance))
wide$holt <- c(kept = "Retained", dropped = "Excluded")[wide$holt]
wide$study13 <- c(merged = "One study", split = "Two samples")[wide$study13]
stopifnot(!anyNA(wide[c("working", "variance", "holt", "study13")]))
qcols <- c("pooled ALL", "pooled IS", "pooled RT", "pooled PE",
           "SNR within ALL", "SNR within IS")
stopifnot(all(qcols %in% names(wide)), !anyNA(wide[qcols]))
body <- apply(wide, 1, function(r) {
  paste(c(r[["working"]], r[["variance"]], r[["holt"]], r[["study13"]],
          r[qcols]), collapse = " & ")
})
body <- paste0(body, " \\\\")
tex_table("si_table_spec.tex",
  "Covariance & Variance & Holt et al., 2020 & Prodi et al., 2019a & Pooled ALL & Pooled IS & Pooled RT & Pooled PE & Slope ALL & Slope IS \\\\",
  body, "llllllllll")

rt <- read.csv(OUTC("analysis_03_r_table", "r_table_headlines.csv"))
rt$cellv <- sprintf("%s (%s)", num(fmt(rt$estimate)), fp(rt$p_CR2))
widr <- reshape(rt[, c("quantity", "r", "cellv")], idvar = "quantity",
                timevar = "r", direction = "wide")
names(widr) <- sub("cellv.", "r = ", names(widr))
ord <- c("empty ALL", "empty IS", "empty RT", "empty PE",
         "SNR slope IS", "SNR slope RT", "SNR slope PE")
widr <- widr[match(ord, widr$quantity), ]
widr$quantity <- sub("^empty", "Pooled", widr$quantity)
rcols <- c("r = delivered", "r = 0.0", "r = 0.3", "r = 0.5", "r = 0.7",
           "r = 0.9")
body <- apply(widr, 1, function(r) paste(c(r[["quantity"]], r[rcols]),
                                         collapse = " & "))
body <- paste0(body, " \\\\")
tex_table("si_table_rtable.tex",
  "Quantity & Independent (primary) & $r=0$ & $r=0.3$ & $r=0.5$ & $r=0.7$ & $r=0.9$ \\\\",
  body, "lllllll")

cwbw <- read.csv(OUTC("analysis_18_inference_layers", "cwb_within_slopes.csv"))
cwbb <- read.csv(OUTC("analysis_20_exit_consolidation",
                      "cwb_between_slopes.csv"))
cellp <- read.csv(OUTC("analysis_18_inference_layers",
                       "cell_level_pooled.csv"))
reuse <- read.csv(OUTC("analysis_19_reuse_effort", "reuse_resamples.csv"))
body <- c(
  sprintf("Within-study SNR slope, %s & %s & CR2 %s; CWB %s \\\\",
          cwbw$family, num(fmt(cwbw$est, 3)), fp(cwbw$p_CR2),
          ifelse(cwbw$p_CWB < 0.001, "$<$0.001", fmt(cwbw$p_CWB, 3))),
  sprintf("Between-study SNR slope, %s & %s & CR2 %s; CWB %s \\\\",
          cwbb$family, num(fmt(cwbb$est, 3)), fp(cwbb$p_CR2),
          fmt(cwbb$p_CWB, 3)),
  sprintf("Pooled $g$ after aggregation, %s & %s [%s, %s] & CR2 %s \\\\",
          cellp$family, num(fmt(cellp$est)), num(fmt(cellp$lo)),
          fmt(cellp$hi), fp(cellp$p_CR2)),
  sprintf("Pooled $g$, one effect per control arm, ALL & %s to %s & Range over 500 resamples \\\\",
          fmt(min(reuse$ALL)), fmt(max(reuse$ALL))),
  sprintf("Pooled $g$, one effect per control arm, IS & %s to %s & Range over 500 resamples \\\\",
          fmt(min(reuse$IS)), fmt(max(reuse$IS))))
tex_table("si_table_inference.tex",
  "Quantity & Estimate & Inference \\\\",
  body, "lll")

snrd <- read.csv(OUTC("analysis_15_pairs_deepened", "snr_decomposition_cwb.csv"))
diffm <- read.csv(OUTC("analysis_15_pairs_deepened",
                       "difficulty_moderation_cwb.csv"))
rowfmt <- function(df, labels) sprintf("%s & %s & %s & %s \\\\", labels,
                                       num(fmt(df$est, 3)), fp(df$p_model),
                                       ifelse(is.na(df$p_CWB), "--",
                                              fmt(df$p_CWB, 3)))
body <- c(rowfmt(snrd, c("Intercept", "SNR, within-study",
                         "SNR, between-study")),
          rowfmt(diffm, c("Intercept", "Baseline accuracy, within-study",
                          "Baseline accuracy, between-study")))
tex_table("si_table_pairs.tex",
  "Term & Estimate & Model $p$ & CWB $p$ \\\\",
  body, "llll")

r12 <- read.csv(OUTC("analysis_14_variability", "family_R12_fdr.csv"))
names(r12) <- c("member", "p", "q_BH")
plv <- read.csv(OUTC("analysis_14_variability", "pooled_lnVR.csv"))
dec <- read.csv(OUTC("analysis_14_variability", "lnVR_SNR_decomposition.csv"))
art <- read.csv(OUTC("analysis_14_variability", "IS_artefact_model.csv"))
disp_row <- function(m) {
  if (grepl("^pooled lnVR ", m)) {
    f <- sub("^pooled lnVR ", "", m); r <- plv[plv$family == f, ]
    return(c(sprintf("Pooled lnVR, %s", f),
             sprintf("%s [%s, %s]", num(fmt(r$est)), num(fmt(r$ci_lo)),
                     num(fmt(r$ci_hi)))))
  }
  if (grepl("^lnVR SNR ", m)) {
    f <- sub("^lnVR SNR (\\w+) .*$", "\\1", m)
    tm <- sub("^.* ", "", m)
    r <- dec[dec$family == f & dec$term == tm, ]
    side <- c(SNR_within = "Within-study", SNR_between_c = "Between-study")[[tm]]
    return(c(sprintf("%s lnVR SNR slope, %s", side, f), num(fmt(r$est, 3))))
  }
  tm <- sub("^IS artefact model: ", "", m)
  r <- art[art$term == tm, ]
  lab <- c(base_c = "IS lnVR model, baseline accuracy",
           drop_c = "IS lnVR model, raw accuracy loss")[[tm]]
  c(lab, num(fmt(r$est, 3)))
}
dl <- t(sapply(r12$member, disp_row))
ord <- order(match(sub(",.*", "", sub(" lnVR SNR slope", "", dl[, 1])),
                   c("Pooled lnVR", "Within-study", "Between-study",
                     "IS lnVR model")),
             match(sub("^.*, ", "", dl[, 1]), c("ALL", "IS", "RT", "PE")),
             dl[, 1])
body <- sprintf("%s & %s & %s & %s \\\\", dl[ord, 1], dl[ord, 2],
                fp(r12$p[ord]), fmt(r12$q_BH[ord], 3))
tex_table("si_table_dispersion.tex",
  "Test & Estimate & $p$ (CR2) & $q$ (BH) \\\\",
  body, "llll")

fam <- read.csv(OUTC("analysis_20_exit_consolidation", "families_at_exit.csv"))
FAMILY_NAMES <- c(R1 = "Pooled effects", R2 = "SNR slopes",
                  R3 = "Moderators of the pooled effect",
                  R4 = "Aligned-pair models", R5 = "Small-study effects",
                  R6 = "Dose--response curvature",
                  R7 = "Control-arm baseline", R8 = "Age",
                  R9 = "Outcome type", R10 = "Raw-scale estimates",
                  R12 = "Dispersion")
stopifnot(setequal(unique(fam$family), names(FAMILY_NAMES)))
MEMBER_RULES <- c(
  "^pooled g (\\w+)$" = "Pooled $g$, \\1",
  "^(within|between) (\\w+)$" = "\\1-study SNR slope, \\2",
  "^age main$" = "Mean age",
  "^outcome (RT|PE) vs IS$" = "Outcome type, \\1 vs IS",
  "^SNR x age$" = "SNR $\\\\times$ age",
  "^pairs SNR (within|between)$" = "\\1-study SNR",
  "^dual vs same task$" = "Cross-task vs same-task pairs",
  "^pairs IS-baseline (within|between)$" = "\\1-study baseline accuracy",
  "^Egger-n/PET (\\w+)$" = "Sample-size Egger test, \\1",
  "^PEESE (\\w+)$" = "PEESE small-study term, \\1",
  "^sample size between-study$" = "Between-study sample size",
  "^SNR\\^2 within (\\w+)$" = "Within-study SNR$^2$, \\1",
  "^IS raw drop$" = "IS raw accuracy loss",
  "^IS control SD \\(artefact term\\)$" = "IS control-arm SD",
  "^(IS|RT) baseline (within|between)$" = "\\1 baseline, \\2-study",
  "^age (within|between) (\\w+)$" = "\\1-study age, \\2",
  "^band (\\d+)-(\\d+) (\\w+)$" = "Age band \\1--\\2 years, \\3",
  "^age\\^2 within (\\w+)$" = "Within-study age$^2$, \\1",
  "^SNR_within x age_within$" = "Within-study SNR $\\\\times$ within-study age",
  "^(RT|PE) vs IS \\(RE\\)$" = "\\1 vs IS, random effects",
  "^(RT|PE) vs IS \\(study FE\\)$" = "\\1 vs IS, study fixed effects",
  "^IS pooled points lost$" = "IS pooled loss (percentage points)",
  "^RT pooled ms$" = "RT pooled increase (ms)",
  "^IS SNR (within|between) \\(pts/dB\\)$" = "IS \\1-study SNR slope (points/dB)",
  "^RT SNR within \\(ms/dB\\)$" = "RT within-study SNR slope (ms/dB)",
  "^lnVR (ALL|IS|RT|PE)$" = "Pooled lnVR, \\1",
  "^lnVR SNR (within|between) (\\w+)$" = "\\1-study lnVR SNR slope, \\2",
  "^lnVR artefact: baseline$" = "IS lnVR model, baseline accuracy",
  "^lnVR artefact: raw drop$" = "IS lnVR model, raw accuracy loss")
member_label <- function(m) {
  hit <- names(MEMBER_RULES)[vapply(names(MEMBER_RULES), grepl, TRUE, x = m)]
  stopifnot(length(hit) == 1)
  lab <- sub(hit, MEMBER_RULES[[hit]], m)
  paste0(toupper(substr(lab, 1, 1)), substr(lab, 2, nchar(lab)))
}
fam$label <- vapply(fam$member, member_label, "")
fam <- fam[order(match(fam$family, names(FAMILY_NAMES)), fam$q, fam$p), ]
body <- unlist(lapply(names(FAMILY_NAMES), function(f) {
  r <- fam[fam$family == f, ]
  c(sprintf("%s\\textit{%s} (%d tests) & & & \\\\",
            if (f == "R1") "" else "\\addlinespace ", FAMILY_NAMES[[f]],
            nrow(r)),
    sprintf("%s & %s & %s & %s \\\\", r$label, r$source,
            ifelse(r$p < 0.001, "$<$0.001", fmt(r$p, 3)),
            ifelse(r$q < 0.001, "$<$0.001", fmt(r$q, 3))))
}))
writeLines(body, file.path(TDIR, "si_table_families.tex"))

sc22 <- read.csv(OUTC("analysis_22_scope_sensitivity", "scope_sensitivity.csv"))
scopes <- c("base", "condition_subset", "no_composite_PE", "both")
cell22 <- function(q, s) {
  r <- sc22[sc22$quantity == q & sc22$scope == s, ]
  if (q == "IS-RT pair r") return(sprintf("$+$%s (%d)", fmt(r$estimate), r$k))
  if (q == "pairs adverse on both") return(sprintf("%d/%d", r$estimate, r$k))
  sprintf("%s (%s)", num(fmt(r$estimate, if (grepl("SNR", q)) 3 else 2)), fp(r$p_CR2))
}
k22 <- sapply(scopes, function(s) {
  r <- sc22[sc22$quantity == "pooled ALL" & sc22$scope == s, ]
  sprintf("%d (%d)", r$k, r$studies)
})
qlab <- c("pooled ALL" = "Pooled $g$, ALL", "pooled IS" = "Pooled $g$, IS",
          "pooled RT" = "Pooled $g$, RT", "pooled PE" = "Pooled $g$, PE",
          "SNR within ALL" = "Within-study SNR slope, ALL",
          "SNR within IS" = "Within-study SNR slope, IS",
          "IS-RT pair r" = "IS--RT pair $r$ (pairs)",
          "pairs adverse on both" = "Pairs adverse on both outcomes")
body <- c(sprintf("Effect sizes (studies) & %s \\\\", paste(k22, collapse = " & ")),
          sapply(names(qlab), function(q)
            sprintf("%s & %s \\\\", qlab[[q]],
                    paste(sapply(scopes, function(s) cell22(q, s)), collapse = " & "))))
tex_table("si_table_scope.tex",
  "Quantity & \\makecell[l]{All\\\\contrasts} & \\makecell[l]{Narrower\\\\condition set} & \\makecell[l]{Without\\\\composite PE} & \\makecell[l]{Both\\\\restrictions} \\\\",
  body, "lllll")

cat("SI tables written\n")

for (f in list.files(TDIR, pattern = "\\.tex$", full.names = TRUE)) {
  txt <- readLines(f, warn = FALSE)
  txt <- gsub("\\bIS\\b", "ACC", txt, perl = TRUE)
  writeLines(txt, f)
}
