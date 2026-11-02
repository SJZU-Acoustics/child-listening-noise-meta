
source(file.path("R", "helpers.R"))

OUTC <- function(a, f) file.path(".cache", "analysis", a, "outcome", f)
TDIR <- file.path("output", "tables")

D <- p26_data()
d <- D$effects

fmt <- function(x, k = 2) formatC(x, format = "f", digits = k)
num <- function(x) gsub("-", "$-$", x)
tex_study <- function(s) {
  s <- gsub("&", "\\\\&", sub(",(\\d)", ", \\1", s))
  s <- gsub("\u00e4", "{\\\"a}", s, fixed = TRUE)
  gsub("\u00f6", "{\\\"o}", s, fixed = TRUE)
}

keys <- read.csv(file.path("output", "numbers", "corpus_bibkeys.csv"))
ohat <- read.csv(file.path("output", "numbers", "ohat_study_summary.csv"))
ohat$study_norm <- gsub(" ", "", ohat$Study)

rows <- lapply(split(d, d$study), function(s) {
  st <- s$study[1]
  cites <- keys$bibkey[keys$study_id == st]
  snr <- sort(unique(s$SNR))
  ages <- range(s$age)
  lr <- ohat$Low_risk_proportion[ohat$study_norm == gsub(" ", "", st)]
  data.frame(
    study = st,
    cite = paste0("\\cite{", paste(cites, collapse = ","), "}"),
    n = max(pmax(s$n1, s$n2)),
    age = if (diff(ages) < 0.05) fmt(ages[1], 1) else
      paste0(fmt(ages[1], 1), "--", fmt(ages[2], 1)),
    snr = {fs <- function(v) ifelse(v == 0, "0", sprintf("%+g", v))
      if (length(snr) == 1) fs(snr) else
      sprintf("%s to %s (%d)", fs(min(snr)), fs(max(snr)), length(snr))},
    noise = paste(sort(unique(s$noise_type)), collapse = ", "),
    masking = paste(sort(unique(as.character(s$masking))), collapse = "/"),
    outcomes = paste(intersect(c("IS", "RT", "PE"), unique(as.character(s$outcome))),
                     collapse = ", "),
    k = nrow(s),
    lowrisk = if (length(lr)) sprintf("%.1f", 100 * lr) else "--")
})
t1 <- do.call(rbind, rows)
t1 <- t1[order(t1$study), ]

n_reports <- nrow(read.csv(file.path(".cache", "data", "corpus_25_reports.csv")))
con <- file(file.path(TDIR, "table1_characteristics.tex"), "w")
writeLines(c(
  "\\begin{sidewaystable}[tbp]",
  sprintf("\\caption{Characteristics of the %d included studies (%d reports). $n$ is the largest analysed group, and $k$ is the number of extracted effect sizes (%d in total). SNR gives the experimental signal-to-noise ratio range, with the number of distinct levels in parentheses. Masking class is energetic (EM) or informational (IM), and EM/IM denotes conditions in both classes. SSN is speech-shaped noise and class denotes classroom noise. The final column gives the proportion of applicable OHAT risk-of-bias items rated low risk.}",
          nrow(t1), n_reports, nrow(d)),
  "\\label{tab:characteristics}",
  "\\footnotesize",
  "\\centering",
  "\\begin{tabular}{lrlllllrr}",
  "\\toprule",
  "Study & $n$ & Age (y) & SNR (dB) & Masker type & Masking class & Outcomes & $k$ & Low risk (\\%) \\\\",
  "\\midrule"), con)
for (i in seq_len(nrow(t1))) {
  r <- t1[i, ]
  lab <- tex_study(r$study)
  writeLines(sprintf("%s~%s & %d & %s & %s & %s & %s & %s & %d & %s \\\\",
                     lab, r$cite, r$n, r$age, num(r$snr), r$noise, r$masking,
                     r$outcomes, r$k, r$lowrisk), con)
}
writeLines(c("\\bottomrule", "\\end{tabular}", "\\end{sidewaystable}"), con)
close(con)

pis <- read.csv(OUTC("analysis_05_heterogeneity", "prediction_intervals.csv"))
rawh <- read.csv(OUTC("analysis_11_raw_scale", "raw_scale_headlines.csv"))
rom <- read.csv(OUTC("analysis_16_log_ratio", "rom_results.csv"))
fam20 <- read.csv(OUTC("analysis_20_exit_consolidation", "families_at_exit.csv"))
het <- read.csv(OUTC("analysis_05_heterogeneity",
                     "heterogeneity_decomposition.csv"))

q1 <- setNames(fam20$q[fam20$family == "R1"],
               sub("pooled g ", "", fam20$member[fam20$family == "R1"]))
kfx <- table(d$outcome)
kst <- sapply(split(d$study, d$outcome), function(x) length(unique(x)))
npart <- function(s) sum(sapply(split(pmax(s$n1, s$n2), s$study), max))
nfx <- sapply(split(d, d$outcome), npart)
i2 <- setNames(het$I2_total, het$family)

fmtq <- function(q) ifelse(q < 0.001, "$<$0.001", fmt(q, 3))
rompct <- function(fam, dir) {
  r <- rom[rom$family == fam, ]
  if (dir == "loss") {
    sprintf("%s [%s, %s]", fmt(100 * (1 - exp(-r$pooled_est)), 1),
            fmt(100 * (1 - exp(-r$pooled_lo)), 1),
            fmt(100 * (1 - exp(-r$pooled_hi)), 1))
  } else {
    sprintf("+%s [+%s, +%s]", fmt(100 * (exp(r$pooled_est) - 1), 1),
            fmt(100 * (exp(r$pooled_lo) - 1), 1),
            fmt(100 * (exp(r$pooled_hi) - 1), 1))
  }
}

t2row <- function(famname, label) {
  p <- pis[pis$family == famname, ]
  raw <- switch(famname,
    IS = sprintf("%s pp [%s, %s]", fmt(rawh$pooled[1], 1),
                 fmt(rawh$ci_lo[1], 1), fmt(rawh$ci_hi[1], 1)),
    RT = sprintf("+%s ms [+%s, +%s]", fmt(rawh$pooled[2], 0),
                 fmt(rawh$ci_lo[2], 0), fmt(rawh$ci_hi[2], 0)),
    "--")
  ks <- if (famname == "ALL") length(unique(d$study)) else kst[[famname]]
  ke <- if (famname == "ALL") nrow(d) else kfx[[famname]]
  nn <- if (famname == "ALL") npart(d) else nfx[[famname]]
  sprintf("%s & %d & %d & %d & %s [%s, %s] & %s & %s & [%s, %s] & %s \\\\",
          label, ks, ke, nn, fmt(p$estimate), fmt(p$ci_lo), fmt(p$ci_hi),
          fmtq(q1[[famname]]), fmt(i2[[famname]], 0),
          num(fmt(p$pi_lo)), fmt(p$pi_hi),
          if (raw == "--") raw else num(raw))
}

con <- file(file.path(TDIR, "table2_pooled.tex"), "w")
writeLines(c(
  "\\begin{table}[tbp]",
  "\\caption{Pooled adverse effects of noise relative to quiet controls. IS is task accuracy, RT response time and PE perceived effort. Positive $g$ indicates an adverse effect. Estimates come from three-level shared-arm models. The 95\\% CIs are model-based, and the $q$ values come from study-clustered CR2 tests corrected across ALL, IS, RT and PE. PIs are standardised 95\\% prediction intervals for a new condition in a new study. Raw changes are percentage points (pp) or milliseconds (ms). Here $m$ counts studies, $k$ effects and $N$ the largest analysed group for that outcome in each study, summed across studies, and $I^2$ includes both heterogeneity levels.}",
  "\\label{tab:pooled}",
  "\\footnotesize",
  "\\setlength{\\tabcolsep}{2.5pt}",
  "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lrrrlllll@{}}",
  "\\toprule",
  "Outcome & $m$ & $k$ & $N$ & \\shortstack{Hedges' $g$\\\\ {[95\\% CI]}} & $q$ & $I^2$ (\\%) & 95\\% PI & Raw scale [95\\% CI] \\\\",
  "\\midrule",
  t2row("ALL", "All"),
  t2row("IS", "IS"),
  t2row("RT", "RT"),
  t2row("PE", "PE"),
  "\\bottomrule", "\\end{tabular*}", "\\end{table}"), con)
close(con)

pp <- read.csv(OUTC("analysis_17_small_study", "pet_peese.csv"))
dd12 <- read.csv(OUTC("analysis_12_design_descriptors", "design_descriptors.csv"))
nb <- dd12[dd12$term == "n_between_c", ]
qnb <- fam20$q[fam20$member == "sample size between-study"]
con <- file(file.path(TDIR, "table3_smallstudy.tex"), "w")
writeLines(c(
  "\\begin{table}[tbp]",
  sprintf("\\caption{Small-study effects and PET--PEESE-adjusted estimates. The detection test is sample-size-based Egger regression, equivalent to the PET slope on the $\\sqrt{n}$ precision scale. This avoids the artefactual correlation between a standardised mean difference and its standard error. PET and PEESE columns report PET--PEESE-adjusted pooled estimates with 95\\%% CIs. A between-study sample-size coefficient was negative (%s per participant, $p = %s$), but not significant after correction within the nine-test small-study family ($q = %s$).}",
          num(fmt(nb$estimate, 3)), fmt(nb$p_CR2, 3), fmt(qnb, 2)),
  "\\label{tab:smallstudy}",
  "\\footnotesize",
  "\\setlength{\\tabcolsep}{3pt}",
  "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lllll@{}}",
  "\\toprule",
  "Outcome & Pooled $g$ & Detection $p$ & PET [95\\% CI] & PEESE [95\\% CI] \\\\",
  "\\midrule"), con)
labs <- c(ALL = "All outcomes", IS = "Task accuracy",
          RT = "Response time", PE = "Perceived effort")
for (f in c("ALL", "IS", "RT", "PE")) {
  r <- pp[pp$family == f, ]
  writeLines(sprintf("%s & %s & %s & %s [%s, %s] & %s [%s, %s] \\\\",
    labs[[f]], fmt(r$pooled), fmt(r$pet_slope_p, 3),
    num(fmt(r$pet_adj)), num(fmt(r$pet_adj_lo)), num(fmt(r$pet_adj_hi)),
    num(fmt(r$peese_adj)), num(fmt(r$peese_adj_lo)), num(fmt(r$peese_adj_hi))),
    con)
}
writeLines(c("\\bottomrule", "\\end{tabular*}", "\\end{table}"), con)
close(con)

wb04 <- read.csv(OUTC("analysis_04_snr_information", "within_between_decomposition.csv"))
wb04 <- wb04[wb04$term == "SNR_within", ]
ploss <- read.csv(OUTC("analysis_11_raw_scale", "IS_per_study_points_lost.csv"))
manifest <- data.frame(quantity = c(
  "studies", "reports", "effects", "IS effects", "RT effects", "PE effects",
  "pairs", "pair studies", "SNR min", "SNR max", "age min", "age max",
  "pooled g ALL", "pooled g IS", "pooled g RT", "pooled g PE",
  "IS points lost", "RT ms added", "IS pct of quiet lost", "RT pct slower",
  "ALL within-slope g/dB", "IS within-slope g/dB", "raw within-slope pp/dB",
  "studies varying SNR", "positive raw-loss studies",
  "children (largest analysed groups)", "children IS", "children RT",
  "children PE"),
  value = c(length(unique(d$study)), n_reports, nrow(d), kfx[["IS"]], kfx[["RT"]],
            kfx[["PE"]], nrow(D$pairs), length(unique(D$pairs$study)),
            min(d$SNR), max(d$SNR), min(d$age), max(d$age),
            fmt(pis$estimate[1]), fmt(pis$estimate[2]), fmt(pis$estimate[3]),
            fmt(pis$estimate[4]), fmt(rawh$pooled[1], 2), fmt(rawh$pooled[2], 0),
            fmt(100 * (1 - exp(-rom$pooled_est[rom$family == "IS"])), 1),
            fmt(100 * (exp(rom$pooled_est[rom$family == "RT"]) - 1), 1),
            fmt(wb04$estimate[wb04$family == "ALL"], 3),
            fmt(wb04$estimate[wb04$family == "IS"], 3), fmt(rawh$snr_within[1], 2),
            sum(tapply(d$SNR, d$study, function(x) length(unique(x))) > 1),
            sum(ploss$points_lost > 0),
            npart(d), nfx[["IS"]], nfx[["RT"]], nfx[["PE"]]))
write.csv(manifest, file.path("output", "numbers", "numbers_manifest.csv"),
          row.names = FALSE)
cat("tables 1-3 + manifest written\n")

for (f in list.files(TDIR, pattern = "\\.tex$", full.names = TRUE)) {
  txt <- readLines(f, warn = FALSE)
  txt <- gsub("\\bIS\\b", "ACC", txt, perl = TRUE)
  writeLines(txt, f)
}
