
source(file.path("R", "helpers.R"))
OUT <- file.path(p26_outdir("analysis_12_design_descriptors"), "outcome")

D <- p26_data()
d <- D$effects
d <- d[order(d$study, d$effect_ID), ]

corpus <- read.csv(file.path(P26_ROOT, ".cache", "data", "corpus_25_reports.csv"),
                   stringsAsFactors = FALSE)
d$report_id <- paste0("R", sprintf("%02d", as.integer(sub("\\.pdf$", "", d$report_source))))
d$year <- setNames(corpus$year, corpus$report_id)[d$report_id]
stopifnot(!any(is.na(d$year)))
d$n_part <- pmin(d$n1, d$n2)

mkV <- function(sub) vcalc(vi = vi, cluster = study, grp1 = experiment_group_id,
                           grp2 = control_group_id, w1 = n1, w2 = n2,
                           data = sub, sparse = FALSE)
report <- function(sub, mods, idx, label) {
  V <- mkV(sub)
  m <- try(rma.mv(yi_al, V, mods = mods, data = sub,
                  random = list(~ 1 | study, ~ 1 | effect_ID),
                  method = "REML", test = "t", dfs = "contain"), silent = TRUE)
  if (inherits(m, "try-error")) { cat(sprintf("  %-34s model failed\n", label)); return(NULL) }
  ct <- clubSandwich::coef_test(m, vcov = "CR2", cluster = sub$study)
  out <- data.frame(term = rownames(m$b)[idx], estimate = as.numeric(m$b)[idx],
                    p_model = m$pval[idx], p_CR2 = ct$p_Satt[idx], df_CR2 = ct$df_Satt[idx])
  for (i in seq_along(idx))
    cat(sprintf("  %-34s %-22s %+10.4f  p_model %-9.3g CR2 p %-8.4f df %.1f\n",
                if (i == 1) label else "", out$term[i], out$estimate[i],
                out$p_model[i], out$p_CR2[i], out$df_CR2[i]))
  out
}

sink(file.path(OUT, "console.txt"), split = TRUE)
cat("A12 - design descriptors\n\n")

cat("=== 1. task (the one A10 left as a gap) ===\n")
print(table(d$task, d$outcome))
nt <- tapply(d$task, d$study, function(x) length(unique(x)))
cat(sprintf("\ntask varies within %d of %d studies (%d effects)\n",
            sum(nt > 1), length(nt), sum(d$study %in% names(nt)[nt > 1])))

dis <- d[d$outcome == "IS" & d$task != "", ]
keep <- names(which(table(dis$task) >= 10))
dis <- dis[dis$task %in% keep, ]
dis$task <- factor(dis$task, levels = c("speech perception", setdiff(keep, "speech perception")))
cat(sprintf("IS analysis on levels with k >= 10: %s (%d effects, %d studies)\n",
            paste(levels(dis$task), collapse = ", "), nrow(dis), length(unique(dis$study))))
cat("\n  random effects (between + within information):\n")
t_re <- report(dis, ~ task, 2:nlevels(dis$task), "task, IS")

nt2 <- tapply(dis$task, dis$study, function(x) length(unique(x)))
dwt <- dis[dis$study %in% names(nt2)[nt2 > 1], ]
cat(sprintf("\n  study fixed effects (purely within-study; %d studies, %d effects):\n",
            length(unique(dwt$study)), nrow(dwt)))
Vw <- mkV(dwt)
mfe <- try(rma.mv(yi_al, Vw, mods = ~ task + factor(study), data = dwt,
                  random = ~ 1 | effect_ID, method = "REML", test = "t", dfs = "contain"),
           silent = TRUE)
if (!inherits(mfe, "try-error")) {
  ctfe <- clubSandwich::coef_test(mfe, vcov = "CR2", cluster = dwt$study)
  k <- grep("^task", rownames(mfe$b))
  t_fe <- data.frame(term = rownames(mfe$b)[k], estimate = as.numeric(mfe$b)[k],
                     p_model = mfe$pval[k], p_CR2 = ctfe$p_Satt[k], df_CR2 = ctfe$df_Satt[k])
  print(t_fe, row.names = FALSE, digits = 4)
} else { t_fe <- NULL; cat("  fixed-effects model failed\n") }

cat("\n=== 2. sample size ===\n")
d$n_between <- ave(d$n_part, d$study, FUN = mean)
d$n_within  <- d$n_part - d$n_between
d$n_between_c <- d$n_between - mean(d$n_between)
cat(sprintf("participants per effect: median %.0f, range %d-%d\n",
            median(d$n_part), min(d$n_part), max(d$n_part)))
g <- ave(d$n_part, d$study, FUN = mean)
cat(sprintf("variation: within-study %.1f%% / between-study %.1f%%\n",
            100 * sum((d$n_part - g)^2) / (sum((d$n_part - g)^2) + sum((g - mean(d$n_part))^2)),
            100 * sum((g - mean(d$n_part))^2) / (sum((d$n_part - g)^2) + sum((g - mean(d$n_part))^2))))
n_res <- report(d, ~ n_within + n_between_c, 2:3, "sample size, all outcomes")

cat("\n=== 3. publication year (between-study by construction) ===\n")
cat(sprintf("range %d-%d; %d distinct years across %d studies\n",
            min(d$year), max(d$year), length(unique(d$year)), length(unique(d$study))))
ny <- tapply(d$year, d$study, function(x) length(unique(x)))
cat(sprintf("studies with >1 publication year: %d (only the two-report study can have one)\n",
            sum(ny > 1)))
d$year_c <- d$year - mean(d$year)
y_res <- report(d, ~ year_c, 2, "year, all outcomes")
y_is  <- report(d[d$outcome == "IS", ], ~ year_c, 2, "year, IS")
cat("\n  Whatever the coefficient, this is a comparison between literatures\n")
cat("  published at different times, confounded with every other way the field\n")
cat("  changed. A04 separates a within-study SNR component from the between-study\n")
cat("  one; year has no within component at all, so nothing can separate it.\n")

res <- rbind(
  if (!is.null(t_re)) cbind(descriptor = "task (random effects)", t_re),
  if (!is.null(t_fe)) cbind(descriptor = "task (study fixed effects)", t_fe),
  if (!is.null(n_res)) cbind(descriptor = "sample size", n_res),
  if (!is.null(y_res)) cbind(descriptor = "year, all", y_res),
  if (!is.null(y_is)) cbind(descriptor = "year, IS", y_is))
write.csv(res, file.path(OUT, "design_descriptors.csv"), row.names = FALSE)
cat("\nwritten to", OUT, "\n")
sink()
