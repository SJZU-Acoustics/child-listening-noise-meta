
source(file.path("R", "helpers.R"))  

p <- read.csv(file.path(".cache", "analysis", "analysis_15_pairs_deepened",
                        "outcome", "pair_level_data.csv"))
stopifnot(nrow(p) == nrow(p26_data(check_freeze = FALSE)$pairs))
NP <- nrow(p); NS <- length(unique(p$study))
BOOT_R <- 1000
OUT <- file.path(".cache", "analysis", "analysis_21_pairs_study_level", "outcome")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUT, "console.txt"), split = TRUE)

boot_ci <- function(stat_fun) {
  studies <- unique(p$study)
  b <- replicate(BOOT_R, {
    ss <- sample(studies, replace = TRUE)
    idx <- unlist(lapply(ss, function(s) which(p$study == s)))
    cl <- rep(seq_along(ss), times = sapply(ss, function(s) sum(p$study == s)))
    stat_fun(p[idx, ], cl)
  })
  b <- b[!is.na(b)]
  c(lo = unname(quantile(b, .025)), hi = unname(quantile(b, .975)),
    draws = length(b))
}

cat("=== A21: unit-of-analysis sensitivity for cor(g_IS, g_RT) ===\n\n")
r_pair <- cor(p$g_IS, p$g_RT)
cat(sprintf("baseline pair-level Pearson r (%d pairs, A15): %+.3f\n\n", NP, r_pair))

sm <- aggregate(cbind(g_IS, g_RT) ~ study, data = p, FUN = mean)
r_sm <- cor(sm$g_IS, sm$g_RT)
set.seed(26)
ci_sm <- boot_ci(function(d, cl) {
  a <- aggregate(cbind(g_IS, g_RT) ~ cl, data = cbind(d, cl = cl), FUN = mean)
  if (nrow(a) < 3 || sd(a$g_IS) == 0 || sd(a$g_RT) == 0) return(NA_real_)
  cor(a$g_IS, a$g_RT)
})
cat(sprintf("1. study-mean r (%d studies): %+.3f, cluster bootstrap 95%% CI [%+.3f, %+.3f] (%d draws)\n",
            NS, r_sm, ci_sm["lo"], ci_sm["hi"], ci_sm["draws"]))

cen <- function(d, g) ave(d[[g]], d$study, FUN = function(x) x - mean(x))
r_wi <- cor(cen(p, "g_IS"), cen(p, "g_RT"))
set.seed(26)
ci_wi <- boot_ci(function(d, cl) {
  ci_ <- ave(d$g_IS, cl, FUN = function(x) x - mean(x))
  cr_ <- ave(d$g_RT, cl, FUN = function(x) x - mean(x))
  if (sd(ci_) == 0 || sd(cr_) == 0) return(NA_real_)
  cor(ci_, cr_)
})
cat(sprintf("2. within-study centred r (%d pairs, %d study means removed): %+.3f, cluster bootstrap 95%% CI [%+.3f, %+.3f]\n",
            NP, NS, r_wi, ci_wi["lo"], ci_wi["hi"]))

loo <- sapply(unique(p$study), function(s) {
  q <- p[p$study != s, ]
  cor(q$g_IS, q$g_RT)
})
cat(sprintf("3. leave-one-study-out pair-level r: range %+.3f to %+.3f (full %+.3f)\n",
            min(loo), max(loo), r_pair))
cat("   by study dropped:\n")
for (s in names(sort(loo))) cat(sprintf("   %-28s %+.3f\n", s, loo[[s]]))

write.csv(data.frame(
  quantity = c("pair-level r", "study-mean r", "within-study centred r",
               "LOSO min", "LOSO max"),
  estimate = c(r_pair, r_sm, r_wi, min(loo), max(loo)),
  boot_lo = c(NA, ci_sm["lo"], ci_wi["lo"], NA, NA),
  boot_hi = c(NA, ci_sm["hi"], ci_wi["hi"], NA, NA)),
  file.path(OUT, "unit_sensitivity.csv"), row.names = FALSE)
write.csv(data.frame(study = names(loo), r_without = unname(loo)),
          file.path(OUT, "leave_one_study_out.csv"), row.names = FALSE)
cat("\nOutputs written to", OUT, "\n")
sink()
