
source(file.path("R", "helpers.R"))
source(file.path("R", "style.R"))
suppressMessages(library(patchwork))

OUTC <- function(a, f) file.path(".cache", "analysis", a, "outcome", f)

per_is <- read.csv(OUTC("analysis_11_raw_scale", "IS_per_study_points_lost.csv"))
head_raw <- read.csv(OUTC("analysis_11_raw_scale", "raw_scale_headlines.csv"))
pe_prof <- read.csv(OUTC("analysis_19_reuse_effort", "PE_study_profile.csv"))
pis <- read.csv(OUTC("analysis_05_heterogeneity", "prediction_intervals.csv"))

D <- p26_data()
d <- D$effects
drt <- d[d$outcome == "RT", ]
e <- escalc(measure = "MD", m1i = m1, sd1i = sd1, n1i = n1,
            m2i = m2, sd2i = sd2, n2i = n2, data = drt)   # noise - control, ms
mkV <- function(sub) vcalc(vi = vi, cluster = study, grp1 = experiment_group_id,
                           grp2 = control_group_id, w1 = n1, w2 = n2,
                           data = sub, sparse = FALSE)
per_rt <- do.call(rbind, lapply(split(e, e$study), function(s) {
  if (nrow(s) == 1) {
    return(data.frame(study = s$study[1], k = 1L, est = s$yi[1],
                      lo = s$yi[1] - 1.96 * sqrt(s$vi[1]),
                      hi = s$yi[1] + 1.96 * sqrt(s$vi[1])))
  }
  Vs <- mkV(s)
  m <- rma.mv(yi, Vs, data = s, random = ~ 1 | effect_ID, method = "REML")
  data.frame(study = s$study[1], k = nrow(s), est = m$b[1],
             lo = m$ci.lb[1], hi = m$ci.ub[1])
}))
write.csv(per_rt, file.path("output", "numbers", "fig2_RT_per_study_ms.csv"),
          row.names = FALSE)

forest_panel <- function(df, pooled, pi, xlab, col, digits = 1) {
  df <- df[order(df$est), ]
  rows <- rbind(
    data.frame(study = "95% prediction interval", est = NA,
               lo = pi[1], hi = pi[2], kind = "pi"),
    data.frame(study = "Pooled", est = pooled[1], lo = pooled[2],
               hi = pooled[3], kind = "pool"),
    data.frame(study = df$study, est = df$est, lo = df$lo, hi = df$hi,
               kind = "study")
  )
  rows$study <- sub(",(\\d)", ", \\1", rows$study)
  rows$study <- factor(rows$study, levels = rows$study)
  ggplot(rows, aes(y = study)) +
    geom_vline(xintercept = 0, linetype = "22", linewidth = 0.3,
               colour = "grey55") +
    geom_errorbarh(aes(xmin = lo, xmax = hi),
                   height = 0, linewidth = 0.35,
                   colour = ifelse(rows$kind == "pi", "grey55", "black")) +
    geom_point(data = rows[rows$kind == "study", ],
               aes(x = est), shape = 15, size = 1.5, colour = col) +
    geom_point(data = rows[rows$kind == "pool", ],
               aes(x = est), shape = 18, size = 3.1, colour = col) +
    labs(x = xlab, y = NULL) +
    theme_pub(base_size = 8, axis_title_size = 9) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
}

pa <- forest_panel(
  transform(per_is, est = points_lost),
  c(head_raw$pooled[1], head_raw$ci_lo[1], head_raw$ci_hi[1]),
  c(head_raw$pi_lo[1], head_raw$pi_hi[1]),
  "Accuracy loss (percentage points)", FAM_COL[["IS"]])

pb <- forest_panel(
  per_rt,
  c(head_raw$pooled[2], head_raw$ci_lo[2], head_raw$ci_hi[2]),
  c(head_raw$pi_lo[2], head_raw$pi_hi[2]),
  "Response-time increase (ms)", FAM_COL[["RT"]])

pe_df <- data.frame(study = pe_prof$study, est = pe_prof$g_PE,
                    lo = pe_prof$g_PE - 1.96 * pe_prof$se,
                    hi = pe_prof$g_PE + 1.96 * pe_prof$se)
pe_row <- pis[pis$family == "PE", ]
pc <- forest_panel(
  pe_df,
  c(pe_row$estimate, pe_row$ci_lo, pe_row$ci_hi),
  c(pe_row$pi_lo, pe_row$pi_hi),
  "Perceived-effort effect (Hedges' g)", FAM_COL[["PE"]])

p <- pa | (pb / pc + plot_layout(heights = c(14, 9)))
p <- p + plot_layout(widths = c(1.08, 1)) + plot_annotation(tag_levels = "a")

save_fig(p, file.path("output", "figures", "fig2_cost.png"),
         "double_column", 6.6)
cat("fig2 written\n")
