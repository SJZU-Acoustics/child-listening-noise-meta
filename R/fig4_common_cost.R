
source(file.path("R", "style.R"))
suppressMessages(library(patchwork))

OUTC <- function(a, f) file.path(".cache", "analysis", a, "outcome", f)
pairs <- read.csv(OUTC("analysis_15_pairs_deepened", "pair_level_data.csv"))
corr <- read.csv(OUTC("analysis_15_pairs_deepened", "pair_correlation.csv"))
conc <- read.csv(OUTC("analysis_19_reuse_effort", "IS_PE_concordance.csv"))

quad <- c(both = sum(pairs$g_IS > 0 & pairs$g_RT > 0),
          acc = sum(pairs$g_IS > 0 & pairs$g_RT <= 0),
          rt = sum(pairs$g_IS <= 0 & pairs$g_RT > 0),
          neither = sum(pairs$g_IS <= 0 & pairs$g_RT <= 0))
lims_a <- list(x = range(pairs$g_IS) + c(-0.25, 0.35),
               y = range(pairs$g_RT) + c(-0.25, 0.35))

pa <- ggplot(pairs, aes(g_IS, g_RT)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "grey92") +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = 0.3,
              colour = "grey55") +
  geom_point(aes(colour = task_group, shape = task_group), size = 1.7) +
  scale_colour_manual(values = c("same-task" = OKABE[1],
                                 "dual-task" = OKABE[4]),
                      labels = c("same-task" = "same-task", "dual-task" = "cross-task")) +
  scale_shape_manual(values = c("same-task" = 16, "dual-task" = 17),
                     labels = c("same-task" = "same-task", "dual-task" = "cross-task")) +
  annotate("text", x = lims_a$x[2] - 0.05, y = lims_a$y[2] - 0.02,
           hjust = 1, vjust = 1, size = 2.9, family = "Helvetica",
           label = sprintf("both adverse: %d of %d", quad[["both"]], nrow(pairs))) +
  annotate("text", x = lims_a$x[2] - 0.05, y = lims_a$y[1] + 0.12,
           hjust = 1, size = 2.9, family = "Helvetica",
           label = sprintf("ACC only: %d", quad[["acc"]])) +
  annotate("label", x = lims_a$x[1] + 0.05, y = lims_a$y[2] - 0.02,
           hjust = 0, vjust = 1, size = 2.9, family = "Helvetica",
           fill = "white", border.colour = NA, label.padding = unit(0.08, "lines"),
           label = sprintf("RT only: %d", quad[["rt"]])) +
  annotate("label", x = lims_a$x[1] + 0.05, y = lims_a$y[1] + 0.12,
           hjust = 0, size = 2.9, family = "Helvetica",
           fill = "white", border.colour = NA, label.padding = unit(0.08, "lines"),
           label = sprintf("neither: %d", quad[["neither"]])) +
  annotate("text", x = lims_a$x[2] - 0.05, y = 0.85, hjust = 1, size = 2.9,
           family = "Helvetica", parse = TRUE,
           label = sprintf("italic(r)*' = +%.2f [%.2f, %.2f]'",
                           corr$pearson, corr$boot_lo, corr$boot_hi)) +
  coord_cartesian(xlim = lims_a$x, ylim = lims_a$y, expand = FALSE) +
  scale_x_continuous(labels = minus_labels) +
  scale_y_continuous(labels = minus_labels) +
  labs(x = expression("ACC effect (Hedges' "*italic(g)*")"),
       y = expression("RT effect (Hedges' "*italic(g)*")")) +
  theme_pub(base_size = 8, axis_title_size = 9) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.98, 0.36),
        legend.justification = c(1, 1))

lims_b <- list(x = c(-0.15, 1.35), y = c(-0.15, 1.35))
pb <- ggplot(conc, aes(g_IS, g_PE)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "grey92") +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_point(size = 2, colour = FAM_COL[["PE"]]) +
  annotate("text", x = lims_b$x[2] - 0.04, y = lims_b$y[1] + 0.05,
           hjust = 1, size = 2.9, family = "Helvetica",
           label = sprintf("both adverse: %d of %d studies",
                           sum(conc$g_IS > 0 & conc$g_PE > 0), nrow(conc))) +
  coord_cartesian(xlim = lims_b$x, ylim = lims_b$y, expand = FALSE) +
  labs(x = expression("ACC effect (Hedges' "*italic(g)*")"),
       y = expression("PE effect (Hedges' "*italic(g)*")")) +
  theme_pub(base_size = 8, axis_title_size = 9)

p <- pa + pb + plot_layout(widths = c(1.35, 1)) + plot_annotation(tag_levels = "a")

save_fig(p, file.path("output", "figures", "fig4_common_cost.png"),
         "double_column", 3.55)
cat("fig4 written\n")
