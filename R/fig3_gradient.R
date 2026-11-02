
source(file.path("R", "helpers.R"))
source(file.path("R", "style.R"))
suppressMessages(library(patchwork))

D <- p26_data()
d <- D$effects

vary <- do.call(rbind, lapply(split(d, list(d$study, d$outcome), drop = TRUE),
  function(s) if (length(unique(s$SNR)) > 1) s else NULL))

study_key <- data.frame(study = sort(unique(vary$study)))
study_key$id <- seq_len(nrow(study_key))
vary$series_id <- study_key$id[match(vary$study, study_key$study)]

series_cols <- setNames(grDevices::hcl.colors(nrow(study_key), "Dark 3"), study_key$study)
traj_panel <- function(fam) {
  s <- vary[vary$outcome == fam, ]
  s <- aggregate(yi_al ~ SNR + study + series_id, data = s, FUN = mean)
  ns <- length(unique(s$study))
  ends <- do.call(rbind, lapply(split(s, s$study), function(z) z[which.max(z$SNR), ]))
  ggplot(s, aes(SNR, yi_al, group = study, colour = study)) +
    geom_hline(yintercept = 0, linetype = "22", linewidth = 0.25, colour = "grey70") +
    geom_line(linewidth = 0.35, alpha = 0.5) +
    geom_point(size = 1.1) +
    ggrepel::geom_text_repel(data = ends, aes(label = series_id),
      size = 2.8, family = "Helvetica", seed = 26, max.overlaps = Inf,
      box.padding = 0.45, point.padding = 0.2, min.segment.length = 0,
      segment.colour = "grey60", segment.size = 0.2, nudge_x = 1.3) +
    scale_colour_manual(values = series_cols) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.4,
             label = sprintf("%s (%d %s)", pub_label(fam), ns,
                             ifelse(ns == 1, "study", "studies")),
             size = 2.9, family = "Helvetica") +
    scale_x_continuous(limits = c(min(vary$SNR), max(vary$SNR) + 3)) +
    scale_y_continuous(limits = range(vary$yi_al) + c(-0.15, 0.15)) +
    labs(x = "SNR (dB)", y = expression("Effect size (Hedges' "*italic(g)*")")) +
    theme_pub(base_size = 8, axis_title_size = 9) + theme(legend.position = "none")
}
pa <- traj_panel("IS")
pb <- traj_panel("RT")
pc <- traj_panel("PE")

dec <- read.csv(file.path(".cache", "analysis", "analysis_04_snr_information",
                          "outcome", "within_between_decomposition.csv"))
dec$component <- ifelse(dec$term == "SNR_within", "Within-study slope",
                        "Between-study slope")
dec$family <- factor(dec$family, levels = c("ALL", "IS", "RT", "PE"))
dec$component <- factor(dec$component,
                        levels = c("Within-study slope", "Between-study slope"))

pd <- ggplot(dec, aes(x = estimate, y = family, colour = component,
                      shape = component)) +
  geom_vline(xintercept = 0, linetype = "22", linewidth = 0.3,
             colour = "grey55") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0,
                 linewidth = 0.4, position = position_dodge(width = 0.55)) +
  geom_point(size = 2.1, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = c("Within-study slope" = "black",
                                 "Between-study slope" = "grey55")) +
  scale_shape_manual(values = c("Within-study slope" = 16,
                                "Between-study slope" = 1)) +
  scale_y_discrete(limits = rev(levels(dec$family)), labels = pub_label) +
  labs(x = expression("SNR slope ("*italic(g)*" per dB)"), y = NULL) +
  theme_pub(base_size = 8, axis_title_size = 9) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.03, 0.97),
        legend.justification = c(0, 1))

key_rows <- ceiling(nrow(study_key) / 2)
study_key$x <- rep(c(0, 0.52), each = key_rows)[seq_len(nrow(study_key))]
study_key$y <- rep(seq(key_rows, 1), 2)[seq_len(nrow(study_key))]
study_key$label <- paste0(study_key$id, "  ", sub(",", ", ", study_key$study, fixed = TRUE))
pkey <- ggplot(study_key, aes(x, y, label = label, colour = study)) +
  geom_text(hjust = 0, size = 2.8, family = "Helvetica") +
  scale_colour_manual(values = series_cols) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.5, key_rows + 0.5), clip = "off") +
  theme_void() + theme(legend.position = "none", plot.margin = margin(0, 8, 0, 12))
plots <- (pa | pb) / (pc | pd) + plot_annotation(tag_levels = "a")
p <- wrap_elements(full = plots) / pkey + plot_layout(heights = c(5, 1.1))
save_fig(p, file.path("output", "figures", "fig3_gradient.png"),
         "double_column", 6.4)
cat("fig3 written\n")
