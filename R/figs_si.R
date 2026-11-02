
source(file.path("R", "helpers.R"))
source(file.path("R", "style.R"))
suppressMessages(library(patchwork))

OUTC <- function(a, f) file.path(".cache", "analysis", a, "outcome", f)
FDIR <- file.path("output", "figures", "si")
dir.create(FDIR, recursive = TRUE, showWarnings = FALSE)

D <- p26_data()
d <- D$effects

fun <- read.csv(OUTC("analysis_06_publication_bias", "funnel_data.csv"))
fun_panel <- function(fam, col) {
  s <- fun[fun$outcome == fam, ]
  ggplot(s, aes(yi, sqrt(n))) +
    geom_vline(xintercept = 0, linetype = "22", linewidth = 0.3,
               colour = "grey55") +
    geom_point(size = 1.3, colour = col, alpha = 0.75) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.4,
             label = sprintf("%s (%d effects)", pub_label(fam), nrow(s)), size = 2.9,
             family = "Helvetica") +
    labs(x = expression("Effect size "*italic(g)*" (aligned)"),
         y = expression(sqrt(italic(n)))) +
    theme_pub(base_size = 8, axis_title_size = 9)
}
p <- fun_panel("IS", FAM_COL[["IS"]]) |
  (fun_panel("RT", FAM_COL[["RT"]]) + theme(axis.title.y = element_blank())) |
  (fun_panel("PE", FAM_COL[["PE"]]) + theme(axis.title.y = element_blank()))
p <- p + plot_annotation(tag_levels = "a")
save_fig(p, file.path(FDIR, "si_fig_funnels.png"), "double_column", 2.6)

oh <- read.csv(file.path("output", "numbers", "ohat_domain_matrix.csv"),
               check.names = FALSE)
doms <- c("Selection bias", "Confounding bias", "Performance bias",
          "Attrition/exclusion bias", "Detection bias",
          "Selective reporting bias", "Other sources of bias")
long <- do.call(rbind, lapply(doms, function(dm)
  data.frame(study = oh$Study, domain = dm, rating = oh[[dm]])))
long$rating <- factor(long$rating, levels = c("++", "+", "-", "--", "NR", "N/A"))
stopifnot(!anyNA(long$rating))
long$domain <- factor(long$domain, levels = doms,
                      labels = c("Selection", "Confounding", "Performance",
                                 "Attrition", "Detection", "Reporting",
                                 "Other"))
long$study <- sub(",(\\d)", ", \\1", long$study)
p <- ggplot(long, aes(domain, study, fill = rating)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = c("++" = "#009E73", "+" = "#8FD3B6",
                               "-" = "#D55E00", "--" = "#7A2E00",
                               "NR" = "#E69F00", "N/A" = "grey85"),
                    labels = c("++" = "++ definitely low", "+" = "+ probably low",
                               "-" = "− probably high", "--" = "−− definitely high",
                               "NR" = "NR not reported", "N/A" = "N/A not applicable"),
                    drop = FALSE) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_pub(base_size = 8, axis_title_size = 9) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom", legend.location = "plot",
        axis.text.x.top = element_text(angle = 30, hjust = 0),
        plot.margin = margin(5, 24, 5, 5))
save_fig(p, file.path(FDIR, "si_fig_ohat.png"), "single_column", 4.6)

sc <- read.csv(OUTC("analysis_13_specification_curve",
                    "specification_curve.csv"))
sc$working <- factor(sc$working,
                     levels = c("shared-arm", "diagonal", "rho=0.5"))
spec_panel <- function(qty, show_leg = FALSE, panel_lab = qty) {
  s <- sc[sc$quantity == qty, ]
  stopifnot(nrow(s) == 48)
  s <- s[order(s$estimate), ]
  s$rank <- seq_len(nrow(s))
  ggplot(s, aes(rank, estimate, colour = working)) +
    geom_hline(yintercept = 0, linetype = "22", linewidth = 0.3,
               colour = "grey55") +
    geom_point(size = 0.9) +
    scale_colour_manual(values = c("shared-arm" = OKABE[1],
                                   "diagonal" = OKABE[3],
                                   "rho=0.5" = OKABE[4]),
                        labels = c("shared-arm V (adopted)", "diagonal V",
                                   expression("constant "*rho*" = 0.5"))) +
    annotate("text", x = 1, y = Inf, hjust = 0, vjust = 2.4, size = 2.7,
             family = "Helvetica", label = pub_label(panel_lab)) +
    labs(x = "Specification (ranked)",
         y = expression("Estimate ("*italic(g)*")")) +
    theme_pub(base_size = 8, axis_title_size = 9) +
    theme(legend.position = if (show_leg) "bottom" else "none")
}
qs <- c("pooled ALL" = "pooled ALL", "pooled IS" = "pooled IS",
        "pooled RT" = "pooled RT", "pooled PE" = "pooled PE",
        "SNR within ALL" = "SNR slope ALL (within-study)",
        "SNR within IS" = "SNR slope IS (within-study)")
ps <- lapply(seq_along(qs), function(i)
  spec_panel(names(qs)[i], i == 5, panel_lab = qs[[i]]))
p <- (ps[[1]] | ps[[2]] | ps[[3]]) / (ps[[4]] | ps[[5]] | ps[[6]]) +
  plot_annotation(tag_levels = "a")
save_fig(p, file.path(FDIR, "si_fig_spec.png"), "double_column", 4.6)

set.seed(26)
pa <- ggplot(d, aes(SNR, outcome, colour = outcome)) +
  geom_jitter(height = 0.18, size = 1.1, alpha = 0.65) +
  scale_colour_manual(values = FAM_COL[c("IS", "RT", "PE")]) +
  scale_y_discrete(limits = c("PE", "RT", "IS"), labels = pub_label) +
  labs(x = "SNR (dB)", y = NULL) +
  theme_pub(base_size = 8, axis_title_size = 9) +
  theme(legend.position = "none")
pb <- ggplot(d, aes(age, outcome, colour = outcome)) +
  geom_jitter(height = 0.18, size = 1.1, alpha = 0.65) +
  scale_colour_manual(values = FAM_COL[c("IS", "RT", "PE")]) +
  scale_y_discrete(limits = c("PE", "RT", "IS"), labels = pub_label) +
  labs(x = "Mean age of group (y)", y = NULL) +
  theme_pub(base_size = 8, axis_title_size = 9) +
  theme(legend.position = "none")
p <- pa | pb
p <- p + plot_annotation(tag_levels = "a")
save_fig(p, file.path(FDIR, "si_fig_distribution.png"), "double_column", 2.5)

cat("SI figures written\n")
