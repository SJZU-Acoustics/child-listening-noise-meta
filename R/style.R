
suppressMessages(library(ggplot2))

WIDTH_MM <- c(single_column = 85, double_column = 178)
mm2in <- function(mm) mm / 25.4

OKABE <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#CC79A7")
FAM_COL <- c(IS = "#0072B2", RT = "#E69F00", PE = "#009E73", ALL = "#000000")

theme_pub <- function(base_size = 9, axis_title_size = 10) {
  theme_classic(base_size = base_size, base_family = "Helvetica") %+replace%
    theme(
      axis.line         = element_line(colour = "black", linewidth = 0.4),
      axis.ticks        = element_line(colour = "black", linewidth = 0.35),
      axis.ticks.length = unit(0.10, "cm"),
      axis.title        = element_text(size = axis_title_size),
      axis.text         = element_text(size = base_size, colour = "black"),
      legend.text       = element_text(size = base_size),
      legend.title      = element_blank(),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.key.size   = unit(0.35, "cm"),
      panel.grid        = element_blank(),
      plot.title        = element_blank(),
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA),
      plot.tag          = element_text(size = 10, face = "bold", family = "Helvetica"),
      complete = FALSE
    )
}

save_fig <- function(p, file, width_mode, height_in) {
  ggsave(file, p, width = mm2in(WIDTH_MM[[width_mode]]), height = height_in,
         dpi = 600, bg = "white")
}

pub_label <- function(x) gsub("\\bIS\\b", "ACC", x, perl = TRUE)

# Negative tick labels take the minus sign (U+2212), never a hyphen.
minus_labels <- function(x) {
  out <- format(x, trim = TRUE, justify = "left")
  out[is.na(x)] <- NA
  gsub("-", "−", out, fixed = TRUE)
}
