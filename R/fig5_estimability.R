
source(file.path("R", "helpers.R"))
source(file.path("R", "style.R"))

D <- p26_data()
d <- D$effects
NS <- length(unique(d$study))

n_vary <- function(v) sum(sapply(split(d[[v]], d$study),
                                 function(x) length(unique(x)) > 1))
tab <- data.frame(
  factor = c("SNR", "Outcome family", "Task category", "Group age",
             "Masker type", "Masking class"),
  varying = c(n_vary("SNR"), n_vary("outcome"), n_vary("task"), n_vary("age"),
              n_vary("noise_type"), n_vary("masking")),
  status = c("within-study estimable", "identified, low power",
             "between-literature only", "weakly identified",
             "between-literature only", "between-literature only")
)
write.csv(tab, file.path("output", "numbers", "fig5_estimability.csv"),
          row.names = FALSE)

tab <- tab[order(tab$varying, tab$factor == "SNR"), ]
tab$factor <- factor(tab$factor, levels = tab$factor)
tab$col <- ifelse(tab$status == "within-study estimable", OKABE[1],
                  ifelse(tab$status == "between-literature only", OKABE[4],
                         OKABE[2]))

p <- ggplot(tab, aes(y = factor)) +
  geom_col(aes(x = NS), fill = "grey90", width = 0.62) +
  geom_col(aes(x = varying), fill = tab$col, width = 0.62) +
  geom_text(aes(x = varying - 0.35, label = varying), hjust = 1, size = 2.9,
            family = "Helvetica", colour = "white", fontface = "bold") +
  geom_text(aes(x = NS - 0.4, label = c(
              "within-study estimable" = "within-study estimate",
              "identified, low power" = "limited power",
              "between-literature only" = "mainly between studies",
              "weakly identified" = "limited variation")[status]), hjust = 1, size = 2.6,
            family = "Helvetica", colour = "grey25") +
  scale_y_discrete(labels = function(x) replace(x, x == "Outcome family", "Outcome type")) +
  scale_x_continuous(limits = c(0, NS), breaks = c(0, NS / 2, NS),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = sprintf("Studies varying the factor within study (of %d)", NS), y = NULL) +
  theme_pub(base_size = 8, axis_title_size = 9) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

save_fig(p, file.path("output", "figures", "fig5_estimability.png"),
         "single_column", 2.5)
cat("fig5 written\n")
