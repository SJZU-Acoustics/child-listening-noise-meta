
source(file.path("R", "style.R"))

box <- function(id, x, y, w, h, label) {
  data.frame(id = id, x = x, y = y, w = w, h = h, label = label)
}

B <- rbind(
  box("db", 24, 95, 38, 14, "Database records\n(n = 3,091)"),
  box("dup", 62, 95, 32, 10, "Duplicates removed\n(n = 798)"),
  box("oth", 100, 95, 32, 14, "Supplementary\nreports (n = 21)"),
  box("scr", 24, 77, 38, 9, "Records screened\n(n = 2,293)"),
  box("scrx", 62, 77, 32, 10, "Records excluded\n(n = 2,221)"),
  box("ft", 24, 59, 38, 10, "Full texts assessed\n(n = 72)"),
  box("ftx", 62, 59, 32, 10, "Reports excluded\n(n = 52)"),
  box("oft", 100, 59, 32, 10, "Full texts assessed\n(n = 21)"),
  box("oftx", 139, 59, 32, 10, "Reports excluded\n(n = 16)"),
  box("inc", 24, 40, 38, 10, "Reports included\n(n = 20)"),
  box("oinc", 100, 40, 32, 10, "Reports included\n(n = 5)"),
  box("rep", 62, 23, 44, 9, "Included reports\n(n = 25)"),
  box("st", 62, 10, 44, 10, "24 studies\n269 effect sizes")
)
flow <- read.csv(".cache/data/search_flow.csv", stringsAsFactors = FALSE)
counts <- setNames(flow$count, flow$quantity)
quantities <- c("database_records", "duplicates", "supplementary_fulltexts", "screened", "screening_exclusions", "database_fulltexts", "database_fulltext_exclusions", "supplementary_fulltexts", "supplementary_exclusions", "database_inclusions", "supplementary_inclusions", "included_reports")
for (i in seq_along(quantities)) B$label[i] <- sub("n = [0-9,]+", paste0("n = ", format(counts[quantities[i]], big.mark = ",", trim = TRUE)), B$label[i])
B$label[B$id == "st"] <- sprintf("%d studies\n%d effect sizes", counts["included_studies"], counts["effect_sizes"])
edge <- function(id, side) {
  b <- B[B$id == id, ]
  switch(side, top = c(b$x, b$y + b$h / 2), bottom = c(b$x, b$y - b$h / 2),
         left = c(b$x - b$w / 2, b$y), right = c(b$x + b$w / 2, b$y))
}
arrow_df <- function(from, fs, to, ts) {
  a <- edge(from, fs); b <- edge(to, ts)
  data.frame(x = a[1], y = a[2], xend = b[1], yend = b[2])
}
arrows <- rbind(
  arrow_df("db", "bottom", "scr", "top"),
  arrow_df("db", "right", "dup", "left"),
  arrow_df("scr", "bottom", "ft", "top"),
  arrow_df("scr", "right", "scrx", "left"),
  arrow_df("ft", "bottom", "inc", "top"),
  arrow_df("ft", "right", "ftx", "left"),
  arrow_df("oth", "bottom", "oft", "top"),
  arrow_df("oft", "right", "oftx", "left"),
  arrow_df("oft", "bottom", "oinc", "top")
)
elbow <- rbind(
  data.frame(x = 24, y = edge("inc", "bottom")[2], xend = 24, yend = 23),
  data.frame(x = 100, y = edge("oinc", "bottom")[2], xend = 100, yend = 23)
)
arrows <- rbind(arrows,
  data.frame(x = 24, y = 23, xend = edge("rep", "left")[1], yend = 23),
  data.frame(x = 100, y = 23, xend = edge("rep", "right")[1], yend = 23),
  arrow_df("rep", "bottom", "st", "top"))

p <- ggplot() +
  geom_rect(data = B, aes(xmin = x - w / 2, xmax = x + w / 2,
                          ymin = y - h / 2, ymax = y + h / 2),
            fill = "white", colour = "black", linewidth = 0.35) +
  geom_text(data = B, aes(x = x, y = y, label = label), size = 2.85,
            family = "Helvetica", lineheight = 0.95) +
  geom_segment(data = elbow, aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.35) +
  geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.35,
               arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  annotate("text", x = 1.8, y = 95, label = "Identification", angle = 90,
           size = 3, family = "Helvetica", fontface = "bold") +
  annotate("text", x = 1.8, y = 68, label = "Screening", angle = 90,
           size = 3, family = "Helvetica", fontface = "bold") +
  annotate("text", x = 1.8, y = 23, label = "Included", angle = 90,
           size = 3, family = "Helvetica", fontface = "bold") +
  coord_cartesian(xlim = c(0, 157), ylim = c(4, 104), expand = FALSE) +
  theme_void(base_family = "Helvetica")

save_fig(p, file.path("output", "figures", "fig1_prisma.png"),
         "double_column", 4.8)
cat("fig1 written\n")
