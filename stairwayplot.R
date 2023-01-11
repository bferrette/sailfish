library(tidyverse)
library(patchwork)

# set working directory
setwd("C:/Users/rapha/Documents/Sync/rcoimbra_phd/project_northern/results/folded_sfs/")

# read input file
tbl <- read_table2("West_African.final.summary")

# generate plot with axis in log scale
p1 <- ggplot(data = tbl, mapping = aes(x = year, y = Ne_median)) +
  geom_ribbon(mapping = aes(ymin = `Ne_2.5%`, ymax = `Ne_97.5%`), fill = "#F0E442", alpha = 0.3, lty = 2) +
  geom_ribbon(mapping = aes(ymin = `Ne_12.5%`, ymax = `Ne_87.5%`), fill = "#F0E442", alpha = 0.35, lty = 2) +
  geom_line(color = "#F0E442", size = 1) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = c(1, 250000)
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = c(1, 20000)
  ) +
  labs(title = "Last 200,000 years (log-scale)", x = "Years", y = expression(italic("N"[e]))) +
  theme_minimal() +
  theme(axis.title = element_text(size = 11),
        axis.text = element_text(size = 9),
        panel.grid.minor = element_blank())

# add log scale ticks to plot
p1.log <- p1 + annotation_logticks(size = 0.2)

# generate zoomed in plot
p2 <- ggplot(data = tbl, mapping = aes(x = year, y = Ne_median)) +
  geom_ribbon(mapping = aes(ymin = `Ne_2.5%`, ymax = `Ne_97.5%`), fill = "#F0E442", alpha = 0.3, lty = 2) +
  geom_ribbon(mapping = aes(ymin = `Ne_12.5%`, ymax = `Ne_87.5%`), fill = "#F0E442", alpha = 0.35, lty = 2) +
  geom_line(color = "#F0E442", size = 1) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(1, 1000)) +
  labs(title = "Last 100 years", x = "Years", y = expression(italic("N"[e]))) +
  theme_minimal() +
  theme(axis.title = element_text(size = 11),
        axis.text = element_text(size = 9))

# generate composed plot
p1.log / p2

# save plot in '.pdf' format
ggsave("stairwayplot.pdf", width = 174, height = 232, units = "mm", dpi = 300)

