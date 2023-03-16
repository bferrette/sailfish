# Load libraries
library(ggplot2)
library(patchwork)

# set working directory
setwd("~/Documents/billfishes/mitogenomes/sailfish/analysis/ebsp")
# Load EBSP tab-delimited data
ATL <- read.table('sailfish_atlantic_lineage_ebsp_combined.tsv', header=TRUE, sep = "\t")
IDWP <- read.table('sailfish_global_lineage_ebsp_combined.tsv', header=TRUE, sep = "\t")

# Make plots
p01 <- ggplot(ATL, aes(x=time)) +
       geom_line(aes(y=median), linetype='solid', size = 1) +
       geom_ribbon(aes(ymin=lower, ymax=upper), fill="#2171B5",
                       size=.2, alpha=.25)+
#       scale_x_continuous(labels = scales::comma) +
       scale_y_continuous(labels = scales::comma) +
       scale_x_continuous(trans='log10', labels = scales::comma) +
#       scale_y_continuous(limits = c(0.01, 10), trans='log10', labels = scales::comma) +
       labs(title = "EBSP Atlantic lineage", x = "Million years ago", y = expression(italic("N"[e]))) +
       theme_bw()
#       theme(axis.title.x = element_blank())
p01

p02 <- ggplot(IDWP, aes(x=time)) +
       geom_line(aes(y=median), linetype='solid', size = 1) +
       geom_ribbon(aes(ymin=lower, ymax=upper), fill="#6A51A3",
                       size=.2, alpha=.25)+
#       scale_x_continuous(labels = scales::comma) +
       scale_y_continuous(labels = scales::comma) +
       scale_x_continuous(trans='log10',labels = scales::comma) +
#       scale_y_continuous(trans='log10', labels = scales::comma) +
       labs(title = "EBSP Global lineage", x = "Million years ago", y = expression(italic("N"[e]))) +
       theme_bw()
p02

# Combine the plots with patchwork
(p01 / p02) + plot_annotation(tag_levels = 'a')

# save plot in '.svg' format
ggsave("sailfish_ebsp.svg", device = "svg", width = 21, height = 15, units = "cm", dpi = 600)
