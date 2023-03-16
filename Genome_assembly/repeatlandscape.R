library(reshape)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(tidyverse)
library(gridExtra)

# set working directory
setwd("")

# load input file RepeatMakser - createRepeatLandscape.pl + the histogram information at the bottom of the output file.
KimuraDistance <- read.delim("sailfish.repeatslandscape", header = TRUE, sep=" ")

#add here the genome size in bp
genomes_size=612904761

kd_melt <- melt(KimuraDistance,id="Div")
kd_melt$norm <- kd_melt$value/genomes_size * 100

p <- ggplot(kd_melt, aes(fill=variable, y=norm, x=Div)) + 
     geom_bar(position="stack", stat="identity",color="black") +
     scale_fill_viridis(option = "turbo", discrete = TRUE) +
     labs(x = "Kimura substitution level", y = "Percent of the genome", fill = "") + 
     coord_cartesian(xlim = c(0, 45)) +
     theme_bw() +
     theme(axis.text=element_text(size=11),
           axis.title =element_text(size=12))
p

# save plot
ggsave("sailfish.repeatlandscape.png", device = "png", width = 35, height = 30, units = "cm", dpi = 600)
ggsave("sailfish.repeatlandscape.svg", device = "svg", width = 35, height = 30, units = "cm", dpi = 600)

