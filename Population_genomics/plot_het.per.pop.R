library(tidyverse)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(patchwork)

# set working directory
setwd("~/sailfish/heterozygosity")

# read list of SFS files
het <- read.delim("sailfish.het.stats.tsv", sep = "\t", header = TRUE)

# create bar plot plot
p1 <- het %>%
      arrange(meanHet) %>%
      mutate(Pop = factor(pop, levels=c("WCA","SWA","ECA","SWI","SEI","WCP"))) %>%
      ggplot(mapping = aes(x = sample, y = meanHet)) +
      geom_col(mapping = aes(fill = Pop),width = 0.8) +
#      geom_errorbar(mapping = aes(ymin = mean-sd, ymax = mean+sd),
#                    position = "dodge",
#                    width = 0.4,
#                    size = 0.5) +
# add matrix of panels defined by two column faceting variables
      facet_grid(cols = vars(Pop), scales = "free_x", space = "free_x") +
#      facet_grid(~fct_inorder(as.character(Pop)),
#             #           switch = "x",
#             scales = 'free',
#             space = 'free') +
#      scale_fill_viridis(option="", discrete=TRUE) +
#      scale_fill_brewer(palette = "Set1", direction=1) + 
      scale_fill_manual(values = c("#2171B5","#6BAED6","#C6DBEF","#6A51A3","#9E9AC8","#DADAEB")) + 
#      scale_y_continuous(labels = scales::number_format(accuracy = 0.0001)) +
#      xlab("Samples") +
#      ylim(0,0.003) +
      ylab("Mean Heterozygosity (%)") +
#      scale_y_continuous(limits = c(0,0.27), labels = scales::percent_format(accuracy = 0.01)) +
      theme_bw() +
      theme(axis.title = element_text(size = 10),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_text(size = 6, angle = 0, hjust = 0.5),
            axis.text.x = element_blank(),
#            axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
#            strip.text.x = element_text(size = 11),
#            strip.text.x = element_blank(),
            legend.position = "none")
p1
# save plot in '.svg' format
ggsave("sailfish.het.ind.svg", device = "svg", width = 21, height = 15, units = "cm", dpi = 600)
ggsave("sailfish.het.ind.png", device = "png", width = 21, height = 15, units = "cm", dpi = 600)

# box plot per population

p2 <- het %>%
      arrange(meanHet) %>%
      mutate(Pop = factor(pop, levels=c("WCA","SWA","ECA","SWI","SEI","WCP"))) %>%
      ggplot( mapping = aes(x = Pop, y = meanHet, fill = Pop)) +
      geom_boxplot() +
      geom_jitter(size = 2, shape=16, position=position_jitter(0.2)) +
      ylab("Mean Heterozygosity (%)") +
      scale_fill_manual(values = c("#2171B5","#6BAED6","#C6DBEF","#6A51A3","#9E9AC8","#DADAEB")) + 
      theme_bw() +
      theme(axis.title = element_text(size = 11),
#            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_text(size = 8, angle = 0, hjust = 0.5),
#            axis.text.x = element_blank(),
            axis.text.x = element_text(size = 8, angle = 0, vjust = 0.5, hjust = 0.5),
#            strip.text.x = element_text(size = 11),
            legend.title = element_blank(),
            legend.position = "none")
p2

# save plot in '.svg' format
ggsave("sailfish.het.pop.svg", device = "svg", width = 21, height = 15, units = "cm", dpi = 600)
ggsave("sailfish.het.pop.png", device = "png", width = 21, height = 15, units = "cm", dpi = 600)

# Missing data

# create bar plot plot
p3 <- het %>%
      arrange(pmiss) %>%
      mutate(Pop = factor(pop, levels=c("WCA","SWA","ECA","SWI","SEI","WCP"))) %>%
      ggplot(mapping = aes(x = sample, y = pmiss)) +
      geom_col(mapping = aes(fill = Pop),width = 0.8, alpha = 0.5) +
# add error bars
#      geom_errorbar(mapping = aes(ymin = mean-sd, ymax = mean+sd),
#                 position = "dodge",
#                 width = 0.4,
#                 size = 0.5) +
# add matrix of panels defined by two column faceting variables
      facet_grid(cols = vars(Pop), scales = "free_x", space = "free_x") +
#      facet_grid(~fct_inorder(as.character(pop)),
#                 switch = "x",
#                 scales = 'free',
#                 space = 'free') +
#      scale_fill_viridis(option="", discrete=TRUE) +
#      scale_fill_brewer(palette = "Set1", direction=1) + 
      scale_fill_manual(values = c("#2171B5","#6BAED6","#C6DBEF","#6A51A3","#9E9AC8","#DADAEB")) + 
#      scale_y_continuous(labels = scales::number_format(accuracy = 0.0001)) +
#      xlab("Samples") +
#      ylim(0,0.003) +
      ylab("Missing data (%)") +
#      scale_y_continuous(limits= c(0,45), labels = scales::percent_format(accuracy = 0.01)) +
      theme_bw() +
      theme(axis.title = element_text(size = 10),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_text(size = 6, angle = 0, hjust = 0.5),
#            axis.text.x = element_blank(),
            axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1),
#            strip.text.x = element_text(size = 11),
            strip.text.x = element_blank(),
            legend.position = "none")
p3

# save plot in '.svg' format
ggsave("sailfish.mis.ind.svg", device = "svg", width = 21, height = 15, units = "cm", dpi = 600)
ggsave("sailfish.mis.ind.png", device = "png", width = 21, height = 15, units = "cm", dpi = 600)

# box plot per population

p4 <- het %>%
      arrange(pmiss) %>%
      mutate(Pop = factor(pop, levels=c("WCA","SWA","ECA","SWI","SEI","WCP"))) %>%
      ggplot(mapping = aes(x = Pop, y = pmiss, fill = Pop)) +
      geom_boxplot(alpha = 0.5) +
      geom_jitter(size = 2, shape=16, position=position_jitter(0.2)) +
      ylab("Missing data (%)") +
      scale_fill_manual(values = c("#2171B5","#6BAED6","#C6DBEF","#6A51A3","#9E9AC8","#DADAEB")) + 
      theme_bw() +
      theme(axis.title = element_text(size = 11),
#            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_text(size = 8, angle = 0, hjust = 0.5),
#            axis.text.x = element_blank(),
            axis.text.x = element_text(size = 8, angle = 0, vjust = 0.5, hjust = 0.5),
#            strip.text.x = element_text(size = 11),
            legend.position = "none")
#            legend.title = element_blank())
p4

# save plot in '.svg' format
ggsave("sailfish.mis.pop.svg", device = "svg", width = 21, height = 15, units = "cm", dpi = 600)
ggsave("sailfish.mis.pop.png", device = "png", width = 21, height = 15, units = "cm", dpi = 600)

# patchwork
p1 / p3 + plot_annotation(tag_levels = 'a')
# save plot in '.svg' format
ggsave("het.mis.ind.svg", device = "svg", width = 21, height = 15, units = "cm", dpi = 600)
ggsave("het.mis.ind.png", device = "png", width = 21, height = 15, units = "cm", dpi = 600)

((p2 / p4) + plot_layout(guides = 'collect')) + plot_annotation(tag_levels = 'a')
# save plot in '.svg' format
ggsave("het.mis.pop.svg", device = "svg", width = 21, height = 15, units = "cm", dpi = 600)
ggsave("het.mis.pop.png", device = "png", width = 21, height = 15, units = "cm", dpi = 600)
