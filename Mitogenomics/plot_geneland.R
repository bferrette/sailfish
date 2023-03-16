library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
# set working directotory
setwd("~/Documents/billfishes/mitogenomes/sailfish/analysis/geneland")
# Posterior probability of cluster membership for sampled individuals
proba.pop.membership.indiv <- read.table("proba.pop.membership.indiv.txt", header = FALSE)

ggplot(proba.pop.membership.indiv) + geom_point(aes(x=V1,y=V2))

# Create color palette
heat_colors_interpolated <- colorRampPalette(paletteer::paletteer_d("RColorBrewer::YlOrRd", n = 9, direction = -1))(15)
colors_interpolated <- colorRampPalette(paletteer::paletteer_d("RColorBrewer::Set1", n = 9, direction = -1))(15)
# Plot Geneland posterior probability of cluster membership for sampled individuals
p0 <- ggplot(proba.pop.membership.indiv)+
      geom_density_2d_filled(aes(x=V1,y=V2)) +
      scale_fill_manual(values = c(heat_colors_interpolated),
                        aesthetics = c("fill", "color")) +
#     scale_alpha_continuous(limits=c(0,0.2),breaks=seq(0,0.2,by=0.025))+
      theme_bw() +
      theme(legend.position="none")
p0

# load data
world <- ne_countries(scale = "large", returnclass = "sf")
# gene world map
ggplot(data = world) + geom_sf()
# plot proba.pop.membership.indiv
tp1 <- p0 + geom_sf(data = world, fill = "lightgray", lty = 1, alpha = 0.2) +
            coord_sf(xlim=c(-75, 140), ylim = c(-30, 30), expand = FALSE) +
            geom_point(data = proba.pop.membership.indiv,
                       aes(x = V1, y = V2), size = 0.5) +
            xlab("Longitute") +
            ylab("Latitude") + 
            theme_bw() +
            theme(legend.position="none",
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank())
tp1

# save plot in '.svg' format
ggsave("proba.pop.membership.indiv.svg", device = "svg", width = 15, height = 8, units = "cm", dpi = 600)
ggsave("proba.pop.membership.indiv.png", device = "png", width = 15, height = 8, units = "cm", dpi = 600)

# Posterior probability of cluster membership for populations
proba.pop.membership <- read.table("proba.pop.membership.txt", header = TRUE)

ggplot(proba.pop.membership) + geom_tile(aes(x=V1,y=V2))

# Create color palette
heat_colors_interpolated <- colorRampPalette(paletteer::paletteer_d("RColorBrewer::YlOrRd", n = 9, direction = -1))(15)

# Plot Geneland posterior probability of cluster membership for sampled individuals
p2 <- ggplot(proba.pop.membership) +
      geom_tile(aes(x=V1,y=V2,fill=V3)) +
      scale_fill_distiller(palette = "YlOrRd", direction = 1)
#      scale_fill_viridis_c(option = "magma", direction = -1)
p2

# load data
world <- ne_countries(scale = "large", returnclass = "sf")
# gene world map
ggplot(data = world) + geom_sf()
# plot proba.pop.membership
tp2 <- p2 + geom_sf(data = world, fill = "lightgray", lty = 1, alpha = 0.2) +
            coord_sf(xlim=c(-75, 140), ylim = c(-30, 30), expand = FALSE) +
            geom_point(data = proba.pop.membership.indiv,
                       aes(x = V1, y = V2), size = 0.5) +
            xlab("Longitute") +
            ylab("Latitude") + 
            theme_bw() +
            theme(legend.position="right",
                  legend.title = element_blank(),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank())
tp2

# save plot in '.svg' format
ggsave("proba.pop.membership.svg", device = "svg", width = 15, height = 8, units = "cm", dpi = 600)
ggsave("proba.pop.membership.pop.png", device = "png", width = 15, height = 8, units = "cm", dpi = 600)

# plot 
# load membership probabilities
k2 <- read.delim("prob.tsv", header = TRUE, sep = "\t")
meta.data <- read.delim("popdata.tsv", header = TRUE, sep = "\t")
# read input Q-matrix
# create data frame
  assign(paste0("k", ".data"),
         k2 %>%
           mutate(pop = meta.data$pop,
                  samples = meta.data$sample
           ) %>%
           pivot_longer(cols = starts_with("C"),
                        names_to = "ancestry",
                        values_to = "probability")
  )

# plot
p3 <- ggplot(k.data, aes(x = samples, y = probability, fill = ancestry)) +
             geom_col(color = "white", size = 0.2, width = 1) +
#             coord_flip() +
             facet_grid(~fct_inorder(pop),
#                        switch = "x",
                        scales = 'free',
                        space = 'free') +
             ylab("K = 2") +
             theme_bw() +
             theme(panel.border = element_blank(),
                   panel.spacing.x = unit(0.01, "lines"),
                   axis.title.x = element_blank(),
                   axis.title.y = element_text(size = 10),
                   axis.text.x = element_text(size = 7, vjust = 0.5, hjust = 1, angle = 90),
                   axis.text.y = element_text(size = 7, hjust = 1),
                   strip.background = element_blank(),
                   strip.text = element_text(size = 10, vjust = 0, hjust = 0.5),
                   panel.grid = element_blank(),
                   strip.placement = NULL) +
            scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
            scale_x_discrete(expand = expansion(add = 0.75)) +
#             scale_fill_manual(values = c("#4DAF4A", "#984EA3"),
#                               guide = "none")
            scale_fill_manual(values = c("#2171B5","#6A51A3"),
                              guide = "none")
p3
ggsave("ancestry.svg", device = "svg", width = 15, height = 8, units = "cm", dpi = 600)
ggsave("ancestry.png", device = "png", width = 15, height = 8, units = "cm", dpi = 600)

# merge plots
( tp1 / tp2 / p3 ) + plot_layout(heights = c(1,1,1)) + plot_annotation(tag_levels = "a")
# save plot
ggsave("geneland.svg", device = "svg", width = 30, height = 25, units = "cm", dpi = 600)
ggsave("geneland.png", device = "png", width = 30, height = 25, units = "cm", dpi = 600)
