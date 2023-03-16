library(tess3r) #https://bcm-uga.github.io/TESS3_encho_sen/index.html
library(rnaturalearth) #https://github.com/ropensci/rnaturalearth
library(rnaturalearthdata)
library(tidyverse) #https://www.tidyverse.org/
library(patchwork) #https://patchwork.data-imaginist.com/index.html
library(data.table) #https://cran.r-project.org/web/packages/data.table/index.html

# set working directory
setwd("~/sailfish_cl/population.analysis/sailfish/tess3r/baq1")
# load saved data
load("~/sailfish_cl/population.analysis/sailfish/tess3r/baq1/tess3r.RData")

# load input files
# missing data MUST be changed to -9 in the structure file
#SFA <- read.table("sailfish.pruned.stru", header = FALSE, sep = "\t") # too slow for large data sets
sfa <- fread("sailfish.pruned.stru", header = FALSE)
coord <- read.table("coords.tsv", header = FALSE, sep = "\t")
# convert to tess3r input file format
snps <- tess2tess3(dataframe = sfa, TESS = FALSE, diploid = TRUE, FORMAT = 2,
                   extra.row = 1, extra.column = 2)
genotype <- snps$X
coordinates <- snps$coord
# run tess3r
sfa.tess <- tess3(X = genotype,
                  coord = coordinates,
                  K=1:6,
                  ploidy=2,
                  lambda = 1,
                  rep = 10,
                  W = NULL,
                  method = "projected.ls",
                  max.iteration = 10000,
                  tolerance = 1e-06,
                  openMP.core.num = 24,
                  Q.init = NULL,
                  mask = 0.2,
                  algo.copy = TRUE,
                  keep = "best",
                  verbose = TRUE)

qobj1 <- qmatrix(sfa.tess, 2, rep = "best")
qobj2 <- qmatrix(sfa.tess, 3, rep = "best")
qobj3 <- qmatrix(sfa.tess, 4, rep = "best")
qobj4 <- qmatrix(sfa.tess, 5, rep = "best")
qobj5 <- qmatrix(sfa.tess, 6, rep = "best")

#bp <- barplot(qobj1, border = NA, space = 0, xlab = "Individuals",
#              ylab = "Ancestry proportions", main = "Ancestry matrix")

#plot(coordinates, pch = 19, cex = .5, 
#     xlab = "Longitude (°E)", ylab = "Latitude (°N)")
#maps::map(add = T, interior = T)

p1 <- ggtess3Q(qobj1, coord, resolution = c(300, 300), window = NULL,
               background = FALSE, map.polygon = FALSE,
               interpolation.model = FieldsKrigModel(10),
               col.palette = CreatePalette(color.vector = c("#377EB8","#4DAF4A"),
                                           palette.length = 15))

p2 <- ggtess3Q(qobj2, coord, resolution = c(300, 300), window = NULL,
               background = FALSE, map.polygon = FALSE,
               interpolation.model = FieldsKrigModel(10),
               col.palette = CreatePalette(color.vector = c("#377EB8","#4DAF4A","#984EA3"),
                                           palette.length = 15))

p3 <- ggtess3Q(qobj3, coord, resolution = c(300, 300), window = NULL,
               background = FALSE, map.polygon = FALSE,
               interpolation.model = FieldsKrigModel(10),
               col.palette = CreatePalette(color.vector = c("#E41A1C","#377EB8","#4DAF4A","#984EA3"),
                                           palette.length = 15))

p4 <- ggtess3Q(qobj4, coord, resolution = c(300, 300), window = NULL,
               background = FALSE, map.polygon = FALSE,
               interpolation.model = FieldsKrigModel(10),
               col.palette = CreatePalette(color.vector = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00"),
                                           palette.length = 15))

p5 <- ggtess3Q(qobj5, coord, resolution = c(300, 300), window = NULL,
               background = FALSE, map.polygon = FALSE,
               interpolation.model = FieldsKrigModel(10),
               col.palette = CreatePalette(color.vector = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33"),
                                           palette.length = 15))

save.image(file = "tess3r.RData")

#map.polygon <- getMap(resolution = "high")
#
#p1 +
#  geom_path(data = map.polygon, aes(x = long, y = lat, group = group)) +
#  xlim(-60, 65) + 
#  ylim(-40, 25) + 
#  coord_equal() + 
#  geom_point(data = as.data.frame(coordinates), aes(x = V1, y = V2), size = 1) + 
#  xlab("Longitute") +
#  ylab("Latitude") + 
#  theme_bw()

# load data
world <- ne_countries(scale = "medium", returnclass = "sf")
# plot world map
#ggplot(data = world) + geom_sf()

# Plot the results
tp1 <- p1 + geom_sf(data = world, fill = "lightgray", lwd = 0.2) +
            coord_sf(xlim=c(-100, 165), ylim = c(-37, 33), expand = FALSE) + 
            geom_point(data = as.data.frame(coordinates), aes(x = V1, y = V2), size = 0.5) +
            labs(title = "K = 2", x = "Longitute", y = "Latitude") +
            theme_bw() + 
            theme(plot.title = element_blank(),
#                  plot.title = element_text(color="black", size=8,
#                                            face = "bold", hjust = 0.5,
#                                            margin = margin(t = 20, r = 10, b = 5, l = 10)),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  axis.text.x = element_blank())
#                  axis.text.x = element_blank(),
#                  axis.title.y = element_blank(),
#                  axis.text.y = element_blank(),
#                  axis.text.x = element_text(color="black", 
#                                             size=8),
#                  axis.text.y = element_text(color="black", 
#                                             size=8))
# save plot
#ggsave("SFA_K2_tess3r.png", device = "png", width = 15, height = 10, units = "cm", dpi = 600)
#ggsave("SFA_K2_tess3r.svg", device = "svg", width = 15, height = 10, units = "cm", dpi = 600)

tp2 <- p2 + geom_sf(data = world, fill = "lightgray", lwd = 0.2) +
            coord_sf(xlim=c(-100, 165), ylim = c(-37, 33), expand = FALSE) + 
            geom_point(data = as.data.frame(coordinates), aes(x = V1, y = V2), size = 0.5) +
            labs(title = "K = 3", x = "Longitute", y = "Latitude") +
            theme_bw() + 
            theme(plot.title = element_blank(),
#                  plot.title = element_text(color = "black", size = 8,
#                                            face = "bold", hjust = 0.5,
#                                            margin = margin(t = 20, r = 10, b = 5, l = 10)),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  axis.text.x = element_blank())
#                  axis.text.y = element_blank())
#                  axis.text.x = element_text(color="black", 
#                                             size=8),
#                  axis.text.y = element_text(color="black", 
#                                             size=8))
# save plot
#ggsave("SFA_K3_tess3r.png", device = "png", width = 15, height = 10, units = "cm", dpi = 600)
#ggsave("SFA_K3_tess3r.svg", device = "svg", width = 15, height = 10, units = "cm", dpi = 600)

tp3 <- p3 + geom_sf(data = world, fill = "lightgray", lwd = 0.2) +
            coord_sf(xlim=c(-100, 165), ylim = c(-37, 33), expand = FALSE) + 
            geom_point(data = as.data.frame(coordinates), aes(x = V1, y = V2), size = 0.5) +
            labs(title = "K = 4", x = "Longitute", y = "Latitude") +
            theme_bw() +
            theme(plot.title = element_blank(),
#                  plot.title = element_text(color = "black", size = 8,
#                                            face = "bold", hjust = 0.5,
#                                            margin = margin(t = 20, r = 10, b = 5, l = 10)),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.title.y = element_blank())
#                  axis.text.x = element_text(color="black", 
#                                             size=8),
#                  axis.text.y = element_text(color="black", 
#                                             size=8))
# save plot
#ggsave("SFA_K4_tess3r.png", device = "png", width = 15, height = 10, units = "cm", dpi = 600)
#ggsave("SFA_K4_tess3r.svg", device = "svg", width = 15, height = 10, units = "cm", dpi = 600)

tp4 <- p4 + geom_sf(data = world, fill = "lightgray", lwd = 0.2) +
            coord_sf(xlim=c(-100, 165), ylim = c(-37, 33), expand = FALSE) + 
            geom_point(data = as.data.frame(coordinates), aes(x = V1, y = V2), size = 0.5) +
            labs(title = "K = 5", x = "Longitute", y = "Latitude") +
            theme_bw() + 
            theme(plot.title = element_blank(),
#                  plot.title = element_text(color = "black", size = 8,
#                                            face = "bold", hjust = 0.5,
#                                            margin = margin(t = 20, r = 10, b = 5, l = 10)),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  axis.text.x = element_blank())
#                  axis.text.x = element_text(color="black", 
#                                             size=8))
#                  axis.text.y = element_text(color="black", 
#                                             size=8))
# save plot
#ggsave("SFA_K5_tess3r.png", device = "png", width = 15, height = 10, units = "cm", dpi = 600)
#ggsave("SFA_K5_tess3r.svg", device = "svg", width = 15, height = 10, units = "cm", dpi = 600)

tp5 <- p5 + geom_sf(data = world, fill = "lightgray", lwd = 0.2) +
            coord_sf(xlim=c(-100, 165), ylim = c(-37, 33), expand = FALSE) + 
            geom_point(data = as.data.frame(coordinates), aes(x = V1, y = V2), size = 0.5) +
            labs(title = "K = 6", x = "Longitute", y = "Latitude") +
            theme_bw() + 
            theme(plot.title = element_blank(),
#                  plot.title = element_text(color = "black", size = 8,
#                                            face = "bold", hjust = 0.5,
#                                            margin = margin(t = 20, r = 10, b = 5, l = 10)),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
#                  axis.text.x = element_blank())
                  axis.text.x = element_text(color="black",
                                             size=8),
                  axis.text.y = element_text(color="black", 
                                             size=8))
# save plot
#ggsave("SFA_K6_tess3r.png", device = "png", width = 15, height = 10, units = "cm", dpi = 600)
#ggsave("SFA_K6_tess3r.svg", device = "svg", width = 15, height = 10, units = "cm", dpi = 600)

# Merge plots
tp1 / tp2 / tp3 / tp4 / tp5 + plot_annotation(tag_levels = 'a')

# save plot
ggsave("SFA_tess3r.png", device = "png", width = 30, height = 15, units = "cm", dpi = 600)
ggsave("SFA_tess3r.svg", device = "svg", width = 30, height = 15, units = "cm", dpi = 600)

