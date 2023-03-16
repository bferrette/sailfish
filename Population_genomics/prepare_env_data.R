# Description:
# Prepare environmental data for redundancy analysis.
#
# Environmental variables:
# Mean sea surface temperature (SST): Present-day (celsius)
# Mean sea surface salinity (SSS): Present-day (practical salinity scale)
# Mean sea surface currents velocity (scvl): Present-day (m-1)
# Mean sea surface primary productivity (sprp): Present-day (g.m-3.day-1)
# Mean pH (PH): Present-day
#  
# Data downloaded from https://www.bio-oracle.org
# A website containing marine data layers for ecological modelling.
# Files are in .asc format.
#
# --------------------------- #

# Load packages
library(adespatial)
library(caret)
library(GGally)
library(ggpubr)
library(patchwork)
library(raster)
library(RColorBrewer)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(sdmpredictors)
library(tidyverse)
library(vegan)
library(viridis)

#set working directory
setwd("/home/bferrette/Documents/billfishes/sailfish/genome/sailfish_cl/populations.analysis/populational/seascape/env_data")
# load environment
load("/home/bferrette/Documents/billfishes/sailfish/genome/sailfish_cl/populations.analysis/populational/seascape/env_data/env.RData")
# Extract marine environmental data from Bio-Oracle v2.2 - https://www.bio-oracle.org/index.php
# List marine data sets
datasets = list_datasets(terrestrial = FALSE, marine = TRUE)
# List layers avaialble in Bio-ORACLE v2
layers.bio2 <- list_layers( datasets="Bio-ORACLE" )
# Prepare data for extraction
bathy <- load_layers("BO_bathymax", datadir = tempdir("/home/bferrette/Documents/billfishes/sailfish/genome/sailfish_cl/populations.analysis/populational/seascape/env_data"))
chlo <- load_layers("BO22_chlomean_ss", datadir = tempdir("/home/bferrette/Documents/billfishes/sailfish/genome/sailfish_cl/populations.analysis/populational/seascape/env_data"))
ph <- load_layers("BO22_ph", datadir = tempdir("/home/bferrette/Documents/billfishes/sailfish/genome/sailfish_cl/populations.analysis/populational/seascape/env_data"))
ptpk <- load_layers("BO22_carbonphytomean_ss", datadir = tempdir("/home/bferrette/Documents/billfishes/sailfish/genome/sailfish_cl/populations.analysis/populational/seascape/env_data"))
scvl <- load_layers("BO22_curvelmean_ss", datadir = tempdir("/home/bferrette/Documents/billfishes/sailfish/genome/sailfish_cl/populations.analysis/populational/seascape/env_data"))
sdO2 <- load_layers("BO22_dissoxltmax_bdmin", datadir = tempdir("/home/bferrette/Documents/billfishes/sailfish/genome/sailfish_cl/populations.analysis/populational/seascape/env_data"))
sprp <- load_layers("BO22_ppmean_ss", datadir = tempdir("/home/bferrette/Documents/billfishes/sailfish/genome/sailfish_cl/populations.analysis/populational/seascape/env_data"))
ssal <- load_layers("BO22_salinitymean_ss", datadir = tempdir("/home/bferrette/Documents/billfishes/sailfish/genome/sailfish_cl/populations.analysis/populational/seascape/env_data"))
sst <- load_layers("BO22_tempmean_ss", datadir = tempdir("/home/bferrette/Documents/billfishes/sailfish/genome/sailfish_cl/populations.analysis/populational/seascape/env_data"))
light <- load_layers("BO22_lightbotltmax_bdmean", datadir = tempdir("/home/bferrette/Documents/billfishes/sailfish/genome/sailfish_cl/populations.analysis/populational/seascape/env_data"))
diff.atn <- load_layers("BO22_damean", datadir = tempdir("/home/bferrette/Documents/billfishes/sailfish/genome/sailfish_cl/populations.analysis/populational/seascape/env_data"))
# Import coordinates of sites
coords <- read.delim("popdata.tsv", header = TRUE, sep = "\t")
names(coords)
# Create SpatialPoints object using coordinates
points <- SpatialPoints(subset(coords, select = c("lon","lat")))
plot(points)
# Extract environmental data for each site and combine into dataframe
df <- data.frame(points,
                 bathy = raster::extract(bathy, points, buffer=100, fun=mean, small=TRUE, exact=TRUE, weights=TRUE, normalizeWeights=TRUE),
                 chlo = raster::extract(chlo, points, buffer=100, fun=mean, small=TRUE, exact=TRUE, weights=TRUE, normalizeWeights=TRUE),
                 scvl = raster::extract(scvl, points, buffer=100, fun=mean, small=TRUE, exact=TRUE, weights=TRUE, normalizeWeights=TRUE),
                 sdO2 = raster::extract(sdO2, points, buffer=100, fun=mean, small=TRUE, exact=TRUE, weights=TRUE, normalizeWeights=TRUE),
                 ph = raster::extract(ph, points, buffer=100, fun=mean, small=TRUE, exact=TRUE, weights=TRUE, normalizeWeights=TRUE),
                 ptpk = raster::extract(ptpk, points, buffer=100, fun=mean, small=TRUE, exact=TRUE, weights=TRUE, normalizeWeights=TRUE),
                 sprp = raster::extract(sprp, points, buffer=100, fun=mean, small=TRUE, exact=TRUE, weights=TRUE, normalizeWeights=TRUE),
                 ssal = raster::extract(ssal, points, buffer=100, fun=mean, small=TRUE, exact=TRUE, weights=TRUE, normalizeWeights=TRUE),
                 light = raster::extract(light, points, buffer=100, fun=mean, small=TRUE, exact=TRUE, weights=TRUE, normalizeWeights=TRUE),
                 diff.atn = raster::extract(diff.atn, points, buffer=100, fun=mean, small=TRUE, exact=TRUE, weights=TRUE, normalizeWeights=TRUE),
                 sst = raster::extract(sst, points, buffer=100, fun=mean, small=TRUE, exact=TRUE, weights=TRUE, normalizeWeights=TRUE))
names(df)
# rename columns names
colnames(df) <- c("lon","lat","bathy","chlo","scvl","sdO2","ph","ptpk","sprp","ssal","light","diff.atn","sst")
# deconstand provides some popular (and effective) standardization methods for community ecologists
env <- decostand(df[3:13], method="standardize", na.rm=TRUE)
env
#env2 <- scale(df[3:13], center = TRUE, scale = TRUE)
#env2
# Export data as a csv file
write.table(env, file="environmental_rawdata.tsv", row.names = FALSE, sep = "\t")
#write.table(df, file="environmental_data_2.tsv", row.names = FALSE, sep = "\t")

# Correlation
# https://stackoverflow.com/questions/35085261/how-to-use-loess-method-in-ggallyggpairs-using-wrap-function/35088740
# Matrix of plots (method "loess" or "lm")
my_fn <- function(data, mapping, method="lm", ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=method, ...)
  p
}
# Default loess curve
ggpairs(env, lower = list(continuous = my_fn)) + theme_bw()
ggsave("colinearity.png", device = "png", width = 35, height = 30, units = "cm", dpi = 600)
ggsave("colinearity.svg", device = "svg", width = 35, height = 30, units = "cm", dpi = 600)
# Correlation matrix plot
ggcorr(env[, 3:13], method = c("pairwise", "pearson"), nbreaks = 10,
       label = TRUE, label_size = 3, label_round = 2, label_alpha = TRUE)
ggsave("colinearity2.png", device = "png", width = 35, height = 30, units = "cm", dpi = 600)
ggsave("colinearity2.svg", device = "svg", width = 35, height = 30, units = "cm", dpi = 600)
# Variable selection (constrained ordination)
# https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
df2 <- cor(env, use = "all.obs", method = "pearson")
#df2 <- cor(env, use = "pairwise.complete.obs", method = "pearson")
hc <- findCorrelation(df2, cutoff=0.3, names = TRUE, exact = TRUE) # putt any value as a "cutoff" 
hc <- sort(hc)
hc
env.data <- subset(env, select = -c(bathy,chlo,diff.atn,ph,ptpk,scvl,sdO2,sst))
env.data
# Variable selection on selected variables
#df3 <- cor(env.data, use = "all.obs", method = "pearson")
df3 <- cor(env.data, use = "pairwise.complete.obs", method = "pearson")
hc2 <- findCorrelation(df3, cutoff=0.3, names = TRUE, exact = TRUE) # putt any value as a "cutoff" 
hc2 <- sort(hc2)
hc2
# Export data as a csv file
write.table(env.data, file="environmental_data.tsv", row.names = FALSE, sep = "\t")

#--------------#
#
# Plot heatmaps
#
#--------------#
# Set map boundary (xmin, xmax, ymin, ymax)
raster::extent(points)
boundary <- raster::extent(-80, 147, -45, 26)
boundary

# Crop rasters to boundary and convert to a dataframe of points
bathy.df <- crop(bathy, y = boundary) %>% rasterToPoints() %>% data.frame()
chlo.df <- crop(chlo, y = boundary) %>% rasterToPoints() %>% data.frame()
scvl.df <- crop(scvl, y = boundary) %>% rasterToPoints() %>% data.frame()
sdO2.df <- crop(sdO2, y = boundary) %>% rasterToPoints() %>% data.frame()
ph.df <- crop(ph, y = boundary) %>% rasterToPoints() %>% data.frame()
ptpk.df <- crop(ptpk, y = boundary) %>% rasterToPoints() %>% data.frame()
sprp.df <- crop(sprp, y = boundary) %>% rasterToPoints() %>% data.frame()
ssal.df <- crop(ssal, y = boundary) %>% rasterToPoints() %>% data.frame()
light.df <- crop(light, y = boundary) %>% rasterToPoints() %>% data.frame()
diff.atn.df <- crop(diff.atn, y = boundary) %>% rasterToPoints() %>% data.frame()
sst.df <- crop(sst, y = boundary) %>% rasterToPoints() %>% data.frame()

# Download a world
world <- ne_countries(scale = "small", returnclass = "sf")

# Sea surface temperature
sst.plot <- ggplot()+
            geom_tile(data = sst.df, aes(x = x, y = y, fill = sst.df[, 3])) +
            geom_sf(data = world, fill = "lightgray", lwd = 0.2) +
            coord_sf(xlim=c(-80, 147), ylim = c(-45, 26), expand = FALSE) +
            xlab("Longitude") +
            ylab("Latitude") +
            ggtitle("Mean sea surface temperature (present-day)") +
            scale_fill_distiller(palette = "Spectral",expression(~degree~C)) +
#            scale_fill_gradientn(expression(~degree~C), colours = temp.cols(20)) +
            theme_bw() +
            theme(axis.text.y = element_text(size = 8, colour = "black"),
                  axis.title.y = element_blank(),
                  axis.text.x = element_blank(),
                  axis.title.x = element_blank(),
                  panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
                  legend.title = element_text(size = 8),
                  legend.text = element_text(size = 8),
                  plot.title = element_text(size = 8, hjust = 0.5),
                  panel.grid = element_blank())
#sst.plot
#ggsave("sst_heatmap.png", device = "png", width = 25, height = 15, units = "cm", dpi = 600)
#ggsave("sst_heatmap.svg", device = "svg", width = 25, height = 15, units = "cm", dpi = 600)

# Sea surface temperature
light.plot <- ggplot()+
              geom_tile(data = light.df, aes(x = x, y = y, fill = light.df[, 3])) +
              geom_sf(data = world, fill = "lightgray", lwd = 0.2) +
              coord_sf(xlim=c(-80, 147), ylim = c(-45, 26), expand = FALSE) +
              xlab("Longitude") +
              ylab("Latitude") +
              ggtitle("Photosynthetically active radiation (PAR) (present-day)") +
              scale_fill_distiller(palette = "Spectral",expression(paste("E.m"^-2~".day"^-1))) +
#              scale_fill_gradientn(expression(~degree~C), colours = temp.cols(20)) +
              theme_bw() +
              theme(axis.text.y = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.x = element_blank(),
                    axis.title.x = element_blank(),
                    panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
                    legend.title = element_text(size = 8),
                    legend.text = element_text(size = 8),
                    plot.title = element_text(size = 8, hjust = 0.5),
                    panel.grid = element_blank())
#light.plot
#ggsave("light_heatmap.png", device = "png", width = 25, height = 15, units = "cm", dpi = 600)
#ggsave("light_heatmap.svg", device = "svg", width = 25, height = 15, units = "cm", dpi = 600)

# Sea surface temperature
diff.atn.plot <- ggplot()+
                 geom_tile(data = diff.atn.df, aes(x = x, y = y, fill = diff.atn.df[, 3])) +
                 geom_sf(data = world, fill = "lightgray", lwd = 0.2) +
                 coord_sf(xlim=c(-80, 147), ylim = c(-45, 26), expand = FALSE) +
                 xlab("Longitude") +
                 ylab("Latitude") +
                 ggtitle("Diffuse attenuation (present-day)") +
                 scale_fill_distiller(palette = "Spectral",expression(paste("m"^-1))) +
#                 scale_fill_gradientn(expression(~degree~C), colours = temp.cols(20)) +
                 theme_bw() +
                 theme(axis.text.y = element_blank(),
                       axis.title.y = element_blank(),
                       axis.text.x = element_blank(),
                       axis.title.x = element_blank(),
                       panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
                       legend.title = element_text(size = 8),
                       legend.text = element_text(size = 8),
                       plot.title = element_text(size = 8, hjust = 0.5),
                       panel.grid = element_blank())
#diff.atn.plot
#ggsave("diff.atn_heatmap.png", device = "png", width = 25, height = 15, units = "cm", dpi = 600)
#ggsave("diff.atn_heatmap.svg", device = "svg", width = 25, height = 15, units = "cm", dpi = 600)

# Sea surface dissolved molecular oxygen
sdO2.plot <- ggplot() +
             geom_tile(data = sdO2.df, aes(x = x, y = y, fill = sdO2.df[, 3])) +
             geom_sf(data = world, fill = "lightgray", lwd = 0.2) +
             coord_sf(xlim=c(-80, 147), ylim = c(-45, 26), expand = FALSE) +
             xlab("Longitude") +
             ylab("Latitude") +
             ggtitle("Sea surface dissolved molecular oxygen (present-day)") +
             scale_fill_distiller(palette = "PuOr", expression(paste("mol.m"^"-3")), direction = 1) +
#             scale_fill_gradientn(expression(paste("mg/m"^"3")), colours = chlor.cols(10)) +
             theme_bw() +
             theme(axis.title = element_blank(),
                   axis.text.x = element_blank(),
                   axis.text.y = element_text(size = 8, colour = "black"),
                   panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
                   legend.title = element_text(size = 8),
                   legend.text = element_text(size = 8),
                   plot.title = element_text(size = 8, hjust = 0.5),
                   panel.grid = element_blank())
#sdO2.plot
#ggsave("sdO2_heatmap.png", device = "png", width = 25, height = 15, units = "cm", dpi = 600)
#ggsave("sdO2_heatmap.svg", device = "svg", width = 25, height = 15, units = "cm", dpi = 600)

# Sea surface current velocity
scvl.plot <- ggplot() +
             geom_tile(data = scvl.df, aes(x = x, y = y, fill = scvl.df[, 3])) +
             geom_sf(data = world, fill = "lightgray", lwd = 0.2) +
             coord_sf(xlim=c(-80, 147), ylim = c(-45, 26), expand = FALSE) +
             xlab("Longitude") +
             ylab("Latitude") +
             ggtitle("Sea surface current velocity (present-day)") +
             scale_fill_distiller(palette = "Spectral","m/s", direction = -1) +
#             scale_fill_gradientn("m-1", colours = scvl.cols(10)) +
             theme_bw() +
             theme(axis.title = element_blank(),
                   axis.text.x = element_blank(),
                   axis.text.y = element_blank(),
                   panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
                   legend.title = element_text(size = 8),
                   legend.text = element_text(size = 8),
                   plot.title = element_text(size = 8, hjust = 0.5),
                   panel.grid = element_blank())
#scvl.plot
#ggsave("scvl_heatmap.png", device = "png", width = 25, height = 15, units = "cm", dpi = 600)
#ggsave("scvl_heatmap.svg", device = "svg", width = 25, height = 15, units = "cm", dpi = 600)

# Bathymetry
bathy.plot <- ggplot() +
              geom_tile(data = bathy.df, aes(x = x, y = y, fill = bathy.df[, 3])) +
              geom_sf(data = world, fill = "lightgray", lwd = 0.2) +
              coord_sf(xlim=c(-80, 147), ylim = c(-45, 26), expand = FALSE) +
              xlab("Longitude") +
              ylab("Latitude") +
              ggtitle("Mean bathymetry (present-day)") +
              scale_fill_distiller(palette = "Blues", expression(paste(m)), direction = -1) +
#              scale_fill_gradientn(expression(~degree~C), colours = temp.cols(20)) +
              theme_bw() +
              theme(axis.text.y = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.x = element_blank(),
                    axis.title.x = element_blank(),
                    panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
                    legend.title = element_text(size = 8),
                    legend.text = element_text(size = 8),
                    plot.title = element_text(size = 8, hjust = 0.5),
                    panel.grid = element_blank())
#bathy.plot
#ggsave("mean.bathymetry.png", device = "png", width = 25, height = 15, units = "cm", dpi = 600)
#ggsave("mean.bathymetry.svg", device = "svg", width = 25, height = 15, units = "cm", dpi = 600)

# Sea surface salinity
ssal.plot <- ggplot() +
             geom_tile(data = ssal.df, aes(x = x, y = y, fill = ssal.df[, 3])) +
             geom_sf(data = world, fill = "lightgray", lwd = 0.2) +
             coord_sf(xlim=c(-80, 147), ylim = c(-45, 26), expand = FALSE) +
             xlab("Longitude") +
             ylab("Latitude") +
             ggtitle("Sea surface salinity (present-day)") +
             scale_fill_distiller("PPS", palette = "PuOr", direction = 1) +
#             scale_fill_gradientn("PPS", colours = sal.cols(10)) +
             theme_bw() +
             theme(axis.title = element_blank(),
                   axis.text.x = element_blank(), 
                   axis.text.y = element_text(size = 8, colour = "black"),
                   panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
                   legend.title = element_text(size = 8),
                   legend.text = element_text(size = 8),
                   plot.title = element_text(size = 8, hjust = 0.5),
                   panel.grid = element_blank())
#ssal.plot
#ggsave("sss_heatmap.png", device = "png", width = 25, height = 15, units = "cm", dpi = 600)
#ggsave("sss_heatmap.svg", device = "svg", width = 25, height = 15, units = "cm", dpi = 600)

# Sea pH
ph.plot <- ggplot() +
           geom_tile(data = ph.df, aes(x = x, y = y, fill = ph.df[, 3])) +
           geom_sf(data = world, fill = "lightgray", lwd = 0.2) +
           coord_sf(xlim=c(-80, 147), ylim = c(-45, 26), expand = FALSE) +
           xlab("Longitude") +
           ylab("Latitude") +
           ggtitle("Sea surface pH (present-day)") +
           scale_fill_distiller(palette = "BrBG", "pH", direction = 1) +
#            scale_fill_gradientn(expression(paste("mol/m"^"3")), colours = calct.cols(10)) +
           theme_bw() +
           theme(axis.title = element_blank(),
                 axis.text.x = element_blank(),
                 axis.text.y = element_blank(),
                 panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
                 legend.title = element_text(size = 8),
                 legend.text = element_text(size = 8),
                 plot.title = element_text(size = 8, hjust = 0.5),
                 panel.grid = element_blank())
#ph.plot
#ggsave("pH_heatmap.png", device = "png", width = 25, height = 15, units = "cm", dpi = 600)
#ggsave("pH_heatmap.svg", device = "svg", width = 25, height = 15, units = "cm", dpi = 600)

# Sea surface Chlorophyll
chlo.plot <- ggplot()+
             geom_tile(data = chlo.df, aes(x = x, y = y, fill = chlo.df[, 3])) +
             geom_sf(data = world, fill = "lightgray", lwd = 0.2 )+
             coord_sf(xlim=c(-80, 147), ylim = c(-45, 26), expand = FALSE) +
             xlab("Longitude") +
             ylab("Latitude") +
             ggtitle("Sea surface Chlorophyll (present-day)") +
             scale_fill_distiller(palette = "RdYlGn", expression(paste(mg.m^-3)), direction = -1) +
#             scale_fill_gradientn(expression(~degree~C), colours = temp.cols(20)) +
             theme_bw() +
             theme(axis.text.y = element_blank(),
                   axis.title.y = element_blank(),
                   axis.text.x = element_blank(),
                   axis.title.x = element_blank(),
                   panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
                   legend.title = element_text(size = 8),
                   legend.text = element_text(size = 8),
                   plot.title = element_text(size = 8, hjust = 0.5),
                   panel.grid = element_blank())
#chlo.plot
#ggsave("clr_heatmap.png", device = "png", width = 25, height = 15, units = "cm", dpi = 600)
#ggsave("clr_heatmap.svg", device = "svg", width = 25, height = 15, units = "cm", dpi = 600)

# Sea surface Phytoplankton
ptpk.plot <- ggplot()+
             geom_tile(data = ptpk.df, aes(x = x, y = y, fill = ptpk.df[, 3])) +
             geom_sf(data = world, fill = "lightgray", lwd = 0.2) +
             coord_sf(xlim=c(-80, 147), ylim = c(-45, 26), expand = FALSE) +
             xlab("Longitude") +
             ylab("Latitude") +
             ggtitle("Sea surface Phytoplankton (present-day)") +
             scale_fill_distiller(palette = "RdYlGn", expression(paste(umol.m^-3)), direction = -1) +
#            scale_fill_gradientn(expression(~degree~C), colours = temp.cols(20)) +
             theme_bw() +
             theme(axis.title = element_blank(),
                   axis.text.x = element_text(size = 8, colour = "black"),
                   axis.text.y = element_text(size = 8, colour = "black"),
                   panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
                   legend.title = element_text(size = 8),
                   legend.text = element_text(size = 8),
                   plot.title = element_text(size = 8, hjust = 0.5),
                   panel.grid = element_blank())
#ptpk.plot
#ggsave("ptpk_heatmap.png", device = "png", width = 25, height = 15, units = "cm", dpi = 600)
#ggsave("ptpk_heatmap.svg", device = "svg", width = 25, height = 15, units = "cm", dpi = 600)

# Sea surface primary productivity
sprp.plot <- ggplot() +
             geom_tile(data = sprp.df, aes(x = x, y = y, fill = sprp.df[, 3])) +
             geom_sf(data = world, fill = "lightgray", lwd = 0.2) +
             coord_sf(xlim=c(-80, 147), ylim = c(-45, 26), expand = FALSE) +
             xlab("Longitude") +
             ylab("Latitude") +
             ggtitle("Sea surface primary productivity (present-day)") +
             scale_fill_distiller(palette = "RdYlGn", expression(paste(mol.m^-3~.day^-1)), direction = -1) +
#             scale_fill_gradientn(expression(paste("mol/m"^"3")), colours = calct.cols(10)) +
             theme_bw() +
             theme(axis.title = element_blank(),
                   axis.text.x = element_text(size = 8, colour = "black"),
                   axis.text.y = element_blank(), 
                   panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
                   legend.title = element_text(size = 8),
                   legend.text = element_text(size = 8),
                   plot.title = element_text(size = 8, hjust = 0.5),
                   panel.grid = element_blank())
#sprp.plot
#ggsave("sprp_heatmap.png", device = "png", width = 25, height = 15, units = "cm", dpi = 600)
#ggsave("sprp_heatmap.svg", device = "svg", width = 25, height = 15, units = "cm", dpi = 600)

# Multiplot of environmental variables
((sst.plot + diff.atn.plot + light.plot) / (sdO2.plot + scvl.plot + bathy.plot) / (ssal.plot + ph.plot + chlo.plot) / (ptpk.plot + sprp.plot + plot_spacer()) +
  plot_layout(nrow = 4))

  plot_annotation(tag_levels = "a")
  
# save multiplot
ggsave("env_data.png", device = "png", width = 45, height = 20, units = "cm", dpi = 600)
#ggsave("env_data.svg", device = "svg", width = 45, height = 20, units = "cm", dpi = 600)
# save R environment
save.image(file = "env.RData")
