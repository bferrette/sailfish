# --------------------------- #
#
# Tutorial:
# Seascape Redundancy Analysis
# 
# Description:
# Prepare environmental data for redundancy analysis.
#
# Environmental variables:
# Mean sea surface temperature (SST): Present-day (celsius)
# Mean sea surface salinity (SSS): Present-day (practical salinity scale)
# Mean sea surface currents velocity (SCV): Present-day (m-1)
# Mean sea surface primary productivity (SPP): Present-day (g.m-3.day-1)
# Mean pH (PH): Present-day
#  
# Data downloaded from https://www.bio-oracle.org
# A website containing marine data layers for ecological modelling.
# Files are in .asc format.
#
# --------------------------- #

# Load packages
library(raster)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(patchwork)

setwd('')

#--------------#
#
# Extract data
#
#--------------#

# Prepare data for extraction
sst.present <- raster("Present.Surface.Temperature.Mean.asc")
sss.present <- raster("Present.Surface.Salinity.Mean.asc")
scv.present <- raster("Present.Surface.Current.Velocity.Mean.asc.BOv2_1.asc")
spp.present <- raster("Present.Surface.Primary.productivity.Mean.asc")
ph.present <- raster("Present.Surface.pH.asc")

# Import coordinates of sites
coords <- read.delim("popdata.tsv", header = TRUE, sep = "\t")
names(coords)

# Create SpatialPoints object using coordinates
points <- SpatialPoints(subset(coords, select = c("lon","lat")))
plot(points)

# Extract environmental data for each site and combine into dataframe
df = data.frame(site = coords,
                sst_mean = extract(sst.present, points),
                sss_mean = extract(sss.present, points),
                scv_mean = extract(scv.present, points),
                spp_mean = extract(spp.present, points),
                ph_mean = extract(ph.present, points)
)

# Export data as a csv file
write.table(df, file="environmental_data.tsv", row.names = FALSE, sep = "\t")

#--------------#
#
# Plot heatmaps
#
#--------------#

# Set map boundary (xmin, xmax, ymin, ymax)
extent(points)
boundary <- extent(-80, 147, -32, 26)
boundary

# Crop rasters to boundary and convert to a dataframe of points
sst.df <- crop(sst.present, y = boundary) %>% rasterToPoints() %>% data.frame()
sss.df <- crop(sss.present, y = boundary) %>% rasterToPoints() %>% data.frame()
scv.df <- crop(scv.present, y = boundary) %>% rasterToPoints() %>% data.frame()
spp.df <- crop(spp.present, y = boundary) %>% rasterToPoints() %>% data.frame()
ph.df <- crop(ph.present, y = boundary) %>% rasterToPoints() %>% data.frame()

# Download a basemap
basemap <- ne_countries(scale = "large", returnclass = "sf")

# Define colour palettes
temp.cols = colorRampPalette(c("blue","white","red"))
sal.cols = colorRampPalette(c("darkred","white"))
scv.cols = colorRampPalette(c("white","green"))
calct.cols = colorRampPalette(c("white","#662506"))

# Sea surface temperature
sst.plt <- ggplot()+
           geom_raster(data = sst.df, aes(x = x, y = y, fill = sst.df[, 3]))+
           geom_sf(data = basemap, fill = "lightgray", lwd = 0.2)+
           coord_sf(xlim=c(-80, 147), ylim = c(-32, 26))+
           xlab("Longitude")+
           ylab("Latitude")+
           ggtitle("Sea surface temperature (present-day)")+
           scale_fill_gradientn(expression(~degree~C), colours = temp.cols(20))+
           theme_bw() +
           theme(axis.title = element_text(size = 12),
                 axis.text = element_text(size = 10, colour = "black"),
                 panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 12),
                 plot.title = element_text(size = 15, hjust = 0.5),
                 panel.grid = element_blank())
sst.plt
ggsave("sst_heatmap.png", device = "png", width = 15, height = 10, units = "cm", dpi = 600)
ggsave("sst_heatmap.svg", device = "svg", width = 15, height = 10, units = "cm", dpi = 600)

# Sea surface salinity
sss.plt = ggplot()+
          geom_raster(data = sss.df, aes(x = x, y = y, fill = sss.df[, 3]))+
          geom_sf(data = basemap, fill = "lightgray", lwd = 0.2)+
          coord_sf(xlim=c(-80, 147), ylim = c(-32, 26))+
          xlab("Longitude")+
          ylab("Latitude")+
          ggtitle("Sea surface salinity (present-day)")+
          scale_fill_gradientn("PPS", colours = sal.cols(10))+
          theme_bw() +
          theme(axis.title = element_text(size = 12),
                axis.text = element_text(size = 10, colour = "black"),
                panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
                legend.title = element_text(size = 13),
                legend.text = element_text(size = 12),
                plot.title = element_text(size = 15, hjust = 0.5),
                panel.grid = element_blank())
sss.plt
ggsave("sss_heatmap.png", device = "png", width = 15, height = 10, units = "cm", dpi = 600)
ggsave("sss_heatmap.svg", device = "svg", width = 15, height = 10, units = "cm", dpi = 600)

# Sea surface current velocity
scv.plt = ggplot()+
          geom_raster(data = scv.df, aes(x = x, y = y, fill = scv.df[, 3]))+
          geom_sf(data = basemap, fill = "lightgray", lwd = 0.2)+
          coord_sf(xlim=c(-80, 147), ylim = c(-32, 26))+
          xlab("Longitude")+
          ylab("Latitude")+
          ggtitle("Sea surface current velocity (present-day)")+
          scale_fill_gradientn("m-1", colours = scv.cols(10))+
          theme_bw() +
          theme(axis.title = element_text(size = 12),
                axis.text = element_text(size = 10, colour = "black"),
                panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
                legend.title = element_text(size = 13),
                legend.text = element_text(size = 12),
                plot.title = element_text(size = 15, hjust = 0.5),
                panel.grid = element_blank())
scv.plt
ggsave("scv_heatmap.png", device = "png", width = 15, height = 10, units = "cm", dpi = 600)
ggsave("scv_heatmap.svg", device = "svg", width = 15, height = 10, units = "cm", dpi = 600)

# Sea surface chlorophyll
ssc.plt = ggplot()+
  geom_tile(data = ssc.df, aes(x = x, y = y, fill = ssc.df[, 3]))+
  geom_polygon(data = basemap, aes(x = long, y = lat, group = group))+
  coord_quickmap(expand = F)+
  xlab("Longitude")+
  ylab("Latitude")+
  ggtitle("Sea surface chlorophyll (present-day)")+
  scale_fill_gradientn(expression(paste("mg/m"^"3")), colours = chlor.cols(10))+
  ggtheme
ssc.plt
ggsave("5.ssc_heatmap.png", width = 10, height = 9, dpi = 600)

# Sea surface calcite
ssca.plt = ggplot()+
  geom_tile(data = ssca.df, aes(x = x, y = y, fill = ssca.df[, 3]))+
  geom_polygon(data = basemap, aes(x = long, y = lat, group = group))+
  coord_quickmap(expand = F)+
  xlab("Longitude")+
  ylab("Latitude")+
  ggtitle("Sea surface calcite (present-day)")+
  scale_fill_gradientn(expression(paste("mol/m"^"3")), colours = calct.cols(10))+
  ggtheme
ssca.plt
ggsave("6.ssca_heatmap.png", width = 10, height = 9, dpi = 600)

# Combine two temperature ggplots
figAB = ggarrange(sst.plt + labs(tag = "A") + ggtheme + theme(axis.title.y = element_blank()),
                  sbt.plt + labs(tag = "B") + ggtheme + theme(axis.title.y = element_blank()),
                  ncol = 2, common.legend = TRUE, legend = "right")
figAB = annotate_figure(figAB,
                        left = text_grob("Latitude", size = 12, rot = 90))

# Combine two salinity ggplots
figCD = ggarrange(sss.plt + labs(tag = "C") + ggtheme + theme(axis.title.y = element_blank()),
                  sbs.plt + labs(tag = "D") + ggtheme + theme(axis.title.y = element_blank()),
                  ncol = 2, common.legend = TRUE, legend = "right")
figCD = annotate_figure(figCD,
                        left = text_grob("Latitude", size = 12, rot = 90))

# Combine temperature and salinity ggplots
fig = ggarrange(figAB, figCD, nrow = 2)
ggsave("7.temp_sal_heatmap.png", width = 10, height = 10, dpi = 600)
# ggsave("6.temp_sal_heatmap.pdf", width = 10, height = 10)
