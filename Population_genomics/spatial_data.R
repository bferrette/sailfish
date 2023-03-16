# --------------------------- #
#
# Tutorial:
# Seascape Redundancy Analysis
# 
# Description:
# Prepare spatial data for redundancy analysis.
#
# Load packages
library(tidyverse)
library(raster)
library(marmap)
library(ade4)
library(adespatial)
library(SoDA)
library(GGally)
library(distances)
library(vegan)

# set working directory
setwd("~/Documents/billfishes/sailfish/genome/sailfish_cl/populations.analysis/populational/seascape/spatial_data")

# Import coordinates of sites
coords <- read.delim("popdata.tsv", header = TRUE, sep = "\t")
coords

#--------------#
#
# Compute dbMEMs
#
#--------------#

# Further info on distance-based Moran's Eigenvector Maps in the vignette
# vignette("tutorial", package = "adespatial")

# Transform geographic coordinates into cartesian coordinates
cart <- geoXY(coords$lat, coords$lon, unit = 1000)
cart
plot(cart[,1], cart[,2], asp = 1)
#plot(coords$lon,coords$lat)
# normalize euclidian distances
#norm.cart <- distances(cart, normalize = "mahalanobize")
#norm.cart
# Calculate euclidian distances
euclidean_distances <- dist(cart, method = "euclidean") 
euclidean_distances
# Compute distance-based Moran's Eigenvector Maps
dbmems <- dbmem(euclidean_distances, thresh = NULL, MEM.autocor="positive") # "positive", "non-null", "all", "negative"
dbmems
#dbmems.dec <- decostand(dbmems, method = "hellinger", MARGIN=2, range.global, na.rm=FALSE)
#dbmems.dec
# Plot dbMEM eigenvalues
# Matrix of plots (method "loess" or "lm")
my_fn <- function(data, mapping, method="loess", ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=method, ...)
  p
}
ggpairs(dbmems, lower = list(continuous = my_fn))
ggsave("dbmem.ggpairs.png", device = "png", width = 35, height = 30, units = "cm", dpi = 600)
ggsave("dbmem.ggpairs.svg", device = "svg", width = 35, height = 30, units = "cm", dpi = 600)

ggcorr(dbmems, label = TRUE, label_round = 2)
ggsave("dbmem.ggcorr.png", device = "png", width = 35, height = 30, units = "cm", dpi = 600)
ggsave("dbmem.ggcorr.svg", device = "svg", width = 35, height = 30, units = "cm", dpi = 600)

# Export dbMEMs
write.table(dbmems, file = "dbmems.tsv", row.names = FALSE, sep = "\t")
