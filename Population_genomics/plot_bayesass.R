# Load packages
library(circlize)
library(tidyverse)
library(RColorBrewer)

# set working directory
setwd("~/Documents/billfishes/sailfish/genome/sailfish_cl/populations.analysis/populational/bayesass.reduced")

# Prepare BayesAss data
baysass <- tibble(Region = c("WCA","SWA","ECA","SWI","SEI","WCP"),
                  WCA = c(0.6875,0.0208,0.2292,0.0208,0.0208,0.0208),
                  SWA = c(0.0209,0.6875,0.2291,0.0209,0.0208,0.0208),
                  ECA = c(0.0196,0.0196,0.9019,0.0196,0.0197,0.0196),
                  SWI = c(0.0185,0.0185,0.0185,0.6852,0.0185,0.2407),
                  SEI = c(0.0417,0.0416,0.0418,0.0417,0.7083,0.1249),
                  WCP = c(0.0158,0.0159,0.0158,0.0158,0.0159,0.9207))
baysass

# Convert to matrix
baysass.mat = as.matrix(baysass[, 2:7])
baysass.mat

# Gene flow matrix
dimnames(baysass.mat) = list(source = baysass$Region, sink = baysass$Region)
baysass.mat

# Define region colours
cols <- c("#4DAF4A","#E41A1C","#FF7F00","#377EB8","#984EA3","#FFFF33")

# Reset the circular layout parameters
circos.clear()

# Parameters for the circular layout
circos.par(start.degree = 30, gap.degree = 3)

# Plot chord diagram
chordDiagram(x = baysass.mat, grid.col = cols, grid.border = "black", transparency = 0.3,
             order = baysass$Region, directional = T, direction.type = "arrows",
             self.link = 1, preAllocateTracks = list(track.height = 0.01),
             annotationTrack = "grid", annotationTrackHeight = c(0.05, 0.05),
             link.border = "NA", link.sort = T, link.decreasing = T,
             link.arr.length = 0.4, link.arr.lty = 3, link.arr.col = "#252525", 
             link.largest.ontop = F)

chordDiagram(x = baysass.mat, grid.col = cols, transparency = 0.4,
             order = baysass$Region, directional = T, direction.type = "arrows", self.link = 1,
             preAllocateTracks = list(track.height = 0.05),
             link.border = "NA", link.sort = T, link.decreasing = T,
             link.arr.length = 0.2, link.arr.lty = 2, link.arr.col = "#252525", 
             link.largest.ontop = F)

chordDiagram(x = baysass.mat, grid.col = cols, transparency = 0.4)

# Add labels to chord diagram
circos.trackPlotRegion(track.index = 1,
                       bg.border = NA,
                       panel.fun = function(x, y) {
                         xlim = get.cell.meta.data("xlim")
                         sector.index = get.cell.meta.data("sector.index")
                         # Text direction
                         theta = circlize(mean(xlim), 1)[1, 1] %% 360
                         dd = ifelse(theta < 180 || theta > 360, "bending.inside", "bending.outside")
                         circos.text(x = mean(xlim), y = 0.8, labels = sector.index, facing = dd,
                                     niceFacing = TRUE, cex = 1, font = 2)
                       }
)

# Add title
mtext("BayesAss", outer = FALSE, cex = 2, font = 2, line = -1)

# save plot
png("sailfish.bayesass.png") 
# 2. Create a plot
plot(pch = 16, frame = FALSE,
     xlab = "wt", ylab = "mpg", col = "#2E9FDF")
# Close the pdf file
dev.off()

# MIGRATE-N
# Prepare migrate data
migrate <- tibble(Region = c("WCA","SWA","ECA","SWI","SEI","WCP"),
                       WCA = c(0,47633.33,0,0,0,0),
                       SWA = c(0,0,47700.00,0,0,0),
                       ECA = c(0,0,0,13833.33,0,0),
                       SWI = c(0,0,0,0,47700.00,0),
                       SEI = c(0,0,0,0,0,47566.67),
                       WCP = c(0,0,0,0,0,0))
migrate

# Convert to matrix
migrate.mat = as.matrix(migrate[, 2:7])
migrate.mat

# Gene flow matrix
dimnames(migrate.mat) = list(source = migrate$Region, sink = migrate$Region)
migrate.mat

# Define region colours
cols <- c("#4DAF4A","#E41A1C","#FF7F00","#377EB8","#984EA3","#FFFF33")

# Reset the circular layout parameters
circos.clear()

# Parameters for the circular layout
circos.par(start.degree = 30, gap.degree = 3)

# Plot chord diagram
chordDiagram(x = migrate.mat, grid.col = cols, grid.border = "black", transparency = 0.3,
             order = migrate$Region, directional = T, direction.type = "arrows",
             self.link = 1, preAllocateTracks = list(track.height = 0.01),
             annotationTrack = "grid", annotationTrackHeight = c(0.05, 0.05),
             link.border = "NA", link.sort = T, link.decreasing = T,
             link.arr.length = 0.4, link.arr.lty = 3, link.arr.col = "#252525", 
             link.largest.ontop = F)

chordDiagram(x = migrate.mat, grid.col = cols, transparency = 0.4,
             order = migrate$Region, directional = T, direction.type = "arrows", self.link = 1,
             preAllocateTracks = list(track.height = 0.05),
             link.border = "NA", link.sort = T, link.decreasing = T,
             link.arr.length = 0.2, link.arr.lty = 2, link.arr.col = "#252525", 
             link.largest.ontop = F)

chordDiagram(x = migrate.mat, grid.col = cols, transparency = 0.4)

# Add labels to chord diagram
circos.trackPlotRegion(track.index = 1,
                       bg.border = NA,
                       panel.fun = function(x, y) {
                         xlim = get.cell.meta.data("xlim")
                         sector.index = get.cell.meta.data("sector.index")
                         # Text direction
                         theta = circlize(mean(xlim), 1)[1, 1] %% 360
                         dd = ifelse(theta < 180 || theta > 360, "bending.inside", "bending.outside")
                         circos.text(x = mean(xlim), y = 0.8, labels = sector.index, facing = dd,
                                     niceFacing = TRUE, cex = 1, font = 2)
                       }
)

# Add title
mtext("migrate", outer = FALSE, cex = 2, font = 2, line = -1)

# save plot
svg("sailfish.migrate.png") 
# 2. Create a plot
plot(pch = 16, frame = FALSE,
     xlab = "wt", ylab = "mpg", col = "#2E9FDF")
# Close the pdf file
dev.off()
