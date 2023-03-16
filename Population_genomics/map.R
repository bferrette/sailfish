# load libraries
library(sf) #https://github.com/r-spatial/sf
library(rnaturalearth) #https://github.com/ropensci/rnaturalearth
library(rnaturalearthdata)
library(tidyverse) #https://www.tidyverse.org/
library(patchwork) #https://patchwork.data-imaginist.com/index.html
library(RColorBrewer)

# set working directory
setwd("~/billfishes/sailfish/genome/sailfish_cl/map")

# geographical coordinates
SFA <- read.table("popdata.tsv", header = TRUE, sep = "\t")
SFA2 <- SFA %>% arrange(sample) %>%
                mutate(Pop = factor(pop, levels=c("WCA","SWA","ECA","SWI","SEI","WCP","BUM","SPF","WHM")))
        
sailfish <- st_read("Istiophorus platypterus.shp")
# load data
world <- ne_countries(scale = "large", returnclass = "sf")
# plot world map
#ggplot(data = world) + geom_sf()
#display.brewer.all()
# 1. Visualize a single RColorBrewer palette 
# by specifying its name
#display.brewer.pal(n, name)
# 2. Return the hexadecimal color code of the palette
#brewer.pal(n, name)
# plot shapefile
#p <- ggplot() + 
#        geom_sf(data = sailfish, size = 1, color = "black", fill = "green") + 
#        coord_sf()
#p

# Plot the results
map <- ggplot() +
       geom_sf(data = sailfish, color = "black", fill = "forestgreen", lwd = 0, alpha = 0.2) +
       geom_sf(data = world, fill = "lightgray", lwd = 0.2) +
       coord_sf(xlim=c(-80, 147), ylim = c(-32, 26), expand = FALSE) +
       geom_point(data = SFA2, aes(x = longitude, y = latitude, colour=Pop), size = 1.5) +
       scale_color_manual(values=c("#4DAF4A","#E41A1C","#FF7F00","#377EB8","#984EA3","#FFFF33","#0000FF","#8DD3C7","#ea07f2"),
                          labels = c("WCA (n = 10)","SWA (n = 10)","ECA (n = 11)","SWI (n = 12)","SEI (n = 2)","WCP (n = 15)",
                                     "blue marlin","longbill spearfish", "white marlin")) +
       theme_bw() + 
       theme(legend.position="bottom",
             legend.title = element_blank(),
             axis.text.x = element_text(color="black", size=8),
             axis.title.x = element_blank(),
             axis.text.y = element_text(color="black", size=8),
             axis.title.y = element_blank()) +
             guides(color = guide_legend(nrow = 2))
map

# save plot
ggsave("sailfish.sampling.svg", device = "svg", width = 20, height = 10, units = "cm", dpi = 600)
ggsave("sailfish.sampling.png", device = "png", width = 20, height = 10, units = "cm", dpi = 600)
