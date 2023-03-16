library(ggplot2) #plotting
library(reshape2) #plotting
library(adegenet) #FST statistics and data storage
library(vcfR) #https://github.com/knausb/vcfR
library(data.table)
library(patchwork)
library(tibble)
library(StAMPP)
library(RColorBrewer)

# set working directory
setwd("/home/bferrette/sailfish_cl/population.analysis/sailfish/stampp")

# load the vcf file
vcf <- read.vcfR("sailfish.pruned.vcf.gz")

# population and geographic coordinates data
pop.data <- read.delim("popdata.tsv", sep = "\t", header = TRUE)
xy <- read.delim("xy.tsv", sep = "\t", header = FALSE)

# check if populations were updated
all(colnames(vcf@gt)[-1] == pop.data$sample)

# converting the vcf to genlight
gl <- vcfR2genlight(vcf, n.cores = 24)
gl@other$xy <- xy
pop(gl) <- pop.data$pop

# Fst paiwise calculations
fst <- stamppFst(gl, nboots = 1000, percent = 95, nclusters = 24)
fst
# plot Fst pairwise
Fst <- data.matrix(fst$Fsts)
Pvalue <- data.matrix(fst$Pvalues)
Fst.melted <- reshape2::melt(Fst, na.rm = TRUE)
Pvalue.melted <- reshape2::melt(Pvalue, na.rm = TRUE)
melted.data <- cbind(Fst.melted, Pvalue.melted)
colnames(melted.data)<- c("Var1","Var2","Fst","Var3","Var4","Pvalue")
melted.data$stars <- cut(melted.data$Pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels

Fst.plot <- ggplot(data = melted.data, aes(Var2, Var1, fill = Fst)) +
#            ggplot(data = melted.data, aes(Var4, Var3, fill = pvalue)) +
            geom_tile(color = "white") +
#            geom_text(aes(label = round(Fst, 3)), size = 3) +
#            geom_text(aes(label = paste(round(Pvalue,1), c(" ","*")[(abs(Pvalue) <= .05)+1]))) +
            geom_text(aes(label=paste0(round(Fst, 3),"\n", stars)), color="black", size=3) +
#               scale_fill_gradient(low = "white", high = "red", name="FST") +
            scale_fill_distiller(palette = "RdYlBu") + 
            ggtitle(expression(atop(Pairwise~F[ST]~WC~(1984), ""))) +
#               labs( x = "Sampling Site", y = "Sampling Site") +
            theme_bw() +
            theme(axis.text.x = element_text(size = 10),
                  axis.text.y = element_text(size = 10),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  legend.title = element_blank()) +
#                     guides(fill = guide_colourbar(barwidth = 2, barheight = 20)) +
            coord_fixed()
Fst.plot
# save plot
ggsave("sailfish.gFst.png", device = "png", width = 20, height = 10, units = "cm", dpi = 600)
ggsave("sailfish.gFst.svg", device = "svg", width = 20, height = 10, units = "cm", dpi = 600)

# save R environment
save.image(file = "stampp.RData")

# mitogenomic FST pairwise
mtFst.data <- read.delim("mtFst.pairwise.tsv", sep = "\t", header = TRUE)
mt.pvalue.data <- read.delim("mtFst.pvalues.tsv", sep = "\t", header = TRUE)

mtFst <- data.matrix(mtFst.data)
mtPvalue <- data.matrix(mt.pvalue.data)
mtFst.melted <- reshape2::melt(mtFst, na.rm = TRUE)
mtPvalue.melted <- reshape2::melt(mtPvalue, na.rm = TRUE)
mtmelted.data <- cbind(mtFst.melted, mtPvalue.melted)
colnames(mtmelted.data)<- c("Var1","Var2","Fst","Var3","Var4","Pvalue")
mtmelted.data$stars <- cut(mtmelted.data$Pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels

mtFst.plot <- ggplot(data = mtmelted.data, aes(Var2, Var1, fill = Fst)) +
#              ggplot(data = melted.data, aes(Var4, Var3, fill = pvalue)) +
              geom_tile(color = "white") +
#              geom_text(aes(label = round(Fst, 3)), size = 3) +
#              geom_text(aes(label = paste(round(Pvalue,1), c(" ","*")[(abs(Pvalue) <= .05)+1]))) +
              geom_text(aes(label=paste0(round(Fst, 3),"\n", stars)), color="black", size=3) +
#               scale_fill_gradient(low = "white", high = "red", name="FST") +
              scale_fill_distiller(palette = "RdYlBu") + 
              ggtitle(expression(atop(Pairwise~Φ[ST]~WC~(1984)~mitogenomic, ""))) +
#               labs( x = "Sampling Site", y = "Sampling Site") +
              theme_bw() +
              theme(axis.text.x = element_text(size = 10),
                    axis.text.y = element_text(size = 10),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    legend.title = element_blank()) +
#                    guides(fill = guide_colourbar(barwidth = 2, barheight = 20)) +
              coord_fixed()
mtFst.plot
ggsave("sailfish.mtFst.png", device = "png", width = 30, height = 20, units = "cm", dpi = 600)
ggsave("sailfish.mtFst.svg", device = "svg", width = 30, height = 20, units = "cm", dpi = 600)

# merge plots
 (Fst.plot | mtFst.plot) + plot_annotation(tag_levels = 'a')
ggsave("sailfish.Fst.pairwise.png", device = "png", width = 30, height = 10, units = "cm", dpi = 600)
ggsave("sailfish.Fst.pairwise.svg", device = "svg", width = 30, height = 10, units = "cm", dpi = 600)

# Isolation-by-distance (IBD)
geo <- read.delim("popdist.tsv", sep = "\t", header = FALSE)

Dgen <- as.dist(fst$Fsts)
Dgen
Dgeo <- as.dist(geo)
Dgeo

# Isolation-by-distance {adegenet}
ibd.pop <- mantel.randtest(Dgen, Dgeo, nrepet = 1000000)
ibd.pop

# plot
ff <- as.vector(Dgen)
gg <- as.vector(Dgeo)
mat <- data.frame(ff,gg)
fit <- tibble(r.squared = 0.7310756)
plot_rsquared <- paste0(
  "R^2 == ",
  fit$r.squared %>% round(3))
ibd.plot <- ggplot(mat, aes(y = ff, x = gg)) +
            annotate("text", x = 15000, y = 0,
                     label = plot_rsquared, angle = 0, vjust = 0,
                     size = 5, colour="black", parse=TRUE) +
            geom_point(size = 2, colour = "black", shape = 19) +
            geom_smooth(method = "lm", colour = "black", se = TRUE, fill = "blue", alpha = 0.3) + 
#                  geom_density_2d_filled(alpha = 0.75, contour_var = "count") +
#                  stat_density_2d(geom = "raster", alpha = 0.8,
#                                  aes(fill = after_stat(density)),
#                                  contour = FALSE) +
#                  scale_fill_viridis_c() +
#                  geom_cor(method = "pearson", xpos = NULL, ypos = NULL, inherit.aes = TRUE) +
            labs(title = "Isolation-By-Distance", x = "Geographic Distance (km)", y = expression(atop(Pairwise~F[ST]~WC~(1984), ""))) + 
            theme_bw()
ibd.plot

# save plot in '.svg' format
ggsave("sailfish.ibd.pop.png", device = "png", width = 15, height = 10, units = "cm", dpi = 600)
ggsave("sailfish.ibd.pop.svg", device = "svg", width = 15, height = 10, units = "cm", dpi = 600)

