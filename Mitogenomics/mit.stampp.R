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
setwd("~/Documents/sailfish/sailfish_cl/stampp")

#loading vcf
vcf <- read.vcfR("sailfish.filtered.vcf.gz")

#populations data and geographic coordinates data
pop.data <- read.delim("popdata.tsv", sep = "\t", header = TRUE)
xy <- read.delim("xy.tsv", sep = "\t", header = FALSE)

#check if populations were updated
all(colnames(vcf@gt)[-1] == pop.data$sample)

#Converting the dataset to a genlight object
gl <- vcfR2genlight(vcf)
gl@other$xy <- xy
pop(gl) <- pop.data$pop

# Fst paiwise
fst <- stamppFst(gl, nboots = 1000, percent = 95, nclusters = 4)
fst

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
            ggtitle(expression(atop(Pairwise~F[ST]~WC~(1984)~genomic, ""))) +
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
ggsave("sailfish.gFst.png", device = "png", width = 20, height = 10, units = "cm", dpi = 600)
ggsave("sailfish.gFst.svg", device = "svg", width = 20, height = 10, units = "cm", dpi = 600)

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
              ggtitle(expression(atop(Pairwise~F[ST]~WC~(1984)~mitogenomic, ""))) +
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
