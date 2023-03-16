library(tidyverse)
library(ggthemes)
library(patchwork)
#library(pophelper)
library(vcfR)
library(adegenet)
#library(starmie)
library(RColorBrewer)

# Set working directory
setwd("~/sailfish/populations.analysis/populational/dapc")
load("./dapc.RData")
# Load files
vcf <- read.vcfR("sailfish.pruned.vcf.gz")
pop.data <- read.delim("popdata.tsv", header = TRUE, sep = "\t")
all(colnames(vcf@gt)[-1] == pop.data$sample)
gl <- vcfR2genlight(vcf, n.cores = 16)
ploidy(gl) <- 2
pop(gl) <- pop.data$pop

#cols <- brewer.pal(n = nPop(gl), name = "Set1")
#display.brewer.all()
# 1. Visualize a single RColorBrewer palette 
# by specifying its name
#display.brewer.pal(n, name)
# 2. Return the hexadecimal color code of the palette
#brewer.pal(n, name)
#display.brewer.all(colorblindFriendly = T)
###############################################################################

# DAPC cross-validation

svg("dapc_xval_step1.svg")
set.seed(123)
cross.val <- xvalDapc(tab(gl, NA.method = "mean"), gl$pop, result = "groupMean",
                      xval.plot = TRUE, n.rep = 100, parallel = "multicore", ncpus = 16) # "multicore" on linux, "snow" on windows 
invisible(dev.off())

svg("dapc_xval_step2.svg")
set.seed(123)
cross.val <- xvalDapc(tab(gl, NA.method = "mean"), gl$pop, result = "groupMean",
                      xval.plot = TRUE, n.pca = 5:20, n.rep = 1000, parallel = "multicore", ncpus = 16)
invisible(dev.off())

cross.val$`Root Mean Squared Error by Number of PCs of PCA`
cross.val$`Number of PCs Achieving Highest Mean Success`
cross.val$`Number of PCs Achieving Lowest MSE`

n.pcs <- as.numeric(cross.val$`Number of PCs Achieving Lowest MSE`)

###############################################################################

# DAPC

# Run a DAPC using site IDs as priors
dapc.eigenval <- tibble(lin_dis = c(seq(1, length(cross.val$DAPC$eig))),
                        explain_var = (cross.val$DAPC$eig/sum(cross.val$DAPC$eig))*100)

dapc.screeplot <- ggplot(dapc.eigenval, aes(x = lin_dis, y = explain_var)) +
                  geom_col(fill = heat.colors(nPop(gl)-1)) +
                  labs(#title = "DA eigenvalues",
                       x = "Linear discriminants")+
#                       y = "Explained variance (%)") +
                  theme_minimal() +
                  theme(title = element_text(size = 10),
                        axis.title = element_text(size = 9),
                        axis.text = element_text(size = 7),
                        axis.title.y = element_blank())

dapc.data <- as_tibble(cross.val$DAPC$ind.coord) %>%
             mutate(sample = row.names(cross.val$DAPC$ind.coord),
                    pop = cross.val$DAPC$grp)

dapc.centroid <- aggregate(cbind(LD1, LD2) ~ pop, data = dapc.data, FUN = mean)

dapc.data <- left_join(dapc.data, dapc.centroid, by = "pop", suffix = c("",".cen"))

dapc.plot <- dapc.data %>%
             arrange(sample) %>%
             mutate(Pop = factor(pop, levels=c("WCA","SWA","ECA","SWI","SEI","WCP"))) %>%
             ggplot(aes(x = LD1, y = LD2, color = factor(Pop))) +
             geom_point(size = 2, alpha = 0.8) +
#             geom_text(size = 3, vjust = 0, nudge_y = 0.05, label = row.names(pca.scores)) +
#             geom_segment(aes(xend = LD1.cen, yend = LD2.cen, color = factor(pop)), show.legend = FALSE) +
             stat_ellipse(lty = 1, lwd = 0.5) +
             geom_hline(yintercept = 0, lty = 2) +
             geom_vline(xintercept = 0, lty = 2) +
             scale_colour_manual(values = c("#4DAF4A","#E41A1C","#FF7F00","#377EB8","#984EA3","#A65628")) +
             labs(title = "DAPC",
                  x = paste0("LD1 (", as.character(round(dapc.eigenval$explain_var[1], 2)), "%)"),
                  y = paste0("LD2 (", as.character(round(dapc.eigenval$explain_var[2], 2)), "%)")) +
             theme_minimal() +
             theme(title = element_text(size = 10),
                   axis.title = element_text(size = 9),
                   axis.text = element_blank(),
                   legend.title = element_blank(),
                   legend.text = element_text(size = 9),
                   legend.position = "bottom",
                   legend.key.size = unit(0.6, "line")) +
                   guides(color=guide_legend(nrow=1))
dapc.plot
ggsave("dapc.png", device = "png", width = 17, height = 17, units = "cm", dpi = 600)
ggsave("dapc.svg", device = "svg", width = 21, height = 12, units = "cm", dpi = 600)
###############################################################################

# PCA
pca <- glPca(gl, nf = n.pcs, parallel = require("parallel"), n.cores = 24)

pca.eigenval <- tibble(prin_comp = c(seq(1, length(pca$eig))),
                       explain_var = (pca$eig/sum(pca$eig))*100)

pca.screeplot <- ggplot(pca.eigenval, aes(x = prin_comp, y = explain_var)) +
                 geom_col(fill = heat.colors(nInd(gl)-1)) +
#                 geom_line(size = 1) +
                 labs(#title = "Eigenvalues",
                      x = "Principal components",
                      y = "Explained variance (%)") +
                 theme_minimal() +
                 theme(title = element_text(size = 10),
                       axis.title = element_text(size = 9),
                       axis.text = element_text(size = 7))

pca.data <- as_tibble(pca$scores) %>%
            mutate(sample = row.names(pca$scores),
                   pop = pop(gl))

pca.centroid <- aggregate(cbind(PC1, PC2) ~ pop, data = pca.data, FUN = mean)

pca.data <- left_join(pca.data, pca.centroid, by = "pop", suffix = c("",".cen"))

pca.plot <- pca.data %>%
            arrange(sample) %>%
            mutate(Pop = factor(pop, levels=c("WCA","SWA","ECA","SWI","SEI","WCP"))) %>%
            ggplot(aes(x = PC1, y = PC2, color = factor(Pop))) +
            geom_point(size = 2, alpha = 0.8) +
#            geom_text(size = 3, vjust = 0, nudge_y = 0.05, label = row.names(pca.scores)) +
#            geom_segment(aes(xend = PC1.cen, yend = PC2.cen, color = factor(pop)), show.legend = FALSE) +
            stat_ellipse(lty = 1, lwd = 0.5) +
            geom_hline(yintercept = 0, lty = 2) +
            geom_vline(xintercept = 0, lty = 2) +
            scale_colour_manual(values = c("#4DAF4A","#E41A1C","#FF7F00","#377EB8","#984EA3","#A65628")) +
            labs(title = "PCA",
                 x = paste0("PC1 (", as.character(round(pca.eigenval$explain_var[1], 2)), "%)"),
                 y = paste0("PC2 (", as.character(round(pca.eigenval$explain_var[2], 2)), "%)")) +
            theme_minimal() +
            theme(title = element_text(size = 10),
                  axis.title = element_text(size = 9),
                  axis.text = element_blank(),
                  legend.title = element_blank(),
                  legend.text = element_text(size = 9),
                  legend.position = "bottom",
                  legend.key.size = unit(0.6, "line")) +
                  guides(color=guide_legend(nrow=1))
pca.plot
ggsave("pca.png", device = "png", width = 170, height = 170, units = "mm", dpi = 600)
ggsave("pca.svg", device = "svg", width = 21, height = 12, units = "cm", dpi = 600)

((pca.plot | dapc.plot) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.title = element_blank())) /
    (pca.screeplot | dapc.screeplot) +
    plot_layout(ncol = 1, nrow = 2, widths = c(6, 10), heights = c(4,2.5,2.5)) + plot_annotation(tag_levels = 'a')
# save plot
ggsave("clustering.png", device = "png", width = 35, height = 30, units = "cm", dpi = 600)
ggsave("clustering.svg", device = "svg", width = 35, height = 30, units = "cm", dpi = 600)

# plot different PCs
((pca.plot | pca.plot2) / (pca.plot3 | pca.screeplot) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.title = element_blank())) +
  plot_annotation(tag_levels = 'a')
# save plot
ggsave("pcas.png", device = "png", width = 35, height = 30, units = "cm", dpi = 600)
ggsave("pcas.svg", device = "svg", width = 35, height = 30, units = "cm", dpi = 600)


save.image(file = "dapc.RData")

###############################################################################

setwd("~/Documents/billfishes/sailfish/genome/sailfish_cl/populations.analysis/populational/ngsadmix")

# Evanno's test for better K value - (http://clumpak.tau.ac.il/index.html)
evanno <- read.delim("evanno.tsv", header = TRUE, sep = "\t")

elpdmean.plot <- ggplot(evanno, aes(x = K, y = mean)) +
                 geom_line(color = "red", alpha = 0.5) +
                 geom_point() +
                 geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 0.1) +
#                 ggtitle("Mean estimated Ln probability") +
#                 xlab("K")+
                 ylab(expression(paste("Mean L(", italic("K"), ") \u00B1 SD"))) +
                 theme_minimal() +
                 theme(axis.title.x = element_blank(),
                       axis.title.y = element_text(size = 8),
                       axis.text.x = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks.x = element_blank())

lnk1.plot <- ggplot(evanno, aes(x = K, y = Ln1)) +
             geom_line(color = "red", alpha = 0.5) +
             geom_point() +
#             geom_errorbar(aes(ymin = lnk1min, ymax = lnk1max), width = 0.1) +
             ylab(expression(paste("L'(", italic("K"), ") \u00B1 SD"))) +
             theme_minimal() +
             theme(axis.title.x = element_blank(),
                   axis.title.y = element_text(size = 8),
                   axis.text.x = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks.x = element_blank())

lnk2.plot <- ggplot(evanno, aes(x = K, y = Ln2)) +
             geom_line(color = "red", alpha = 0.5) +
             geom_point() +
#             geom_errorbar(aes(ymin = lnk2min, ymax = lnk2max), width = 0.1) +
             ylab(expression(paste("L''(", italic("K"), ") \u00B1 SD"))) +
             theme_minimal() +
             theme(legend.title=element_blank(),
                   axis.title.x = element_text(size = 8),
                   axis.title.y = element_text(size = 8),
                   axis.text.x = element_text(size = 8),
                   axis.text.y = element_blank(),
                   axis.ticks.x = element_blank(),
                   legend.position = "none")

bestK <- ggplot(data=evanno, aes(x=K, y=DeltaK, group=1)) +
         geom_line(size=0.5, linetype=1, color = "red", alpha = 0.5) +
#          geom_path()+
         geom_point(size=2)+
#          ggtitle("Evanno's Method")+
         xlab("K")+
         ylab("ΔK")+
         theme_minimal()+
         theme(legend.title=element_blank(),
               axis.title.x = element_text(size = 8),
               axis.title.y = element_text(size = 8),
               axis.text.x = element_text(size = 8),
               axis.text.y = element_blank(),
               axis.ticks.x = element_blank(),
               legend.position = "none")
########################################################################################
# admixture
# get input file names
files <- dir("./", pattern = ".qopt")

# get metadata
meta.data <- read.delim("popdata.tsv", header = TRUE, sep = "\t")

# iterate over files
for(i in 1:length(files)){
# read input Q-matrix
q.matrix <- read.table(paste0("./", files[i]))
# create data frame
assign(paste0("k", i+1, ".data"),
       q.matrix %>%
       mutate(pop = meta.data$pop,
              samples = meta.data$sample
              ) %>%
              pivot_longer(cols = starts_with("V"),
                           names_to = "ancestry",
                           values_to = "probability")
  )
}

k2.plot <- ggplot(k2.data, aes(x = samples, y = probability, fill = ancestry)) +
           geom_col(color = "white", size = 0.2, width = 1) +
#           coord_flip() +
           facet_grid(~fct_inorder(pop),
#           switch = "x",
           scales = 'free',
           space = 'free') +
           ylab("K = 2") +
           theme_minimal() +
           theme(panel.spacing.x = unit(0.001, "lines"),
                 axis.title.x = element_blank(),
                 axis.title.y = element_text(size = 8),
                 axis.text.x = element_blank(),
                 axis.text.y = element_text(size = 7, hjust = 1),
                 strip.background = element_blank(),
                 strip.text = element_text(size = 10, vjust = 0.5, hjust = 0.5),
                 panel.grid = element_blank(),
                 strip.placement = NULL) +
           scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
           scale_x_discrete(expand = expansion(add = 0.75)) +
#           scale_fill_manual(values = c("#4DAF4A", "#984EA3"),
#                             guide = "none")
           scale_fill_manual(values = c("#2171B5","#6A51A3"),
                             guide = "none")

k3.plot <- ggplot(k3.data, aes(x = samples, y = probability, fill = ancestry)) +
           geom_col(color = "white", size = 0.2, width = 1) +
#           coord_flip() +
           facet_grid(~fct_inorder(pop),
#           switch = "x",
           scales = 'free',
           space = 'free') +
           ylab("K = 3") +
           theme_minimal() +
           theme(panel.spacing.x = unit(0.01, "lines"),
                 axis.title.x = element_blank(),
                 axis.title.y = element_text(size = 8),
                 axis.text.x = element_blank(),
                 axis.text.y = element_text(size = 7, hjust = 1),
                 strip.background = element_blank(),
                 strip.text = element_blank(),
                 panel.grid = element_blank(),
                 strip.placement = NULL) +
           scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
           scale_x_discrete(expand = expansion(add = 0.75)) +
#           scale_fill_manual(values = c("#4DAF4A", "#984EA3", "#E41A1C"),
#                             guide = "none")
           scale_fill_manual(values = c("#6A51A3","#6BAED6","#2171B5"),
                             guide = "none")

k4.plot <- ggplot(k4.data, aes(x = samples, y = probability, fill = ancestry)) +
           geom_col(color = "white", size = 0.2, width = 1) +
#           coord_flip() +
           facet_grid(~fct_inorder(pop),
#           switch = "x",
           scales = 'free',
           space = 'free') +
           ylab("K = 4") +
           theme_minimal() +
           theme(panel.spacing.x = unit(0.01, "lines"),
                 axis.title.x = element_blank(),
                 axis.title.y = element_text(size = 8),
                 axis.text.x = element_blank(),
#                 axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
                 axis.text.y = element_text(size = 7, hjust = 1),
                 strip.background = element_blank(),
                 strip.text = element_blank(),
                 panel.grid = element_blank(),
                 strip.placement = NULL) +
            scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
            scale_x_discrete(expand = expansion(add = 0.75)) +
#            scale_fill_manual(values = c("#4DAF4A","#E41A1C","#FF7F00","#984EA3"),
#                              guide = "none")
            scale_fill_manual(values = c("#4292C6","#6A51A3","#2171B5","#6BAED6"),
                              guide = "none")

k5.plot <- ggplot(k5.data, aes(x = samples, y = probability, fill = ancestry)) +
           geom_col(color = "white", size = 0.2, width = 1) +
#           coord_flip() +
           facet_grid(~fct_inorder(pop),
#           switch = "x",
           scales = 'free',
           space = 'free') +
           ylab("K = 5") +
           theme_minimal() +
           theme(panel.spacing.x = unit(0.01, "lines"),
                 axis.title.x = element_blank(),
                 axis.title.y = element_text(size = 8),
                 axis.text.x = element_blank(),
#                 axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
                 axis.text.y = element_text(size = 7, hjust = 1),
                 strip.background = element_blank(),
                 strip.text = element_blank(),
                 panel.grid = element_blank(),
                 strip.placement = NULL) +
           scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
           scale_x_discrete(expand = expansion(add = 0.75)) +
#           scale_fill_manual(values = c("#0022ff","#984EA3","#4DAF4A","#E41A1C","#FF7F00"),
#                             guide = "none")
           scale_fill_manual(values = c("#6A51A3","#4292C6","#2171B5","#6BAED6","#9ECAE1"),
                             guide = "none")

k6.plot <- ggplot(k6.data, aes(x = samples, y = probability, fill = ancestry)) +
           geom_col(color = "white", size = 0.2, width = 1) +
#           coord_flip() +
           facet_grid(~fct_inorder(pop),
#           switch = "x",
           scales = 'free',
           space = 'free') +
           ylab("K = 6") +
           theme_minimal() +
           theme(panel.spacing.x = unit(0.01, "lines"),
                 axis.title.x = element_blank(),
                 axis.title.y = element_text(size = 8),
#                 axis.text.x = element_blank(),
                 axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
#                 axis.text.y = element_blank(),
                 axis.text.y = element_text(size = 7, hjust = 1),
                 strip.background = element_blank(),
                 strip.text = element_blank(),
                 panel.grid = element_blank(),
                 strip.placement = NULL) +
           scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
           scale_x_discrete(expand = expansion(add = 0.75)) +
#           scale_fill_manual(values = c("#0022ff","#4DAF4A","#E41A1C","#984EA3","#FF7F00","#f6ff00"),
#                             guide = "none")
           scale_fill_manual(values = c("#C6DBEF","#6A51A3","#4292C6","#2171B5","#6BAED6","#9ECAE1"),
                             guide = "none")

(k2.plot / k3.plot / k4.plot / k5.plot / k6.plot) & theme(strip.placement = NULL)
# save plot
ggsave("sailfish.admixture.png", device = "png", width = 30, height = 25, units = "cm", dpi = 600)
ggsave("sailfish.admixture.svg", device = "svg", width = 30, height = 25, units = "cm", dpi = 600)

(((pca.plot | dapc.plot) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.title = element_blank())) /
  (pca.screeplot | dapc.screeplot) / (elpdmean.plot | lnk1.plot) / (lnk2.plot | bestK) +
   plot_layout(ncol = 1, nrow = 4, widths = c(6, 10), heights = c(4,2.5,2.5,2.5)) |
  (k2.plot / k3.plot / k4.plot / k5.plot / k6.plot)) + plot_annotation(tag_levels = 'a')

# save plot
ggsave("sailfish.clustering.png", device = "png", width = 45, height = 35, units = "cm", dpi = 600)
ggsave("sailfish.clustering.svg", device = "svg", width = 45, height = 35, units = "cm", dpi = 600)

# plot all the clustering and admixtures
  (((pca.plot | dapc.plot) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.title = element_blank())) /
    (pca.screeplot | dapc.screeplot) / (elpdmean.plot | lnk1.plot) / (lnk2.plot | bestK) + plot_layout(nrow = 4, widths = c(15), heights = c(10,5,5,5)) |
    (k2.plot / k3.plot / k4.plot / k5.plot / k6.plot) |
    (tp1 / tp2 / tp3 / tp4 / tp5) + plot_layout(heights = c(25))) + plot_annotation(tag_levels = 'a')

# save plot
ggsave("sailfish.clustering.all.png", device = "png", width = 50, height = 35, units = "cm", dpi = 600)
ggsave("sailfish.clustering.all.svg", device = "svg", width = 50, height = 35, units = "cm", dpi = 600)
ggsave("sailfish.clustering.all.pdf", device = "pdf", width = 50, height = 35, units = "cm", dpi = 600)
