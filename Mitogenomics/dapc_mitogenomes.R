# load libraries
library(adegenet)
library(ape)
library(dartR)
library(tidyverse)
library(patchwork)

# set working directory
setwd('~/Documents/billfishes/mitogenomes/sailfish/analysis/dapc')
# load .RData
load("~/Documents/billfishes/mitogenomes/sailfish/analysis/dapc/mit.dapc.RData")

# read.fasta: read FASTA formatted files
fasta <- read.FASTA("mitogenomes.fasta")

# SequencesToGenind: Importing data from an alignement of sequences to a genind object
gen <- DNAbin2genind(fasta, pop=NULL, exp.char=c("a","t","g","c"), polyThres=1/100)
gen
# population data
pop.data <- read.delim("popdata.tsv", header = TRUE, sep = "\t")
pop(gen) <- pop.data$pop

# DAPC cross-validation
# first round
svg("dapc_xval_step1.svg")
set.seed(123)
cross.val <- xvalDapc(tab(gen, NA.method = "mean"), gen$pop, result = "groupMean",
                      xval.plot = TRUE, n.rep = 100, parallel = "multicore", ncpus = 8) # "multicore" on linux, "snow" on windows 
invisible(dev.off())
# second round
svg("dapc_xval_step2.svg")
set.seed(123)
cross.val <- xvalDapc(tab(gen, NA.method = "mean"), gen$pop, result = "groupMean",
                      xval.plot = TRUE, n.pca = 5:20, n.rep = 1000, parallel = "multicore", ncpus = 8)
invisible(dev.off())
# cross validation
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
                  geom_col(fill = heat.colors(nPop(gen)-1)) +
                  labs(#title = "DA eigenvalues",
                       x = "Linear discriminants",
                       y = "Explained variance (%)") +
                  theme_minimal() +
                  theme(title = element_text(size = 10),
                        axis.title = element_text(size = 9),
                        axis.text = element_text(size = 7))
#                        axis.title.y = element_blank())
dapc.screeplot
ggsave("mit.dapc.screeplot.svg", device = "svg", width = 21, height = 12, units = "cm", dpi = 600)

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

# merge plots
p1 <- dapc.plot + inset_element(dapc.screeplot, left = 0.8, bottom = 0.5, right = 0.98, top = 0.98)
p1
# save plot
ggsave("mit.dapc.png", device = "png", width = 21, height = 12, units = "cm", dpi = 600)
ggsave("mit.dapc.svg", device = "svg", width = 21, height = 12, units = "cm", dpi = 600)

# PCA
# Converts a genind object into a genlight object
gl <- gi2gl(gen, parallel = TRUE, verbose = 5)
gl

pca <- glPca(gl, nf = n.pcs, parallel = require("parallel"), n.cores = 8)

pca.eigenval <- tibble(prin_comp = c(seq(1, length(pca$eig))),
                       explain_var = (pca$eig/sum(pca$eig))*100)

pca.screeplot <- ggplot(pca.eigenval, aes(x = prin_comp, y = explain_var)) +
                 geom_col(fill = heat.colors(59)) +
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
ggsave("mit.pca.png", device = "png", width = 170, height = 170, units = "mm", dpi = 600)
ggsave("mit.pca.svg", device = "svg", width = 21, height = 12, units = "cm", dpi = 600)

p2 <- pca.plot + inset_element(pca.screeplot, left = 0.8, bottom = 0.5, right = 0.98, top = 0.98)
p2

# merge plots
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
# save the analyses
save.image(file = "mit.dapc.RData")
