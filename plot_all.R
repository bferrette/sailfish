library(tidyverse)
library(ggthemes)
library(patchwork)
library(pophelper)
library(vcfR)
library(adegenet)
library(RColorBrewer)

# Set working directory
setwd("/home/bferrette/sailfish_cl/dapc")
# Load files
vcf <- read.vcfR("sailfish.filtered.vcf.gz")
pop.data <- read.delim("popdata.tsv", header = TRUE, sep = "\t")
all(colnames(vcf@gt)[-1] == pop.data$sample)
gl <- vcfR2genlight(vcf)
ploidy(gl) <- 2
pop(gl) <- pop.data$pop

cols <- brewer.pal(n = nPop(gl), name = "Set1")

###############################################################################

# DAPC cross-validation

pdf("dapc_xval_step1.pdf")
set.seed(123)
cross.val <- xvalDapc(tab(gl, NA.method = "mean"), gl$pop, result = "groupMean", xval.plot = TRUE, n.rep = 100, parallel = "multicore", ncpus = 12)
invisible(dev.off())

pdf("dapc_xval_step2.pdf")
set.seed(123)
cross.val <- xvalDapc(tab(gl, NA.method = "mean"), gl$pop, result = "groupMean", xval.plot = TRUE, n.pca = 5:20, n.rep = 1000, parallel = "multicore", ncpus = 12)
invisible(dev.off())

cross.val$`Root Mean Squared Error by Number of PCs of PCA`
cross.val$`Number of PCs Achieving Highest Mean Success`
cross.val$`Number of PCs Achieving Lowest MSE`

n.pcs <- as.numeric(cross.val$`Number of PCs Achieving Lowest MSE`)

###############################################################################

# DAPC

dapc.eigenval <- tibble(lin_dis = c(seq(1, length(cross.val$DAPC$eig))),
                        explain_var = (cross.val$DAPC$eig/sum(cross.val$DAPC$eig))*100)

dapc.screeplot <- ggplot(dapc.eigenval, aes(x = lin_dis, y = explain_var)) +
  geom_col(fill = heat.colors(nPop(gl)-1)) +
  labs(title = "DA eigenvalues",
       x = "Linear discriminants",
       y = "Explained variance (%)") +
  theme_minimal() +
  theme(title = element_text(size = 10),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 7))

dapc.data <- as_tibble(cross.val$DAPC$ind.coord) %>%
  mutate(sample = row.names(cross.val$DAPC$ind.coord),
         pop = cross.val$DAPC$grp)

dapc.centroid <- aggregate(cbind(LD1, LD2) ~ pop, data = dapc.data, FUN = mean)

dapc.data <- left_join(dapc.data, dapc.centroid, by = "pop", suffix = c("",".cen"))

dapc.plot <- ggplot(dapc.data, aes(x = LD1, y = LD2, color = factor(pop))) +
             geom_point(size = 2, alpha = 0.5) +
#             geom_text(size = 3, vjust = 0, nudge_y = 0.05, label = row.names(pca.scores)) +
             geom_segment(aes(xend = LD1.cen, yend = LD2.cen, color = factor(pop)), show.legend = FALSE) +
             stat_ellipse(lty = 1, lwd = 0.5) +
             geom_hline(yintercept = 0, lty = 2) +
             geom_vline(xintercept = 0, lty = 2) +
             scale_colour_manual(values = cols) +
             labs(title = "DAPC",
                  x = paste0("LD1 (", as.character(round(dapc.eigenval$explain_var[1], 2)), "%)"),
                  y = paste0("LD2 (", as.character(round(dapc.eigenval$explain_var[2], 2)), "%)")) +
             theme_minimal() +
             theme(title = element_text(size = 10),
                   axis.title = element_text(size = 9),
                   axis.text = element_blank(),
                   legend.title = element_blank(),
                   legend.text = element_text(size = 9),
                   legend.key.size = unit(0.6, "line"))

ggsave("dapc.png", device = "png", width = 17, height = 17, units = "cm", dpi = 600)
ggsave("dapc.svg", device = "svg", width = 21, height = 12, units = "cm", dpi = 600)
###############################################################################

# PCA

pca <- glPca(gl, nf = n.pcs)

pca.eigenval <- tibble(prin_comp = c(seq(1, length(pca$eig))),
                       explain_var = (pca$eig/sum(pca$eig))*100)

pca.screeplot <- ggplot(pca.eigenval, aes(x = prin_comp, y = explain_var)) +
                 geom_col(fill = heat.colors(nInd(gl)-1)) +
                 geom_line(size = 1) +
                 labs(title = "Eigenvalues",
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

pca.plot <- ggplot(pca.data, aes(x = PC1, y = PC2, color = factor(pop))) +
            geom_point(size = 2, alpha = 0.5) +
#            geom_text(size = 3, vjust = 0, nudge_y = 0.05, label = row.names(pca.scores)) +
            geom_segment(aes(xend = PC1.cen, yend = PC2.cen, color = factor(pop)), show.legend = FALSE) +
            stat_ellipse(lty = 1, lwd = 0.5) +
            geom_hline(yintercept = 0, lty = 2) +
            geom_vline(xintercept = 0, lty = 2) +
            scale_colour_manual(values = cols) +
            labs(title = "PCA",
                 x = paste0("PC1 (", as.character(round(pca.eigenval$explain_var[1], 2)), "%)"),
                 y = paste0("PC2 (", as.character(round(pca.eigenval$explain_var[2], 2)), "%)")) +
            theme_minimal() +
            theme(title = element_text(size = 10),
                  axis.title = element_text(size = 9),
                  axis.text = element_blank(),
                  legend.title = element_blank(),
                  legend.text = element_text(size = 9),
                  legend.key.size = unit(0.6, "line"))

ggsave("pca.png", device = "png", width = 170, height = 170, units = "mm", dpi = 600)
ggsave("pca.svg", device = "svg", width = 21, height = 12, units = "cm", dpi = 600)

###############################################################################

# Evanno
#trace(evannoMethodStructure, edit=TRUE)
slist <- readQ(files=choose.files(multi=TRUE))
tr1 <- tabulateQ(slist)
sr1 <- summariseQ(tr1)
evanno <- evannoMethodStructure(
  data = sr1,
  writetable = TRUE,
  exportplot = TRUE,
  returnplot = FALSE,
  returndata = TRUE,
  pointsize = NA,
  pointtype = 20,
  pointcol = "black",
  linesize = NA,
  linecol = "black",
  ebwidth = 0.2,
  ebcol = "grey30",
  textcol = "grey30",
  xaxisbreaks = waiver(),
  xaxislabels = waiver(),
  basesize = 6,
  gridsize = 0.18,
  imgtype = "tiff",
  height = NA,
  width = NA,
  dpi = 300,
  units = "cm",
  theme = "theme_gray",
  font = "",
  na.rm = TRUE,
  outputfilename = "evannoMethodStructure",
  exportpath = getwd()
)  

tbl <- read_tsv("evannoMethodStructure.txt")

elpdmean.plot <- ggplot(tbl, aes(x = k, y = elpdmean/1000)) +
                 geom_line(color = "red", alpha = 0.5) +
                 geom_point() +
                 geom_errorbar(aes(ymin = elpdmin/1000, ymax = elpdmax/1000), width = 0.1) +
                 ylab(expression(paste("Mean L(", italic("K"), ") \u00B1 SD (", 10^3, ")"))) +
                 theme_minimal() +
                 theme(axis.title.x = element_blank(),
                 axis.title.y = element_text(size = 9),
                 axis.text.x = element_blank(),
                 axis.text.y = element_text(size = 7),
                 axis.ticks.x = element_blank())

lnk1.plot <- ggplot(tbl, aes(x = k, y = lnk1)) +
             geom_line(color = "red", alpha = 0.5) +
             geom_point() +
             geom_errorbar(aes(ymin = lnk1min, ymax = lnk1max), width = 0.1) +
             ylab(expression(paste("L'(", italic("K"), ") \u00B1 SD"))) +
             theme_minimal() +
             theme(axis.title.x = element_blank(),
                   axis.title.y = element_text(size = 9),
                   axis.text.x = element_blank(),
                   axis.text.y = element_text(size = 7),
                   axis.ticks.x = element_blank())

lnk2.plot <- ggplot(tbl, aes(x = k, y = lnk2)) +
             geom_line(color = "red", alpha = 0.5) +
             geom_point() +
             geom_errorbar(aes(ymin = lnk2min, ymax = lnk2max), width = 0.1) +
             ylab(expression(paste("L''(", italic("K"), ") \u00B1 SD"))) +
             theme_minimal() +
             theme(axis.title.x = element_blank(),
                   axis.title.y = element_text(size = 9),
                   axis.text.x = element_blank(),
                   axis.text.y = element_text(size = 7),
                   axis.ticks.x = element_blank())

deltaK.plot <- ggplot(tbl, aes(x = k, y = deltaK)) +
               geom_line(color = "red", alpha = 0.5) +
               geom_point() +
               labs(x = expression(italic("K")),
               y = expression(paste(Delta, italic("K")))) +
               scale_x_continuous(breaks = c(seq(1, 6))) +
               theme_minimal() +
               theme(axis.title = element_text(size = 9),
                     axis.text = element_text(size = 7))

###############################################################################

# Structure

files <- dir("./", pattern = ".tsv")

for(i in 1:length(files)){
  fin <- read_tsv(files[i], col_names = FALSE)
  assign(paste0("k", i+1, ".data"),
         data <- fin %>%
           rename(sample = X1) %>%
           left_join(select(pop.data, 1:2)) %>%
           pivot_longer(cols = starts_with("X"),
                        names_to = "ancestry",
                        values_to = "probability") %>%
           group_by(sample) %>%
           as_tibble()
  )
}

k2.plot <- ggplot(k2.data, aes(x = sample, y = probability, fill = ancestry)) +
           geom_col(color = "grey", size = 0.025) +
           coord_flip() +
           facet_grid(fct_inorder(as.character(pop)) ~ .,
#             switch = "x",
             scales = 'free',
             space = 'free') +
           labs(title = "K = 2") +
           theme_minimal() +
           theme(panel.spacing.x = unit(0.01, "lines"),
           title = element_text(size = 10),
           axis.title = element_blank(),
           axis.text.x = element_blank(),
           axis.text.y = element_blank(),
           strip.background = element_blank(),
           strip.text = element_blank(),
           panel.grid = element_blank()) +
           scale_y_continuous(expand = c(0, 0)) +
           scale_x_discrete(expand = expansion(add = 1)) +
           scale_fill_manual(values = c("#E41A1C","#377EB8"),
                             guide = FALSE)

k3.plot <- ggplot(k3.data, aes(x = sample, y = probability, fill = ancestry)) +
  geom_col(color = "grey", size = 0.025) +
  coord_flip() +
  facet_grid(fct_inorder(as.character(pop)) ~ .,
             #switch = "x",
             scales = 'free',
             space = 'free') +
  labs(title = "K = 3") +
  theme_minimal() +
  theme(panel.spacing.x = unit(0.01, "lines"),
        title = element_text(size = 10),
        axis.title = element_blank(),
        axis.text = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  scale_fill_manual(values = c("#4DAF4A","#E41A1C","#377EB8"),
                    guide = FALSE)

k4.plot <- ggplot(k4.data, aes(x = sample, y = probability, fill = ancestry)) +
  geom_col(color = "grey", size = 0.025) +
  coord_flip() +
  facet_grid(fct_inorder(as.character(pop)) ~ .,
             #switch = "x",
             scales = 'free',
             space = 'free') +
  labs(title = "K = 4") +
  theme_minimal() +
  theme(panel.spacing.x = unit(0.01, "lines"),
        title = element_text(size = 10),
        axis.title = element_blank(),
        axis.text = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, hjust = 0.5),
        panel.grid = element_blank()) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  scale_fill_manual(values = c("#E41A1C","#984EA3","#377EB8","#4DAF4A"),
                    guide = FALSE)

###############################################################################

(elpdmean.plot + lnk1.plot) / (lnk2.plot + deltaK.plot)
ggsave("evanno.pdf", device = "pdf", width = 170, height = 170, units = "mm", dpi = 600)
ggsave("evanno.png", device = "png", width = 170, height = 170, units = "mm", dpi = 600)
ggsave("evanno.svg", device = "svg", width = 21, height = 12, units = "cm", dpi = 600)


( ( ( ( ( pca.plot + dapc.plot ) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.title = element_blank()) ) / 
        ( pca.screeplot + dapc.screeplot ) / ( elpdmean.plot / deltaK.plot ) ) ) + plot_layout(nrow = 3, heights = c(1.5, 0.75, 1.5)) | 
    ( k2.plot + k3.plot + k4.plot ) )
ggsave("spz.atlantic.pdf", device = "pdf", width = 210, height = 210, units = "mm", dpi = 600)
ggsave("spz.atlantic.png", device = "png", width = 170, height = 170, units = "mm", dpi = 600)
ggsave("spz.atlantic.svg", device = "svg", width = 35, height = 20, units = "cm", dpi = 600)
