# --------------------------- #
#
# Seascape Redundancy Analysis
# 
# Description:
# Perform a redundancy analysis (RDA) - https://pmassicotte.github.io/stats-denmark-2019/07_rda.html#/ , https://mb3is.megx.net/gustame/constrained-analyses/rda
# https://pmassicotte.github.io/stats-denmark-2019/07_rda.html#/
# --------------------------- #

# Load packages
library(adespatial)
library(GGally)
library(ggrepel)
library(gtools)
library(patchwork)
library(RColorBrewer)
library(tidyverse)
library(vegan)

# set working directora
setwd("/home/bferrette/sailfish_cl/population.analysis/sailfish/seascape/rda")

# Import genetic data
allele_freqs <- read.delim("allele_freqs.tsv", header = TRUE, sep = "\t")
allele_freqs
#outliers <- decostand(allele_freqs, method = "hellinger")
#outliers
#allele_freqs <- read.delim("allele_freqs.gwds.tsv", header = TRUE, sep = "\t")
# Import environmental data
env <- read_delim("environmental_rawdata.tsv", delim = "\t")
# Import spatial data
dbmem <- read_delim("dbmems.tsv", delim = "\t")
# Set seed
set.seed(123)

#--------------#
#
# Identify significant variables
#
#--------------#
#
# Global model with all environmental variables
my_rda <- rda(allele_freqs ~ ., data = env, scale=TRUE) # Predicting allele frequency using all variables contained in env
anova.cca(my_rda, permutations = how(nperm=9999), by = NULL, model = "full",  parallel = getOption("mc.cores")) # Global RDA significance
# Since p = 0.001 we reject H0 meaning that the RDA model is significant.
plot(my_rda) # Show the triplot

# Multicollinearity before variables selection using the square root of Variable Inflation Factors (VIF)
sqrt(vif.cca(my_rda)) 
# subset only significant independent environmental variables to include in the RDA
env.vif <- subset(env, select = -c(chlo,ptpk,diff.atn,sdO2))
str(env.vif)
my_rda <- rda(allele_freqs ~ ., data = env.vif, scale = TRUE)
anova.cca(my_rda, permutations = how(nperm=9999), by = NULL, model = "full",  parallel = getOption("mc.cores")) # Global RDA significance
sqrt(vif.cca(my_rda))  # Multicollinearity before variables selection
plot(my_rda)
# global R2
# Attention: R2 as the relative contribution of each eigenvectors are unadjusted and are therefore biased.
# For a proper computation of unbiased, adjusted R2 one should use the RsquareAdj() function.
adj.r2 <- RsquareAdj(my_rda, permutations = 9999)$adj.r.squared # Total variance explained by the RDA
adj.r2
# Use forward selection to identify significant environmental variables - https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
env.fs <- forward.sel(Y = allele_freqs,
                      X = env.vif,
                      adjR2thresh = adj.r2,
                      nperm = 9999,
                      R2more = 0.05,
                      alpha = 0.05)
env.fs
# P-values from the forward selection above are in the object sel.fs$pval, and the adjustement can be done in the following way (Holm's correction):
n.tests <- ncol(env.fs)  # number of tests equals to number of all variables from which is being selected
pval.adj <- p.adjust(env.fs$pval, method = 'holm', n = n.tests)
pval.adj
env.fs$pval.adj <- pval.adj
env.fs
# subset only significant independent environmental variables to include in the RDA
env.sig <- subset(env, select = env.fs$variables)
str(env.sig)

# create another RDA using the selected variables
my_rda2 <- rda(allele_freqs ~ ssal + sst + sprp, data = env.sig, scale = TRUE)
anova.cca(my_rda, permutations = how(nperm=9999), by = NULL, model = "full",  parallel = getOption("mc.cores")) # Global RDA significance
sqrt(vif.cca(my_rda2))  # Multicollinearity before variables selection
# global R2
adj.r2.2 <- RsquareAdj(my_rda2, permutations = 9999)$adj.r.squared
adj.r2.2

# Use forward selection to identify significant dbmems
dbmem.for <- forward.sel(Y = allele_freqs, X = dbmem,
                         adjR2thresh = 0.95,
                         nperm = 9999,
                         R2more = 0.05,
                         alpha = 0.05)
dbmem.for
# Subset only significant independent variables to include in the RDA
dbmem.sig <- subset(dbmem, select = dbmem.for$variables)
str(dbmem.sig)
# Combine environmental variables and dbmems
data.sig <- cbind(env.vif, dbmem.sig)
str(data.sig)

#--------------#
#
# Redundancy analysis
#
#--------------#

# Perform RDA
rda <- rda(allele_freqs ~ ., data = data.sig, scale = TRUE)
rda
anova.cca(rda, by = NULL, permutations = how(nperm=9999), model = "full",  parallel = getOption("mc.cores")) # Global RDA significance
anova.cca(rda, by='axis', permutations = how(nperm=9999), model = "full",  parallel = getOption("mc.cores")) # Test which axis are significant
anova.cca(rda, by='terms', permutations = how(nperm=9999), model = "full",  parallel = getOption("mc.cores")) # Test which terms are significant
coef(rda)
sqrt(vif.cca(rda)) # variance inflation factor (<2 OK)
adj.R2 <- round(RsquareAdj(rda, permutations = 9999)$adj.r.squared, 3)
adj.R2

# Variance explained by each canonical axis
summary(eigenvals(rda, model = "constrained"))
#screeplot(rda)
eig <- eigenvals(rda, model = "constrained")
eig.tbl <- tibble(RDA = c(seq(1, length(eig))), explain_var = (eig))
eig.plot <- ggplot(eig.tbl, aes(RDA,explain_var)) +
            geom_col(fill = heat.colors(length(eig.tbl$RDA)))+
#            geom_text(aes(label = explain_var, size = 3))+
#            scale_y_continuous(labels = scales::percent_format(accuracy = 0.001))+
            labs(title = "Eigenvalues",
#                 x = "RDA",
                 y = "Explained variance (%)")+
            theme_bw()
eig.plot  
ggsave("eigenvalues.svg", device = "svg", width = 15, height = 10, units = "cm", dpi = 600)
# RDA triplot - https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html
# Scaling 2 shows the effects of explanatory variables.
# Longer arrows mean this variable strongly drives the variation in the community matrix.
# Arrows pointing in opposite directions have a negative relationship.
# Arrows pointing in the same direction have a positive relationship.
rda.scaling <- summary(rda, scaling = 1)
rda.site <- data.frame(rda.scaling$sites)[1:2]
rda.species <- data.frame(rda.scaling$species)[1:2]
rda.env <- data.frame(rda.scaling$biplot)[1:2]

# Population data
popdata <- read.delim("popdata.tsv", header = TRUE, sep = "\t")

#cols <- brewer.pal(n = nPop(gl), name = "Set1")
#display.brewer.all()
# 1. Visualize a single RColorBrewer palette 
# by specifying its name
#display.brewer.pal(n, name)
# 2. Return the hexadecimal color code of the palette
#brewer.pal(n, name)
#display.brewer.all(colorblindFriendly = T)

# RDA ggplot2
rda.site <- cbind(rda.site, popdata[1:2])

rda.plot <- rda.site %>%
            arrange(sample) %>%
            mutate(Pop = factor(pop, levels=c("WCA","SWA","ECA","SWI","SEI","WCP"))) %>%
            ggplot(aes(RDA1, RDA2)) +
            geom_point(aes(color = Pop), size = 3, alpha = 0.5) +
            stat_ellipse(aes(color = Pop), level = 0.95,
                         show.legend = FALSE, lty = "solid") +
#                         scale_color_brewer(palette = "Set1", 5) +
            scale_color_manual(values = c("#4DAF4A","#E41A1C","#FF7F00","#377EB8","#984EA3","#A65628")) +
            labs(#title = "Atlantic X Indo-Pacific",
                 x = paste0("RDA1 (", as.character(round(eig.tbl$explain_var[1], 2)), "%)"),
                 y = paste0("RDA2 (", as.character(round(eig.tbl$explain_var[2], 2)), "%)")) +
            geom_vline(xintercept = 0, color = 'gray', size = 0.5, lty = "dashed") + 
            geom_hline(yintercept = 0, color = 'gray', size = 0.5, lty = "dashed") +
            geom_segment(data = rda.env,
                         aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
                         arrow = arrow(length = unit(0.2, 'cm')),
                         size = 1, color = c("#000000","#000000","#000000","#000000","#2171B5","#E41A1C","#66A61E","#000000"), alpha = 0.5)+
            geom_text(data = rda.env,
                      aes(RDA1 * 1.1, RDA2 * 1.1,
                      label = rownames(rda.env)),
                      color = c("#2171B5","#E41A1C","#66A61E","#000000"), size = 4,
                      position=position_jitter(width=0.03,height=0.05))+
                      annotate("text", x = -0.5, y = 0.25,
                      label = paste0("adj.R^2== ", as.character(round(adj.R2[1], 3)), ""), angle = 0, vjust = 0,
                      size = 5, colour="black", parse=TRUE) +
            theme_bw() +
            theme(legend.title = element_blank(),
                  legend.position = "bottom") +
            guides(color=guide_legend(nrow=1))
rda.plot
ggsave("rda.png", device = "png", width = 20, height = 12, units = "cm", dpi = 600)
ggsave("rda.svg", device = "svg", width = 25, height = 20, units = "cm", dpi = 600)

# plot with eigenvals
p1 <- rda.plot + inset_element(eig.plot, left = 0.6, bottom = 0.6, right = 0.98, top = 0.98)
p1

#--------------#
#
# Redundancy analysis 2
#
#--------------#

#--------------#
#
# Identify significant variables
#
#--------------#
#
# Global model with all environmental variables
my_rda2 <- rda(allele_freqs ~ ., data = env, scale = TRUE)
anova.cca(my_rda2, permutations = how(nperm=9999), by = NULL, model = "full",  parallel = getOption("mc.cores")) # Global RDA significance
sqrt(vif.cca(my_rda2))  # Multicollinearity before variables selection
# subset only significant independent environmental variables to include in the RDA
env.vif2 <- subset(env, select = -c(chlo,ptpk,diff.atn,sst))
str(env.vif2)
my_rda2 <- rda(allele_freqs ~ ., data = env.vif2, scale = TRUE)
anova.cca(my_rda2, permutations = how(nperm=9999), by = NULL, model = "full",  parallel = getOption("mc.cores")) # Global RDA significance
sqrt(vif.cca(my_rda2))  # Multicollinearity before variables selection
# global R2
adj.r2.3 <- RsquareAdj(my_rda2, permutations = 9999)$adj.r.squared
adj.r2.3
# Use forward selection to identify significant environmental variables - https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
env.fs2 <- forward.sel(Y = allele_freqs,
                      X = env.vif2,
                      adjR2thresh = adj.r2.3,
                      nperm = 9999,
                      R2more = 0.05,
                      alpha = 0.05)
env.fs2
# P-values from the forward selection above are in the object sel.fs$pval, and the adjustement can be done in the following way (Holm's correction):
n.tests2 <- ncol(env.fs2)  # number of tests equals to number of all variables from which is being selected
pval.adj2 <- p.adjust(env.fs2$pval, method = 'holm', n = n.tests2)
pval.adj2
env.fs2$pval.adj2 <- pval.adj2
env.fs2
# subset only significant independent environmental variables to include in the RDA
env.sig2 <- subset(env, select = env.fs2$variables)
str(env.sig2)
# Combine environmental variables and dbmems
data.sig2 <- cbind(env.sig2, dbmem.sig)
str(data.sig2)

# Perform RDA with sea surface temperature (SST)
rda2 <- rda(allele_freqs ~ ., data = data.sig2, scale = TRUE)
rda2
anova.cca(rda2, by = NULL, permutations = how(nperm=9999), model = "full",  parallel = getOption("mc.cores")) # Global RDA significance
anova.cca(rda2, by='axis', permutations = how(nperm=9999), model = "full",  parallel = getOption("mc.cores")) # Test which axis are significant
anova.cca(rda2, by='terms', permutations = how(nperm=9999), model = "full",  parallel = getOption("mc.cores")) # Test which terms are significant
coef(rda2)
sqrt(vif.cca(rda2)) # variance inflation factor (<2 OK)
adj.R2.2 <- round(RsquareAdj(rda2, permutations = 9999)$adj.r.squared, 3)
adj.R2.2

# Variance explained by each canonical axis
summary(eigenvals(rda2, model = "constrained"))
#screeplot(rda)
eig2 <- eigenvals(rda2, model = "constrained")
eig.tbl2 <- tibble(RDA = c(seq(1, length(eig2))), explain_var = (eig2))
eig.plot2 <- ggplot(eig.tbl2, aes(RDA,explain_var)) +
             geom_col(fill = heat.colors(length(eig.tbl2$RDA)))+
#             geom_text(aes(label = explain_var, size = 3))+
#             scale_y_continuous(labels = scales::percent_format(accuracy = 0.001))+
             labs(title = "Eigenvalues",
#                  x = "RDA",
                  y = "Explained variance (%)")+
             theme_bw()
eig.plot2
ggsave("eigenvalues2.svg", device = "svg", width = 15, height = 10, units = "cm", dpi = 600)

rda.scaling2 <- summary(rda2, scaling = 1)
rda.site2 <- data.frame(rda.scaling2$sites)[1:2]
rda.species2 <- data.frame(rda.scaling2$species)[1:2]
rda.env2 <- data.frame(rda.scaling2$biplot)[1:2]

# RDA ggplot2
rda.site2 <- cbind(rda.site2, popdata[1:2])

rda.plot2 <- rda.site2 %>%
             arrange(sample) %>%
             mutate(Pop = factor(pop, levels=c("WCA","SWA","ECA","SWI","SEI","WCP"))) %>%
             ggplot(aes(RDA1, RDA2)) +
             geom_point(aes(color = Pop), size = 3, alpha = 0.5) +
             stat_ellipse(aes(color = Pop), level = 0.95,
                          show.legend = FALSE, lty = "solid") +
#             scale_color_brewer(palette = "Set1", 5) +
             scale_color_manual(values = c("#4DAF4A","#E41A1C","#FF7F00","#377EB8","#984EA3","#A65628")) +
             labs(#title = "Atlantic X Indo-Pacific",
                  x = paste0("RDA1 (", as.character(round(eig.tbl2$explain_var[1], 2)), "%)"),
                  y = paste0("RDA2 (", as.character(round(eig.tbl2$explain_var[2], 2)), "%)")) +
             geom_vline(xintercept = 0, color = 'gray', size = 0.5, lty = "dashed") + 
             geom_hline(yintercept = 0, color = 'gray', size = 0.5, lty = "dashed") +
             geom_segment(data = rda.env2,
                          aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
                          arrow = arrow(length = unit(0.2, 'cm')),
                          size = 1, color = c("#2171B5","#7570B3","#66A61E","#000000"), alpha = 0.5) +
             geom_text(data = rda.env2,
                       aes(RDA1 * 1.1, RDA2 * 1.1,
                       label = rownames(rda.env2)),
                       color = c("#2171B5","#7570B3","#66A61E","#000000"), size = 4,
                       position=position_jitter(width=0.03,height=0.05)) +
             annotate("text", x = -0.5, y = 0.25,
                       label = paste0("adj.R^2== ", as.character(round(adj.R2.2[1], 3)), ""), angle = 0, vjust = 0,
                       size = 5, colour="black", parse=TRUE) +
             theme_bw() +
             theme(legend.title = element_blank(),
                   legend.position = "bottom") +
             guides(color=guide_legend(nrow=1))
rda.plot2
ggsave("rda2.png", device = "png", width = 25, height = 20, units = "cm", dpi = 600)
ggsave("rda2.svg", device = "svg", width = 25, height = 20, units = "cm", dpi = 600)

# plot with eigenvals
p2 <- rda.plot2 + inset_element(eig.plot2, left = 0.6, bottom = 0.6, right = 0.98, top = 0.98)
p2
# Merge plot
((rda.plot / rda.plot2) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.title = element_blank())) + plot_annotation(tag_levels = 'a')

((p1 / p2) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.title = element_blank())) + plot_annotation(tag_levels = 'a')
#ggsave("rdas.png", device = "png", width = 20, height = 30, units = "cm", dpi = 600)
ggsave("rdas.svg", device = "svg", width = 20, height = 30, units = "cm", dpi = 600)
save.image(file = "rda.RData")
