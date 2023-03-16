library(tidyverse)
library(patchwork)
library(viridis)
library(RColorBrewer)

# set working directory
setwd("")

# get input file names
files <- dir("./", pattern = ".tsv")

# meta data from samples
meta.data <- read.table("popdata", header = TRUE, sep = "\t")

# create a data frame with all samples
data <- files %>%
  map_dfr(read_tsv, col_names = FALSE) %>%
  select(X2:X7) %>%
  rename(sample_id = X2,
         contig = X3,
         start = X4,
         end = X5,
         length = X6,
         snps = X7
  ) %>%
  as_tibble()

# calculates the number and the accumulated lengths of ROH, and the F_ROH
data.sum <- data %>%
  group_by(sample_id) %>%
  summarise(n_roh = n(),
            sum_roh = sum(length)
  ) %>%
  mutate(f_roh = sum_roh/619037273) %>%
  add_column(sample = meta.data$sample, pop = meta.data$pop)

# save table
write.table(data.sum,file="sfa.roh.csv",quote = FALSE,sep = ",",row.names = FALSE)

# create basic plot for F_ROH
p2 <- ggplot(data = data.sum, mapping = aes(x = sample, y = f_roh)) +
# add barplot
      geom_col(mapping = aes(fill = pop), width = 0.9) +
# add matrix of panels defined by two column faceting variables
      facet_grid(cols = vars(pop), scales = "free_x", space = "free_x") +
#  facet_grid(fct_inorder(as.character(pop)) ~ .,
#                        # switch = "x",
#                         scales = 'free_x',
#                         space = 'free_x') +
# change color palette
#  scale_fill_viridis(option="", discrete=TRUE) +
      scale_fill_brewer(palette = "Set1", direction=1) +
      scale_y_continuous(expand = c(0, 0),labels = scales::percent_format(accuracy = 0.01)) +
#  scale_fill_manual(values = cbPalette) +
# rename Y label
      ylab(expression(F[ROH])) +
      theme_bw() +
# adjust plot appearance
      theme(axis.title.x = element_blank(),
            axis.text.y = element_text(size = 8, angle = 0, hjust = 0.5),
            axis.text.x = element_blank(),
#            axis.text.x = element_text(size = 8,angle = 90, vjust = 0.5,hjust = 1),
#            axis.ticks.x = element_blank(),
            strip.background.x = element_blank(),
            strip.text.x = element_blank(),
            legend.position = "none")

p3 <- ggplot(data = data.sum, mapping = aes(x = sum_roh/1000000, y = n_roh)) +
  # add barplot
  geom_point(mapping = aes(colour = pop)) +
  # change color palette
  #scale_fill_manual(values = cbPalette) +
  #scale_color_viridis(option="", discrete=TRUE) +
  scale_color_brewer(palette = "Set1") +
  # rename Y label
  labs(x = "Total length of RoH (Mbp)", y = "Number of RoH") +
  theme_bw() +
  # adjust plot appearance
  theme(legend.title = element_blank(),
        axis.title = element_text(size = 11),
#        axis.title.x = element_blank(),
        axis.text = element_text(size = 9, hjust = 0.5))

# create a data frame with all samples
data <- files %>%
  map_dfr(read_tsv, col_names = FALSE) %>%
  select(X2:X7) %>%
  rename(sample = X2,
         contig = X3,
         start = X4,
         end = X5,
         length = X6,
         snps = X7
  ) %>%
  as_tibble()

# join dataframes
metadata <- as_tibble(meta.data)
roh <-full_join(data, metadata, by = "sample")

# split data
d <- roh %>% 
  group_by(sample) %>%
  mutate(category=cut(length,
                 breaks=c(0, 100, 500, 1000, 5000, 10000, 15000, 20000, 25000, 30000,
                          35000, 40000, 45000, 50000, 100000, 150000, 200000, 250000,
                          300000, 350000, 400000, 450000, 500000, 550000, 1000000,
                          1500000, 2000000, 2500000),
                 #            breaks=c(0, 100, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 75000, 100000, 150000, 200000, 250000, 300000, 350000, 400000, 450000, 500000, 550000, 1000000, 1500000, 2000000)))
                 labels=c("0-0.1","0.1-0.5","0.5-1","1-5","5-10","10-15","15-20","20-25","25-30","30-35","35-40","40-45","45-50","50-100","100-150","150-200","200-250","250-300","300-350","350-400","400-450","450-500","500-550","550-1000","1000-1500","1500-2000","2000-2500")))

# plot
p4 <- ggplot(d, aes(x=sample, y=category, color=pop)) + 
      geom_jitter() +
      facet_grid(cols = vars(pop), scales = "free_x", space = "free_x") +
      scale_color_brewer(palette = "Set1", direction=1) +
      ylab("ROH length (Kb)") +
#      scale_x_continuous(breaks = waiver(), n.breaks = NULL) +
      theme_bw() +
      theme(legend.position = "none",
            strip.text.x = element_blank(),
            axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
            legend.title = element_blank(),
            axis.title.x = element_blank(),
            axis.title = element_text(size = 11),
            axis.text = element_text(size = 9, hjust = 0.5))

# assemble all the plots (patchwork)
p1 / p2 / p3 / p4 + plot_layout(ncol = 1, nrow = 4,
                           widths = c(5, 5, 5, 5), heights = unit(c(5, 5, 5, 10),'cm')) +
  plot_annotation(tag_levels = 'a')

# save  all the plots together in '.svg' format
ggsave("sailfish.genomic_summary.png", device = "png", width = 21, height = 35, units = "cm", dpi = 600)
ggsave("sailfish.genomic_summary.svg", device = "svg", width = 21, height = 35, units = "cm", dpi = 600)

# save  all the plots together in '.tiff' format
ggsave("sailfish.genomic_summary.tiff", device = "tiff", width = 28, height = 21, units = "cm", dpi = 500)
