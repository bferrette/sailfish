library(tidyverse)
library(reshape2)

# set working directory
setwd("")

# meta data of samples
meta.data <- read_csv("popfile.csv", col_names = TRUE)

# read list of SFS files
sfs.list <- read.table("sfs.list")[,1]

# extract input file names from SFS list
fin <- sub("", "", sfs.list)

# create an empty data frame
df <- data.frame(sample_id = character(),
                 sample = character(),
                 heterozygosity = double())

# iterate over input files
for (f in fin) {
  # get sample name
  sample_id <- sub("*.sfs", "", f)
  # read input file
  sfs <- read.table(f)
  # calculate heterozygosity
  het <- sfs[2]/rowSums(sfs)
  # add header
  colnames(het) <- "heterozygosity"
  # create new data frame entry
  new.df <- data.frame(sample_id = rep(sample_id, length(het)),
                       heterozygosity = het)
  # append new entries to data frame
  df <- rbind(df, new.df)
}

# calculates mean, standard deviation (SD), standard error (SE),
# and confidence interval (CI)
df.sum <- df %>%
  group_by(sample_id) %>%
  summarise(n = n(),
            mean = mean(heterozygosity),
            sd = sd(heterozygosity)
  ) %>%
  mutate(se = sd/sqrt(n)) %>%
  mutate(ci = se*qt((1-0.05)/2+0.5, n-1)) %>%
  add_column(sample = meta.data$sample, pop = meta.data$pop)
# save table
write.table(df.sum,file="sfa.heterozygosity.csv",quote = FALSE,sep = ",",row.names = FALSE)

# create basic plot
p1 <- ggplot(data = df.sum,
       mapping = aes(x = sample, y = mean)) +
# add barplot
  geom_col(mapping = aes(fill = pop),
           width = 0.8) +
# add error bars
  geom_errorbar(mapping = aes(ymin = mean-sd, ymax = mean+sd),
                position = "dodge",
                width = 0.4,
                size = 0.5) +
# add matrix of panels defined by two column faceting variables
  facet_grid(cols = vars(pop), scales = "free_x", space = "free_x") +
# change color palette
#  scale_fill_viridis(option="", discrete=TRUE) +
  scale_fill_brewer(palette = "Set1", direction=1) +
# rename X and Y labels
#  xlab("Samples") +
  ylab("Heterozygosity") +
# change title and axis font size and remove legend
  theme(axis.title = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8, angle = 90, hjust = 0.5),
        axis.text.x = element_blank(),
#        axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        strip.text.x = element_text(size = 8),
        legend.position = "none")

# save plot in '.svg' format
ggsave("sailfish.heterozygsity.svg", device = "svg", width = 21, height = 15, units = "cm", dpi = 500)
