library(tidyverse)
library(patchwork)

# set working directory
setwd("~/Sync/rcoimbra_phd/project_northern/results/heterozygosity/")

# find input files
files <- list.files(pattern = "*.sfs")

# set factor levels
fct.lvls <- c("West African", "Kordofan", "Nubian")

# create a tibble with the number of homozygous, heterozygous, and total sites
# for all samples and calculate heterozygosity
het <- files %>%
  set_names(str_remove(files, ".sfs")) %>%
  map_dfr(read_table, .id = "sample_id", col_names = FALSE) %>%
  select(!X3) %>%
  rename(hom_sites = X1,
         het_sites = X2
  ) %>%
  mutate(total_sites = rowSums(across(where(is.numeric))),
         heterozygosity = het_sites/total_sites,
         subspecies = case_when(
           str_detect(sample_id, regex("WA|Niger"))        ~ "West African",
           str_detect(sample_id, regex("GNP|SNR|ZAK|ZNP")) ~ "Kordofan",
           str_detect(sample_id, regex("ETH|MF"))          ~ "Nubian",
           TRUE                                            ~ NA_character_
           )
  ) %>%
  mutate(subspecies = fct_relevel(subspecies, fct.lvls))

# calculates mean, standard deviation (SD), standard error (SE),
# and confidence interval (CI)
het.sum <- het %>%
  group_by(sample_id, subspecies) %>%
  summarise(n = n(),
            mean = mean(heterozygosity) * 100,
            sd = sd(heterozygosity)
  )

# set color palette
cbPalette <- c("West African" = "#F0E442",
               "Kordofan"     = "#E69F00",
               "Nubian"       = "#D55E00")

# create basic plot for heterozygosity
p1 <- ggplot(data = het.sum,
             mapping = aes(x = sample_id,
                           y = mean)) +
  # add barplot
  geom_col(mapping = aes(fill = subspecies),
           width = 0.8) +
  # add error bars
  geom_errorbar(mapping = aes(ymin = mean-sd,
                              ymax = mean+sd),
                position = "dodge",
                width = 0.4,
                size = 0.25) +
  # add matrix of panels defined by two column faceting variables
  facet_grid(cols = vars(subspecies),
             scales = "free_x",
             space = "free_x") +
  # change color palette
  scale_fill_manual(values = cbPalette) +
  # change Y axis limit
  scale_y_continuous(limits = c(0, 5e-02)) +
  # rename Y label
  ylab("Heterozygosity (%)") +
  # adjust plot appearance
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 11),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 11),
        legend.position = "none")

###############################################################################

# set working directory
setwd("~/Sync/rcoimbra_phd/project_northern/results/roh/")

# read bamlist used with ANGSD
bam.list <- read.table("bamlist")[,1]

# extract sample IDs from bamlist
sample.id <- bam.list %>%
  str_replace("/.*/", "") %>%
  str_replace(".clean.bam", "") %>%
  as_tibble_col(column_name = "sample_id")

# read tsv of realized inbreeding coefficients (F_ROH)
f.roh <- tibble(read_tsv("realized_autozygosity.txt"))

# create a tibble with F_ROH for all samples
f.roh <- f.roh %>%
  add_column(sample.id,
             subspecies = case_when(
               str_detect(sample_id, regex("WA|Niger"))        ~ "West African",
               str_detect(sample_id, regex("GNP|SNR|ZAK|ZNP")) ~ "Kordofan",
               str_detect(sample_id, regex("ETH|MF"))          ~ "Nubian",
               TRUE                                            ~ NA_character_
             ),
             .before = "R_2"
  ) %>%
  select(!NonHBD) %>%
  pivot_longer(cols = starts_with("R_"),
               names_to = "hbd_class",
               values_to = "f_roh"
  ) %>%
  mutate(hbd_class = case_when(
           hbd_class == "R_2"     ~ "2",
           hbd_class == "R_4"     ~ "4",
           hbd_class == "R_8"     ~ "8",
           hbd_class == "R_16"    ~ "16",
           hbd_class == "R_32"    ~ "32",
           hbd_class == "R_64"    ~ "64",
           hbd_class == "R_128"   ~ "128",
           hbd_class == "R_256"   ~ "256",
           hbd_class == "R_512"   ~ "512",
           hbd_class == "R_1024"  ~ ">= 1024",
           hbd_class == "R_2048"  ~ ">= 1024",
           hbd_class == "R_4096"  ~ ">= 1024",
           hbd_class == "R_8192"  ~ ">= 1024",
           hbd_class == "R_16384" ~ ">= 1024",
           hbd_class == "R_32768" ~ ">= 1024"
         ),
         hbd_class = fct_relevel(
           hbd_class,
           ">= 1024", "512", "256", "128", "64", "32", "16", "8", "4", "2"
         ),
         subspecies = fct_relevel(subspecies, fct.lvls)
  )

# create basic plot for F_ROH
p2 <- ggplot(data = f.roh,
             mapping = aes(x = sample_id,
                           y = f_roh,
                           fill = hbd_class)) +
  # add barplot
  geom_bar(stat = "identity",
           width = 0.8) +
  # add matrix of panels defined by two column faceting variables
  facet_grid(cols = vars(subspecies),
             scales = "free_x",
             space = "free_x") +
  # change color palette
  scale_fill_viridis_d(name = "HBD classes",
                       option = "C") +
  # change Y axis limit
  scale_y_continuous(limits = c(0, 0.5)) +
  # rename Y label
  ylab(expression("F"[ROH])) +
  # adjust plot appearance
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 11),
        axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 9),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.title = element_text(face = "bold", size = 9),
        legend.text = element_text(size = 9),
        legend.key.size = unit(1, "line"),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, reverse = TRUE))

###############################################################################

# get input file names
files <- list.files(pattern = "hbd_segments_")

# create a tibble with HBD segments for all samples
hbd.segs <- files %>%
  map_dfr(read_tsv) %>%
  select(!start_snp:end_snp) %>%
  filter(length >= 100000) %>%
  add_column(subspecies = NA,
             .before = "chrom"
  ) %>%
  mutate(subspecies = case_when(
           id <= 10            ~ "West African",
           id >= 11 & id <= 20 ~ "Kordofan",
           id >= 21            ~ "Nubian",
           TRUE                ~ NA_character_
         )
  ) %>%
  mutate(subspecies = fct_relevel(subspecies, fct.lvls))

# calculate the number of ROH segments (N_ROH) and
# the total length of ROH (S_ROH) per sample
s.roh_n.roh <- hbd.segs %>%
  select(id, subspecies, length) %>%
  group_by(id, subspecies) %>%
  summarise(n_roh = n(),
            s_roh = sum(length)/10^6
  )

# create basic plot for S_ROH x N_ROH
p3 <- ggplot(data = s.roh_n.roh,
             mapping = aes(x = s_roh,
                           y = n_roh)) +
  # add scatterplot
  geom_point(mapping = aes(color = subspecies,
                           shape = subspecies),
             size = 2.5) +
  # change color palette
  scale_color_manual(name = "Subspecies",
                     labels = fct.lvls,
                     values = cbPalette) +
  # change data point shapes
  scale_shape_manual(name = "Subspecies",
                     labels = fct.lvls,
                     values = c(18, 18, 18, 15, 17, 17, 16, 16)) +
  # rename X and Y labels
  labs(x = "Sum of ROH (Mbp)",
       y = "Number of ROH") +
  # adjust plot appearance
  theme_minimal() +
  theme(axis.title = element_text(size = 11),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9),
        legend.text = element_text(size = 9),
        legend.key = element_blank(),
        legend.key.size = unit(2, "line"),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# generate composed plot
p1 / p2 / p3 +
  plot_annotation(tag_levels = "A")+
  plot_layout(nrow = 3, heights = c(1/3, 1/3, 1/3))

# save plot in '.tiff' format
ggsave("het_froh.tiff", path = "~/Sync/rcoimbra_phd/project_northern/figures/",
       width = 180, height = 232, units = "mm", dpi = 300)
