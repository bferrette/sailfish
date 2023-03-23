# Geneland in Parallel: https://github.com/carlopacioni/GenelandHandleR
# Geneland 4.9.2 (https://i-pri.org/special/Biostatistics/Software/Geneland/)
library(GenelandHandleR)
# Set working directory
setwd('~/billfishes/mitogenomes/geneland')
# load Data
hap <- read.table("sailfish_geneland.txt", header = FALSE, sep = " ")
coord <- read.table("sailfish_coord.txt", header = FALSE, sep = " ")
path.mcmc <- "~/billfishes/mitogenomes/geneland"
burnin <- 100000
nthin <- 1000
# Run Geneland in parallel
# Description
# Run nruns in parallel with either Geneland correlated or uncorrelated model with non-admixture.
# Usage
run_paral_geneland(nrun = 10,
                   ncores = NULL,
                   model = "Correlated",
                   main.dir = "output_corr",
                   spatial = TRUE,
                   path = path.mcmc,
                   gen = hap,
                   coords = coord,
                   jitter = 0.01,
                   nxdom = "auto",
                   nydom = "auto",
                   burnin = burnin,
                   npopmax = 4,
                   npopinit = 3,
                   npopmin = 2,
                   niter = 100000000,
                   nthin = nthin,
                   nullMatrix = NULL
)

# Make a plot of the sampled number of clusters during the chain
make.npop.plots(d, bur = burnin, th = nthin)

# Plot samples to map colour-coded as for cluster asignment
PlotSample2Map(
  dirIn = "./output_corrR4",
  map = "world",
  txt = "sailfish_coord.txt",
  pal = "Set1",
  UTM = FALSE,
  w = NULL,
  h = NULL,
  a = 0.2
)

# Generate a bar plot with probability of membership for each sample
path.mcmc <- "~/billfishes/mitogenomes/geneland/output_corrR1"
popdata <- read.delim("popdata.tsv", header = TRUE, sep = "\t")
PlotProbMembership(path.mcmc, pal = "Set1", popdata, "pop")
