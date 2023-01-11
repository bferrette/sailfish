# Geneland: Detection of Structure from Multilocus Genetic Data
# Geneland 4.9.2 (https://i-pri.org/special/Biostatistics/Software/Geneland/)
# CRAN - https://mran.microsoft.com/snapshot/2017-04-22/web/packages/Geneland/index.html
# Set working directory
setwd('/home/bferrette/sailfish/mitogenomes/geneland')
# load Data
hap <- read.table("sailfish_geneland.txt", header = FALSE, sep = " ")
coord <- read.table("sailfish_coord.txt", header = FALSE, sep = " ")
## Loop for multiple runs
path.mcmc <- "/home/bferrette/sailfish/mitogenomes/geneland"
nrun <- 10
burnin <- 1000
for(irun in 1:nrun)
{
## define path to MCMC directory
path.mcmc <- paste("./",irun,"/",sep="")
system(paste("mkdir ",path.mcmc))
MCMC(coordinates = coord,
     geno.hap = hap,
     qtc = NULL,
     qtd = NULL,
     ql = NULL,
     path.mcmc = path.mcmc,
     rate.max = 63,
     delta.coord = 0.01,
     shape1 = 2,
     shape2 = 20,
     npopmin = 2,
     npopinit = 4,
     npopmax = 6,
     nb.nuclei.max = 189,
     nit = 10000000,
     thinning = 1000,
     freq.model = "Correlated",
     varnpop = TRUE,
     spatial = TRUE,
     jcf = TRUE,
     filter.null.alleles = FALSE,
     prop.update.cell = 0.1,
     write.rate.Poisson.process = FALSE,
     write.number.nuclei = TRUE,
     write.number.pop = TRUE,
     write.coord.nuclei = TRUE,
     write.color.nuclei = TRUE,
     write.freq = TRUE,
     write.ancestral.freq = TRUE,
     write.drifts = TRUE,
     write.logposterior = TRUE,
     write.loglikelihood = TRUE,
     write.true.coord = TRUE,
     write.size.pop = FALSE,
     write.mean.quanti = TRUE,
     write.sd.quanti = TRUE,
     write.betaqtc = FALSE,
     miss.loc = NULL)
## MCMC postprocessing
PostProcessChain(coordinates=coord,
                 path.mcmc=path.mcmc,
                 nxdom = 600,
                 nydom = 600,
                 burnin = burnin)
}
## Computing average posterior probability
## with a burnin of 200 (* 100) iterations
lpd <- rep(NA,nrun)
for(irun in 1:nrun)
{
path.mcmc <- paste("./",irun,"/",sep="")
path.lpd <- paste(path.mcmc,"log.posterior.density.txt",sep="")
lpd[irun] <- mean(scan(path.lpd)[-(1:burnin)])
}
## Runs sorted by decreasing average posterior probability:
order(lpd,decreasing=TRUE)
# Plot of number of populations along the MCMC run
Plotnpop(path.mcmc = "/home/bferrette/sailfish/mitogenomes/geneland/1/",
         burnin,
         printit = TRUE,
         file = "Number_populations.ps",
         format = "ps")
# Maps of posterior probability of membership
PlotTessellation(coordinates=coord,
                 path.mcmc,
                 printit = TRUE,
                 path = "./")
# Map of mode of posterior distribution of population membership
PosteriorMode(coordinates=coord,
              path.mcmc,
              plotit = FALSE,
              format = "ps",
              new.dev = TRUE,
              printit = TRUE,
              file = "Map_mode_posterior_distribution_population_membership.ps",
              main.title = "Map of mode of posterior distribution of population membership")
