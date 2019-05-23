library(R.utils) 
library(foreach)
library(phytools)
library(diversitree)
library(geiger)
# install.packages("doMC", repos="http://R-Forge.R-project.org")
# installs doMC mirror for windows machines
library(doMC)
library(profvis)
registerDoMC(30)
source("fnx.R")


statistic = "mean"
# statistic = "min"
model = "ind"
# model = "dep" # not supported
h1 = "greater"
# h1 = "lesser"
# h1 = 'twotail'
nul.iter = 100 # number of null data points desired
stoc.iter = 10 # number of stochastic maps per data point desired
message = T
plot = F
TO <- NULL

# setting the rate that is varied. .1 == no correlation
# higher values == stronger correlation
oddrate <- .55
# different tree sizes
ntaxa <- 50

smp.size <- 30

opts <- list(preschedule=FALSE)
tests <- foreach(h = 1:smp.size, .options.multicore=opts, .combine = 'c') %dopar% {
  # creating a random tree
  set.seed(h)
  tree <- trees(pars=c(1,.1), type="bd", max.taxa=ntaxa)[[1]]
  tree$edge.length <- tree$edge.length/max(branching.times(tree))
  # creating simulated tip state and organizing into a data frame
  par <- rbind(c(0, .1, .5, 0), c(.1, 0, 0 ,.5), c(.5, 0, 0, oddrate), c(0, .5, .1, 0))
  # tests if each state occurs in the stochastic mapping
  good.data2 <- F
  while(good.data2 == F){ 
    # tests if the simulated tip states have an acceptable number of each state
    good.data1 <- F
    while(good.data1 == F){
      data <- as.integer(rTraitDisc(phy = tree, 
                                    model = par, 
                                    k = 4,
                                    states = 1:4,
                                    root.value = sample(1:4, 1)))
      # changes the 1 to 4 from rTraitDisc to 11,12,21,22 for entering testcorr function
      data1 <- data2 <- c()
      for(g in 1:length(data)){
        if(data[g] == 1){
          data1 <- c(data1, 1)
          data2 <- c(data2, 1)
        }else if(data[g] == 2){
          data1 <- c(data1, 1)
          data2 <- c(data2, 2)
        }else if(data[g] == 3){
          data1 <- c(data1, 2)
          data2 <- c(data2, 1)
        }else if(data[g] == 4){
          data1 <- c(data1, 2)
          data2 <- c(data2, 2)
        }
      }
      if(sum(data1 == 1)>.3*ntaxa && sum(data1 == 1)<.7*ntaxa &&
         sum(data2 == 1)>.3*ntaxa && sum(data2 == 1)<.7*ntaxa) good.data1 <- T
    }
    tip.dat <- data.frame(tree$tip.label,data1,data2)
    tester <- VecCom(tip.dat)
    if('12' %in% tester && '11' %in% tester && '21' %in% tester && '22' %in% tester) good.data2 <- T
  }
  # performing the tempcorr test on the simulated tree and data
  set.seed(Sys.time())
  rap1 <- TestCorr(tree, 
                  tip.dat, 
                  statistic, 
                  h1, 
                  model, 
                  nul.iter,
                  stoc.iter,
                  message,
                  plot)
  set.seed(Sys.time())
  rap2 <- TestCorr(tree, 
                   tip.dat, 
                   statistic, 
                   h1, 
                   model, 
                   nul.iter,
                   stoc.iter,
                   message,
                   plot)
  pval <- rbind(rap1$`p-values`,rap2$`p-values`)
  obs <- rbind(rap1$observed,rap2$observed)
  sim <- rbind(rap1$null,rap2$null)
  list(pval, obs, sim)

}
tests