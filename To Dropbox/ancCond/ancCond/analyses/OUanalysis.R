# install.packages("phytools")
# install.packages("diversitree")
# install.packages("geiger")
library(R.utils)
library(phytools)
library(diversitree)
library(geiger)
library(doSNOW)
library(foreach)
library(OUwie)
# cl<-makeCluster(3, type="SOCK")
# on.exit(stopCluster(cl))
# opts <- list(preschedule = FALSE)
# registerDoSNOW(cl)

n.trees <- 100
n.taxa <- 200
message <- T
rate <- .6
source('AncCond2.R')



# We begin with a single tree and test it at every scaling factor then move to the next tree
# first the tree
tree <- trees(pars = c(3,1),
              type = "bd",
              n = 1,
              max.taxa = n.taxa,
              include.extinct = F)[[1]]
tree$edge.length <- tree$edge.length / max(branching.times(tree))

good.sim <- F
while(good.sim == F){
  disc.trait <- sim.char(phy = tree,
                         par = matrix(c(-rate, rate, rate, -rate), 2),
                         model = 'discrete',
                         root = sample(c(1,2),1))[,1,1]
  if((0.05 * n.taxa) < sum(disc.trait == min(disc.trait)) &&
     sum(disc.trait == min(disc.trait)) < (.95 * n.taxa)){
    good.sim <- T
  }
}

anc.state.dt <- make.simmap(tree, disc.trait,
                            model = matrix(c(0,2,1,0),2),
                            nsim = 1,
                            pi = 'equal',
                            message = message)

cont.trait <- OUwie.sim(phy = anc.state.dt, 
                        simmap.tree = T, 
                        alpha = c(1,.5), 
                        sigma.sq = c(.45,.9),
                        theta0 = 1,
                        theta = c(1,2))

dat <- data.frame(cont.trait, disc.trait)

res <- AncCond(tree,
               dat,
               drop.state=NULL,
               mat=c(0,2,1,0),
               pi="estimated",
               n.tails = 1,
               nsim = 12,
               iter = 101,
               message = T,
               make.plot = F)

