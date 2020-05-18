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
              max.taxa = 100,
              include.extinct = F)[[1]]
tree$edge.length <- tree$edge.length / max(branching.times(tree))



good.sim <- F
while(good.sim == F){

  anc.state.dt <- sim.history(tree=tree,
                              Q=matrix(c(-rate, rate, rate, -rate), 2,2),
                              nsim=1)

  if((0.05 * n.taxa) < sum(anc.state.dt$states == min(anc.state.dt$states)) &&
     sum(anc.state.dt$states == min(anc.state.dt$states)) < (.95 * n.taxa)){
    good.sim <- T
  }
}





#
# good.sim <- F
# while(good.sim == F){
#   disc.trait <- sim.char(phy = tree,
#                          par = matrix(c(-rate, rate, rate, -rate), 2),
#                          model = 'discrete',
#                          root = sample(c(1,2),1))[,1,1]
#   if((0.05 * n.taxa) < sum(disc.trait == min(disc.trait)) &&
#      sum(disc.trait == min(disc.trait)) < (.95 * n.taxa)){
#     good.sim <- T
#   }
# }
#
# anc.state.dt <- make.simmap(tree, disc.trait,
#                             model = matrix(c(0,2,1,0),2),
#                             nsim = 1,
#                             pi = 'equal',
#                             message = message)
plot(anc.state.dt)
cont.trait <- OUwie.sim(phy = anc.state.dt,
                        simmap.tree = T,
                        alpha = c(1,1),
                        sigma.sq = c(.75,.75),
                        theta0 = 2.5,
                        theta = c(1,5))
tip.data <- cont.trait$X
names(tip.data) <- cont.trait$Genus_species
contMap(tree, x=tip.data)
dat <- data.frame(cont.trait, anc.state.dt$states)

res.drop1 <- AncCond(tree,
               dat,
               drop.state=NULL,
               mat=c(0,2,1,0),
               pi="estimated",
               n.tails = 1,
               nsim = 12,
               iter = 101,
               message = T)

