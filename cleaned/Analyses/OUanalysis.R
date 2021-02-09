# install.packages("phytools")
# # install.packages("diversitree")
# # install.packages("geiger")
# library(R.utils)
# library(phytools)
# library(diversitree)
# library(geiger)
# library(doMC)
# library(foreach)
# library(OUwie)
# #cl<-makeCluster(3, type="SOCK")
# #on.exit(stopCluster(cl))
# opts <- list(preschedule = FALSE)
# registerDoMC(3)

n.trees <- 100
n.taxa <- 200
message <- T
rate <- .6
source('./PackageFunctions/AncCond.R')


p.vals <- array(dim = c(n.trees,3))
colnames(p.vals) <- c('NoDrop','Drop1','Drop2')
p.vals <- foreach(t = 1:n.trees, .options.multicore=opts, .combine = 'rbind', 
                  .packages=c("phytools","diversitree","geiger",'R.utils','OUwie')) %dopar%{
                    results <- vector(length = 3)
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
                    
                    results[1] <- AncCond(tree,
                                         dat,
                                         drop.state=NULL,
                                         model = 'ARD',
                                         pi="estimated",
                                         n.tails = 1,
                                         nsim = 100,
                                         iter = 100,
                                         message = T)
                    results[2] <- AncCond(tree,
                                          dat,
                                          drop.state=1,
                                          model = 'ARD',
                                          pi="estimated",
                                          n.tails = 1,
                                          nsim = 100,
                                          iter = 100,
                                          message = T)
                    results[3] <- AncCond(tree,
                                          dat,
                                          drop.state=2,
                                          model = 'ARD',
                                          pi="estimated",
                                          n.tails = 1,
                                          nsim = 100,
                                          iter = 100,
                                          message = T)
                    results
                    
                  }
save(p.vsls, file = 'OUresults.RData')