# library(R.utils)
# library(phytools)
# library(diversitree)
# library(geiger)
# library(doMC)
# library(foreach)
# # cl<-makeCluster(3, type="SOCK")
# # on.exit(stopCluster(cl))
# opts <- list(preschedule = FALSE)
# registerDoMC(3)

tree <- read.tree(file = 'Data/whales.tre')
tree$edge.length <- tree$edge.length / max(branching.times(tree))
sizes <- read.csv('Data/whale_sizes.csv')
source('./PackageFunctions/AncCond.R')



rate <- .2
reordered.sizes <- rep(NA, length = length(tree$tip.label))
names(reordered.sizes) <- tree$tip.label
for(i in 1:length(sizes$Australophocaena.dioptrica)){
  reordered.sizes[match(gsub( ' ', '_', as.character(sizes$Australophocaena.dioptrica[[i]])),
                        names(reordered.sizes))] <- sizes$X1.86[[i]]
  # if(is.na(match(gsub( ' ', '_', as.character(sizes$Australophocaena.dioptrica[[i]])),
  #                 names(reordered.sizes)))){stop(cat('i = ',i))}
}
reordered.sizes <- reordered.sizes[!is.na(reordered.sizes)]
tree <- keep.tip(tree, names(reordered.sizes))
sig.vector <- foreach(i = 1:300, .options.multicore=opts, .combine = 'c', 
                      .packages=c("phytools","diversitree","geiger")) %dopar% {
                        good.sim <- F
                        while(good.sim == F){
                          disc.trait <- sim.char(phy = tree,
                                                 par = matrix(c(-rate, 0, rate, 0), 2),
                                                 model = 'discrete',
                                                 root = 1)
                          if(10 < sum(disc.trait == min(disc.trait)) && 
                             sum(disc.trait == min(disc.trait)) < (length(tree$tip.label) - 10)){
                            good.sim <- T
                          }
                        }
                        
                        dat <- data.frame(tree$tip.label, reordered.sizes, disc.trait)
                        rslt <- AncCond(trees = tree, 
                                        data = dat, 
                                        drop.state = 2, 
                                        model = 'UNI', 
                                        pi = c(1,0), 
                                        message = T)
                        rslt$pval < .05
                      }
paste('With no relation the AncCond reported significant correlation ', sum(sig.vector)/3, '% of the time.')
whale.percent.sig <- sum(sig.vector) / 3
save(whale.percent.sig, file = 'whaleFP.RData')
