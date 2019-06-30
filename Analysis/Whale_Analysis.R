tree <- read.tree(file = 'Data/whales.tre')
sizes <- read.csv('Data/whale_sizes.csv')
source('AncCond.R')

library(phytools)
library(diversitree)
library(geiger)


reordered.sizes <- rep(NA, length = length(tree$tip.label))
names(reordered.sizes) <- tree$tip.label

for(i in 1:length(sizes$Australophocaena.dioptrica)){
  reordered.sizes[match(gsub( ' ', '_', as.character(sizes$Australophocaena.dioptrica[[i]])),
                        names(reordered.sizes))] <- sizes$X1.86[[i]]
  # if(is.na(match(gsub( ' ', '_', as.character(sizes$Australophocaena.dioptrica[[i]])),
  #                 names(reordered.sizes)))){stop(cat('i = ',i))}
}
good.sim <- F
rate <- .1
while(good.sim == F){
  disc.trait <- sim.char(phy = tree,
                         par = matrix(c(-rate, 0, rate, 0), 2),
                         model = 'discrete',
                         root = 1)
  if(5 < sum(disc.trait == min(disc.trait)) && 
     sum(disc.trait == min(disc.trait)) < (length(tree$tip.label) - 5)){
    good.sim <- T
  }
}

dat <- cbind(tree$tip.label, reordered.sizes, disc.trait)
rslt <- AncCond(trees = tree, 
                data = dat, 
                drop.state = 2, 
                mat = c(0,0,1,0), 
                pi = c(1,0), 
                message = T)

