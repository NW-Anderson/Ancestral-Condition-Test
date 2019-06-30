tree <- read.tree(file = 'Data/whales.tre')
sizes <- read.csv('Data/whale_sizes.csv')
source('AncCond.R')

library(R.utils)
library(phytools)
library(diversitree)
library(geiger)
library(doSNOW)
library(foreach)
cl<-makeCluster(3, type="SOCK")
on.exit(stopCluster(cl))
opts <- list(preschedule = FALSE)
registerDoSNOW(cl)


reordered.sizes <- rep(NA, length = length(tree$tip.label))
names(reordered.sizes) <- tree$tip.label

for(i in 1:length(sizes$Australophocaena.dioptrica)){
  reordered.sizes[match(gsub( ' ', '_', as.character(sizes$Australophocaena.dioptrica[[i]])),
                        names(reordered.sizes))] <- sizes$X1.86[[i]]
  # if(is.na(match(gsub( ' ', '_', as.character(sizes$Australophocaena.dioptrica[[i]])),
  #                 names(reordered.sizes)))){stop(cat('i = ',i))}
}
sig.vector <- foreach(i = 1:100, .options.multicore=opts, .combine = 'c', 
                      .packages=c("phytools","diversitree","geiger")) %dopar% {
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
                        rslt$pval < .05
                      }
paste('With no relation the AncCond reportedsignificant correlation ', sum(sig.vector), '% of the time.')
ct.vec <- reordered.sizes[!is.na(reordered.sizes)]
##### scaling tree #####
cont.trait.AC <- anc.ML(tree, ct.vec, model = "BM")
# this will hold all of the branch means in the same order they are given in trees
branch.means <- c()
# branch names is essentially paste(rootward node, tipward node)
branch.names <- c()
# then for each branch we go through and calculate the name and mean
for(j in 1:nrow(tree$edge)){
  # we first find the cont trait value at the rootward node
  node.o.int <- tree$edge[j,1]
  # we have to look in two different places for cont trait values, either in the cont.trait vector 
  # (if the node is a tip) or in the ASR if it is an interior node
  if(node.o.int <= 87){
    one <- reordered.sizes[node.o.int]
  }else{
    one <- cont.trait.AC$ace[names(cont.trait.AC$ace) == as.character(node.o.int)]
  }
  # we do the same for the tipward node
  node.o.int <- tree$edge[j,2]
  if(node.o.int <= 87){
    two <- reordered.sizes[node.o.int]
  }else{
    two <- cont.trait.AC$ace[names(cont.trait.AC$ace) == as.character(node.o.int)]
  }
  # to find the mean we avg the rootward and the tipward cont trait values
  if(sum(c(is.na(one), is.na(two))) == 0){
    branch.means <- c(branch.means, mean(one, two, na.rm = T))
  }else if(is.na(one)){
    branch.means <- two
  }else if(is.na(two)){
    branch.means <- one
  }
  # we create branch names by pasting the rootwward and tipward node labels together
  branch.names <- c(branch.names, paste(as.character(tree$edge[j,1]),as.character(tree$edge[j,2])))
}
# we name the branch names for nice bookkeeping
# names(branch.means) <- branch.names
rm(branch.names)
# finding upper and lower quartiles
upper <- summary(branch.means)[[5]]
lower <- summary(branch.means)[[2]]
scale.factor <- 5
# we leave the original trees un altered 
alt.tree <- tree 

# we then manipulate the branch lengths of those branches whose cont trait means are in the upper or lower quartiles
for(j in 1:length(branch.means)){
  if(branch.means[j] < lower){alt.tree$edge.length[j] <- alt.tree$edge.length[j] / scale.factor}
  if(branch.means[j] > upper){alt.tree$edge.length[j] <- alt.tree$edge.length[j] * scale.factor}
}
sig.vector <- foreach(i = 1:100, .options.multicore=opts, .combine = 'c', 
                      .packages=c("phytools","diversitree","geiger")) %dopar% {
                        good.sim <- F
                        rate <- .1
                        while(good.sim == F){
                          disc.trait <- sim.char(phy = alt.tree,
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
                        rslt$pval < .05
                      }
paste('With a scaling factor of 5, the AncCond test found a significant relationship correlation', sum(sig.vector),
      '% of the time.')