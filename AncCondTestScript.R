library(R.utils)
# library(foreach)
library(phytools)
library(diversitree)
library(geiger)
# *****------------------------------------------- *****
# install.packages("doMC", repos="http://R-Forge.R-project.org")
# installs doMC mirror for windows machines
# ***** does not work in parallel for testing only *****
# library(doMC)
# library(profvis)




## trees: a phylo or multiPhylo object

## data: a dataframe with three collumns (tips, cont trait, disc trait)
## -tips should match taxa labels in the phylogeny
## -continuous trait should be numeric values
## -discrete trait 

# size <- 30
# trees <- trees(pars=c(1,.1), type="bd", max.taxa=size)[[1]]
# # plot(tree)
# 
# tips <- seq(1, size)
# tree$tip.label <- tips
# cont.trait <- runif(size, 0, 10)
# disc.trait <- sample(2, size, replace = T)
# 
# data <- data.frame(cbind(tips, cont.trait, disc.trait))
# rm(tips, cont.trait, disc.trait, size)
# 
# drop.stae <- 2
# mc = 1000 

# ------------------------------------- #


source('AncCond.R')

sig.array <- array(dim = c(200,10))

for(i in 1:20){
  size <- 20
  good.tree <- F
  while(good.tree == F){
    tryCatch({
      #### figure out how to prune extinct taxa ####
      trees <- sim.bdtree(d = 0, stop = 'taxa', n = size, extinct = F)
      # trees <- trees(pars=c(1,.1), type="bd", max.taxa=size)[[1]]
      cont.trait <- sim.char(trees, 0.2, model = 'BM')
      names(cont.trait) <- trees$tip.label
      # cont.trait.AC <- contMap(trees, cont.trait, plot = T)
      cont.trait.AC <- anc.ML(trees, cont.trait, model = "BM")
      branch.means <- c()
      branch.names <- c()
      for(j in 1:nrow(trees$edge)){
        node.o.int <- trees$edge[j,1]
        if(node.o.int <= size){
          one <- cont.trait[node.o.int]
        }else{
          one <- cont.trait.AC$ace[names(cont.trait.AC$ace) == as.character(node.o.int)]
        }
        node.o.int <- trees$edge[j,2]
        if(node.o.int <= size){
          two <- cont.trait[node.o.int]
        }else{
          two <- cont.trait.AC$ace[names(cont.trait.AC$ace) == as.character(node.o.int)]
        }
        branch.means <- c(branch.means, mean(one, two))
        branch.names <- c(branch.names, paste(as.character(trees$edge[j,1]),as.character(trees$edge[j,2])))
      }
      names(branch.means) <- branch.names
      rm(branch.names)
      #### or put scale factor thing here ####
      for(k in 1:10){
        scale.factor <- k
        for(j in 1:length(branch.means)){
          if(branch.means[j] < summary(branch.means)[2]){trees$edge.length[j] <- trees$edge.length[j] / scale.factor}
          if(branch.means[j] > summary(branch.means)[4]){trees$edge.length[j] <- trees$edge.length[j] * scale.factor}
        }
        sig.results <- 0
        good.tree <- T
        for(l in 1:1000){
          good <- F
          ######### ADD TIMEOUT ##########
          
          suppressWarnings(withTimeout({while(good == F){
            disc.trait <- as.integer(rTraitDisc(phy = trees, 
                                                model = matrix(c(0, 0, 1, 0), 2), 
                                                k = 2,
                                                states = 1:2,
                                                root.value = 1))
            if(min(disc.trait) != max(disc.trait)){good == T}
          }}, timeout = 100, onTimeout = "error"))
          dat <- cbind(trees$tip.label, cont.trait, disc.trait)
          rslt <- AncCond(trees, dat, drop.state = 2)
          if(rslt$pval < .05){sig.results <- sig.results + 1}
          cat('l')
        }
        sig.array[i,k] <- sig.results / 1000
      }  
      cat('k ')
    }, error = function(e){good.tree <- F})
    if(good.tree == F){cat('Timeout Occured')}}
  cat('i')
} 

# disc.trait <- sim.char(trees, list(matrix(c(0, 0, 10^6, 0), 2)), model="discrete")

