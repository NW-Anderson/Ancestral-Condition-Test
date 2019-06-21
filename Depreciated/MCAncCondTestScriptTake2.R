library(R.utils)
library(phytools)
library(diversitree)
library(geiger)
library(doSNOW)
library(foreach)
cl<-makeCluster(3)
on.exit(stopCluster(cl))


source('AncCond.R')

n.trees <- 200
n.taxa <- 200


##### Making fig2 ######
# we do the following for each of 200 trees
# this will hold the p.val for each of 200 tests for the 10 scaling factors
# dont need this for mc
# p.val.array <- array(dim = c(n.trees, 10))
message <- T

opts <- list(preschedule = FALSE)
registerDoSNOW(cl)
p.val.array <-foreach(t = 1:n.trees, .options.multicore=opts, .combine = 'rbind') %dopar%{
  p.val.vec <- c()
  good.tree <- F
  while(good.tree == F){
    # some trees take a real long time for simulating discrete traits due to incredibly short branch lengths where a transition 
    #   must occur. To get around this we added a timeout and trycatch. In the event of a tree that is taking an unacceptable time,
    #   we begin back here by resimulating a tree
    tryCatch({
      # We begin with a single tree and test it at every scaling factor then move to the next tree
      # first the tree
      trees <- trees(pars = c(3,1),
                     type = "bd",
                     n = 1,
                     max.taxa = n.taxa,
                     include.extinct = F)[[1]]
      trees$edge.length <- trees$edge.length / max(branching.times(trees))

      # we then simulate the continious character
      cont.trait <- sim.char(trees, 0.2, model = 'BM')
      names(cont.trait) <- trees$tip.label # this line somehow makes anc.ML work????
      
      # identifying which branch had a mean cont trait value in the upper and lower quartiles
      # we do this by 1st doing an ASR for the continious trait
      cont.trait.AC <- anc.ML(trees, cont.trait, model = "BM")
      # this will hold all of the branch means in the same order they are given in trees
      branch.means <- c()
      # branch names is essentially paste(rootward node, tipward node)
      branch.names <- c()
      # then for each branch we go through and calculate the name and mean
      for(j in 1:nrow(trees$edge)){
        # we first find the cont trait value at the rootward node
        node.o.int <- trees$edge[j,1]
        # we have to look in two different places for cont trait values, either in the cont.trait vector 
        # (if the node is a tip) or in the ASR if it is an interior node
        if(node.o.int <= n.taxa){
          one <- cont.trait[node.o.int]
        }else{
          one <- cont.trait.AC$ace[names(cont.trait.AC$ace) == as.character(node.o.int)]
        }
        # we do the same for the tipward node
        node.o.int <- trees$edge[j,2]
        if(node.o.int <= n.taxa){
          two <- cont.trait[node.o.int]
        }else{
          two <- cont.trait.AC$ace[names(cont.trait.AC$ace) == as.character(node.o.int)]
        }
        # to find the mean we avg the rootward and the tipward cont trait values
        branch.means <- c(branch.means, mean(one, two))
        # we create branch names by pasting the rootwward and tipward node labels together
        branch.names <- c(branch.names, paste(as.character(trees$edge[j,1]),as.character(trees$edge[j,2])))
      }
      # we name the branch names for nice bookkeeping
      names(branch.means) <- branch.names
      rm(branch.names)
      # finding upper and lower quartiles
      upper <- summary(branch.means)[[4]]
      lower <- summary(branch.means)[[2]]
      
      # next we perform the following analysis on this tree for each of the scaling factors
      for(s in 1:10){
        scale.factor <- s
        # we leave the original trees un altered 
        alt.tree <- trees 
        
        # we then manipulate the branch lengths of those branches whose cont trait means are in the upper or lower quartiles
        for(j in 1:length(branch.means)){
          if(branch.means[j] < lower){alt.tree$edge.length[j] <- alt.tree$edge.length[j] / scale.factor}
          if(branch.means[j] > upper){alt.tree$edge.length[j] <- alt.tree$edge.length[j] * scale.factor}
        }
        # next we simulated a discrete trait on this altered tree
        ####### Make sure par is set up correctly ######
        # while loop is set up to make sure sufficient transitions occur on the tree
        ####### should i not do this it is making things take quite a while ???????? ########
        good.sim <- F
        count <- 0
        rate <- .02
        ####### i think the timeouts were caused by not unit izing the branch lengths ########
        withTimeout({while(good.sim == F){
          disc.trait <- sim.char(phy = alt.tree, 
                                 par = matrix(c(-rate, 0, rate, 0), 2), 
                                 model = 'discrete', 
                                 root = 1)
          if((0.1 * n.taxa) < sum(disc.trait == min(disc.trait)) && 
             sum(disc.trait == min(disc.trait)) < (0.9 * n.taxa)){
            good.sim <- T
            if(message == T){cat(min(disc.trait), max(disc.trait), ' good sim ')}
          }
          if(message == T && count %% 50 == 0){cat(min(disc.trait), 
                                                   max(disc.trait), 
                                                   '    ', 
                                                   sum(disc.trait == min(disc.trait)), 
                                                   '      ')}
          count <- count + 1
        }}, timeout = 600, onTimeout = "error")
        if(message == T){cat('\n')}
        # we now apply the AncCond test to our simulated data and record its result
        dat <- cbind(alt.tree$tip.label, cont.trait, disc.trait)
        ####### do we still want to run the anccond test with ancestral and derived traits or run the general case??? ########
        withTimeout({rslt <- AncCond(trees = trees, 
                                     data = dat, 
                                     drop.state = 2, 
                                     mat = c(0,0,1,0), 
                                     pi = c(1,0), 
                                     message = F)}, 
                    timeout = 600, onTimeout = "error")
        # saving results in arrays
        #p.val.array[t,s] <- rslt$pval
        p.val.vec[s] <- rslt$pval
        if(message == T){cat(' s = ', s)}
      }
      if(message == T){cat('\n')}
      # closting the while loop if all goes well
      good.tree <- T
    }, error = function(e){good.tree <- F})
  }
  p.val.vec
  if(message == T){
    cat('\n')
    cat(' t = ', t)
  }
}
fig2.data <- p.val.array
save(fig2.data, file = 'fig2Data.RData')
##### END FIGURE 2 #####

n.trees <- 100
scaling.factor <- 5
n.taxa <- seq(0, 180, length.out = 10)


##### Making fig3 ######

## this will hold the p.val for each of 100 tests for the 10 tree sizes
# p.val.array <- array(dim = c(n.trees, 10))
message <- T
# this time we vary the size of the tree
opts <- list(preschedule = FALSE)
registerDoSNOW(cl)
p.val.array <-foreach(s = 1:n.taxa, .options.multicore=opts, .combine = 'cbind') %dopar%{
  p.val.vec <- c()
  # for each tree size we repeat the same analysis for the same number of trees
  for(t in 1:n.trees){
    good.tree <- F
    while(good.tree == F){
      # some trees take a real long time for simulating discrete traits due to incredibly short branch lengths where a transition 
      #   must occur. To get around this we added a timeout and trycatch. In the event of a tree that is taking an unacceptable time,
      #   we begin back here by resimulating a tree
      tryCatch({
        # We begin with a single tree and test it at every scaling factor then move to the next tree
        # first the tree
        trees <- trees(pars = c(3,1),
                       type = "bd",
                       n = 1,
                       max.taxa = n.taxa,
                       include.extinct = F)[[1]]
        trees$edge.length <- trees$edge.length / max(branching.times(trees))
        
        
        
        
        # we then simulate the continious character
        cont.trait <- sim.char(trees, 0.2, model = 'BM')
        names(cont.trait) <- trees$tip.label # this line somehow makes anc.ML work????
        
        # identifying which branch had a mean cont trait value in the upper and lower quartiles
        # we do this by 1st doing an ASR for the continious trait
        cont.trait.AC <- anc.ML(trees, cont.trait, model = "BM")
        # this will hold all of the branch means in the same order they are given in trees
        branch.means <- c()
        # branch names is essentially paste(rootward node, tipward node)
        branch.names <- c()
        # then for each branch we go through and calculate the name and mean
        for(j in 1:nrow(trees$edge)){
          # we first find the cont trait value at the rootward node
          node.o.int <- trees$edge[j,1]
          # we have to look in two different places for cont trait values, either in the cont.trait vector 
          # (if the node is a tip) or in the ASR if it is an interior node
          if(node.o.int <= n.taxa){
            one <- cont.trait[node.o.int]
          }else{
            one <- cont.trait.AC$ace[names(cont.trait.AC$ace) == as.character(node.o.int)]
          }
          # we do the same for the tipward node
          node.o.int <- trees$edge[j,2]
          if(node.o.int <= n.taxa){
            two <- cont.trait[node.o.int]
          }else{
            two <- cont.trait.AC$ace[names(cont.trait.AC$ace) == as.character(node.o.int)]
          }
          # to find the mean we avg the rootward and the tipward cont trait values
          branch.means <- c(branch.means, mean(one, two))
          # we create branch names by pasting the rootwward and tipward node labels together
          branch.names <- c(branch.names, paste(as.character(trees$edge[j,1]),as.character(trees$edge[j,2])))
        }
        # we name the branch names for nice bookkeeping
        names(branch.means) <- branch.names
        rm(branch.names)
        # finding upper and lower quartiles
        upper <- summary(branch.means)[[4]]
        lower <- summary(branch.means)[[2]]
        
          # we then manipulate the branch lengths of those branches whose cont trait means are in the upper or lower quartiles
          for(j in 1:length(branch.means)){
            if(branch.means[j] < lower){alt.tree$edge.length[j] <- alt.tree$edge.length[j] / scale.factor}
            if(branch.means[j] > upper){alt.tree$edge.length[j] <- alt.tree$edge.length[j] * scale.factor}
          }
          # next we simulated a discrete trait on this altered tree
          ####### Make sure par is set up correctly ######
          # while loop is set up to make sure sufficient transitions occur on the tree
          ####### should i not do this it is making things take quite a while ???????? ########
          good.sim <- F
          count <- 0
          rate <- .02
          ####### i think the timeouts were caused by not unit izing the branch lengths ########
          withTimeout({while(good.sim == F){
            disc.trait <- sim.char(phy = alt.tree, 
                                   par = matrix(c(-rate, 0, rate, 0), 2), 
                                   model = 'discrete', 
                                   root = 1)
            if((0.1 * n.taxa) < sum(disc.trait == min(disc.trait)) && 
               sum(disc.trait == min(disc.trait)) < (0.9 * n.taxa)){
              good.sim <- T
              if(message == T){cat(min(disc.trait), max(disc.trait), ' good sim ')}
            }
            if(message == T && count %% 50 == 0){cat(min(disc.trait), 
                                                     max(disc.trait), 
                                                     '    ', 
                                                     sum(disc.trait == min(disc.trait)), 
                                                     '      ')}
            count <- count + 1
          }}, timeout = 360, onTimeout = "error")
          if(message == T){cat('\n')}
          # we now apply the AncCond test to our simulated data and record its result
          dat <- cbind(alt.tree$tip.label, cont.trait, disc.trait)
          ####### do we still want to run the anccond test with ancestral and derived traits or run the general case??? ########
          withTimeout({rslt <- AncCond(trees = trees, 
                                       data = dat, 
                                       drop.state = 2, 
                                       mat = c(0,0,1,0), 
                                       pi = c(1,0), 
                                       message = F)}, 
                      timeout = 600, onTimeout = "error")
          # saving results in arrays
          # p.val.array[t,s] <- rslt$pval
        
        if(message == T){cat('\n')}
        # closting the while loop if all goes well
        good.tree <- T
      }, error = function(e){good.tree <- F})
    }
    if(message == T){
      cat('\n')
      cat(' t = ', t)
    }
    p.val.vec[t] <- rslt$pval
  }
  if(message == T){
    cat('\n')
    cat(' s = ', s)
  }
}
fig3.data <- p.val.array
save(fig3.data, file = 'fig3Data.RData')
#######