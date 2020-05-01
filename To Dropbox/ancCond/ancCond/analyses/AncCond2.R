##### Creating Testing Data #####
library(phytools)
library(geiger)
library(diversitree)
tree <- trees(pars = c(3,1),
              type = "bd",
              n = 1,
              max.taxa = 100,
              include.extinct = F)[[1]]
tree$edge.length <- tree$edge.length / max(branching.times(tree))

# we then simulate the continious character
cont.trait <- sim.char(tree, 0.2, model = 'BM')[,,1]

# identifying which branch had a mean cont trait value in the upper and lower quartiles
# we do this by 1st doing an ASR for the continious trait
cont.trait.AC <- anc.ML(tree, cont.trait, model = "BM")
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
  if(node.o.int <= 100){
    one <- cont.trait[node.o.int]
  }else{
    one <- cont.trait.AC$ace[names(cont.trait.AC$ace) == as.character(node.o.int)]
  }
  # we do the same for the tipward node
  node.o.int <- tree$edge[j,2]
  if(node.o.int <= 100){
    two <- cont.trait[node.o.int]
  }else{
    two <- cont.trait.AC$ace[names(cont.trait.AC$ace) == as.character(node.o.int)]
  }
  # to find the mean we avg the rootward and the tipward cont trait values
  branch.means <- c(branch.means, mean(one, two))
  # we create branch names by pasting the rootwward and tipward node labels together
  branch.names <- c(branch.names, paste(as.character(tree$edge[j,1]),as.character(tree$edge[j,2])))
}
# we name the branch names for nice bookkeeping
names(branch.means) <- branch.names
# finding upper and lower quartiles
upper <- summary(branch.means)[[5]]
lower <- summary(branch.means)[[2]]

# next we perform the following analysis on this tree for each of the scaling factors

scale.factor <- 5
# we leave the original trees un altered
alt.tree <- tree

# we then manipulate the branch lengths of those branches whose cont trait means are in the upper or lower quartiles
for(j in 1:length(branch.means)){
  if(branch.means[j] < lower){alt.tree$edge.length[j] <- alt.tree$edge.length[j] / scale.factor}
  if(branch.means[j] > upper){alt.tree$edge.length[j] <- alt.tree$edge.length[j] * scale.factor}
}
# next we simulated a discrete trait on this altered tree
# while loop is set up to make sure sufficient transitions occur on the tree
good.sim <- F
rate <- .3
while(good.sim == F){
  disc.trait <- sim.char(phy = alt.tree,
                         par = matrix(c(-rate, rate, rate, -rate), 2),
                         model = 'discrete',
                         root = sample(c(1,2),1))
  if((0.05 * 100) < sum(disc.trait == min(disc.trait)) &&
     sum(disc.trait == min(disc.trait)) < (.95 * 100)){
    good.sim <- T
  }
}
# we now apply the AncCond test to our simulated data and record its result
data <- data.frame(alt.tree$tip.label, cont.trait, disc.trait)

mc <- 1000
drop.state <- NULL
mat <- c(0,2,1,0)
pi <- 'estimated'
n.tails <- 1
message <- T


rm(list=ls()[-c(21,6,13,8,12,18,15,14)])  
trees <- tree
rm(tree)
###############################

##################################################################
#                                                                #
#  Nathan W Anderson Heath Blackmon & Richard Adams              #
#  Continuous value at nodes producing a                         #
#  derived state: August 10 2015                                 #
#                                                                #
##################################################################
##
## INPUT DATA
## trees: a phylo object
## - for multiphylo analysis. Loop through every tree and compile the the
##   of each loop to create a distribution of observed values and a multiphylo
##   null distribution

## data: a dataframe with three columns (Labels, cont trait, disc trait)
## - Labels should match taxa labels in the phylogeny
## - continuous trait should be numeric values
## - discrete trait must be coded as 1 and 2 if one is ancestral then it must be coded as 1
## There should be no missing data. If a taxa does not have available cont
## and discrete data, please prune it from the tree

## mc: the number of Monte Carlo simulations per simulated dataset
## used to calc p-value

## drop.state: should be NULL unless working under the assumption
## that one state is ancestral and the other derived and back
## transitions are not possible. Using this assumption will
## ignore continuous data from taxa in the derived state

## mat: transition matrix. Should contain the rate
## matrix for evolution of the discrete trait. Acceptable matrices are
## c(0,0,1,0), c(0,1,1,0), c(0,2,1,0)

## pi: The probabilities the root of the tree are either of the
## discrete character states same values possible as make.simmap:
## "equal", "estimated", or vector length 2 with probabilities
## for each state

## n.tails: either 1 or 2 depending on whether user has apriori hypothesis about a certain state
##################################################################




AncCond <- function(trees, data, mc = 1000, drop.state=NULL, mat=c(0,2,1,0), pi="equal", n.tails = 1, message = T) {
  ##### testing inputs #####
  if(class(trees) != 'phylo') {stop('trees must be class phylo')}
  if(!is.data.frame(data) & ncol(data) == 3){stop('data should be a dataframe with 3 columns\n(tip labels, cont data, discrete data)')}
  if(class(mc) != 'numeric' | round(mc) != mc | mc < 1){stop('mc should be a numeric positive integer integer')}
  if(!is.null(drop.state)) if(!drop.state %in% c(1,2)){stop('drop.state must be NULL, or numeric 1 or 2')}
  if(!sum(mat == c(0,0,1,0)) == 4 & !sum(mat == c(0,1,1,0)) == 4 & !sum(mat == c(0,2,1,0)) == 4){
    stop('mat must be a vector of the form c(0,0,1,0), c(0,1,1,0), or c(0,2,1,0)')
  }
  
  if((!pi %in% c('equal', 'estimated'))[1]){
    if(!is.numeric(pi)) stop('pi must be equal, estimated or a vector of length 2\nwith probabilities for the state of the discrete character at the root')
    if(length(pi) != 2 | sum(pi) != 1) stop('pi must be equal, estimated or a vector of length 2\nwith probabilities for the state of the discrete character at the root')
  }
  
  if(n.tails != 1 & n.tails != 2){stop('n.tails should be numeric 1 or 2')}
  #####
  #####
  
  ##### create named vector for disc trait for all taxa #####
  dt.vec <- data[, 3]
  names(dt.vec) <- data[, 1]
  
  
  ##### create named vector for cont trait taxa not in derived state #####
  if(!is.null(drop.state)){
    ct.data <- data[(data[, 3] != drop.state),]
    ct.vec <- as.numeric(ct.data[, 2])
    names(ct.vec) <- ct.data[, 1]
  }else{
    ct.data <- data
    ct.vec <- as.numeric(ct.data[, 2])
    names(ct.vec) <- ct.data[, 1]
  }
  if(sum(is.na(ct.vec)) > 0 | sum(is.na(dt.vec)) > 0){
    stop('There exists missing trait data for some species in the phylogeny.\n
         Please remove such taxa from the tree.')
  }
  
  #### ASR for the continuous trait  ####
  anc.states.cont.trait <- anc.ML(trees, ct.vec, model = "BM")
  
  #### Stochastic map for discrete trait  ####
  ## using stochastic mappings to nail down specific transition points
  ### TODO convert to q=MCMC to take into account uncertainty in rates of
  ### evolution for discrete trait like 100 for nsim
  # anc.state.dt <- make.simmap(trees, dt.vec,
  #                             model = matrix(mat, 2),
  #                             nsim = 1,
  #                             pi = pi,
  #                             message = F)
  
  
  anc.state.dt <- make.simmap(trees, dt.vec, 
                              model = matrix(mat,2),
                              nsim = 100,
                              pi = pi,
                              Q = 'mcmc',
                              message = T)
  
  observed.anc.cond <- list()
  null.anc.cond <- list()
  for(j in 1:100){
    current.map <-  anc.state.dt[[j]]
    #### Parse simmap to get producing nodes ####
    # the mapped edge object has time spent in a state in
    # two columns so only branches with a change have an entry
    # in both columns
    ss_nodes <- current.map$mapped.edge[, 1] > 0 &
      current.map$mapped.edge[, 2] > 0
    
    # this returns the node pairs describing a branch with origins
    wanted_branches <- ss_nodes[ss_nodes == T]
    wanted_nodes <- names(wanted_branches)
    if(sum(mat) > 1){
      # for the general model we partition the producing nodes for 1->2 and 1<-2 transitions
      producing.nodes12 <- c()
      producing.nodes21 <-c()
      trans.maps <- current.map$maps[ss_nodes == T]
      # now we take the rootward node of each branch and get rid of duplicates
      wanted_nodes <- gsub(",.*", "", wanted_nodes)
      ##### Just realized we can do this with describe.simmap :(
      ##### But i dont want to change it, it would require match function
      for(i in 1:length(wanted_nodes)){
        if(names(trans.maps[[i]])[1] == '1'){
          producing.nodes12 <- c(producing.nodes12, wanted_nodes[i])
        }else if(names(trans.maps[[i]])[1] == '2'){
          producing.nodes21 <- c(producing.nodes21, wanted_nodes[i])
        }
      }
      
      producing.nodes12 <- unique(producing.nodes12)
      producing.nodes21 <- unique(producing.nodes21)
      
      
      ##### get estimated ancestral conditions ######
      observed.anc.cond[[j]] <- list('12' = anc.states.cont.trait$ace[names(anc.states.cont.trait$ace) %in%
                                                                        producing.nodes12],
                                     '21' = anc.states.cont.trait$ace[names(anc.states.cont.trait$ace) %in%
                                                                        producing.nodes21])
      null.anc.cond[[j]] <- CreateNull(tree = trees,                     # a tree type phylo
                                       iter = 10,                        # number of simulations for null
                                       anc.state.dt = current.map,       # for Q-matrix
                                       mat = mat,                        # model specification matrix
                                       anc.states.cont.trait = anc.states.cont.trait)    # ancestral state reconstruction for continuous
      
      
      
      
      ##### DEPRECIEATED ######
      # ## get the mean ancestral value for the cont trait
      # ## at nodes producing the derived state marginalizing across trees
      # orig.val12 <- mean(anc.states.cont.trait$ace[names(anc.states.cont.trait$ace) %in%
      #                                                producing.nodes12])
      # orig.val21 <- mean(anc.states.cont.trait$ace[names(anc.states.cont.trait$ace) %in%
      #                                                producing.nodes21])
      # ## Produce the null distribution of nodes in ancestral cond
      # null.orig.val12 <- vector(length = mc)
      # null.orig.val21 <- vector(length = mc)
      # number.of.trans12 <- length(producing.nodes12)
      # number.of.trans21 <- length(producing.nodes21)
      # anc.dt <- anc.state.dt
      # anc.ct <- anc.states.cont.trait
      # node.states <- describe.simmap(anc.dt)$states
      # anc.cond.nodes12 <- anc.ct$ace[names(anc.ct$ace) %in%
      #                                  names(node.states)[node.states != '2']]
      # anc.cond.nodes21 <- anc.ct$ace[names(anc.ct$ace) %in%
      #                                  names(node.states)[node.states != '1']]
      # ##### creating null distribution #
      # ##### this is where the create null function will go
      # for (j in 1:mc){
      #   # set.seed(j)
      #   null.orig.val12[j] <- mean(sample(anc.cond.nodes12,
      #                                     length(producing.nodes12)))
      #   null.orig.val21[j] <- mean(sample(anc.cond.nodes21,
      #                                     length(producing.nodes21)))
      # }
      # ## how many more extreme
      ###### DEPRECIATED ######
      
      ###### MOVED TO END #####
      # bigger12 <- (sum(null.orig.val12 >= orig.val12) / mc)
      # smaller12 <- (sum(null.orig.val12 <= orig.val12) / mc)
      # if(!is.null(producing.nodes12)){
      #   if (bigger12 <= smaller12){pval12 <- bigger12}
      #   if (smaller12 < bigger12){pval12 <- smaller12}
      #   if (n.tails == 2){pval12 <- 2 * pval12}
      # }else{
      #   pval12 <- NA
      # }
      # 
      # bigger21 <- (sum(null.orig.val21 >= orig.val21) / mc)
      # smaller21 <- (sum(null.orig.val21 <= orig.val21) / mc)
      # if(!is.null(producing.nodes21)){
      #   if (bigger21 <= smaller21){pval21 <- bigger21}
      #   if (smaller21 < bigger21){pval21 <- smaller21}
      #   if (n.tails == 2){pval21 <- 2 * pval21}
      # }else{
      #   pval21 <- NA
      # }
      # 
      # ## print results to terminal
      # if (message == T){
      #   cat(paste(
      #     "Mean value for the continuous trait at origin oftrait 2:",
      #     round(orig.val12, digits = 4),
      #     "\n"
      #   ))
      #   cat(paste(
      #     "Mean value for the continuous trait at origin of trait 1:",
      #     round(orig.val21, digits = 4),
      #     "\n"
      #   ))
      #   cat(paste("Number of producing nodes 1->2:", round(mean(number.of.trans12),
      #                                                      digits = 4), "\n"))
      #   cat(paste("Number of producing nodes 2->1:", round(mean(number.of.trans21),
      #                                                      digits = 4), "\n"))
      #   cat(paste("Mean of null dist 1->2:", round(mean(null.orig.val12),
      #                                              digits = 4), "\n"))
      #   cat(paste("Mean of null dist 2->1:", round(mean(null.orig.val21),
      #                                              digits = 4), "\n"))
      #   cat(paste("SD of null dist 1->2:", round(sd(null.orig.val12), digits = 4), "\n"))
      #   cat(paste("SD of null dist 2->1:", round(sd(null.orig.val21), digits = 4), "\n"))
      #   
      #   cat(paste("pvalue 1->2:", round(pval12, digits = 4), "\n"))
      #   cat(paste("pvalue 2->1:", round(pval21, digits = 4), "\n\n"))
      #   if(is.null(producing.nodes12)){cat('No 1 -> 2 transitions occured NA and NaN values produced.')}
      #   if(is.null(producing.nodes21)){cat('No 2 -> 1 transitions occured NA and NaN values produced.')}
      # }
      # 
      # ## return results to user
      # results <- list()
      # results[[1]] <- orig.val12
      # results[[2]] <- number.of.trans12
      # results[[3]] <- null.orig.val12
      # results[[4]] <- pval12
      # results[[5]] <- orig.val21
      # results[[6]] <- number.of.trans21
      # results[[7]] <- null.orig.val21
      # results[[8]] <- pval21
      # names(results) <- c("OriginatingVals1->2", "NTrans1->2",
      #                     "NullDist1->2", "pval1->2","OriginatingVals2->1", "NTrans2->1",
      #                     "NullDist2->1", "pval2->1")
      ##### MOVED TO END #####
    }else{
      # now we take the rootward node of each branch and get rid of duplicates
      wanted_nodes <- gsub(",.*", "", wanted_nodes)
      producing.nodes <- unique(wanted_nodes)
      ## get the mean ancestral value for the cont trait
      ## at nodes producing the derived state marginalizing across trees
      observed.anc.cond[[j]] <-  list(anc.states.cont.trait$ace[names(anc.states.cont.trait$ace) %in%
                                                       producing.nodes])
      null.anc.cond[[j]] <- CreateNull(tree = trees,                     # a tree type phylo
                                       iter = 10,                        # number of simulations for null
                                       anc.state.dt = current.map,       # for Q-matrix
                                       mat = mat,                        # model specification matrix
                                       anc.states.cont.trait = anc.states.cont.trait)    # ancestral state reconstruction for continuous
      ##### DEPRECIATED #####
      # anc.states <- anc.states.cont.trait
      # orig.val <- mean(anc.states$ace[names(anc.states$ace) %in%
      #                                   producing.nodes])
      # 
      # ## Produce the null distribution of nodes in ancestral cond
      # null.orig.val <- vector(length = mc)
      # number.of.trans <- length(producing.nodes)
      # anc.dt <- anc.state.dt
      # anc.ct <- anc.states.cont.trait
      # node.states <- describe.simmap(anc.dt)$states
      # anc.cond.nodes <- anc.ct$ace[names(anc.ct$ace) %in%
      #                                names(node.states)[node.states != '2']]
      # 
      # for (j in 1:mc){
      #   # set.seed(j)
      #   null.orig.val[j] <- mean(sample(anc.cond.nodes,
      #                                   length(producing.nodes)))
      # }
      ##### DEPRECIATED ######
      
      ##### MOVED TO END #####
      # ## how many more extreme
      # 
      # bigger <- (sum(null.orig.val >= orig.val) / mc)
      # smaller <- (sum(null.orig.val <= orig.val) / mc)
      # if (bigger <= smaller){pval <- bigger}
      # if (smaller < bigger){pval <- smaller}
      # if (n.tails == 2){
      #   pval <- 2 * pval
      # }
      # ## print results to terminal
      # if(message == T){cat(paste(
      #   "Mean value for the continuous trait at origin of derived trait:",
      #   round(orig.val, digits = 4),
      #   "\n"
      # ))
      #   cat(paste("Number of producing nodes:", round(mean(number.of.trans),
      #                                                 digits = 4), "\n"))
      #   cat(paste("Mean of null dist:", round(mean(null.orig.val),
      #                                         digits = 4), "\n"))
      #   cat(paste("SD of null dist:", round(sd(null.orig.val), digits = 4), "\n"))
      #   
      #   cat(paste("pvalue:", round(pval, digits = 4), "\n\n"))
      # }
      # 
      # ## return results to user
      # results <- list()
      # results[[1]] <- orig.val
      # results[[2]] <- number.of.trans
      # results[[3]] <- null.orig.val
      # results[[4]] <- pval
      # names(results) <- c("OriginatingVals", "NTrans",
      #                     "NullDist", "pval")
      ##### MOVED TO END ######
    }
  }
  ###### MOVED TO END #####
  # bigger12 <- (sum(null.orig.val12 >= orig.val12) / mc)
  # smaller12 <- (sum(null.orig.val12 <= orig.val12) / mc)
  # if(!is.null(producing.nodes12)){
  #   if (bigger12 <= smaller12){pval12 <- bigger12}
  #   if (smaller12 < bigger12){pval12 <- smaller12}
  #   if (n.tails == 2){pval12 <- 2 * pval12}
  # }else{
  #   pval12 <- NA
  # }
  # 
  # bigger21 <- (sum(null.orig.val21 >= orig.val21) / mc)
  # smaller21 <- (sum(null.orig.val21 <= orig.val21) / mc)
  # if(!is.null(producing.nodes21)){
  #   if (bigger21 <= smaller21){pval21 <- bigger21}
  #   if (smaller21 < bigger21){pval21 <- smaller21}
  #   if (n.tails == 2){pval21 <- 2 * pval21}
  # }else{
  #   pval21 <- NA
  # }
  # 
  # ## print results to terminal
  # if (message == T){
  #   cat(paste(
  #     "Mean value for the continuous trait at origin oftrait 2:",
  #     round(orig.val12, digits = 4),
  #     "\n"
  #   ))
  #   cat(paste(
  #     "Mean value for the continuous trait at origin of trait 1:",
  #     round(orig.val21, digits = 4),
  #     "\n"
  #   ))
  #   cat(paste("Number of producing nodes 1->2:", round(mean(number.of.trans12),
  #                                                      digits = 4), "\n"))
  #   cat(paste("Number of producing nodes 2->1:", round(mean(number.of.trans21),
  #                                                      digits = 4), "\n"))
  #   cat(paste("Mean of null dist 1->2:", round(mean(null.orig.val12),
  #                                              digits = 4), "\n"))
  #   cat(paste("Mean of null dist 2->1:", round(mean(null.orig.val21),
  #                                              digits = 4), "\n"))
  #   cat(paste("SD of null dist 1->2:", round(sd(null.orig.val12), digits = 4), "\n"))
  #   cat(paste("SD of null dist 2->1:", round(sd(null.orig.val21), digits = 4), "\n"))
  #   
  #   cat(paste("pvalue 1->2:", round(pval12, digits = 4), "\n"))
  #   cat(paste("pvalue 2->1:", round(pval21, digits = 4), "\n\n"))
  #   if(is.null(producing.nodes12)){cat('No 1 -> 2 transitions occured NA and NaN values produced.')}
  #   if(is.null(producing.nodes21)){cat('No 2 -> 1 transitions occured NA and NaN values produced.')}
  # }
  # 
  # ## return results to user
  # results <- list()
  # results[[1]] <- orig.val12
  # results[[2]] <- number.of.trans12
  # results[[3]] <- null.orig.val12
  # results[[4]] <- pval12
  # results[[5]] <- orig.val21
  # results[[6]] <- number.of.trans21
  # results[[7]] <- null.orig.val21
  # results[[8]] <- pval21
  # names(results) <- c("OriginatingVals1->2", "NTrans1->2",
  #                     "NullDist1->2", "pval1->2","OriginatingVals2->1", "NTrans2->1",
  #                     "NullDist2->1", "pval2->1")
  ##### MOVED TO END #####
  ##### MOVED TO END #####
  # ## how many more extreme
  # 
  # bigger <- (sum(null.orig.val >= orig.val) / mc)
  # smaller <- (sum(null.orig.val <= orig.val) / mc)
  # if (bigger <= smaller){pval <- bigger}
  # if (smaller < bigger){pval <- smaller}
  # if (n.tails == 2){
  #   pval <- 2 * pval
  # }
  # ## print results to terminal
  # if(message == T){cat(paste(
  #   "Mean value for the continuous trait at origin of derived trait:",
  #   round(orig.val, digits = 4),
  #   "\n"
  # ))
  #   cat(paste("Number of producing nodes:", round(mean(number.of.trans),
  #                                                 digits = 4), "\n"))
  #   cat(paste("Mean of null dist:", round(mean(null.orig.val),
  #                                         digits = 4), "\n"))
  #   cat(paste("SD of null dist:", round(sd(null.orig.val), digits = 4), "\n"))
  #   
  #   cat(paste("pvalue:", round(pval, digits = 4), "\n\n"))
  # }
  # 
  # ## return results to user
  # results <- list()
  # results[[1]] <- orig.val
  # results[[2]] <- number.of.trans
  # results[[3]] <- null.orig.val
  # results[[4]] <- pval
  # names(results) <- c("OriginatingVals", "NTrans",
  #                     "NullDist", "pval")
  ##### MOVED TO END ######
  return(results)
}



























