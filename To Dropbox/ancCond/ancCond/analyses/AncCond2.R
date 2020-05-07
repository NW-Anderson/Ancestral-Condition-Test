source('internal.functions.R')
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
## tree: a phylo object
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




AncCond <- function(tree,
                    data,
                    drop.state=NULL,
                    mat=c(0,2,1,0),
                    pi="equal",
                    n.tails = 1,
                    nsim = 100,
                    iter = 100,
                    message = F,
                    make.plot = F){
  
  
  # This will generate warnings and stop if we violate
  # any basic assumptions on incoming data
  InputTesting(tree,
               data,
               drop.state,
               mat,
               pi,
               n.tails,
               nsim,
               iter)
  # should we unit length the tree???
  # prepare data format
  unpackeddata <- UnpackData(data, drop.state)
  dt.vec <- unpackeddata[[1]]
  ct.vec <- unpackeddata[[2]]
  rm(unpackeddata)
  
  # ASR for the continuous trait
  anc.states.cont.trait <- anc.ML(tree, ct.vec, model = "BM")
  
  # Stochastic map for discrete trait using stochastic mappings to nail
  # down specific transition points
  anc.state.dt <- make.simmap(tree, dt.vec,
                              model = matrix(mat,2),
                              nsim = as.numeric(nsim),
                              pi = pi,
                              Q = 'mcmc',
                              message = message)
  
  # processing our stoch maps to extract the ancestral condition and
  # construct the null
  observed.anc.cond <- list()
  null.anc.cond <- list()
  for(j in 1:nsim){
    if(message){
      cat('\014')
      cat('Analyzing map: ',j,' of ', nsim)
    }
    observed.anc.cond[[j]] <- exctractAncestral(current.map = anc.state.dt[[j]],
                                                anc.states.cont.trait)
    
    # creating the null
    current.Q <- anc.state.dt[[j]]$Q
    current.Q[current.Q == 0] <- c(10^(-25),-10^(-25))
    null.anc.cond[[j]] <- CreateNull(tree = tree,
                                     iter = iter,
                                     current.Q = current.Q,
                                     anc.states.cont.trait = anc.states.cont.trait)
  }
  
  obs.dist <- ProcessObserved(observed.anc.cond)
  null.dist <- ProcessNull(null.anc.cond)
  
  results <- list(obs.dist, null.dist)
  names(results) <- c("observed","null")
  pvals <- CalcPVal(results)
  
  results <- list(obs.dist, null.dist,pvals)
  names(results) <- c("observed","null",'pvals')
  class(results) <- "AncCond"
  
  if(message) summary(results)
  if(make.plot) plot(results)
  
  return(results)
}

