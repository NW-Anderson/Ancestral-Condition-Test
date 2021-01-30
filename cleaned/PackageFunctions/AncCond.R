###############################

##################################################################
#                                                                #
#  Nathan W Anderson Heath Blackmon & Richard Adams              #
#  Continuous value at nodes producing a                         #
#  derived state: August 10 2015                                 #
#                                                                #
##################################################################



source('internal.functions.R')

AncCond <- function(tree, # phylo object
                    data, # dataframe 3 columns species names, continuous, discrete
                    drop.state = NULL, # state to drop
                    #TODO fix to be the same as ACE
                    model = "ARD",
                    # mat = c(0, 2, 1, 0), # numeric vector default description of model for discrete trait
                    pi="estimated", # text string
                    n.tails = 1, # calc p-value on high or low values
                    #TODO evaluate for a better naming and usage structure
                    nsim = 100, # number of qmats 1 stochmap per Q
                    iter = 100, # point estimates for the null
                    message = F){


  # This will generate warnings and stop if we violate
  # any basic assumptions on incoming data
  InputTesting(tree, data, drop.state, mat,
               pi, n.tails, nsim, iter)

  # convert input tree to unit length
  tree$edge.length <- tree$edge.length / max(branching.times(tree))

  # prepare data format
  unpackeddata <- UnpackData(data, drop.state)
  dt.vec <- unpackeddata[[1]]
  ct.vec <- unpackeddata[[2]]
  rm(unpackeddata)

  if(message) cat("Estimating ancestral states for the continuous trait\n")
  # ASR for the continuous trait
  anc.states.cont.trait <- anc.ML(tree, ct.vec, model = "BM")

  if(message) cat('Simulating stochastic mappings:\n')
  # Stochastic maps for discrete trait using mcmc to account
  # for uncertainty in rates and history
  anc.state.dt <- make.simmap(tree, dt.vec,
                              model = make.model.mat(model),
                              nsim = nsim,
                              pi = pi,
                              Q = 'mcmc',
                              message = message)

  observed <- ProcessObserved(anc.state.dt = anc.state.dt,
                              nsim = nsim,
                              anc.states.cont.trait = anc.states.cont.trait,
                              message = message)

  null <- CreateNull(tree,                     # a tree type phylo
                     iter,                     # number of simulations for null
                     anc.states.cont.trait = anc.states.cont.trait,   # ancestral state reconstruction for continuous
                     anc.state.dt = anc.state.dt,
                     dt.vec = dt.vec,
                     message = message,
                     nsim = nsim)

  results <- list(observed, null)
  names(results) <- c("observed","null")
  pvals <- CalcPVal(results = results,
                    n.tails)

  results <- list(observed, null, pvals)
  names(results) <- c("observed","null",'pvals')
  class(results) <- "AncCond"
  ### TODO fix summary
  if(message) summary(results)
  return(results)
}

