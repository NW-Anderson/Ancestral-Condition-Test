#' Ancestral Condition Test
#'
#' Tests whether transitions in a discrete trait are associated with extreme
#' values of a continuous trait.
#'
#' If the rate identity matrix is returned rates are indicated by the numbers:
#' \cr rate13 demipolyploidy for state2 - odd \cr \cr \cr The argument
#' \code{constrain} can contain any of the items listed below.  The default
#' state is shown. \cr \cr \code{drop.poly=F} \cr Sets polyploidy rate to zero
#' \cr
#'
#' @param data A data.frame with 3 columns. The first is the species name.
#' 2nd is continious trait values.
#' Third is the discrete trait values.
#' @param tree A phylo object. For multiphyo analysis, Loop through every tree
#' and compile the output of each
#' @return constrained likelihood function is returned
#' @note %% ~~further notes~~
#' @author Heath Blackmon
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @examples
#'
#' @export AncCond
#' @export plot.AncCond
#' @export summary.AncCond
#' @import phytools
#' @import geiger
#' @import diversitree
#'
#'
#'
#'
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

summary.AncCond <- function(results){
  ## print results to terminal
  cat('\n')
  cat(paste(
    "Mean value for the continuous trait at 1 - > 2 transitions:",
    round(results$observed[1], digits = 4),
    "\n"
  ))
  cat(paste(
    "Mean value for the continuous trait at 2 - > 1 transitions:",
    round(results$observed[2], digits = 4),
    "\n\n"
  ))
  # cat(paste('Mean number of 1 -> 2 transitions:', round(results$`mean n trans`[1], digits = 4), '\n'))
  # cat(paste('Mean number of 2 -> 1 transitions:', round(results$`mean n trans`[2], digits = 4), '\n\n'))
  # # cat(paste("Number of producing nodes 1->2:", round(mean(number.of.trans12),
  # #                                                    digits = 4), "\n"))
  # # cat(paste("Number of producing nodes 2->1:", round(mean(number.of.trans21),
  # #                                                    digits = 4), "\n"))
  cat(paste("Mean of null dist 1->2:", round(mean(results$null$`12`, na.rm = T),
                                             digits = 4), "\n"))
  cat(paste("Mean of null dist 2->1:", round(mean(results$null$`21`, na.rm = T),
                                             digits = 4), "\n\n"))
  cat(paste("SD of null dist 1->2:", round(sd(results$null$`12`, na.rm = T), digits = 4), "\n"))
  cat(paste("SD of null dist 2->1:", round(sd(results$null$`21`, na.rm = T), digits = 4), "\n\n"))

  cat(paste("pvalue 1->2:", round(results$pvals[1], digits = 4), "\n"))
  cat(paste("pvalue 2->1:", round(results$pvals[2], digits = 4), "\n\n\n"))

  if(is.na(results$pvals[1])){
    cat('NA and NaN values are produced when no transitions of a type have occured. \n\n')
  }
  if(is.na(results$pvals[2])){
    cat('NA and NaN values are produced when no transitions of a type have occured. \n\n')
  }
}
plot.AncCond <- function(results){
  if(!is.na(results$pvals[1])){
    plot(density(results$null$`12`, na.rm = T),
         main = '1 -> 2',
         xlim= c(min(c(results$null$`12`,
                       results$observed[1])),
                 max(c(results$null$`12`,
                       results$observed[1]))),
         xlab = 'Ancestral Condition',
         ylab = 'Frequency')
    abline(v=results$observed[1], col = 'red')
    legend(x = 'topright', legend = c('Observed', 'Null'), col  = c('red', 'black'), lwd = 2)
  }
  if(!is.na(results$pvals[2])){
    plot(density(results$null$`21`, na.rm = T),
         main = '2 -> 1',
         xlim= c(min(c(results$null$`21`,
                       results$observed[2]), na.rm = T),
                 max(c(results$null$`21`,
                       results$observed[2]), na.rm = T)),
         xlab = 'Ancestral Condition',
         ylab = 'Frequency')
    abline(v=results$observed[2], col = 'red')
    legend(x = 'topright', legend = c('Observed', 'Null'), col  = c('red', 'black'), lwd = 2)
  }
}

