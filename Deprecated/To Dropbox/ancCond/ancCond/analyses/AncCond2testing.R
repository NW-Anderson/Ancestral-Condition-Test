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



source('internal.functions.R')

AncCond <- function(tree,
                    data,
                    drop.state=NULL,
                    mat=c(0,2,1,0),
                    pi="estimated",
                    n.tails = 1,
                    nsim = 100,
                    iter = 100,
                    message = F){


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
  if(message) cat("Estimating ancestral states for the continuous trait\n")
  anc.states.cont.trait <- anc.ML(tree, ct.vec, model = "BM")

  # Stochastic map for discrete trait using stochastic mappings to nail
  # down specific transition points
  if(message) cat('Simulating stochastic mappings:\n')
  anc.state.dt <- make.simmap(tree, dt.vec,
                              model = matrix(mat,2),
                              nsim = nsim,
                              pi = pi,
                              Q = 'mcmc',
                              message = message)

  # processing our stoch maps to extract the ancestral condition and
  # construct the null
  observed.anc.cond <- list()
  null.anc.cond <- list()
  meantrans <- vector(length = 2)
  names(meantrans) <- c('12','21')
  for(j in 1:nsim){
    if(nsim == 1){
      current.map <- anc.state.dt
    }else{current.map <- anc.state.dt[[j]]}
    if(message){
      cat('\014')
      cat('Analyzing map: ',j,' of ', nsim)
    }
    observed.anc.cond[[j]] <- exctractAncestral(current.map = current.map,
                                                anc.states.cont.trait, count = T)
    trans12 <- observed.anc.cond[[j]]$ntrans[1]
    trans21 <- observed.anc.cond[[j]]$ntrans[2]
    meantrans[1] <- meantrans[1] + trans12
    meantrans[2] <- meantrans[2] + trans21
    observed.anc.cond[[j]]$ntrans <- NULL
    
    # creating the null
    null.anc.cond[[j]] <- CreateNull(tree = tree,
                                     iter = iter,
                                     current.map = current.map,
                                     anc.states.cont.trait = anc.states.cont.trait,
                                     dt.vec = dt.vec,
                                     message = message,
                                     j = j,
                                     nsim = nsim,
                                     trans12 = trans12,
                                     trans21 = trans21)
  }
  meantrans <- meantrans / nsim
  
  #########
  
  obs.dist <- ProcessObserved(observed.anc.cond)
  null.dist <- ProcessNull(null.anc.cond, iter)
  plot(density(null.dist[[1]], na.rm = T),
       main = paste('nsim = ', nsim),
       xlim= c(min(c(null.dist[[1]],
                     obs.dist[1]), na.rm = T),
               max(c(null.dist[[1]],
                     obs.dist[1]), na.rm = T)),
       xlab = 'Ancestral Condition',
       ylab = 'Frequency')
  abline(v=obs.dist[1], col = 'red')
  legend(x = 'topright', legend = c('Observed', 'Null'), col  = c('red', 'black'), lwd = 2)
  results10 <- list(obs.dist, null.dist)
  names(results10) <- c("observed","null")
  pvals <- CalcPVal(results10, n.tails)

  results10 <- list(obs.dist, null.dist,pvals, meantrans)
  names(results10) <- c("observed","null",'pvals', 'mean n trans')
  class(results10) <- "AncCond"
  
  #######
  
  observed.anc.cond5 <- observed.anc.cond[1:5]
  null.anc.cond5 <- null.anc.cond[1:5]
  
  obs.dist <- ProcessObserved(observed.anc.cond)
  null.dist <- ProcessNull(null.anc.cond, iter)
  plot(density(null.dist[[1]], na.rm = T),
       main = paste('nsim = ', 5),
       xlim= c(min(c(null.dist[[1]],
                     obs.dist[1]), na.rm = T),
               max(c(null.dist[[1]],
                     obs.dist[1]), na.rm = T)),
       xlab = 'Ancestral Condition',
       ylab = 'Frequency')
  abline(v=obs.dist[1], col = 'red')
  legend(x = 'topright', legend = c('Observed', 'Null'), col  = c('red', 'black'), lwd = 2)
  results5 <- list(obs.dist, null.dist)
  names(results5) <- c("observed","null")
  pvals <- CalcPVal(results5, n.tails)
  
  results5 <- list(obs.dist, null.dist,pvals, meantrans)
  names(results5) <- c("observed","null",'pvals', 'mean n trans')
  class(results5) <- "AncCond"
  
  ########
  
  observed.anc.cond2 <- observed.anc.cond[1:2]
  null.anc.cond2 <- null.anc.cond[1:2]
  
  obs.dist <- ProcessObserved(observed.anc.cond)
  null.dist <- ProcessNull(null.anc.cond, iter)
  plot(density(null.dist[[1]], na.rm = T),
       main = paste('nsim = ', 2),
       xlim= c(min(c(null.dist[[1]],
                     obs.dist[1]), na.rm = T),
               max(c(null.dist[[1]],
                     obs.dist[1]), na.rm = T)),
       xlab = 'Ancestral Condition',
       ylab = 'Frequency')
  abline(v=obs.dist[1], col = 'red')
  legend(x = 'topright', legend = c('Observed', 'Null'), col  = c('red', 'black'), lwd = 2)
  results2 <- list(obs.dist, null.dist)
  names(results2) <- c("observed","null")
  pvals <- CalcPVal(results2, n.tails)
  
  results2 <- list(obs.dist, null.dist,pvals, meantrans)
  names(results2) <- c("observed","null",'pvals', 'mean n trans')
  class(results2) <- "AncCond"
  
  #########
  
  observed.anc.cond <- list(observed.anc.cond[[1]])
  null.anc.cond <- list(null.anc.cond[[1]])
  
  obs.dist <- ProcessObserved(observed.anc.cond)
  null.dist <- ProcessNull(null.anc.cond, iter)
  plot(density(null.dist[[1]], na.rm = T),
       main = paste('nsim = ', 1),
       xlim= c(min(c(null.dist[[1]],
                     obs.dist[1]), na.rm = T),
               max(c(null.dist[[1]],
                     obs.dist[1]), na.rm = T)),
       xlab = 'Ancestral Condition',
       ylab = 'Frequency')
  abline(v=obs.dist[1], col = 'red')
  legend(x = 'topright', legend = c('Observed', 'Null'), col  = c('red', 'black'), lwd = 2)
  results1 <- list(obs.dist, null.dist)
  names(results1) <- c("observed","null")
  pvals <- CalcPVal(results1, n.tails)
  
  results1 <- list(obs.dist, null.dist,pvals, meantrans)
  names(results1) <- c("observed","null",'pvals', 'mean n trans')
  class(results1) <- "AncCond"
  
  # summary(results10)
  summary(results5)
  summary(results2)
  summary(results1)
  
  # obs.dist2 <- ProcessObservedOneMean(observed.anc.cond)
  # null.dist2 <- ProcessNullOneMean(null.anc.cond,iter)
  # plot(density(null.dist2$`12`, na.rm = T),
  #      main = 'One Mean',
  #      xlim= c(min(c(null.dist2$`12`,
  #                    obs.dist2$`12`),na.rm = T),
  #              max(c(null.dist2$`12`,
  #                    obs.dist2$`12`),na.rm = T)),
  #      xlab = 'Ancestral Condition',
  #      ylab = 'Frequency')
  # lines(density(obs.dist2$`12`,na.rm = T), col = 'red')
  # abline(v = obs.dist2$`12`, col = 'red')
  # legend(x = 'topright', legend = c('Observed', 'Null'), col  = c('red', 'black'), lwd = 2)
  # 
  # obs.dist3 <- ProcessObservedNoMean(observed.anc.cond)
  # null.dist3 <- ProcessNullNoMean(null.anc.cond,iter)
  # plot(hist(null.dist3$`12`, na.rm = T),
  #      main = 'No Mean',
  #      xlim= c(min(c(null.dist3$`12`,
  #                    obs.dist3$`12`)),
  #              max(c(null.dist3$`12`,
  #                    obs.dist3$`12`))),
  #      xlab = 'Ancestral Condition',
  #      ylab = 'Frequency')
  # lines(density(obs.dist3$`12`,na.rm = T), col = 'red')
  # legend(x = 'topright', legend = c('Observed', 'Null'), col  = c('red', 'black'), lwd = 2)
  

  
  
  # 
  # observed.anc.cond <- list(observed.anc.cond[[1]])
  # null.anc.cond <- list(null.anc.cond[[1]])
  # 
  # obs.dist <- ProcessObserved(observed.anc.cond)
  # null.dist <- ProcessNull(null.anc.cond, iter)
  # plot(density(null.dist[[1]], na.rm = T),
  #      main = 'Mean of Means',
  #      xlim= c(min(c(null.dist[[1]],
  #                    obs.dist[1]), na.rm = T),
  #              max(c(null.dist[[1]],
  #                    obs.dist[1]), na.rm = T)),
  #      xlab = 'Ancestral Condition',
  #      ylab = 'Frequency')
  # abline(v=obs.dist[1], col = 'red')
  # legend(x = 'topright', legend = c('Observed', 'Null'), col  = c('red', 'black'), lwd = 2)
  # results1 <- list(obs.dist, null.dist)
  # names(results1) <- c("observed","null")
  # pvals <- CalcPVal(results1, n.tails)
  # 
  # results1 <- list(obs.dist, null.dist,pvals, meantrans)
  # names(results1) <- c("observed","null",'pvals', 'mean n trans')
  # class(results1) <- "AncCond"
  # summary(results)
  # summary(results1)
  # 
  # 
  # 
  # 
  # 
  # 
  # return(results)
  # 
  # sim.anc.state.dt <- sim.history(tree=tree, Q=anc.state.dt[[1]]$Q,
  #                                 nsim=1, message = F,
  #                                 anc = 1)
}

