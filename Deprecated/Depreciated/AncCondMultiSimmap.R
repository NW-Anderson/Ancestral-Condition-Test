################################################
#                                              #
#  Nathan Anderson & Heath Blackmon & Richard Adams              #
#  Continuous value at nodes producing a       #
#  derived state: July 8th 2019              #
#                                              #
################################################
##
## INPUT DATA
## trees: a phylo or multiPhylo object

## data: a dataframe with three collumns (tips, cont trait, disc trait)
## -tips should match taxa labels in the phylogeny
## -continuous trait should be numeric values
## -discrete trait 

## mc: the number of Monte Carlo simulations per simulated dataset
## used to calc p-value

## drop.state: should be NULL unless working under the assumption 
## that one state is ancestral and the other derived and back 
## transitions are not possible. Using this assumption will
## fix the base of the tree in the ancestral condition and will
## ignore continuous data from taxa in the derived state

## mat transition matrix for make.simmap

## pi discrete state probabilities at the root
## pi same values possible as make.simmap: "equal", "estimated", 
## vector length 2 with probabilities for each state

## n.tails whether the p value will be calculated with 1 or 2 tails



AncCond <- function(trees, data, mc = 1000, drop.state=NULL, mat=c(0,2,1,0), pi="equal", n.tails = 1, message = T) {
  ## create named vector for disc trait for all extant taxa from input data
  dt.vec <- data[, 3]
  names(dt.vec) <- data[, 1]
  
  ## create named vector for cont trait from input data
  ## taxa not in derived state optionally not included from drop.state param
  if(!is.null(drop.state)){
    ct.data <- data[(data[, 3] != drop.state),]
    ct.vec <- as.numeric(ct.data[, 2])
    names(ct.vec) <- ct.data[, 1]
    ct.vec <- ct.vec[!is.na(ct.vec)]
  }else{
    ct.data <- data
    ct.vec <- as.numeric(ct.data[, 2])
    names(ct.vec) <- ct.data[, 1]
    ct.vec <- ct.vec[!is.na(ct.vec)]
  }
  
  ## Ancestral state reconstruction for the continuous trait
  anc.states.cont.trait <- anc.ML(trees, ct.vec, model = "BM")
  
  ## ASR for discrete trait
  ## using stochastic mappings to nail down specific transition points
  
  ## Here we differentiate between uni and bi directional discrete evolution
  if(sum(mat == c(0,2,1,0)) == 4){
    ## total null will hold all 1000 null data points for each of the 10 stoch simulations in a 1000x 10 array
    ## null data points are calculated as the mean cont trait value of randomly selected nodes in the ancestral state
    ## we differentiate between transitions in the 1 -> 2 and 2 -> 1 directions
    total.null12 <- array(dim = c(1000,10))
    total.null21 <- array(dim = c(1000,10))
    ## total orig will be a vector of length 10 holding the mean of the cont trait value at the producing nodes 
    ## in each simulation
    total.orig12 <- c()
    total.orig21 <- c()
    # looping for 10 simulations
    for(nsim in 1:10){
      ## creating the stochastic simulations
      anc.state.dt <- make.simmap(trees, dt.vec,
                                  model = matrix(mat, 2),
                                  nsim = 1,
                                  pi = pi,
                                  message = F)
      ## Parse simmap to get producing nodes
      ## the mapped edge object has time spent in a state in
      ## two columns so only branches with a change have an entry
      ## in both columns
      ss_nodes <- anc.state.dt$mapped.edge[, 1] > 0 &
        anc.state.dt$mapped.edge[, 2] > 0
      
      ## this returns the node pairs describing a branch with origins
      wanted_branches <- ss_nodes[ss_nodes == T]
      wanted_nodes <- names(wanted_branches)
      
      
      ## for this scenario we partition the producing nodes for 1->2 and 1<-2 transitions
      producing.nodes12 <- c()
      producing.nodes21 <-c()
      ## transpams is useed to determine if the rootward node from a transition in=s in state 1 or 2
      trans.maps <- anc.state.dt$maps[ss_nodes == T]
      # now we take the rootward node of each branch and get rid of duplicates
      wanted_nodes <- gsub(",.*", "", wanted_nodes)
      ## here we determine if the rootward node is in state 1 or 2
      ##### Just realized we can do this with describe.simmap :( 
      ##### But i dont want to change it, it would require match function
      for(i in 1:length(wanted_nodes)){
        if(names(trans.maps[[i]])[1] == '1'){
          producing.nodes12 <- c(producing.nodes12, wanted_nodes[i])
        }else if(names(trans.maps[[i]])[1] == '2'){
          producing.nodes21 <- c(producing.nodes21, wanted_nodes[i])
        }
      }
      ## getting rid of duplicates
      producing.nodes12 <- unique(producing.nodes12)
      producing.nodes21 <- unique(producing.nodes21)
      ## get the mean ancestral value for the cont trait
      ## at nodes producing the derived state marginalizing across trees
      anc.states <- anc.states.cont.trait
      orig.val12 <- mean(anc.states$ace[names(anc.states$ace) %in%
                                          producing.nodes12])
      orig.val21 <- mean(anc.states$ace[names(anc.states$ace) %in%
                                          producing.nodes21])
      ## creating the null data points as the mean of the ancestral nodes with a sample 
      ## equal to that of the number of producing nodes
      null.orig.val12 <- vector(length = mc)
      null.orig.val21 <- vector(length = mc)
      ## finding number of transitions in either direction
      number.of.trans12 <- length(producing.nodes12)
      number.of.trans21 <- length(producing.nodes21)
      anc.dt <- anc.state.dt
      anc.ct <- anc.states.cont.trait
      node.states <- describe.simmap(anc.dt)$states
      anc.cond.nodes12 <- anc.ct$ace[names(anc.ct$ace) %in%
                                       names(node.states)[node.states != '2']]
      anc.cond.nodes21 <- anc.ct$ace[names(anc.ct$ace) %in%
                                       names(node.states)[node.states != '1']]
      
      for (j in 1:mc){
        # set.seed(j)
        null.orig.val12[j] <- mean(sample(anc.cond.nodes12,
                                          length(producing.nodes12)))
        null.orig.val21[j] <- mean(sample(anc.cond.nodes21,
                                          length(producing.nodes21)))
      }
      
      
      ## putting this sumulations data into the overarching data frames
      total.null12[,nsim] <- null.orig.val12
      total.null21[,nsim] <- null.orig.val21
      total.orig12[nsim] <- orig.val12
      total.orig21[nsim] <- orig.val21
    }
    ## to caculate the final emperical data and null data points we find the 
    ## mean of the orig val and a null data point from each of the sumulations
    ## im starting to doubt this was the correct way of doing things
    final.null12 <- vector(length = 1000)
    final.null21 <- vector(length = 1000)
    final.orig12 <- mean(total.orig12)
    final.orig21 <- mean(total.orig21)
    for(r in 1:1000){
      final.null12[r] <- mean(total.null12[r,])
      final.null21[r] <- mean(total.null21[r,])
    }
    ## using final results to calculate a pvalue
    bigger12 <- (sum(final.null12 >= final.orig12) / mc)
    smaller12 <- (sum(final.null12 <= final.orig12) / mc)
    if(!is.null(producing.nodes12)){ 
      if (bigger12 <= smaller12){pval12 <- bigger12}
      if (smaller12 < bigger12){pval12 <- smaller12}
      if (n.tails == 2){pval12 <- 2 * pval12}
    }else{
      pval12 <- NA
    }
    
    bigger21 <- (sum(final.null21 >= final.orig21) / mc)
    smaller21 <- (sum(final.null21 <= final.orig21) / mc)
    if(!is.null(producing.nodes21)){ 
      if (bigger21 <= smaller21){pval21 <- bigger21}
      if (smaller21 < bigger21){pval21 <- smaller21}
      if (n.tails == 2){pval21 <- 2 * pval21}
    }else{
      pval21 <- NA
    }
    ## print results to terminal
    if (message == T){
      cat(paste(
        "Mean value for the continuous trait at origin oftrait 2:",
        round(orig.val12, digits = 4),
        "\n"
      ))
      cat(paste(
        "Mean value for the continuous trait at origin of trait 1:",
        round(orig.val21, digits = 4),
        "\n"
      ))
      cat(paste("Number of producing nodes 1->2:", round(mean(number.of.trans12), 
                                                         digits = 4), "\n"))
      cat(paste("Number of producing nodes 2->1:", round(mean(number.of.trans21), 
                                                         digits = 4), "\n"))
      cat(paste("Mean of null dist 1->2:", round(mean(null.orig.val12), 
                                                 digits = 4), "\n"))
      cat(paste("Mean of null dist 2->1:", round(mean(null.orig.val21), 
                                                 digits = 4), "\n"))
      cat(paste("SD of null dist 1->2:", round(sd(null.orig.val12), digits = 4), "\n"))
      cat(paste("SD of null dist 2->1:", round(sd(null.orig.val21), digits = 4), "\n"))
      
      cat(paste("pvalue 1->2:", round(pval12, digits = 4), "\n"))
      cat(paste("pvalue 2->1:", round(pval21, digits = 4), "\n\n"))
      if(is.null(producing.nodes12)){cat('No 1 -> 2 transitions occured NA and NaN values produced.')}
      if(is.null(producing.nodes21)){cat('No 2 -> 1 transitions occured NA and NaN values produced.')}
    }
    
    ## return results to user
    results <- list()
    results[[1]] <- orig.val12
    results[[2]] <- number.of.trans12
    results[[3]] <- null.orig.val12
    results[[4]] <- pval12
    results[[5]] <- orig.val21
    results[[6]] <- number.of.trans21
    results[[7]] <- null.orig.val21
    results[[8]] <- pval21
    names(results) <- c("OriginatingVals1->2", "NTrans1->2",
                        "NullDist1->2", "pval1->2","OriginatingVals2->1", "NTrans2->1",
                        "NullDist2->1", "pval2->1")
  ## this is the unidirectional case
  }else{
    ## total null will hold all 1000 null data points for each of the 10 stoch simulations in a 1000x 10 array
    ## null data points are calculated as the mean cont trait value of randomly selected nodes in the ancestral state
    total.null <- array(dim = c(1000,10))
    ## total orig will be a vector of length 10 holding the mean of the cont trait value at the producing nodes 
    ## in each simulation
    total.orig <- c()
    ## looping through 10 stochastic mappings
    for(nsim in 1:10){
      ## creating thye first stochastic mapping
      anc.state.dt <- make.simmap(trees, dt.vec,
                                  model = matrix(mat, 2),
                                  nsim = 1,
                                  pi = pi,
                                  message = F)
      
      ## Parse simmap to get producing nodes
      ## the mapped edge object has time spent in a state in
      ## two columns so only branches with a change have an entry
      ## in both columns
      ss_nodes <- anc.state.dt$mapped.edge[, 1] > 0 &
        anc.state.dt$mapped.edge[, 2] > 0
      
      ## this returns the node pairs describing a branch with origins
      wanted_branches <- ss_nodes[ss_nodes == T]
      wanted_nodes <- names(wanted_branches)
      ## now we take the rootward node of each branch and get rid of duplicates
      wanted_nodes <- gsub(",.*", "", wanted_nodes)
      ## getting rid of duplicates
      producing.nodes <- unique(wanted_nodes)
      ## get the mean ancestral value for the cont trait
      ## at nodes producing the derived state marginalizing across trees
      anc.states <- anc.states.cont.trait
      ## calculating the mean of the producitng nodes
      orig.val <- mean(anc.states$ace[names(anc.states$ace) %in%
                                        producing.nodes])
      
      ## Produce the null distribution of nodes in ancestral cond
      ## this vector will contain the null data points
      null.orig.val <- vector(length = mc)
      
      ## creating the null data points as the mean of the ancestral nodes with a sample 
      ## equal to that of the number of producing nodes
      number.of.trans <- length(producing.nodes)
      anc.dt <- anc.state.dt
      anc.ct <- anc.states.cont.trait
      node.states <- describe.simmap(anc.dt)$states
      anc.cond.nodes <- anc.ct$ace[names(anc.ct$ace) %in%
                                     names(node.states)[node.states != '2']]
      
      for (j in 1:mc){
        # set.seed(j)
        null.orig.val[j] <- mean(sample(anc.cond.nodes,
                                        length(producing.nodes)))
      }
      
      ## putting this sumulations data into the overarching data frames
      total.null[,nsim] <- null.orig.val
      total.orig[nsim] <- orig.val
    }
    ## to caculate the final emperical data and null data points we find the 
    ## mean of the orig val and a null data point from each of the sumulations
    ## im starting to doubt this was the correct way of doing things
    final.null <- vector(length = 1000)
    final.orig <- mean(total.orig)
    for(r in 1:1000){
      final.null[r] <- mean(total.null[r,])
    }
    ## using the final values to calculate a pvalue
    bigger <- (sum(final.null >= final.orig) / mc)
    smaller <- (sum(final.null <= final.orig) / mc)
    if(!is.null(producing.nodes)){ 
      if (bigger <= smaller){pval <- bigger}
      if (smaller < bigger){pval <- smaller}
      if (n.tails == 2){pval <- 2 * pval}
    }else{
      pval <- NA
    }
    ## print results to terminal
    if(message == T){cat(paste(
      "Mean value for the continuous trait at origin of derived trait:",
      round(orig.val, digits = 4),
      "\n"
    ))
      cat(paste("Number of producing nodes:", round(mean(number.of.trans), 
                                                    digits = 4), "\n"))
      cat(paste("Mean of null dist:", round(mean(null.orig.val), 
                                            digits = 4), "\n"))
      cat(paste("SD of null dist:", round(sd(null.orig.val), digits = 4), "\n"))
      
      cat(paste("pvalue:", round(pval, digits = 4), "\n\n"))
    }
    
    ## return results to user
    results <- list()
    results[[1]] <- orig.val
    results[[2]] <- number.of.trans
    results[[3]] <- null.orig.val
    results[[4]] <- pval
    names(results) <- c("OriginatingVals", "NTrans",
                        "NullDist", "pval")
  }
  return(results)
}



























