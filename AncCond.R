################################################
#                                              #
#  Heath Blackmon & Richard Adams              #
#  Continuous value at nodes producing a       #
#  derived state: August 10 2015               #
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

## pi same values possible as make.simmap: "equal", "estimated" vector length
## 2 with probabilities for each state
# pi <- c(1,0)
# drop.state <- 2
AncCond <- function(trees, data, mc = 1000, drop.state=NULL, mat=c(0,2,1,0), pi="equal", message = T) {
  ## create named vector for disc trait for all taxa
  dt.vec <- data[, 3]
  names(dt.vec) <- data[, 1]
  
  ## create named vector for cont trait taxa not in derived state
  if(!is.null(drop.state)){
    ct.data <- data[data[, 3] != drop.state,]
    ct.vec <- as.numeric(ct.data[, 2])
    names(ct.vec) <- ct.data[, 1]
  }else{
    ct.data <- data
    ct.vec <- as.numeric(ct.data[, 2])
    names(ct.vec) <- ct.data[, 1]
  }
  
  ## ASR for the continuous trait
  anc.states.cont.trait <- anc.ML(trees, ct.vec, model = "BM")
  
  ## ASR for discrete trait
  ## using stochastic mappings to nail down specific transition points
  ## if(!is.null(drop.state)){
  ###### should i add a timeout????? #######
  ###### change this to sim.char???  #######
  anc.state.dt <- make.simmap(trees, dt.vec,
                              model = matrix(mat, 2),
                              nsim = 1,
                              pi = pi,
                              message = F
  )
  
  ## Parse simmap to get producing nodes
  # the mapped edge object has time spent in a state in
  # two columns so only branches with a change have an entry
  # in both columns
  ss_nodes <- anc.state.dt$mapped.edge[, 1] > 0 &
    anc.state.dt$mapped.edge[, 2] > 0
  
  # this returns the node pairs describing a branch with origins
  wanted_branches <- ss_nodes[ss_nodes == T]
  wanted_nodes <- names(wanted_branches)
  
  ##### HOW TO DO THIS #####
  if(sum(mat == c(0,2,1,0)) == 4){
    # for the general model we partition the producing nodes for 1->2 and 1<-2 transitions
    producing.nodes12 <- c()
    producing.nodes21 <-c()
    trans.maps <- anc.state.dt$maps[ss_nodes == T]
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
    ## get the mean ancestral value for the cont trait
    ## at nodes producing the derived state marginalizing across trees
    anc.states <- anc.states.cont.trait
    orig.val12 <- mean(anc.states$ace[names(anc.states$ace) %in%
                                        producing.nodes12])
    orig.val21 <- mean(anc.states$ace[names(anc.states$ace) %in%
                                        producing.nodes21])
    ## Produce the null distribution of nodes in ancestral cond
    null.orig.val12 <- vector(length = mc)
    null.orig.val21 <- vector(length = mc)
    number.of.trans12 <- length(producing.nodes12)
    number.of.trans21 <- length(producing.nodes21)
    anc.dt <- anc.state.dt
    anc.ct <- anc.states.cont.trait
    node.states <- describe.simmap(anc.dt)$states
    if(!is.null(drop.state)){
      anc.cond.nodes <- anc.ct$ace[names(anc.ct$ace) %in%
                                     names(node.states)[node.states != drop.state]]
    }
    if(is.null(drop.state)){
      anc.cond.nodes <- anc.ct$ace[names(anc.ct$ace) %in% names(node.states)]
    }
    
    for (j in 1:mc){
      # set.seed(j)
      null.orig.val12[j] <- mean(sample(anc.cond.nodes,
                                        length(producing.nodes12)))
      null.orig.val21[j] <- mean(sample(anc.cond.nodes,
                                        length(producing.nodes21)))
    }
    ## how many more extreme
    
    bigger12 <- (sum(null.orig.val12 >= orig.val12) / mc)
    smaller12 <- (sum(null.orig.val12 <= orig.val12) / mc)
    if (bigger12 <= smaller12){pval12 <- bigger12}
    if (smaller12 < bigger12){pval12 <- smaller12}
    
    bigger21 <- (sum(null.orig.val21 >= orig.val21) / mc)
    smaller21 <- (sum(null.orig.val21 <= orig.val21) / mc)
    if (bigger21 <= smaller21){pval21 <- bigger21}
    if (smaller21 < bigger21){pval21 <- smaller21}
    ## print results to terminal
    if(message == T){
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
      
      cat(paste("pvalue 1->2:", round(pval12, digits = 4), "\n\n"))
      cat(paste("pvalue 2->1:", round(pval21, digits = 4), "\n\n"))
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
  }else{
    # now we take the rootward node of each branch and get rid of duplicates
    wanted_nodes <- gsub(",.*", "", wanted_nodes)
    producing.nodes <- unique(wanted_nodes)
    ## get the mean ancestral value for the cont trait
    ## at nodes producing the derived state marginalizing across trees
    anc.states <- anc.states.cont.trait
    orig.val <- mean(anc.states$ace[names(anc.states$ace) %in%
                                      producing.nodes])
    
    ## Produce the null distribution of nodes in ancestral cond
    null.orig.val <- vector(length = mc)
    number.of.trans <- length(producing.nodes)
    anc.dt <- anc.state.dt
    anc.ct <- anc.states.cont.trait
    node.states <- describe.simmap(anc.dt)$states
    if(!is.null(drop.state)){
      anc.cond.nodes <- anc.ct$ace[names(anc.ct$ace) %in%
                                     names(node.states)[node.states != drop.state]]
    }
    if(is.null(drop.state)){
      anc.cond.nodes <- anc.ct$ace[names(anc.ct$ace) %in% names(node.states)]
    }
    
    for (j in 1:mc){
      # set.seed(j)
      null.orig.val[j] <- mean(sample(anc.cond.nodes,
                                      length(producing.nodes)))
    }
    ## how many more extreme
    
    bigger <- (sum(null.orig.val >= orig.val) / mc)
    smaller <- (sum(null.orig.val <= orig.val) / mc)
    if (bigger <= smaller){pval <- bigger}
    if (smaller < bigger){pval <- smaller}
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



























