
library(coda)
library(diversitree)
library(geiger)
library(phytools)
sig.array <- array(dim = c(100, 3))
source('AncCond.R')
l = 1
sig.array <- array(dim = c(100, 3))
while(l < 100){
  tryCatch({
    # for(l in 1:100){
    # make data
    tree <- trees(pars= c(3, 1), max.taxa = 200, type="bd", include.extinct = F)[[1]]
    tree$edge.length <- tree$edge.length/max(branching.times(tree))
    cont.trait <- sim.char(tree, 0.2)[,,1]
    # reconstruct values
    anc.states <- anc.ML(tree, cont.trait)
    node.vals <- anc.states$ace
    
    # order the reconstructions
    node.vals <- node.vals[order(node.vals)]
    searching <- T
    i <- 0
    while(searching){
      i <- i + 1
      curr.count <- sum(getDescendants(tree, node=names(node.vals)[i])<=200)
      if(curr.count>=20){
        searching <- F
      }
    }
    
    # this is the key node 
    # names(node.vals)[i]
    
    # get the tips that belong to it
    node.plus.tips <- getDescendants(tree, node=names(node.vals)[i])
    tips <- node.plus.tips[node.plus.tips<=200]
    disc.char <- rep(1, 200)
    names(disc.char) <- tree$tip.label
    
    disc.char[names(disc.char) %in% tree$tip.label[tips]] <- 2
    
    # fitDiscrete(tree, disc.char, model = matrix(c(0,0,1,0), 2))
    anc.state.dt <- make.simmap(tree, disc.char,
                                model = matrix(c(0,0,1,0), 2),
                                nsim = 1,
                                pi = c(1,0),
                                message = F)
    plot(anc.state.dt)
    
    
    # creating the discretized cont trait for pagels test
    mdn <- summary(cont.trait)[2]
    disc.trait <- disc.char
    disc.cont.trait <- cont.trait < mdn
    disc.cont.trait <- as.character(as.vector(disc.cont.trait) + 1)
    disc.trait <- as.vector(as.character(disc.trait))
    names(disc.cont.trait) <- names(disc.trait) <- tree$tip.label
    # doing pagels test
    pagel <- fitPagel(tree, disc.trait, disc.cont.trait, method = 'fitDiscrete')$P
    
    # threshold test
    X <- cbind((as.numeric(disc.trait) - 1),as.vector(cont.trait))
    colnames(X) <- c('disc.trait', 'cont.trait')
    row.names(X) <- tree$tip.label
    X <- as.matrix(X)
    sample <- 1000 # sample every 1000 steps
    ngen <- 50000 # chain length, > 2 million is suggested
    burnin <- 0.2 * ngen # 20% of all data is discarded as burnin
    thresh <- threshBayes(tree, X, ngen = ngen,
                          control = list(sample = sample, quiet=TRUE))
    thresh1 <- thresh$par[(burnin/sample + 1):nrow(thresh$par), "r"]
    class(thresh1) <- 'mcmc'
    thresh2 <- HPDinterval(thresh1)
    if(sign(thresh2[1,1]) == sign(thresh2[1,2])){thresh3 <- T}
    if(sign(thresh2[1,1]) != sign(thresh2[1,2])){thresh3 <- F}
    results <- c((pagel < 0.05), thresh3)
    
    # anccond test
    dat <- data.frame(tree$tip.label, cont.trait, disc.char)
    
    rslt <- AncCond(trees = tree, 
                    data = dat, 
                    drop.state = 2, 
                    mat = c(0,0,1,0), 
                    pi = c(1,0), 
                    message = T)
    #if(is.na(rslt$`pval2->1`)){
    results <- c(results, (rslt$`pval` < .05))
    #}else {
    #  results <- c(results, (rslt$`pval1->2` < .025 | rslt$`pval2->1` < .025))
    #}
    sig.array[l,] <- results
    cat(paste(l,sum(sig.array[,1],na.rm = T),sum(sig.array[,2],na.rm = T),
              sum(sig.array[,1],na.rm = T)),'\n')
    l <- l + 1
  }, error = function(e){cat(' ERROR ')})
  
}
