
CreateNull <- function(tree,                     # a tree type phylo
                       iter,                     # number of simulations for null
                       current.map,             # for Q-matrix
                       anc.states.cont.trait){   # ancestral state reconstruction for continuous
  # if(sum(mat) > 1){
     nulldist <- list()
  # }else{
  #   nulldist <- vector(length = iter)
  # }
  for(n in 1:iter){
    # while loop is set up to make sure sufficient transitions occur on the tree
    good.sim <- F
    while(good.sim == F){
        null.disc.trait <- sim.char(phy = tree,
                                    par = current.map$Q,
                                    model = 'discrete',
                                    root = sample(c(1,2),1))[,,1]
        if(length(unique(null.disc.trait))>1) good.sim <- T
    }
    # nullnames <- names(null.disc.trait)
    # null.disc.trait <- as.factor(null.disc.trait)
    # names(null.disc.trait) <- nullnames
    sim.anc.state.dt <- make.simmap(tree = tree,
                                    null.disc.trait,
                                    Q = current.map$Q,
                                    nsim = 1,
                                    pi = pi,
                                    message = F)
    nulldist[[n]] <- exctractAncestral(current.map = sim.anc.state.dt,
                                       anc.states.cont.trait = anc.states.cont.trait)
  }
  return(nulldist)
}



# this takes a stochastic map and continuous trait
exctractAncestral <- function(current.map,
                              anc.states.cont.trait){
#### Parse simmap to get producing nodes ####
# the mapped edge object has time spent in a state in
# two columns so only branches with a change have an entry
# in both columns
#######
# gets branches with transitions
ss_nodes <- current.map$mapped.edge[, 1] > 0 &
  current.map$mapped.edge[, 2] > 0

# this returns the node pairs describing a branch with transitions
wanted_branches <- ss_nodes[ss_nodes == T]
wanted_nodes <- names(wanted_branches)

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
  observed.anc.cond <- list('12' = anc.states.cont.trait$ace[names(anc.states.cont.trait$ace) %in%
                                                                    producing.nodes12],
                                 '21' = anc.states.cont.trait$ace[names(anc.states.cont.trait$ace) %in%
                                                                    producing.nodes21])
  return(observed.anc.cond)
}
