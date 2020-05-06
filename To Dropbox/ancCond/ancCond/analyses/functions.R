quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
}

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
    quiet(sim.anc.state.dt <- make.simmap(tree = tree,
                                    null.disc.trait,
                                    Q = current.map$Q,
                                    nsim = 1,
                                    pi = pi,
                                    message = F))
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

InputTesting <- function(tree,
             data,
             mc,
             drop.state,
             mat,
             pi,
             n.tails,
             message){
  ##### testing inputs #####
  
  # TODO SUBFUNCTION THIS TESTING STUFF
  
  if(class(tree) != 'phylo') {stop('tree must be class phylo')}
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
}

UnpackData <- function(data, drop.state){
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
  return(list(dt.vec, ct.vec))
}

ProcessObserved <- function(observed.anc.cond){
  vals12 <- c()
  vals21 <- c()
  for(i in 1:length(observed.anc.cond)){
    vals12 <- c(vals12, observed.anc.cond[[i]]$'12')
    vals21 <- c(vals21, observed.anc.cond[[i]]$'21')
  }
  return(list('12' = vals12,
              '21' = vals21))
}

ProcessNull <- function(null.anc.cond){
  vals12 <- c()
  vals21 <- c()
  for(i in 1:length(null.anc.cond)){
    for(j in 1:length(null.anc.cond[[1]])){
      vals12 <- c(vals12, null.anc.cond[[i]][[j]]$'12')
      vals21 <- c(vals21, null.anc.cond[[i]][[j]]$'21')
    }
  }
  return(list('12' = vals12,
              '21' = vals21))
}

ProcessObserved2 <- function(observed.anc.cond){
  vals12 <- vector(length = length(observed.anc.cond))
  vals21 <- vector(length = length(observed.anc.cond))
  for(i in 1:length(observed.anc.cond)){
    vals12[i] <- mean(observed.anc.cond[[i]]$'12')
    vals21[i] <- mean(observed.anc.cond[[i]]$'21')
  }
  return(list('12' = vals12,
              '21' = vals21))
}

ProcessNull2 <- function(null.anc.cond){
  vals12 <- vector(length = 10 * length(null.anc.cond))
  vals21 <- vector(length = 10 * length(null.anc.cond))
  for(i in 1:length(null.anc.cond)){
    for(j in 1:length(null.anc.cond[[i]])){
      vals12[10*(i-1)+j] <- mean(null.anc.cond[[i]][[j]]$'12')
      vals21[10*(i-1)+j] <- mean(null.anc.cond[[i]][[j]]$'21')
    }
  }
  return(list('12' = vals12,
              '21' = vals21))
}