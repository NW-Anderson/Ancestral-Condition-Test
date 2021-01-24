##### testing functions not necessary for the actual package
quickplot <- function(obs.dist, null.dist, lab){
  plot(density(null.dist[[1]], na.rm = T),
       main = paste('nsim = ', lab),
       xlim= c(min(c(null.dist[[1]],
                     obs.dist[1]), na.rm = T),
               max(c(null.dist[[1]],
                     obs.dist[1]), na.rm = T)),
       xlab = 'Ancestral Condition',
       ylab = 'Frequency')
  abline(v=obs.dist[1], col = 'red')
  legend(x = 'topright', legend = c('Observed', 'Null'), col  = c('red', 'black'), lwd = 2)
}

ProcessObservedOneMean <- function(observed.anc.cond){
  vals12 <- vector(length = length(observed.anc.cond))
  vals21 <- vector(length = length(observed.anc.cond))
  for(i in 1:length(observed.anc.cond)){
    vals12[i] <- mean(observed.anc.cond[[i]]$'12', na.rm = T)
    vals21[i] <- mean(observed.anc.cond[[i]]$'21', na.rm = T)
  }

  res <- list('12' = vals12,
              '21' = vals21)
  return(res)
}

ProcessNullOneMean <- function(null.anc.cond, iter){
  vals12 <- vector(length = length(null.anc.cond))
  vals21 <- vector(length = length(null.anc.cond))
  for(j in 1:iter){
    cur.sim12 <- cur.sim21 <- c()
    for(i in 1:length(null.anc.cond)){
      cur.sim12[i] <- mean(null.anc.cond[[i]][[j]]$'12', na.rm = T)
      cur.sim21[i] <- mean(null.anc.cond[[i]][[j]]$'21', na.rm = T)
    }
    vals12 <- c(vals12, cur.sim12)
    vals21 <- c(vals21, cur.sim21)
  }
  return(list('12' = vals12, '21' = vals21))
}

ProcessObservedNoMean <- function(observed.anc.cond){
  vals12 <- vector(length = length(observed.anc.cond))
  vals21 <- vector(length = length(observed.anc.cond))
  for(i in 1:length(observed.anc.cond)){
    vals12 <- c(vals12,observed.anc.cond[[i]]$'12')
    vals21 <- c(vals21,observed.anc.cond[[i]]$'21')
  }

  res <- list('12' = vals12,
              '21' = vals21)
  return(res)
}

ProcessNullNoMean <- function(null.anc.cond, iter){
  vals12 <- vector(length = length(null.anc.cond))
  vals21 <- vector(length = length(null.anc.cond))
  for(j in 1:iter){
    cur.sim12 <- cur.sim21 <- c()
    for(i in 1:length(null.anc.cond)){
      cur.sim12 <- c(cur.sim12, null.anc.cond[[i]][[j]]$'12')
      cur.sim21 <- c(cur.sim21, null.anc.cond[[i]][[j]]$'21')
    }
    vals12 <- c(vals12, cur.sim12)
    vals21 <- c(vals21, cur.sim21)
  }
  return(list('12' = vals12, '21' = vals21))
}
