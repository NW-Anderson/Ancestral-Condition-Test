fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
library(plotfunctions)
library(R.utils)
library(phytools)
library(diversitree)
library(geiger)
library(doSNOW)
library(foreach)
##### Fig 1 #####
par(mfrow = c(2,2), mar = c(4,4,0,0) + .1)
# trees <- trees(pars = c(3,1),
#                        type = "bd",
#                        n = 1,
#                        max.taxa = 30,
#                        include.extinct = F)[[1]]
# trees$edge.length <- trees$edge.length / max(branching.times(trees))
# cont.trait <- sim.char(trees, 0.2, model = 'BM')
# names(cont.trait) <- trees$tip.label
load('Data/Fig1Tree.RData')
load('Data/Fig1ContTrait.RData')

smp <- contMap(trees,cont.trait, ftype = 'off', legend = F, lims = c(.24,2), plot = F)
n<-length(smp$cols)
## change to grey scale
smp$cols[1:n]<-rainbow(n, end = 4/6)
plot(smp, legend = F,ftype = 'off')
gradientLegend(depth = .03, valRange = c(.24,2), side = 1, pos = .17, color = rainbow(n, end = 4/6))
legend(x = 'bottomleft', legend = '', title = '         Cont Trait Value', bg="transparent", bty = 'n')
fig_label('A',cex = 2.5)

# cont.trait.AC <- anc.ML(trees, cont.trait, model = "BM")
# branch.means <- c()
# branch.names <- c()
# for(j in 1:nrow(trees$edge)){
#   node.o.int <- trees$edge[j,1]
#   if(node.o.int <= 30){
#     one <- cont.trait[node.o.int]
#   }else{
#     one <- cont.trait.AC$ace[names(cont.trait.AC$ace) == as.character(node.o.int)]
#   }
#   node.o.int <- trees$edge[j,2]
#   if(node.o.int <= 30){
#     two <- cont.trait[node.o.int]
#   }else{
#     two <- cont.trait.AC$ace[names(cont.trait.AC$ace) == as.character(node.o.int)]
#   }
#   branch.means <- c(branch.means, mean(one, two))
#   branch.names <- c(branch.names, paste(as.character(trees$edge[j,1]),as.character(trees$edge[j,2])))
# }
# names(branch.means) <- branch.names
# rm(branch.names)
# upper <- summary(branch.means)[[5]]
# lower <- summary(branch.means)[[2]]
# scale.factor <- 50
# alt.tree <- trees
# for(j in 1:length(branch.means)){
#   if(branch.means[j] < lower){alt.tree$edge.length[j] <- alt.tree$edge.length[j] / scale.factor}
#   if(branch.means[j] > upper){alt.tree$edge.length[j] <- alt.tree$edge.length[j] * scale.factor}
# }
# good.sim <- F
# rate <- .2
# while(good.sim == F){
#   disc.trait <- sim.char(phy = trees,
#                          par = matrix(c(-rate, 0, rate, 0), 2),
#                          model = 'discrete',
#                          root = 1)
#   if(4 < sum(disc.trait == min(disc.trait)) &&
#      sum(disc.trait == min(disc.trait)) < 26){
#     good.sim <- T
#   }
# }
# names(disc.trait) <- trees$tip.label
# anc.state.dt <- make.simmap(trees, disc.trait,
#                             model = matrix(c(0,0,1,0), 2),
#                             nsim = 1,
#                             pi = c(1,0),
#                             message = F)

load('Data/Fig1DiscSimMap.RData')
plotSimmap(anc.state.dt, lwd = 4, ftype = 'off')
legend(x = 'bottomleft', legend = c('Ancestral','Derived'), col = c('black', 'red'), pch = 15, bty = 'n')
fig_label('B',cex = 2.5)


pies <- array(dim = c(anc.state.dt$Nnode, 3))
pies[1:4,] <- rep.row(c(1,0,0),4)
pies[5,] <- t(c(0,0,1))
pies[6:7,] <- rep.row(c(0,1,0),2)
pies[8:11,] <- rep.row(c(1,0,0),4)
pies[12,] <- t(c(0,0,1))
pies[13:20,] <- rep.row(c(1,0,0),8)
pies[21,] <- t(c(0,0,1))
pies[22:23,] <- rep.row(c(0,1,0),2)
pies[24:29,] <- rep.row(c(1,0,0),6)
# plot(trees, tip.color = 'transparent', edge.width = 3)
plotSimmap(anc.state.dt, lwd = 4, ftype = 'off')
nodelabels(pie = pies, piecol = c('blue','green', 'red'),cex = .8)
legend(x = 'bottomleft', legend = c('Ancestral','Producing (Ancestral)','Derived'),
       col = c('blue', 'red','green'), pch = 16, bg="transparent", bty = 'n')
fig_label('C',cex = 2.5)

ss_nodes <- anc.state.dt$mapped.edge[, 1] > 0 &
  anc.state.dt$mapped.edge[, 2] > 0
wanted_branches <- ss_nodes[ss_nodes == T]
wanted_nodes <- names(wanted_branches)
wanted_nodes <- gsub(",.*", "", wanted_nodes)
producing.nodes <- unique(wanted_nodes)
anc.states <- anc.ML(trees, cont.trait, model = "BM")
orig.val <- mean(anc.states$ace[names(anc.states$ace) %in% producing.nodes])
null.orig.val <- vector(length = 1000)
number.of.trans <- length(producing.nodes)
anc.dt <- anc.state.dt
anc.ct <- anc.states
node.states <- describe.simmap(anc.dt)$states
anc.cond.nodes <- anc.ct$ace[names(anc.ct$ace) %in% names(node.states)[node.states != '2']]
for (j in 1:1000){
  # set.seed(j)
  null.orig.val[j] <- mean(sample(anc.cond.nodes,
                                  length(producing.nodes)))
}
par(mar = c(4,4,0,0) + .1)
plot(density(null.orig.val, bw = .025), ylab = 'Frequency', xlab = 'Mean Cont Trait', main = '')
abline(v = orig.val, col = 'red')
legend(x = 'topright', legend = c('Producing','Ancestral'), col = c('red', 'black'), pch = 15, bty = 'n')
fig_label('D',cex = 2.5)

##### Fig 2 #####
par(mfrow = c(1,1), mar = c(5,4,4,2) + .1)
load('Data/Fig2UniData.RData')


probs <- vector()
for(i in 1:10){
  probs[i] <- sum(fig2.data[1:100, i] <= .05) / 100
}
x <- rep(1:10)
load('Data/Fig2BiData.RData')


y <- vector()
for(j in 1:10){
  for(i in 1:100){
    y[2 * (100 * (j - 1) + i - 1) + 1] <- as.numeric(gsub(",.*", "",fig4.data[i,j]))
    y[2 * (100 * (j - 1) + i - 1) + 2] <- as.numeric(substr(fig4.data[i,j], (nchar(gsub(",.*", "",fig4.data[i,j])) + 2), 
                                                            nchar(fig4.data[i,j])))
  }
}
probs2 <- vector()
for(i in 1:10){
  probs2[i] <- round(sum(y[(200 * (i - 1) + 1):(200 * i)] <= .025, na.rm = T) / 
                      sum(!is.na(y[(200 * (i - 1) + 1):(200 * i)])),digits = 2)
}

plot(x, probs, type = 'b', xaxt="n",xlab="", ylab= "",
     pch='o',cex=1.5, 
     main = 'Strength of Correlation vs Power and False Positive', 
     adj = 0, ylim = c(0,.8), col = 'blue', lwd = 4)
points(x = 1, probs[1], pch = 'o', cex = 1.5, col = 'red')
lines(x, probs2, type = 'b', pch = '+', cex = 1.5, col = 'blue', lwd = 4, lty = 2)
points(x = 1, probs2[1], pch = '+', cex = 1.5, col = 'red')
mtext(1:10, side=1, at=1:10, cex=.85)
mtext("Scaling Factor", side=1, line=1)
mtext("Percent Significant", side=2, line=2.2)
abline(h = .05, lty = 2)
legend(x = 'topleft', 
       legend = c('Unidirectional Power','Bidirectional Power',
                  'Unidirectional False Positive', 'Bidirectional False Positive'), 
       col = c('blue', 'blue','red','red'), pch = c(NA,NA,'o','+'), bty = 'n',
       lwd = 2, lty = c(1,2,NA,NA))


##### Fig 3 #####

load('Data/Fig3UniFPData.RData')
fig3pt5.data <- fig3.data
load('Data/Fig3UniPowerData.RData')

x <- seq(from=20, to=200, by=20)
y <- vector()
for(i in 1:10){
  y <- c(y, fig3.data[1:100, i])
}
probs <- vector()
for(i in 1:10){
  probs[i] <- sum(fig3.data[1:100, i] <= .05)
}
probsfp <- vector()
for(i in 1:10){
  probsfp[i] <- sum(fig3pt5.data[1:100, i] <= .05)
}

load('Data/Fig3BiPowerData.RData')
load('Data/Fig3BiFPData.RData')

y <- vector()
for(i in 1:10){
  for(j in 1:10){
    for(i in 1:100){
      y[2 * (100 * (j - 1) + i - 1) + 1] <- as.numeric(gsub(",.*", "",fig5.data[i,j]))
      y[2 * (100 * (j - 1) + i - 1) + 2] <- as.numeric(substr(fig5.data[i,j], (nchar(gsub(",.*", "",fig5.data[i,j])) + 2), 
                                                              nchar(fig5.data[i,j])))
    }
  }
}
biprobs <- vector()
for(i in 1:10){
  biprobs[i] <- 
      round(100 * sum(y[(200 * (i - 1) + 1):(200 * i)] <= .025, na.rm = T) / sum(!is.na(y[(200 * (i - 1) + 1):(200 * i)])),
            digits = 0)
}

y <- vector()
for(i in 1:10){
  for(j in 1:10){
    for(i in 1:100){
      y[2 * (100 * (j - 1) + i - 1) + 1] <- as.numeric(gsub(",.*", "",fig5pt5.data[i,j]))
      y[2 * (100 * (j - 1) + i - 1) + 2] <- as.numeric(substr(fig5pt5.data[i,j], (nchar(gsub(",.*", "",fig5pt5.data[i,j])) + 2), 
                                                              nchar(fig5pt5.data[i,j])))
    }
  }
}
biprobsfp <- vector()
for(i in 1:10){
  biprobsfp[i] <- round(100 * sum(y[(200 * (i - 1) + 1):(200 * i)] <= .025, na.rm = T) / 
                          sum(!is.na(y[(200 * (i - 1) + 1):(200 * i)])),
            digits = 0)
  
}

plot(x, (probs/100), type = 'b', xaxt="n",xlab="", ylab= "",
     pch='o',cex=1.1, 
     main = 'Taxa Number vs Power and False Positive', 
     adj = 0, ylim = c(0,.8), col = 'blue', lwd = 4)
lines(x, (probsfp/100), type = 'b', pch = 'o', cex = 1.1, col = 'red', lwd = 4, lty = 1)
lines(x, (biprobs/100), type = 'b', pch = '+', cex = 1.5, col = 'blue', lwd = 4, lty = 2)
lines(x, (biprobsfp/100), type = 'b', pch = '+', cex = 1.5, col = 'red', lwd = 4, lty = 2)
abline(h = .05, lty = 2)
mtext(c(20,40,60,80,100,120,140,160,180,200), side=1, 
      at=c(20,40,60,80,100,120,140,160,180,200), cex=.85)
mtext("Taxa", side=1, line=1)
mtext("Percent Significant", side=2, line=2.2)
legend(x = 'topleft', 
       legend = c('Unidirectional Power','Bidirectional Power',
                  'Unidirectional False Positive', 'Bidirectional False Positive'), 
       col = c('blue', 'blue','red','red'), bty = 'n',
       lwd = 2, lty = c(1,2,1,2))

