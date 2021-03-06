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
par(mfrow = c(2,2))
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
fig_label('A:',cex = 2.5)

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
plotSimmap(anc.state.dt, lwd = 3, ftype = 'off')
legend(x = 'bottomleft', legend = c('Ancestral','Derived'), col = c('black', 'red'), pch = 15, bty = 'n')
fig_label('B:',cex = 2.5)


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
plotSimmap(anc.state.dt, lwd = 3, ftype = 'off')
nodelabels(pie = pies, piecol = c('blue','green', 'red'),cex = .8)
legend(x = 'bottomleft', legend = c('Ancestral','Producing (Ancestral)','Derived'),
       col = c('blue', 'red','green'), pch = 16, bg="transparent", bty = 'n')
fig_label('C:',cex = 2.5)

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
legend(x = 'topright', legend = c('Producing','All'), col = c('red', 'black'), pch = 15, bty = 'n')
fig_label('D:',cex = 2.5)

##### Fig 2 #####
par(mfrow = c(1,1), mar = c(5,4,4,2) + .1)
load('Data/Fig2Data.RData')



x <- rep(1:10, each=100)
x <- jitter(x, factor=1.5)
y <- vector()
for(i in 1:10){
  y <- c(y, fig2.data[1:100, i])
}
probs <- vector()
for(i in 1:10){
  probs[i] <- paste(as.character(sum(fig2.data[1:100, i] <= .05)),'%')
}
plot(x = x, y = y, xaxt="n",xlab="", ylab= "",
     pch=16,cex=.49, 
     main = 'Unidirectional Evolution', 
     adj = 0, ylim = c(0,.57))
points(x = x[1:100], y = y[1:100], pch = 16,
       cex = .49, col = 'red')
points(x = x[101:1000], y = y[101:1000], pch = 16,
       cex = .49, col = 'blue')
legend(x = 'topright', 
       legend = c('False Positive (No Correlation)',
                  'Power (Correlation)'), 
       col = c('red', 'blue'), pch = 16, bty = 'n')
mtext(probs[1], side=3, at=1, cex=.7, col = 'red')
mtext(probs[2:10], side=3, at=2:10, cex=.7, col = 'blue')
mtext(1:10, side=1, at=1:10, cex=.85)
mtext("Scaling Factor", side=1, line=1)
mtext("p-value", side=2, line=2.2)
abline(h = .05, lty = 2, lwd = .7)

##### Origins Figure ##### ????????

##### Fig 3 #####
load('Data/Fig3pt5Data.RData')
fig3pt5.data <- fig3.data
load('Data/Fig3Data.RData')

# data <- cbind(rep(1:10, each = 100), 
#               as.vector(fig2.data))
# colnames(data) <- c('Scale.Factor','Pval')


x <- rep(seq(from=20, to=200, by=20), each=100)
x <- jitter(x, factor=.9)
y <- vector()
for(i in 1:10){
  y <- c(y, fig3.data[1:100, i])
}
probs <- vector()
for(i in 1:10){
  probs[i] <- paste(as.character(sum(fig3.data[1:100, i] <= .05)),'%')
}
probsfp <- vector()
for(i in 1:10){
  probsfp[i] <- paste(as.character(sum(fig3pt5.data[1:100, i] <= .05)), '%')
}
plot(x = (x+3), y = y, xlab="", ylab= "", xaxt="n", 
     pch=16,cex=.49, xlim=c(10, 210), 
     main = 'Unidirectional Evolution', adj = 0,
     col = 'blue', ylim = c(0,.57))

x <- rep(seq(from=20, to=200, by=20), each=100)
x <- jitter(x, factor=.9)
y <- vector()
for(i in 1:10){
  y <- c(y, fig3pt5.data[1:100, i])
}
points(x = (x-3), y = y, col = 'red', pch = 16, cex = .49)
mtext(probs, 
      side=3, at=seq(from=23, to=203, by=20), 
      cex=.7, col = 'blue')
mtext(probsfp, side = 3, at=seq(from=17, to=197, by=20), 
      cex = .7, col = 'red')
mtext(c(20,40,60,80,100,120,140,160,180,200), side=1, 
      at=c(20,40,60,80,100,120,140,160,180,200), cex=.85)
mtext("Taxa", side=1, line=1)
mtext("p-value", side=2, line=2.2)
abline(h = .05, lty = 2, lwd = .7)
legend(x = 'topright', 
       legend = c('False Positive (No Correlation)',
                  'Power (Correlation)'), 
       col = c('red', 'blue'), pch = 16, bty = 'n')

##### Fig 4 #####
load('Data/Fig4Data.RData')

x <- rep(1:10, each=200)
x <- jitter(x, factor=1.5)
y <- vector()
for(j in 1:10){
  for(i in 1:100){
    y[2 * (100 * (j - 1) + i - 1) + 1] <- as.numeric(gsub(",.*", "",fig4.data[i,j]))
    y[2 * (100 * (j - 1) + i - 1) + 2] <- as.numeric(substr(fig4.data[i,j], (nchar(gsub(",.*", "",fig4.data[i,j])) + 2), 
                                                            nchar(fig4.data[i,j])))
  }
}
probs <- vector()
probs2 <- c()
for(i in 1:10){
  probs[i] <- paste(
    as.character(
      round(100 * sum(y[(200 * (i - 1) + 1):(200 * i)] <= .025, na.rm = T) / sum(!is.na(y[(200 * (i - 1) + 1):(200 * i)])),
            digits = 0)),'%')
  probs2[i] <- paste(
    as.character(
      round(100 * sum(y[(200 * (i - 1) + 1):(200 * i)] <= .05, na.rm = T) / sum(!is.na(y[(200 * (i - 1) + 1):(200 * i)])),
            digits = 0)),'%')
}
plot(x = x, y = y, xaxt="n",xlab="",
     ylab= "", pch=16,cex=.49, 
     main = 'Bidirectional Evolution', adj = 0,
     ylim = c(0,.57))
points(x = x[1:200], y = y[1:200], pch = 16,
       cex = .49, col = 'red')
points(x = x[201:2000], y = y[201:2000], pch = 16,
       cex = .49, col = 'blue')
mtext(probs[1], side=3, at=1, cex=.7, col = 'red')
mtext(probs[2:10], side=3, at=2:10, cex=.7, col = 'blue')
# mtext(probs2, side = 3, at = 1:10, cex = .7, line = .6, col = 'red')
mtext(1:10, side=1, at=1:10, cex=.85)
mtext("Scaling Factor", side=1, line=1)
mtext("p-value", side=2, line=2.2)
abline(h = .025, lty = 2, lwd = .7)
# abline(h = .05, lty = 2, lwd = .7)
legend(x = 'topright', 
       legend = c('False Positive (No Correlation)',
                  'Power (Correlation)'), 
       col = c('red', 'blue'), pch = 16, bty = 'n')

##### Fig 5 #####
load('Data/Fig5Data.RData')
load('MC_BiNtaxa_fp.RData')
# with rate = 3 there is still 14%NA values. Should I go higher??? #
x <- rep(seq(from=20, to=200, by=20), each=200)
x <- jitter(x, factor=.9)
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
probs <- vector()
probs2 <- vector()
for(i in 1:10){
  probs[i] <- paste(
    as.character(
      round(100 * sum(y[(200 * (i - 1) + 1):(200 * i)] <= .025, na.rm = T) / sum(!is.na(y[(200 * (i - 1) + 1):(200 * i)])),
            digits = 0)),'%')
  probs2[i] <- paste(
    as.character(
      round(100 * sum(y[(200 * (i - 1) + 1):(200 * i)] <= .05, na.rm = T) / sum(!is.na(y[(200 * (i - 1) + 1):(200 * i)])),
            digits = 0)),'%')
}
plot(x = (x+3), y = y, xlab="", ylab= "", xaxt="n", 
     pch=16,cex=.49, xlim=c(10, 210),
     main = 'Bidirectional Evolution', adj = 0,
     ylim = c(0,.57), col = 'blue')

x <- rep(seq(from=20, to=200, by=20), each=200)
x <- jitter(x, factor=.9)
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
probsfp <- vector()
for(i in 1:10){
  probsfp[i] <- paste(
    as.character(
      round(100 * sum(y[(200 * (i - 1) + 1):(200 * i)] <= .025, na.rm = T) / sum(!is.na(y[(200 * (i - 1) + 1):(200 * i)])),
            digits = 0)),'%')

}
points(x = (x-3), y = y, col = 'red', pch = 16, cex = .49)
mtext(probs, 
      side=3, at=seq(from=23, to=203, by=20), 
      cex=.7, col = 'blue')
mtext(probsfp, side = 3, at=seq(from=17, to=197, by=20), 
      cex = .7, col = 'red')
# mtext(probs2, 
#       side = 3, at=seq(from=20, to=200, by=20), cex = .7, line = .6)
mtext(c(20,40,60,80,100,120,140,160,180,200), side=1, 
      at=c(20,40,60,80,100,120,140,160,180,200), cex=.85)
mtext("Taxa", side=1, line=1)
mtext("p-value", side=2, line=2.2)
abline(h = .025, lty = 2, lwd = .7)
# abline(h = .05, lty = 2, lwd = .7)
legend(x = 'topright', 
       legend = c('False Positive (No Correlation)',
                  'Power (Correlation)'), 
       col = c('red', 'blue'), pch = 16, bty = 'n')

