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
load('Fig1Tree.RData')
load('Fig1ContTrait.RData')
smp <- contMap(trees,cont.trait, ftype = 'off', legend = F)
fig_label('A:',cex = 2.5)

# good.sim <- F
# rate <- .1
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
load('Fig1DiscSimMap.RData')
plotSimmap(anc.state.dt, lwd = 3, ftype = 'off')
legend(x = 'bottomleft', legend = c('Ancestral','Derived'), col = c('black', 'red'), pch = 15)
fig_label('B:',cex = 2.5)

##### Fig 2 #####
load('AncCondFig2DataPostBlackmon.RData')
# data <- cbind(rep(1:10, each = 100), 
#               as.vector(fig2.data))
# colnames(data) <- c('Scale.Factor','Pval')
# 
# data <- as.data.frame(data)
# ggraptR(data)


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
plot(x = x, y = y, xaxt="n",xlab="", ylab= "", pch=16,cex=.7)
mtext(probs, side=3, at=1:10, cex=.7)
mtext(1:10, side=1, at=1:10, cex=.85)
mtext("Scaling Factor", side=1, line=1)
mtext("p-value", side=2, line=2.2)
abline(h = .05, lty = 2, lwd = .7)

# sigres <- c()
# for(i in 1:10){
#   sigres[i] <- sum(fig2.data[,i] < .05)
# }
# sigres


##### Origins Figure ##### ????????

##### Fig 3 #####
load('AncCondFig3DataPostBlackmon.RData')
# data <- cbind(rep(1:10, each = 100), 
#               as.vector(fig2.data))
# colnames(data) <- c('Scale.Factor','Pval')


x <- rep(seq(from=20, to=200, by=20), each=100)
x <- jitter(x, factor=1.5)
y <- vector()
for(i in 1:10){
  y <- c(y, fig3.data[1:100, i])
}
probs <- vector()
for(i in 1:10){
  probs[i] <- paste(as.character(sum(fig3.data[1:100, i] <= .05)),'%')
}
plot(x = x, y = y, xlab="", ylab= "", xaxt="n", 
     pch=16,cex=.7, xlim=c(10, 210))
mtext(probs, 
      side=3, at=seq(from=20, to=200, by=20), cex=.7)
mtext(c(20,40,60,80,100,120,140,160,180,200), side=1, 
      at=c(20,40,60,80,100,120,140,160,180,200), cex=.85)
mtext("Taxa", side=1, line=1)
mtext("p-value", side=2, line=2.2)
abline(h = .05, lty = 2, lwd = .7)


##### Fig 4 #####
load('AncCondFig4DataPostBlackmon.RData')

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
            digits = 1)),'%')
  probs2[i] <- paste(
    as.character(
      round(100 * sum(y[(200 * (i - 1) + 1):(200 * i)] <= .05, na.rm = T) / sum(!is.na(y[(200 * (i - 1) + 1):(200 * i)])),
            digits = 1)),'%')
}
plot(x = x, y = y, xaxt="n",xlab="", ylab= "", pch=16,cex=.6)
mtext(probs, side=3, at=1:10, cex=.7, col = 'red')
mtext(probs2, side = 3, at = 1:10, cex = .7, line = .6)
mtext(1:10, side=1, at=1:10, cex=.85)
mtext("Scaling Factor", side=1, line=1)
mtext("p-value", side=2, line=2.2)
abline(h = .025, lty = 2, lwd = .7, col = 'red')
abline(h = .05, lty = 2, lwd = .7)


##### Fig 5 #####
load('AncCondFig5DataPostBlackmon.RData')
#### with rate = 3 there is still 14%NA values. Should I go higher??? ####
x <- rep(seq(from=20, to=200, by=20), each=200)
x <- jitter(x, factor=1.5)
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
            digits = 1)),'%')
  probs2[i] <- paste(
    as.character(
      round(100 * sum(y[(200 * (i - 1) + 1):(200 * i)] <= .05, na.rm = T) / sum(!is.na(y[(200 * (i - 1) + 1):(200 * i)])),
            digits = 1)),'%')
}
plot(x = x, y = y, xlab="", ylab= "", xaxt="n", 
     pch=16,cex=.6, xlim=c(10, 210))
mtext(probs, 
      side=3, at=seq(from=20, to=200, by=20), cex=.7, col = 'red')
mtext(probs2, 
      side = 3, at=seq(from=20, to=200, by=20), cex = .7, line = .6)
mtext(c(20,40,60,80,100,120,140,160,180,200), side=1, 
      at=c(20,40,60,80,100,120,140,160,180,200), cex=.85)
mtext("Taxa", side=1, line=1)
mtext("p-value", side=2, line=2.2)
abline(h = .025, lty = 2, lwd = .7, col = 'red')
abline(h = .05, lty = 2, lwd = .7)
