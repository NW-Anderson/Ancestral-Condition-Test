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
