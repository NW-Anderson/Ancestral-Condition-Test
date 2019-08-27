
##### Fig 2 #####
par(mfrow = c(1,1), mar = c(5,4,4,2) + .1)
load('../results/Fig2UniData.RData')


probs <- vector()
for(i in 1:10){
  probs[i] <- sum(fig2.data[1:100, i] <= .05) / 100
}
x <- rep(1:10)
load('../results/Fig2BiData.RData')


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

