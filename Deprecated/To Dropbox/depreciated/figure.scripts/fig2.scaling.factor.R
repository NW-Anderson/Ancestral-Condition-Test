
##### Fig 2 #####
par(mfrow = c(1,1), mar = c(5,4,4,2) + .1)
load('../../results/UnidirectionalScalingAnalysisResults.RData')


probs <- vector()
for(i in 1:10){
  probs[i] <- sum(scaling.uni.results[1:100, i] <= .05) / 100
}
x <- rep(1:10)
load('../../results/BidirectionalScalingAnalysisResults.RData')


y <- vector()
for(j in 1:10){
  for(i in 1:100){
    y[2 * (100 * (j - 1) + i - 1) + 1] <- as.numeric(gsub(",.*", "",scaling.bi.results[i,j]))
    y[2 * (100 * (j - 1) + i - 1) + 2] <- as.numeric(substr(scaling.bi.results[i,j], (nchar(gsub(",.*", "",scaling.bi.results[i,j])) + 2), 
                                                            nchar(scaling.bi.results[i,j])))
  }
}
probs2 <- vector()
for(i in 1:10){
  probs2[i] <- round(sum(y[(200 * (i - 1) + 1):(200 * i)] <= .025, na.rm = T) / 
                       sum(!is.na(y[(200 * (i - 1) + 1):(200 * i)])),digits = 2)
}
library(viridis)
cols <- viridis(10)
plot(x[2:10], probs[2:10], type = 'b', xaxt="n",xlab="", ylab= "",
     pch=16,cex=1, 
     main = 'Strength of Correlation vs Power and False Positive', 
     adj = 0, ylim = c(0,.8), col = cols[1], lwd = 1, xlim=c(.75,10))
lines(x[2:10], probs2[2:10], type = 'b', pch = 16, cex = 1, col = cols[5], lwd = 1, lty = 3)
points(x = 1, probs[1], pch = 'o', cex = 1.5, col = 'red')
points(x = 1, probs2[1], pch = '+', cex = 1.5, col = 'red')
mtext(1:10, side=1, at=1:10, cex=.85)
mtext("Scaling Factor", side=1, line=1)
mtext("Percent Significant", side=2, line=2.2)
abline(h = .05, lty = 2)
legend(x = 'topleft', 
       legend = c('Unidirectional Power','Bidirectional Power',
                  'Unidirectional False Positive', 'Bidirectional False Positive'), 
       col = c(cols[1], cols[5],'red','red'), pch = c(NA,NA,'o','+'), bty = 'n',
       lwd = 1, lty = c(1,3,NA,NA))







results <- matrix(,20,3)
colnames(results) <- 


