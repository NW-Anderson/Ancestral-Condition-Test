
##### Fig 3 #####

load('../results/UnidirectionalTaxaFPResults.RData')
load('../results/UnidirectionalTaxaPowerResults.RData')

x <- seq(from=20, to=200, by=20)
y <- vector()
for(i in 1:10){
  y <- c(y, taxa.uni.power.results[1:100, i])
}
probs <- vector()
for(i in 1:10){
  probs[i] <- sum(taxa.uni.power.results[1:100, i] <= .05)
}
probsfp <- vector()
for(i in 1:10){
  probsfp[i] <- sum(taxa.uni.fp.results[1:100, i] <= .05)
}

load('../results/BidirectionalTaxaPowerResults.RData')
load('../results/BidirectionalTaxaFPResults.RData')

## THIS LOOP IS GENERATING WARNINGS CAN YOU COMMENT SO I CAN FIGURE OUT WHATS UP
y <- vector()
for(i in 1:10){
  for(j in 1:10){
    for(i in 1:100){
      y[2 * (100 * (j - 1) + i - 1) + 1] <- as.numeric(gsub(",.*", "",taxa.bi.power.results[i,j]))
      y[2 * (100 * (j - 1) + i - 1) + 2] <- as.numeric(substr(taxa.bi.power.results[i,j], (nchar(gsub(",.*", "",taxa.bi.power.results[i,j])) + 2), 
                                                              nchar(taxa.bi.power.results[i,j])))
    }
  }
}
biprobs <- vector()
for(i in 1:10){
  biprobs[i] <- 
      round(100 * sum(y[(200 * (i - 1) + 1):(200 * i)] <= .025, na.rm = T) / sum(!is.na(y[(200 * (i - 1) + 1):(200 * i)])),
            digits = 0)
}

### SAME FOR THIS LOOP
y <- vector()
for(i in 1:10){
  for(j in 1:10){
    for(i in 1:100){
      y[2 * (100 * (j - 1) + i - 1) + 1] <- as.numeric(gsub(",.*", "",taxa.bi.fp.results[i,j]))
      y[2 * (100 * (j - 1) + i - 1) + 2] <- as.numeric(substr(taxa.bi.fp.results[i,j], (nchar(gsub(",.*", "",taxa.bi.fp.results[i,j])) + 2), 
                                                              nchar(taxa.bi.fp.results[i,j])))
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

