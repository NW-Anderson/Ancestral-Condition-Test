library(ggraptR)
##### Fig 2 #####
load('AncCondFig2Data.RData')
data <- cbind(rep(1:10, each = 100), 
              as.vector(fig2.data))
colnames(data) <- c('Scale.Factor','Pval')

data <- as.data.frame(data)
ggraptR(data)


ggplot(data)

sigres <- c()
for(i in 1:10){
  sigres[i] <- sum(fig2.data[,i] < .05)
}
sigres

##### Fig 3 #####
load('AncCondFig3Data.RData')
data <- cbind(rep(1:10, each = 100), 
              as.vector(fig2.data))
colnames(data) <- c('Scale.Factor','Pval')