results <- as.data.frame(matrix(,2,3))
row.names(results) <- c("power","false.positive")
colnames(results) <- c("pagel","thresh","ancond")
load("results/PagelThreshPower.RData")
pval.array[,2] <- !pval.array[,2]
results[1,1:3] <- c(colSums(pval.array), 50)

load("results/PagelThreshFP.RData")
results[2,1:3] <- c(colSums(pval.array), 6)

res <- as.data.frame(cbind(unlist(results),
      rep(c("power","false positive"), times=3),
      rep(c("Pagel's","Threshold","Ancestral condition"), each=2)))
colnames(res) <- c("pos","type","test")
res$pos <- as.numeric(as.character(res$pos))
library(ggraptR)
#ggraptR(res)

ggplot(res, aes(y=pos, x=as.factor(test))) + 
  geom_bar(aes(fill=as.factor(test)), stat="identity", position="stack", alpha=0.9) + 
  facet_grid(. ~ type) + 
  guides(fill=guide_legend(title="Method")) + 
  scale_fill_brewer(palette="Set1")+
  xlab("") + 
  ylab("rate")
