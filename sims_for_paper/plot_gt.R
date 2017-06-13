library(ggplot2)
library(reshape2)
library(plyr)

gt = read.table("./gt_comparisons2.txt", header=TRUE, fill=T)
#nsample segsites cov          method precision recall
gt$nsample<-factor(gt$nsample)
gt$segsites<-factor(gt$segsites)

#get mean/var
gt_summary<-ddply(gt, c("nsample", "segsites", "cov","method"), summarise,
    mean_prec = mean(precision), sd_precision = sd(precision), sem_precision = sd(precision)/sqrt(length(precision)),
    mean_recall = mean(recall), sd_recall = sd(recall), sem_recall = sd(recall)/sqrt(length(recall)))
    
samplelabels <- c("5" = 'sample size = 5', "10" = 'sample size = 10', "20" = 'sample size = 20')
sitelabels <- c("100" = '100 sites', "500" = '500 sites', "1000" = '1000 sites')

p <- ggplot(gt_summary, aes(x=cov, y=mean_prec, colour=method)) + 
  facet_grid(segsites ~ nsample, labeller=labeller(nsample = samplelabels, segsites=sitelabels)) +
  geom_point(aes(shape=method)) + geom_line(aes(shape=method)) +  ylim(0, 1) + 
  theme(legend.position="top", legend.title=element_blank(), legend.key = element_blank()) +
  xlab("Coverage Depth") + ylab("Mean Precision")
  
ggsave("./precision.pdf", plot = p, width = 6.5, height = 9, dpi = 300)

p2 <- ggplot(gt_summary, aes(x=cov, y=mean_recall, colour=method)) + 
  facet_grid(segsites ~ nsample, labeller=labeller(nsample = samplelabels, segsites=sitelabels)) +
  geom_point(aes(shape=method)) + geom_line(aes(shape=method)) +  ylim(0, 1) + 
  theme(legend.position="top", legend.title=element_blank(), legend.key = element_blank()) +
  xlab("Coverage Depth") + ylab("Mean Recall")
  
ggsave("./recall.pdf", plot = p2, width = 6.5, height = 9, dpi = 300)