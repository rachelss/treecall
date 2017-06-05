library(ggplot2)
library(reshape2)
library(plyr)

rf = read.delim("./treecomp2.txt",sep = " ", header=FALSE)
colnames(rf)<-c("num_samples","seg_sites","rep","cov","method","rf_dist","V7","V8","V9","V10")
#rf_tall <- melt(rf, id.vars=c("coverage","samples"))
#p <- ggplot(rf_tall, aes(x=coverage, y=value, colour=variable))+ facet_grid(~ samples) +
#  geom_point(aes(shape=variable)) + geom_line()+ ylim(0, 1) +
#  theme(legend.position="top", legend.title=element_blank(), legend.key = element_blank())

#need to get mean/var
rf_summary<-ddply(rf, c("num_samples", "seg_sites", "cov","method"), summarise,
      mean = mean(rf_dist), sd = sd(rf_dist),
      sem = sd(rf_dist)/sqrt(length(rf_dist)))

colnames(rf_summary)[5]<-"mean_rf"

p <- ggplot(rf_summary, aes(x=cov, y=mean_rf, colour=method)) +
  facet_grid(seg_sites ~ num_samples) +
  geom_point(aes(shape=method)) +  ylim(0, .6) +
  theme(legend.position="top", legend.title=element_blank(), legend.key = element_blank())

ggsave("./treecomp_rf.pdf", plot = p, width = 6.5, height = 9, dpi = 300)

p <- ggplot(subset(rf_summary,rf_summary$seg_sites==100), aes(x=cov, y=mean_rf, colour=method)) +
  facet_grid(seg_sites ~ num_samples) +
  geom_point(aes(shape=method)) + geom_line(aes(shape=method)) + ylim(0, .6) +
  theme(legend.position="top", legend.title=element_blank(), legend.key = element_blank())

ggsave("./treecomp_rf.pdf", plot = p, width = 6.5, height = 3, dpi = 300)
