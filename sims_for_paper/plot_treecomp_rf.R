library(ggplot2)
library(reshape2)
library(plyr)

args <- commandArgs(trailingOnly = TRUE)
rf = read.delim("./treecomp2.txt",sep = " ", header=FALSE)
colnames(rf)<-c("num_samples","seg_sites","rep","cov","method","rf_dist","ref_edge_in_tree","tree_edge_in_ref","diff_adj","ratio_adj")
rf$num_samples<-factor(rf$num_samples)
rf$seg_sites<-factor(rf$seg_sites)

#rf_tall <- melt(rf, id.vars=c("coverage","samples"))
#p <- ggplot(rf_tall, aes(x=coverage, y=value, colour=variable))+ facet_grid(~ samples) +
#  geom_point(aes(shape=variable)) + geom_line()+ ylim(0, 1) +
#  theme(legend.position="top", legend.title=element_blank(), legend.key = element_blank())

#need to get mean/var
rf_summary<-ddply(rf, c("num_samples", "seg_sites", "cov","method"), summarise,
    mean_rf = mean(rf_dist), sd_rf = sd(rf_dist), sem_rf = sd(rf_dist)/sqrt(length(rf_dist)),
    mean_ra = mean(ratio_adj), sd_ra = sd(ratio_adj), sem_ra = sd(ratio_adj)/sqrt(length(ratio_adj)))

samplelabels <- c("5" = 'sample size = 5', "10" = 'sample size = 10', "20" = 'sample size = 20')
sitelabels <- c("100" = '100 sites', "500" = '500 sites', "1000" = '1000 sites')

p <- ggplot(rf_summary, aes(x=cov, y=mean_rf, colour=method)) +
  facet_grid(seg_sites ~ num_samples, labeller=labeller(num_samples = samplelabels, seg_sites=sitelabels)) +
  geom_point(aes(shape=method)) + geom_line(aes(shape=method)) +  ylim(0, .6) +
  theme(legend.position="top", legend.title=element_blank(), legend.key = element_blank()) +
  xlab("Coverage Depth") + ylab("Mean RF Distance")

ggsave(paste0("./treecomp_rf",args[1],".pdf"), plot = p, width = 6.5, height = 9, dpi = 300)

p2 <- ggplot(rf_summary, aes(x=cov, y=mean_ra, colour=method)) +
  facet_grid(seg_sites ~ num_samples, labeller=labeller(num_samples = samplelabels, seg_sites=sitelabels)) +
  geom_point(aes(shape=method)) + geom_line(aes(shape=method)) +  ylim(0.5, 1) +
  theme(legend.position="top", legend.title=element_blank(), legend.key = element_blank()) +
  xlab("Coverage Depth") + ylab("Mean RA Distance")

ggsave(paste0("./treecomp_ra",args[1],".pdf"), plot = p2, width = 6.5, height = 9, dpi = 300)

#p <- ggplot(subset(rf_summary,rf_summary$seg_sites==100), aes(x=cov, y=mean_rf, colour=method)) +
#  facet_grid(seg_sites ~ num_samples) +
#  geom_point(aes(shape=method)) + geom_line(aes(shape=method)) + ylim(0, .6) +
#  theme(legend.position="top", legend.title=element_blank(), legend.key = element_blank())

#ggsave("./treecomp_rf.pdf", plot = p, width = 6.5, height = 3, dpi = 300)
