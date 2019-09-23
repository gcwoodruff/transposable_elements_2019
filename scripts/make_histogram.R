###

#this will print some branch length histograms for the clusters

library(ggplot2)

args <- commandArgs(TRUE)

depfile <- args[1]
depfile_chr <- as.character(args[1]) 

dep <- read.table(depfile, sep="\t", header=F)

title <- depfile_chr

plot_1 <- ggplot(dep,aes(x=V6)) + geom_histogram() +geom_vline(xintercept=mean(dep$V6),colour="red") + geom_vline(xintercept=median(dep$V6),colour="blue") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=12,colour = "black"), axis.title=element_text(size=14,colour = "black"),legend.position="none",axis.ticks = element_line(colour = "black")) + ggtitle(title) + labs(x="Length(bp)",y="Frequency")

pdftitle_1 <- paste("/projects/phillipslab/gavincw/repeats_12-18-18/54_cluster_distributions/04_size_distributions/02_make_histogram/briggsae/",title, ".pdf", sep = "", collapse = NULL)
	#make sure to change output directory depending on species
ggsave(pdftitle_1, width=9, height=7,units = "in", useDingbats=FALSE)  
print(plot_1)
dev.off
