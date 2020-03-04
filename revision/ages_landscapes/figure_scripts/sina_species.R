library(ggplot2)
library(lemon)
library(ggforce)

args <- commandArgs(TRUE)

rep_file <- args[1]
rep_file_type <- as.character(args[1]) 

title <- rep_file_type

rep_data <- read.table(rep_file, sep="\t", header=TRUE)

rep_data$species <- factor(rep_data$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata"))

levels(rep_data$species)[levels(rep_data$species)=="briggsae"] <- "C. briggsae"
levels(rep_data$species)[levels(rep_data$species)=="nigoni"] <- "C. nigoni"
levels(rep_data$species)[levels(rep_data$species)=="elegans"] <- "C. elegans"
levels(rep_data$species)[levels(rep_data$species)=="inopinata"] <- "C. inopinata"
levels(rep_data$species)[levels(rep_data$species)=="remanei"] <- "C. remanei"

rep_plot <-  ggplot(rep_data, aes(x = species, y = kimura_distance)) + geom_sina(size=0.15, alpha=0.13,scale="width")  + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + scale_colour_brewer(palette="Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=13, angle = 45, hjust = 1,face="italic"), axis.text.y = element_text(colour="black", size=13),legend.text=element_text(colour="black", size=13),legend.title=element_text(colour="black", size=13),axis.title=element_text(size=15), strip.background = element_rect(colour="white", fill="white"), strip.text = element_text(size=13, colour="black"),legend.key=element_blank()) + xlab("Species") + ylab("Kimura distance") + ggtitle(title) 




pdftitle_1 <- paste("/projects/phillipslab/gavincw/repeats_12-18-18/66_kimura_distance_genomic_landscapes/01_class/08_genomic_landscape_rep_type_pdf/",title, ".pdf", sep = "", collapse = NULL)
ggsave(pdftitle_1, width=9, height=9,units = "in", useDingbats=FALSE)  
print(rep_plot)
dev.off
