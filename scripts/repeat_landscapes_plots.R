#this plots the genomic landscape of a repeat taxon in five caenorhabditis species

library(ggplot2)
library(lemon)
#set up command line arguments
args <- commandArgs(TRUE)

rep_file <- args[1]
rep_file_type <- as.character(args[1]) 

title <- rep_file_type
#get data in there
rep_data <- read.table(rep_file, sep="\t", header=TRUE)
#turn bp into MB
rep_data$MB <- rep_data$BP/1000000
#turn number of bp into percentage
rep_data$perc_N <- ((rep_data$bp_rep/10000)*100)
#get species in phylogenetic order
rep_data$species <- factor(rep_data$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata"))
#make the plot
rep_plot <- ggplot(rep_data, aes(x = MB, y = perc_N)) + geom_point(alpha=0.25, size=1,colour="black")  + facet_rep_grid(species ~ Chr) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=28), axis.text.x=element_text(colour="black", size=22),axis.text.y=element_text(colour="black", size=24),strip.text.x = element_text(size=24),strip.text.y = element_text(size=24), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("Percent repetitive elements") + ggtitle(title) + theme(plot.title = element_text(size=24))

#set up output name
	#remember to change this to the appropriate output directory
pdftitle_1 <- paste("/projects/phillipslab/gavincw/repeats_12-18-18/46_new_taxonomy_num_bp/04_superfamily/07_perc_N_rep_class_pdf/",title, ".pdf", sep = "", collapse = NULL)

#save the plot
ggsave(pdftitle_1, width=18, height=18,units = "in", useDingbats=FALSE)  
print(rep_plot)
dev.off

