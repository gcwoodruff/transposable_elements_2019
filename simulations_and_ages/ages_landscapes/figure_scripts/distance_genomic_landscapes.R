library(ggplot2)
library(lemon)

args <- commandArgs(TRUE)

rep_file <- args[1]
rep_file_type <- as.character(args[1]) 

title <- rep_file_type

rep_data <- read.table(rep_file, sep="\t", header=TRUE)

rep_data$MB <- rep_data$BP/1000000


rep_plot <- ggplot(rep_data, aes(x = species, y = kimura_distance)) + geom_point(alpha=0.25, size=1,colour="black")  + facet_rep_grid(species ~ Chr) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=28), axis.text.x=element_text(colour="black", size=22),axis.text.y=element_text(colour="black", size=24),strip.text.x = element_text(size=24),strip.text.y = element_text(size=24), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("Kimura distance") + ggtitle(title) + theme(plot.title = element_text(size=24))

pdftitle_1 <- paste("/projects/phillipslab/gavincw/repeats_12-18-18/66_kimura_distance_genomic_landscapes/01_class/08_genomic_landscape_rep_type_pdf/",title, ".pdf", sep = "", collapse = NULL)

ggsave(pdftitle_1, width=18, height=18,units = "in", useDingbats=FALSE)  
print(rep_plot)
dev.off
