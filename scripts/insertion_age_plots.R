#this plots insertion age by chromosome position

library(ggplot2)
library(lemon)

args <- commandArgs(TRUE)

rep_file <- args[1]
rep_file_type <- as.character(args[1]) 

rep_data <- read.table(rep_file, sep="\t", header=TRUE)

title <- paste(rep_file_type,as.character(unique(levels(rep_data$repeat_superfamily))), sep = " ", collapse = NULL)


rep_data$MB <- rep_data$BP_start/1000000

rep_plot <- ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.25, size=0.5,colour="black") + facet_rep_wrap( ~ Chr,nrow=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("Insertion age (substitutions/site)") + ggtitle(title) + theme(plot.title = element_text(size=12))

pdftitle_1 <- paste("/projects/phillipslab/gavincw/repeats_12-18-18/55_inopinata_insertion_ages/15_merge_age_insertion_ch_position_taxon/03_R_plot_pdf/",rep_file_type, ".pdf", sep = "", collapse = NULL)

ggsave(pdftitle_1, width=7, height=4,units = "in", useDingbats=FALSE)  
print(rep_plot)
dev.off
