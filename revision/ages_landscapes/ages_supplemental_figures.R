
#############
#Revision: Supplemental Figure 21
#############


library(ggplot2)
library(lemon)

#get data in there (line 1494 of age_pipeline.sh)
rep_data <- read.table("kimura_distances_global_genomic_landscape.tsv", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

#re-order species levels
rep_data$species <- factor(rep_data$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata"))

big_plot_data <- rep_data

levels(big_plot_data$species)[levels(big_plot_data$species)=="briggsae"] <- "C. briggsae"
levels(big_plot_data$species)[levels(big_plot_data$species)=="nigoni"] <- "C. nigoni"
levels(big_plot_data$species)[levels(big_plot_data$species)=="elegans"] <- "C. elegans"
levels(big_plot_data$species)[levels(big_plot_data$species)=="inopinata"] <- "C. inopinata"
levels(big_plot_data$species)[levels(big_plot_data$species)=="remanei"] <- "C. remanei"


ggplot(big_plot_data, aes(x = MB, y = kimura_distance)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5, se=FALSE) + facet_rep_grid(species ~ Chr) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("Repetitive element age (Kimura distance)") + theme(plot.title = element_text(size=12))



#############
#Revision: Supplemental Figure 22
#############

rep_data <- read.table("kimura_distances_global_landscape_norm_dist_center.tsv.tsv", sep="\t", header=TRUE)

library(ggplot2)
library(lemon)
library(ggforce)

rep_data$MB <- rep_data$BP/1000000

rep_data$species <- factor(rep_data$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata"))

#define chromosome arms and centers
rep_data$chr_str_type <- ifelse(rep_data$norm_dist_center >= 0.25,"arms", "centers")

rep_data$chr_str_type <- factor(rep_data$chr_str_type, levels = c("centers","arms"))

remanei <- expression(paste(italic("C. remanei")))
briggsae <- expression(paste(italic("C. briggsae")))
nigoni <- expression(paste(italic("C. nigoni")))
inopinata <- expression(paste(italic("C. inopinata")))
elegans <- expression(paste(italic("C. elegans")))

species_labels <- c(briggsae,nigoni,remanei,elegans,inopinata)


ggplot(rep_data, aes(x = species, y = kimura_distance)) + geom_sina(aes(colour=chr_str_type),size=0.15, alpha=0.13,scale="width") + stat_summary(aes(group=chr_str_type),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + scale_colour_brewer(palette="Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=13, angle = 45, hjust = 1), axis.text.y = element_text(colour="black", size=13),legend.text=element_text(colour="black", size=13),legend.title=element_text(colour="black", size=13),axis.title=element_text(size=15), strip.background = element_rect(colour="white", fill="white"), strip.text = element_text(size=13, colour="black"),legend.key=element_blank()) + xlab("Species") + ylab("Kimura distance") + scale_x_discrete(labels = species_labels) + labs(colour = "Chromosome\nregion",size=12) + guides(colour = guide_legend(override.aes = list(size=1,alpha=1))) + scale_y_continuous(limits=c(0,60),breaks=c(0,10,20,30,40,50,60))





#############
#Revision: Supplemental Figure 23
#############

library(ggforce)


rep_data <- read.table("kimura_distances_global_genomic_landscape.tsv", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

#re-order species levels
rep_data$species <- factor(rep_data$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata"))

big_plot_data <- rep_data

levels(big_plot_data$species)[levels(big_plot_data$species)=="briggsae"] <- "C. briggsae"
levels(big_plot_data$species)[levels(big_plot_data$species)=="nigoni"] <- "C. nigoni"
levels(big_plot_data$species)[levels(big_plot_data$species)=="elegans"] <- "C. elegans"
levels(big_plot_data$species)[levels(big_plot_data$species)=="inopinata"] <- "C. inopinata"
levels(big_plot_data$species)[levels(big_plot_data$species)=="remanei"] <- "C. remanei"

ggplot(big_plot_data, aes(x=species,y=kimura_distance)) + geom_sina(colour="black",size=0.25,alpha=0.1) + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="yellow") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12,face="italic"),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12), strip.background = element_blank(), axis.ticks = element_line(colour = "black"),legend.position = "none") + xlab("Species") + ylab("Repetitive element age (Kimura distance)") + scale_y_continuous(limits=c(0,60),breaks=c(0,5,10,15,20,25,30,35,40,45,50,55,60))


#############
#Revision: Supplemental Figure 24
#############


library(ggplot2)
library(lemon)
library(ggforce)


#get data in there (see ages.sh to see how this is generated)
rep_data <- read.table("kimura_distances.tsv", sep="\t", header=TRUE)

#get species factors right
rep_data$species <- factor(rep_data$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata"))


#set aside species

cb <- rep_data[rep_data$species == "briggsae",]
ce <- rep_data[rep_data$species == 'elegans',]
ci <- rep_data[rep_data$species == 'inopinata',]
cn <- rep_data[rep_data$species == 'nigoni',]
cr <- rep_data[rep_data$species == 'remanei',]


#set aside superfamily

cb_superfamily <- cb[cb$rep_taxononmic_rank == 'superfamily',]
ce_superfamily <- ce[ce$rep_taxononmic_rank == 'superfamily',]
ci_superfamily <- ci[ci$rep_taxononmic_rank == 'superfamily',]
cn_superfamily <- cn[cn$rep_taxononmic_rank == 'superfamily',]
cr_superfamily <- cr[cr$rep_taxononmic_rank == 'superfamily',]

#remove unused factor levels

cb_superfamily$rep_type <- droplevels(cb_superfamily$rep_type)
ce_superfamily$rep_type <- droplevels(ce_superfamily$rep_type)
ci_superfamily$rep_type <- droplevels(ci_superfamily$rep_type)
cn_superfamily$rep_type <- droplevels(cn_superfamily$rep_type)
cr_superfamily$rep_type <- droplevels(cr_superfamily$rep_type)


#make inopinata superfamily thing rank-ordered by superfamily abundance


#ok, get some colors going, add the previously-determined data to this df ; see supplemental_tables.xls

ci_superfamily$arm_cen_effect_size <- ifelse(ci_superfamily$rep_type == "hAT",0.419398812841078, 0)

ci_superfamily[ci_superfamily$rep_type == "hAT",]$arm_cen_effect_size <- 0.419398812841078
ci_superfamily[ci_superfamily$rep_type == "Bel-Pao",]$arm_cen_effect_size <- -0.182923071309108
ci_superfamily[ci_superfamily$rep_type == "Tc1-Mariner",]$arm_cen_effect_size <- -0.131119051032375
ci_superfamily[ci_superfamily$rep_type == "Jockey",]$arm_cen_effect_size <- -0.121445137273785
ci_superfamily[ci_superfamily$rep_type == "Helitron",]$arm_cen_effect_size <- 0.161883531213795
ci_superfamily[ci_superfamily$rep_type == "PiggyBac",]$arm_cen_effect_size <- 0.088579826520128
ci_superfamily[ci_superfamily$rep_type == "Mutator",]$arm_cen_effect_size <- 0.114910059036189
ci_superfamily[ci_superfamily$rep_type == "Gypsy",]$arm_cen_effect_size <- -0.052004494014452
ci_superfamily[ci_superfamily$rep_type == "Sola",]$arm_cen_effect_size <- -0.038733559177296
ci_superfamily[ci_superfamily$rep_type == "R2",]$arm_cen_effect_size <- -0.015892649575379
ci_superfamily[ci_superfamily$rep_type == "DIRS",]$arm_cen_effect_size <- -0.022518655439846
ci_superfamily[ci_superfamily$rep_type == "Maverick",]$arm_cen_effect_size <- -0.032115572248369
ci_superfamily[ci_superfamily$rep_type == "Merlin",]$arm_cen_effect_size <- -0.018914330512782
ci_superfamily[ci_superfamily$rep_type == "Penelope",]$arm_cen_effect_size <- 0.016369694606429
ci_superfamily[ci_superfamily$rep_type == "RTE",]$arm_cen_effect_size <- 0.046987988040966
ci_superfamily[ci_superfamily$rep_type == "L1",]$arm_cen_effect_size <- -0.031738341786068
ci_superfamily[ci_superfamily$rep_type == "tRNA",]$arm_cen_effect_size <- 0.02056728608969
ci_superfamily[ci_superfamily$rep_type == "PIF-Harbinger",]$arm_cen_effect_size <- 0.005763991770564
ci_superfamily[ci_superfamily$rep_type == "I",]$arm_cen_effect_size <- 0.003264985757307
ci_superfamily[ci_superfamily$rep_type == "CACTA",]$arm_cen_effect_size <- 0.007496096555879

#perc genome repeat 


ci_superfamily$percent_genome_repeat <- ifelse(ci_superfamily$rep_type == "hAT",2.94881513702492, 0)

ci_superfamily[ci_superfamily$rep_type == "hAT",]$percent_genome_repeat <- 2.94881513702492
ci_superfamily[ci_superfamily$rep_type == "Bel-Pao",]$percent_genome_repeat <- 2.60998040841266
ci_superfamily[ci_superfamily$rep_type == "Tc1-Mariner",]$percent_genome_repeat <- 12.8711552836572
ci_superfamily[ci_superfamily$rep_type == "Jockey",]$percent_genome_repeat <- 1.00351606829511
ci_superfamily[ci_superfamily$rep_type == "Helitron",]$percent_genome_repeat <- 1.04580202747984
ci_superfamily[ci_superfamily$rep_type == "PiggyBac",]$percent_genome_repeat <- 0.093782082504237
ci_superfamily[ci_superfamily$rep_type == "Mutator",]$percent_genome_repeat <- 0.464091651653432
ci_superfamily[ci_superfamily$rep_type == "Gypsy",]$percent_genome_repeat <- 2.5486270061669
ci_superfamily[ci_superfamily$rep_type == "Sola",]$percent_genome_repeat <- 0.162578239988973
ci_superfamily[ci_superfamily$rep_type == "R2",]$percent_genome_repeat <- 0.021324863735346
ci_superfamily[ci_superfamily$rep_type == "DIRS",]$percent_genome_repeat <- 0.174248749523139
ci_superfamily[ci_superfamily$rep_type == "Maverick",]$percent_genome_repeat <- 0.352705615464847
ci_superfamily[ci_superfamily$rep_type == "Merlin",]$percent_genome_repeat <- 0.102525613285042
ci_superfamily[ci_superfamily$rep_type == "Penelope",]$percent_genome_repeat <- 0.001631814033423
ci_superfamily[ci_superfamily$rep_type == "RTE",]$percent_genome_repeat <- 7.6468629123317
ci_superfamily[ci_superfamily$rep_type == "L1",]$percent_genome_repeat <- 0.002566951015327
ci_superfamily[ci_superfamily$rep_type == "tRNA",]$percent_genome_repeat <- 0.002004933689203
ci_superfamily[ci_superfamily$rep_type == "PIF-Harbinger",]$percent_genome_repeat <- 0.111310290093051
ci_superfamily[ci_superfamily$rep_type == "I",]$percent_genome_repeat <- 0.001139931980941
ci_superfamily[ci_superfamily$rep_type == "CACTA",]$percent_genome_repeat <- 0.062005192722147



ggplot(ci_superfamily, aes(x=rep_type,y=kimura_distance)) + geom_sina(aes(colour=arm_cen_effect_size),size=0.25) + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12, angle = 45, hjust = 1),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12), strip.background = element_blank(), axis.ticks = element_line(colour = "black"),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=14)) + xlab("Repeat superfamily") + ylab("Repetitive element age (Kimura distance)") + scale_y_continuous(limits=c(0,60),breaks=c(0,5,10,15,20,25,30,35,40,45,50,55,60)) + scale_colour_gradient2(midpoint=0, low="blue", mid="gray", high="red", name = "Arm-center\ndifference", breaks=c(-0.6,-0.3,0,0.3,0.6), limits=c(-0.6,0.6))


#############
#Revision: Supplemental Figures 25-27
#############


library(ggplot2)

rep_data <- read.table("all_cluster_insertion_age_abundance.tsv", sep="\t", header=TRUE)

rep_data$repeat_subclass <- as.factor(rep_data$repeat_subclass)

#just look at inopinata

rep_inopinata <- rep_data[rep_data$species == "inopinata",]

# change NA's to "Unknown" where appropriate ; ie if RM classification is "Unknown" all other repeat taxa are "Unknown"

#add factor level function from stackoverflow ; https://stackoverflow.com/questions/23316815/add-extra-level-to-factors-in-dataframe

addLevel <- function(x, newlevel=NULL) {
  if(is.factor(x)) {
    if (is.na(match(newlevel, levels(x))))
      return(factor(x, levels=c(levels(x), newlevel)))
  }
  return(x)
}


rep_inopinata$repeat_class <- addLevel(rep_inopinata$repeat_class, "Unknown")

rep_inopinata$repeat_class[rep_inopinata$rm_classification %in% "Unknown"] <- "Unknown"

rep_inopinata$repeat_subclass <- addLevel(rep_inopinata$repeat_subclass, "Unknown")

rep_inopinata$repeat_subclass[rep_inopinata$rm_classification %in% "Unknown"] <- "Unknown"

rep_inopinata$repeat_order <- addLevel(rep_inopinata$repeat_order, "Unknown")

rep_inopinata$repeat_order[rep_inopinata$rm_classification %in% "Unknown"] <- "Unknown"

rep_inopinata$repeat_superfamily <- addLevel(rep_inopinata$repeat_superfamily, "Unknown")

rep_inopinata$repeat_superfamily[rep_inopinata$rm_classification %in% "Unknown"] <- "Unknown"

rep_inopinata$repeat_family <- addLevel(rep_inopinata$repeat_family, "Unknown")

rep_inopinata$repeat_family[rep_inopinata$rm_classification %in% "Unknown"] <- "Unknown"

#cool we did it

#now, subset superfamily, not including NA

rep_inopinata_no_na_superfamily <- rep_inopinata[!is.na(rep_inopinata$repeat_superfamily), ]

#drop unused factor levels for superfamily

rep_inopinata_no_na_superfamily$repeat_superfamily <- droplevels(rep_inopinata_no_na_superfamily$repeat_superfamily) 

#alright, are we ready for some plotting.... ah but first, what I want is some extra columns ; is Tc1-Mariner, is_Gypsy, is_Bel-Pao, is_RTE

rep_inopinata_no_na_superfamily$is_Tc1.Mariner <- ifelse(rep_inopinata_no_na_superfamily$repeat_superfamily == "Tc1-Mariner","yes", "no")
rep_inopinata_no_na_superfamily$is_Gypsy <- ifelse(rep_inopinata_no_na_superfamily$repeat_superfamily == "Gypsy","yes", "no")
rep_inopinata_no_na_superfamily$is_Bel.Pao <- ifelse(rep_inopinata_no_na_superfamily$repeat_superfamily == "Bel-Pao","yes", "no")
rep_inopinata_no_na_superfamily$is_RTE <- ifelse(rep_inopinata_no_na_superfamily$repeat_superfamily == "RTE","yes", "no")

#ok, now we can plot I think

library(RColorBrewer)

colourCount = length(unique(rep_inopinata_no_na_superfamily$repeat_superfamily))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


#exclude clusters with less than 30 insertions


ino_greater_than_30_insertions <- rep_inopinata_no_na_superfamily[rep_inopinata_no_na_superfamily$insertion_count > 29,]

#Supplemental Figure 25
ggplot(ino_greater_than_30_insertions, aes(x = average, y = insertion_count)) + geom_point(aes(colour=repeat_superfamily),size=1) + geom_smooth(size=1,method="lm",se=FALSE,linetype="dotted") + theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title=element_text(size=15), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + xlab("Mean Kimura distance") + ylab("Number of insertions") + scale_colour_manual(values=getPalette(colourCount)) + scale_y_continuous(limits=c(0,9000),breaks=c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000))


#Supplemental Figure 26
ggplot(ino_greater_than_30_insertions, aes(x = log(average+1), y = log(insertion_count+1))) + geom_point(aes(colour=repeat_superfamily),size=1) + geom_smooth(size=1,method="lm",se=FALSE,linetype="dotted") + theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title=element_text(size=15), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + xlab("ln(Mean Kimura distance+1)") + ylab("ln(Number of insertions+1)") + scale_colour_manual(values=getPalette(colourCount)) 

#Supplemental Figure 27


ino_greater_than_30_insertions$is_Tc1.Mariner <- as.factor(ino_greater_than_30_insertions$is_Tc1.Mariner)
ino_greater_than_30_insertions$is_Gypsy <- as.factor(ino_greater_than_30_insertions$is_Gypsy)
ino_greater_than_30_insertions$is_Bel.Pao <- as.factor(ino_greater_than_30_insertions$is_Bel.Pao)
ino_greater_than_30_insertions$is_RTE <- as.factor(ino_greater_than_30_insertions$is_RTE)


ggplot(ino_greater_than_30_insertions, aes(x = log(average+1), y = log(insertion_count+1))) + geom_point(aes(colour=is_Tc1.Mariner),size=1) + geom_smooth(size=1,method="lm",se=FALSE,linetype="dotted") + theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title=element_text(size=15), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + xlab("ln(Mean Kimura distance+1)") + ylab("ln(Number of insertions+1)") + scale_colour_manual(values=c("grey","red")) 

ggplot(ino_greater_than_30_insertions, aes(x = log(average+1), y = log(insertion_count+1))) + geom_point(aes(colour=is_Gypsy),size=1) + geom_smooth(size=1,method="lm",se=FALSE,linetype="dotted") + theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title=element_text(size=15), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + xlab("ln(Mean Kimura distance+1)") + ylab("ln(Number of insertions+1)") + scale_colour_manual(values=c("grey","red")) 

ggplot(ino_greater_than_30_insertions, aes(x = log(average+1), y = log(insertion_count+1))) + geom_point(aes(colour=is_Bel.Pao),size=1) + geom_smooth(size=1,method="lm",se=FALSE,linetype="dotted") + theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title=element_text(size=15), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + xlab("ln(Mean Kimura distance+1)") + ylab("ln(Number of insertions+1)") + scale_colour_manual(values=c("grey","red")) 
ggplot(ino_greater_than_30_insertions, aes(x = log(average+1), y = log(insertion_count+1))) + geom_point(aes(colour=is_RTE),size=1) + geom_smooth(size=1,method="lm",se=FALSE,linetype="dotted") + theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title=element_text(size=15), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + xlab("ln(Mean Kimura distance+1)") + ylab("ln(Number of insertions+1)") + scale_colour_manual(values=c("grey","red")) 



