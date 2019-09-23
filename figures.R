#figures for Woodruff and Teterina 2019.

#############
#############
#############
#############
#section "Repeat density covaries with chromosomal position in all species but C. inopinata"
#############
#############
#############
#############

#figure 1 -- phylogeny


library(tidyverse)
library(ggplot2)
library(ggtree)
library(ape)

	#this is the bayesian tree from Stevens et al. 2019 ; retrieved from (https://zenodo.org/record/1402254)
tree <- read.tree("PhyloBayes_species_tree.nwk")

tree_1 <- keep.tip(tree,c("CKAMA","CELEG","CSP34","CREMA","CBRIG","CNIGO"))


#root the tree

tree_2 <- root(tree_1, outgroup= "CKAMA")

#new tip labels


tip_labels <- c("C. kamaaina","C. inopinata","C. elegans","C. remanei","C. nigoni","C. briggsae")

tree_2$tip.label <- tip_labels

ggtree(tree_2) + geom_tiplab() + geom_treescale()

	#cleaned it up  in illustrator



#figure 2 -- all the repetitve regions for all species

library(lemon)

#get data in there (see lines 706-1095 of repeats.sh to see how this is generated)
rep_data <- read.table("global_repeat_density_10kb_win_norm_dist_cent.tsv", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000
#bp repetitive to % repetitive
rep_data$perc_N <- (rep_data$bp_rep/10000)*100

#re-order species levels
rep_data$species <- factor(rep_data$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata"))

big_plot_data <- rep_data

levels(big_plot_data$species)[levels(big_plot_data$species)=="briggsae"] <- "C. briggsae"
levels(big_plot_data$species)[levels(big_plot_data$species)=="nigoni"] <- "C. nigoni"
levels(big_plot_data$species)[levels(big_plot_data$species)=="elegans"] <- "C. elegans"
levels(big_plot_data$species)[levels(big_plot_data$species)=="inopinata"] <- "C. inopinata"
levels(big_plot_data$species)[levels(big_plot_data$species)=="remanei"] <- "C. remanei"


ggplot(big_plot_data, aes(x = MB, y = perc_N)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5, se=FALSE) + facet_rep_grid(species ~ Chr) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("Percent repetitive region") + theme(plot.title = element_text(size=12))

#Supplemental figure 1. Sina plot of repetitive content in arms and centers in all species ; same data as figure 2


rep_data <- read.table("global_repeat_density_10kb_win_norm_dist_cent.tsv", sep="\t", header=TRUE)


library(ggplot2)
library(lemon)
library(ggforce)

rep_data$MB <- rep_data$BP/1000000

rep_data$perc_N <- (rep_data$num_rep/10000)*100

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


ggplot(rep_data, aes(x = species, y = perc_N)) + geom_sina(aes(colour=chr_str_type),size=0.15, alpha=0.13,scale="width") + stat_summary(aes(group=chr_str_type),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + scale_colour_brewer(palette="Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=13, angle = 45, hjust = 1), axis.text.y = element_text(colour="black", size=13),legend.text=element_text(colour="black", size=13),legend.title=element_text(colour="black", size=13),axis.title=element_text(size=15), strip.background = element_rect(colour="white", fill="white"), strip.text = element_text(size=13, colour="black"),legend.key=element_blank()) + xlab("Species") + ylab("Percent repetitive region") + scale_x_discrete(labels = species_labels) + labs(colour = "Chromosome\nregion",size=12) + guides(colour = guide_legend(override.aes = list(size=1,alpha=1)))




#supplemental figure 2 -- repeat density by normalized chromosomal position for all species ; same data as figure 2


rep_data <- read.table("global_repeat_density_10kb_win_norm_dist_cent.tsv", sep="\t", header=TRUE)

rep_data$MB <- rep_data$BP/1000000

rep_data$perc_N <- (rep_data$num_rep/10000)*100

rep_data$species <- factor(rep_data$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata"))

levels(rep_data$species)[levels(rep_data$species)=="briggsae"] <- "C. briggsae"
levels(rep_data$species)[levels(rep_data$species)=="nigoni"] <- "C. nigoni"
levels(rep_data$species)[levels(rep_data$species)=="elegans"] <- "C. elegans"
levels(rep_data$species)[levels(rep_data$species)=="inopinata"] <- "C. inopinata"
levels(rep_data$species)[levels(rep_data$species)=="remanei"] <- "C. remanei"



ggplot(rep_data, aes(x = norm_dist_center, y = perc_N)) + geom_point(alpha=0.05, size=0.5,colour="black") + geom_smooth(size=0.5) + facet_rep_wrap( ~ species,ncol=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=12), axis.text.x=element_text(colour="black", size=11),axis.text.y=element_text(colour="black", size=11),strip.text.x = element_text(size=11,face = "italic"),strip.text.y = element_text(size=11), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Normalized distance\nfrom chromosome midpoint") + ylab("Percent repetitive region")

#supplemental figure 2 -- linear model coefficients for the relationship between repeat density by normalized chromosomal position for all species


beta_data <- read.table("global_repeat_density_chr_pos_lm_coefficients.tsv", sep="\t", header=TRUE)
beta_data$species <- factor(beta_data$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata"))

levels(beta_data$species)[levels(beta_data$species)=="briggsae"] <- "C. briggsae"
levels(beta_data$species)[levels(beta_data$species)=="nigoni"] <- "C. nigoni"
levels(beta_data$species)[levels(beta_data$species)=="elegans"] <- "C. elegans"
levels(beta_data$species)[levels(beta_data$species)=="inopinata"] <- "C. inopinata"
levels(beta_data$species)[levels(beta_data$species)=="remanei"] <- "C. remanei"


ggplot(beta_data,aes(x=species, y=beta)) + geom_bar(stat="identity",aes(fill=rep_mode)) + geom_errorbar(aes(ymin=beta-error, ymax=beta+error), colour="black", width=.1, position=position_dodge(0.1)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black",size=13,face="italic"), axis.text.y = element_text(colour="black", size=15),legend.text=element_text(colour="black", size=15),legend.title=element_text(colour="black", size=16),axis.title=element_text(size=20), strip.background = element_rect(colour="white", fill="white")) + xlab("Species") + ylab("Slope") + scale_fill_brewer(name = "Reproductive mode")



#############
#############
#############
#############
#section "Repeat density covaries with chromosomal position in all species but C. inopinata"
#############
#############
#############
#############



#Figure 3 

rep_dat <- read.table("all_repeat_taxa_density_RTE_briggsae_zeros.tsv", sep="\t", header=TRUE)
	#to get an identical figure, add "0"'s for every briggsae genomic window and label it RTE; briggsae does not have this repeat superfamily

# rep_dat <- read.table("all_repeat_taxa_density_RTE_briggsae_zeros.tsv", sep="\t", header=TRUE)

rep_dat$MB <- rep_dat$BP/1000000

rep_dat$perc_N <- (rep_dat$num_bp_rep/10000)*100

rep_dat$species <- factor(rep_dat$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata"))


levels(rep_dat$species)[levels(rep_dat$species)=="briggsae"] <- "C. briggsae"
levels(rep_dat$species)[levels(rep_dat$species)=="nigoni"] <- "C. nigoni"
levels(rep_dat$species)[levels(rep_dat$species)=="elegans"] <- "C. elegans"
levels(rep_dat$species)[levels(rep_dat$species)=="inopinata"] <- "C. inopinata"
levels(rep_dat$species)[levels(rep_dat$species)=="remanei"] <- "C. remanei"


#break up by repeat taxonomic rank

rep_superfamily <- rep_dat[rep_dat$taxonomic_rank == "superfamily",]

cool_superfamilies <- rep_superfamily[rep_superfamily$rep_class == "hAT" | rep_superfamily$rep_class == "Mutator" | rep_superfamily$rep_class == "PiggyBac" | rep_superfamily$rep_class == "Bel-Pao" | rep_superfamily$rep_class == "Tc1-Mariner" | rep_superfamily$rep_class == "RTE", ]

chr_iii <-  cool_superfamilies[cool_superfamilies$Chr == "III" ,]

#remove unused factor levels

chr_iii$rep_class <- droplevels(chr_iii$rep_class)
chr_iii$species <- droplevels(chr_iii$species)
chr_iii$Chr <- droplevels(chr_iii$Chr)

chr_iii$rep_class <- factor(chr_iii$rep_class, levels = c("hAT","Mutator","PiggyBac","Bel-Pao","Tc1-Mariner","RTE"))


ggplot(chr_iii, aes(x = MB, y = perc_N)) + geom_point(alpha=0.25, size=0.5,aes(colour=rep_class)) + facet_rep_grid(species ~ rep_class) + scale_colour_brewer(palette="Dark2") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=11),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12,face="italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black"),legend.position="none") + xlab("Position (MB)") + ylab("Percent repetitive element")


#Figure 4a

#get data in there and prepare for plotting
rep_data <- read.table("global_repeat_density_10kb_win_norm_dist_cent.tsv", sep="\t", header=TRUE)

rep_data$MB <- rep_data$BP/1000000

rep_data$perc_N <- (rep_data$num_rep/10000)*100

rep_data$species <- factor(rep_data$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata"))

ce_and_ci <- rep_data[rep_data$species == "elegans" | rep_data$species == "inopinata",]

levels(ce_and_ci$species)[levels(ce_and_ci$species)=="elegans"] <- "C. elegans"
levels(ce_and_ci$species)[levels(ce_and_ci$species)=="inopinata"] <- "C. inopinata"


ggplot(ce_and_ci, aes(x = norm_dist_center, y = perc_N)) + geom_point(alpha=0.05, size=0.5,colour="black") + geom_smooth(size=0.5) + geom_vline(xintercept=0.25, linetype="dashed", color = "red",size=0.25) + facet_rep_wrap( ~ species,ncol=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=12), axis.text.x=element_text(colour="black", size=11),axis.text.y=element_text(colour="black", size=11),strip.text.x = element_text(size=11, face="italic"),strip.text.y = element_text(size=11), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Normalized distance\nfrom chromosome center") + ylab("Percent repetitive region")


#Figure 4B


#get data in there and prepare for plotting
rep_dat <- read.table("all_repeat_taxa_density.tsv", sep="\t", header=TRUE)

rep_dat$MB <- rep_dat$BP/1000000

rep_dat$perc_N <- (rep_dat$num_bp_rep/10000)*100

rep_dat$species <- factor(rep_dat$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata"))


levels(rep_dat$species)[levels(rep_dat$species)=="briggsae"] <- "C. briggsae"
levels(rep_dat$species)[levels(rep_dat$species)=="nigoni"] <- "C. nigoni"
levels(rep_dat$species)[levels(rep_dat$species)=="elegans"] <- "C. elegans"
levels(rep_dat$species)[levels(rep_dat$species)=="inopinata"] <- "C. inopinata"
levels(rep_dat$species)[levels(rep_dat$species)=="remanei"] <- "C. remanei"


#break up by repeat taxonomic rank

low_ace <- rep_dat[rep_dat$species == "C. inopinata" & rep_dat$rep_class == "Bel-Pao",]
mid_ace <- rep_dat[rep_dat$species == "C. inopinata" & rep_dat$rep_class == "RTE",]
high_ace <- rep_dat[rep_dat$species == "C. elegans" & rep_dat$rep_class == "Helitron",]

ace <- rbind(low_ace,mid_ace,high_ace)

#remove unused factor levels

ace$rep_class <- droplevels(ace$rep_class)
ace$species <- droplevels(ace$species)
ace$Chr <- droplevels(ace$Chr)

ace$rep_class <- factor(ace$rep_class, levels = c("Bel-Pao","RTE","Helitron"))

ace_na_omit <- na.omit(ace)



levels(ace_na_omit$rep_class)[levels(ace_na_omit$rep_class)=="Bel-Pao"] <- "C. inopinata Bel-Pao"
levels(ace_na_omit$rep_class)[levels(ace_na_omit$rep_class)=="RTE"] <- "C. inopinata RTE"
levels(ace_na_omit$rep_class)[levels(ace_na_omit$rep_class)=="Helitron"] <- "C. elegans Helitron"

ace_na_omit$rep_class <- droplevels(ace_na_omit$rep_class)

anno <- data.frame(x=c(0.425,0.425,0.425),y=c(99,93,95),lab=c("Arm-center difference = -0.18","Arm-center difference = 0.047","Arm-center difference = 0.54"),rep_class=c("C. inopinata Bel-Pao","C. inopinata RTE","C. elegans Helitron"))

ggplot(ace_na_omit, aes(x = norm_dist_chr_center, y = perc_N)) + geom_point(size=0.5, alpha=0.5,aes(colour=rep_class)) + geom_smooth(se=FALSE,size=0.5) + facet_rep_wrap( ~ rep_class, ncol=1) + scale_colour_manual(values=c("#A78DD4","#C9B5B0","#FD351B")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=10), axis.text.x=element_text(colour="black", size=10),axis.text.y=element_text(colour="black", size=10),strip.text.x = element_text(size=9),strip.text.y = element_text(size=9), strip.background = element_blank(), axis.ticks = element_line(colour = "black"),legend.position="none") + xlab("Normalized distance\nfrom chromosome midpoint") + ylab("Percent repetitive element") + geom_text(data = anno, aes(x = x,  y = y, label = lab),size=2) + geom_vline(xintercept=0.25, linetype="dashed", color = "red",size=0.25)


#Figure 4C


stat_dat <- read.table("all_superfamilies_inopinata_labels.tsv", sep="\t", header=TRUE)
	#this  is the same as all_models_superfamily.tsv generated in lines 394-425 in statistics.R but with a "labels" column added to point out notable superfamilies in C. inopinata


stat_dat$species <- factor(stat_dat$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata"))

levels(stat_dat$species)[levels(stat_dat$species)=="briggsae"] <- "C. briggsae"
levels(stat_dat$species)[levels(stat_dat$species)=="nigoni"] <- "C. nigoni"
levels(stat_dat$species)[levels(stat_dat$species)=="elegans"] <- "C. elegans"
levels(stat_dat$species)[levels(stat_dat$species)=="inopinata"] <- "C. inopinata"
levels(stat_dat$species)[levels(stat_dat$species)=="remanei"] <- "C. remanei"

library(ggrepel)

ggplot(stat_dat, aes(x = log(arm_cen_effect_size+1), y = log(percent_genome_repeat+1))) + geom_point(aes(colour=species),size=2) + geom_smooth(aes(colour=species),size=1,method="lm",se=FALSE,linetype="dotted") + geom_text_repel(aes(label = label),color = "black",segment.colour = "black") + theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12,face="italic"),legend.title=element_text(colour="black", size=15),axis.title=element_text(size=15), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + xlab("Arm-center difference (log-transformed)") + ylab("Percentage of genome (log-transformed)") + scale_colour_manual(values=c("#C0C0C0","#808080","#404040","#000000","red")) 

#Supplemental Figure 11 (same thing but not transformed)

ggplot(stat_dat, aes(x = arm_cen_effect_size, y = percent_genome_repeat)) + geom_point(aes(colour=species),size=2) + geom_smooth(aes(colour=species),size=1,method="lm",se=FALSE,linetype="dotted") + geom_text_repel(aes(label = label),color = "black",segment.colour = "black") + theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12,face="italic"),legend.title=element_text(colour="black", size=15),axis.title=element_text(size=15), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + xlab("Arm-center difference") + ylab("Percentage of genome") + scale_colour_manual(values=c("#C0C0C0","#808080","#404040","#000000","red")) 

#Figure 4D
stat_dat <- read.table("all_models_superfamily.tsv", sep="\t", header=TRUE)

ci <- stat_dat[stat_dat$species == "inopinata", ]

ggplot(ci, aes(x = reorder(rep_class,-percent_genome_repeat), y = percent_genome_repeat)) + geom_col(aes(fill=arm_cen_effect_size),size=0.25) + theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=10, angle = 45, hjust =1), axis.text.y = element_text(colour="black", size=11),legend.text=element_text(colour="black", size=10),legend.title=element_text(colour="black", size=11),axis.title=element_text(size=12), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.position = c(0.85, 0.7),legend.title.align=0.5) + xlab("Repeat superfamily") + ylab("Percentage of genome") + scale_fill_gradient2(midpoint=0, low="blue", mid="gray", high="red", name = "Arm-center\ndifference", breaks=c(-0.6,-0.3,0,0.3,0.6), limits=c(-0.6,0.6)) + scale_y_continuous(limits=c(0,12.5), breaks=c(0,2.5,5,7.5,10,12.5))

#Supplemental Figures heat maps and bar plots

#classes

stat_dat <- read.table("all_models_class.tsv", sep="\t", header=TRUE)


stat_dat$species <- factor(stat_dat$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata"))

levels(stat_dat$rep_class)[levels(stat_dat$rep_class)=="Low_complexity"] <- "Low complexity"
levels(stat_dat$rep_class)[levels(stat_dat$rep_class)=="Simple_repeat"] <- "Simple repeat"

remanei <- expression(paste(italic("C. remanei")))
briggsae <- expression(paste(italic("C. briggsae")))
nigoni <- expression(paste(italic("C. nigoni")))
inopinata <- expression(paste(italic("C. inopinata")))
elegans <- expression(paste(italic("C. elegans")))

species_labels <- c(briggsae,nigoni,remanei,elegans,inopinata)

levels(stat_dat$rep_class)[levels(stat_dat$rep_class)=="I"] <- "Class I retrotransposon"
levels(stat_dat$rep_class)[levels(stat_dat$rep_class)=="II"] <- "Class II DNA transposon"


levels(stat_dat$species)[levels(stat_dat$species)=="briggsae"] <- "C. briggsae"
levels(stat_dat$species)[levels(stat_dat$species)=="nigoni"] <- "C. nigoni"
levels(stat_dat$species)[levels(stat_dat$species)=="elegans"] <- "C. elegans"
levels(stat_dat$species)[levels(stat_dat$species)=="inopinata"] <- "C. inopinata"
levels(stat_dat$species)[levels(stat_dat$species)=="remanei"] <- "C. remanei"

#Supplemental Figure 4


stat_dat_b <- stat_dat[stat_dat$rep_class != "ARTEFACT" & stat_dat$rep_class != "rRNA" & stat_dat$rep_class != "snRNA",]

stat_dat_b$rep_class <- droplevels(stat_dat_b$rep_class)

ggplot(stat_dat_b, aes(x = species, y = percent_genome_repeat)) + geom_col(aes(fill=arm_cen_effect_size)) + facet_rep_wrap( ~ rep_class) +  theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=12, angle = 45, hjust =1), axis.text.y = element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title=element_text(size=15), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black")) + xlab("Species") + ylab("Percentage of genome") + scale_fill_gradient2(midpoint=0, low="blue", mid="gray", high="red", name = "Arm-center\ndifference", breaks=c(-1,-0.5,0,0.5,1), limits=c(-1,1)) + scale_x_discrete(labels = species_labels) 


#Supplemental Figure 5
ggplot(stat_dat, aes(x = species, y = rep_class, fill = arm_cen_effect_size)) + geom_tile() + geom_raster() + scale_fill_gradient2(midpoint=0, low="blue", mid="gray", high="red", name = "Arm-center\ndifference", breaks=c(-1,-0.5,0,0.5,1), limits=c(-1,1)) + geom_text(aes(label = paste(ifelse(stat_dat$percent_genome_repeat > 0.01, signif(stat_dat$percent_genome_repeat,2), "<0.01"),"%")), size=3) +  theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=11, angle = 45, hjust =1,face="italic"), axis.text.y = element_text(colour="black", size=11),legend.text=element_text(colour="black", size=11),legend.title=element_text(colour="black", size=11),axis.title=element_text(size=12), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 11, colour = "black")) + xlab("Species") + ylab("Repeat class/category")

#orders

stat_dat <- read.table("all_models_order.tsv", sep="\t", header=TRUE)


stat_dat$species <- factor(stat_dat$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata"))

levels(stat_dat$rep_class)[levels(stat_dat$rep_class)=="Low_complexity"] <- "Low complexity"
levels(stat_dat$rep_class)[levels(stat_dat$rep_class)=="Simple_repeat"] <- "Simple repeat"

remanei <- expression(paste(italic("C. remanei")))
briggsae <- expression(paste(italic("C. briggsae")))
nigoni <- expression(paste(italic("C. nigoni")))
inopinata <- expression(paste(italic("C. inopinata")))
elegans <- expression(paste(italic("C. elegans")))

species_labels <- c(briggsae,nigoni,remanei,elegans,inopinata)

levels(stat_dat$species)[levels(stat_dat$species)=="briggsae"] <- "C. briggsae"
levels(stat_dat$species)[levels(stat_dat$species)=="nigoni"] <- "C. nigoni"
levels(stat_dat$species)[levels(stat_dat$species)=="elegans"] <- "C. elegans"
levels(stat_dat$species)[levels(stat_dat$species)=="inopinata"] <- "C. inopinata"
levels(stat_dat$species)[levels(stat_dat$species)=="remanei"] <- "C. remanei"

stat_dat_omit <- na.omit(stat_dat)

#Supplemental Figure 6

ggplot(stat_dat_omit, aes(x = species, y = percent_genome_repeat)) + geom_col(aes(fill=arm_cen_effect_size)) + facet_rep_wrap( ~ rep_class) +  theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=12, angle = 45, hjust =1), axis.text.y = element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title=element_text(size=15), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black")) + xlab("Species") + ylab("Percentage of genome") + scale_fill_gradient2(midpoint=0, low="blue", mid="gray", high="red", name = "Arm-center\ndifference", breaks=c(-1,-0.5,0,0.5,1), limits=c(-1,1)) + scale_x_discrete(labels = species_labels) 


#Supplemental Figure 7

ggplot(stat_dat_omit, aes(x = species, y = rep_class, fill = arm_cen_effect_size)) + geom_tile() + geom_raster() + scale_fill_gradient2(midpoint=0, low="blue", mid="gray", high="red", name = "Arm-center\ndifference", breaks=c(-0.8,-0.4,0,0.4,0.8), limits=c(-0.8,0.8)) + geom_text(aes(label = paste(ifelse(stat_dat_omit$percent_genome_repeat > 0.01, signif(stat_dat_omit$percent_genome_repeat,2), "<0.01"),"%")), size=3) +  theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=11, angle = 45, hjust =1,face="italic"), axis.text.y = element_text(colour="black", size=11),legend.text=element_text(colour="black", size=11),legend.title=element_text(colour="black", size=11),axis.title=element_text(size=12), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 11, colour = "black")) + xlab("Species") + ylab("Repeat order")

#superfamilies

stat_dat <- read.table("all_models_superfamily.tsv", sep="\t", header=TRUE)


stat_dat$species <- factor(stat_dat$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata"))

levels(stat_dat$rep_class)[levels(stat_dat$rep_class)=="Low_complexity"] <- "Low complexity"
levels(stat_dat$rep_class)[levels(stat_dat$rep_class)=="Simple_repeat"] <- "Simple repeat"

remanei <- expression(paste(italic("C. remanei")))
briggsae <- expression(paste(italic("C. briggsae")))
nigoni <- expression(paste(italic("C. nigoni")))
inopinata <- expression(paste(italic("C. inopinata")))
elegans <- expression(paste(italic("C. elegans")))

species_labels <- c(briggsae,nigoni,remanei,elegans,inopinata)

levels(stat_dat$species)[levels(stat_dat$species)=="briggsae"] <- "C. briggsae"
levels(stat_dat$species)[levels(stat_dat$species)=="nigoni"] <- "C. nigoni"
levels(stat_dat$species)[levels(stat_dat$species)=="elegans"] <- "C. elegans"
levels(stat_dat$species)[levels(stat_dat$species)=="inopinata"] <- "C. inopinata"
levels(stat_dat$species)[levels(stat_dat$species)=="remanei"] <- "C. remanei"

#order by more inclusive taxa (orders, classes)
stat_dat$rep_class <- factor(stat_dat$rep_class, levels = c("Tc1-Mariner","Sola","PiggyBac","PIF-Harbinger","P","Mutator","Merlin","Kolobok","hAT","Dada","CACTA","Maverick","Helitron","Crypton","tRNA","Penelope","Gypsy","ERV","Copia","Bel-Pao","RTE","R2","L1","Jockey","I","DIRS"))
	

stat_dat_omit <- na.omit(stat_dat)

#Supplemental Figure 8

ggplot(stat_dat_omit, aes(x = species, y = percent_genome_repeat)) + geom_col(aes(fill=arm_cen_effect_size)) + facet_rep_wrap( ~ rep_class) +  theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=12, angle = 45, hjust =1), axis.text.y = element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title=element_text(size=15), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black")) + xlab("Species") + ylab("Percentage of genome") + scale_fill_gradient2(midpoint=0, low="blue", mid="gray", high="red", name = "Arm-center\ndifference", breaks=c(-1,-0.5,0,0.5,1), limits=c(-1,1)) + scale_x_discrete(labels = species_labels) 


#Supplemental Figure 9


ggplot(stat_dat_omit, aes(x = species, y = rep_class, fill = arm_cen_effect_size)) + geom_tile() + geom_raster() + scale_fill_gradient2(midpoint=0, low="blue", mid="gray", high="red", name = "Arm-center\ndifference", breaks=c(-0.8,-0.4,0,0.4,0.8), limits=c(-0.8,0.8)) + geom_text(aes(label = paste(ifelse(stat_dat_omit$percent_genome_repeat > 0.01, signif(stat_dat_omit$percent_genome_repeat,2), "<0.01"),"%")), size=3) +  theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=11, angle = 45, hjust =1,face="italic"), axis.text.y = element_text(colour="black", size=11),legend.text=element_text(colour="black", size=11),legend.title=element_text(colour="black", size=11),axis.title=element_text(size=12), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 11, colour = "black")) + xlab("Species") + ylab("Repeat order")
	#colored and labeled by higher taxa in illustrator 



#families

stat_dat <- read.table("all_models_family.tsv", sep="\t", header=TRUE)


stat_dat$species <- factor(stat_dat$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata"))

remanei <- expression(paste(italic("C. remanei")))
briggsae <- expression(paste(italic("C. briggsae")))
nigoni <- expression(paste(italic("C. nigoni")))
inopinata <- expression(paste(italic("C. inopinata")))
elegans <- expression(paste(italic("C. elegans")))

species_labels <- c(briggsae,nigoni,remanei,elegans,inopinata)

levels(stat_dat$species)[levels(stat_dat$species)=="briggsae"] <- "C. briggsae"
levels(stat_dat$species)[levels(stat_dat$species)=="nigoni"] <- "C. nigoni"
levels(stat_dat$species)[levels(stat_dat$species)=="elegans"] <- "C. elegans"
levels(stat_dat$species)[levels(stat_dat$species)=="inopinata"] <- "C. inopinata"
levels(stat_dat$species)[levels(stat_dat$species)=="remanei"] <- "C. remanei"

stat_dat_omit <- na.omit(stat_dat)


#Supplemental Figure 10


ggplot(stat_dat_omit, aes(x = species, y = rep_class, fill = arm_cen_effect_size)) + geom_tile() + geom_raster() + scale_fill_gradient2(midpoint=0, low="blue", mid="gray", high="red", name = "Arm-center\ndifference", breaks=c(-0.8,-0.4,0,0.4,0.8), limits=c(-0.8,0.8)) + geom_text(aes(label = paste(ifelse(stat_dat_omit$percent_genome_repeat > 0.01, signif(stat_dat_omit$percent_genome_repeat,2), "<0.01"),"%")), size=3) +  theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=11, angle = 45, hjust =1,face="italic"), axis.text.y = element_text(colour="black", size=11),legend.text=element_text(colour="black", size=11),legend.title=element_text(colour="black", size=11),axis.title=element_text(size=12), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 11, colour = "black")) + xlab("Species") + ylab("Repeat order")

#Supplemental Figure 12 ; insertion ages untrimmed alignments

library(ggforce)

dat <- read.table("inopinata_branch_length_summaries_untrimmed_alignments.tsv", sep="\t", header=TRUE)
	#same as data "branch_lengths_summaries" line 6166 in repeats.sh


dat$exceptional_superfamily <- ifelse(dat$repeat_superfamily == "Tc1-Mariner" | dat$repeat_superfamily == "RTE" | dat$repeat_superfamily == "Bel-Pao" | dat$repeat_superfamily == "Gypsy","yes", "no")

ggplot(dat, aes(x=exceptional_superfamily,y=mean)) + geom_sina(aes(colour=exceptional_superfamily),size=1) + scale_colour_brewer(palette="Set1") + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12), strip.background = element_blank(), axis.ticks = element_line(colour = "black"),legend.position = "none") + xlab("Exceptional superfamily") + ylab("average branch length") 

#Supplemental Figure 13 ; insertion ages untrimmed alignments

library(RColorBrewer)
colourCount = length(unique(dat$repeat_superfamily))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


ggplot(dat, aes(x=repeat_superfamily,y=mean)) + geom_sina(aes(colour=repeat_superfamily),size=1) + scale_colour_manual(values=getPalette(colourCount)) + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12, angle = 45, hjust = 1),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12), strip.background = element_blank(), axis.ticks = element_line(colour = "black"),legend.position = "none") + xlab("Repeat superfamily") + ylab("average branch length") 

#Supplemental Figure 14 ; insertion ages trimmed alignments


dat <- read.table("inopinata_branch_length_summaries_trimmed_alignments.tsv", sep="\t", header=TRUE)
	#same as data "branch_lengths_summaries" line 7204 in repeats.sh

dat$exceptional_superfamily <- ifelse(dat$repeat_superfamily == "Tc1-Mariner" | dat$repeat_superfamily == "RTE" | dat$repeat_superfamily == "Bel-Pao" | dat$repeat_superfamily == "Gypsy","yes", "no")


ggplot(dat, aes(x=exceptional_superfamily,y=mean)) + geom_sina(aes(colour=exceptional_superfamily),size=1) + scale_colour_brewer(palette="Set1") + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12), strip.background = element_blank(), axis.ticks = element_line(colour = "black"),legend.position = "none") + xlab("Exceptional superfamily") + ylab("average branch length") 

#Supplemental Figure 15 ; insertion ages untrimmed alignments

library(RColorBrewer)
colourCount = length(unique(dat$repeat_superfamily))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


ggplot(dat, aes(x=repeat_superfamily,y=mean)) + geom_sina(aes(colour=repeat_superfamily),size=1) + scale_colour_manual(values=getPalette(colourCount)) + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12, angle = 45, hjust = 1),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12), strip.background = element_blank(), axis.ticks = element_line(colour = "black"),legend.position = "none") + xlab("Repeat superfamily") + ylab("average branch length") 

#Figure 5, remove four superfamiliies from C. inopinata to see change in global repetitive genomic landscape


rep_data <- read.table("global_repeat_density_remove_four_superfamilies.tsv", sep="\t", header=TRUE)

rep_data$MB <- rep_data$BP/1000000

rep_data$perc_N <- (rep_data$num_rep/10000)*100

rep_data$species <- factor(rep_data$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata"))


ce_and_ci <- rep_data[rep_data$species == "elegans" | rep_data$species == "inopinata",]

levels(ce_and_ci$species)[levels(ce_and_ci$species)=="elegans"] <- "C. elegans"
levels(ce_and_ci$species)[levels(ce_and_ci$species)=="inopinata"] <- "C. inopinata"


ggplot(ce_and_ci, aes(x = MB, y = perc_N)) + geom_point(alpha=0.12, size=0.25,colour="black") + geom_smooth(size=0.5,se = FALSE) + facet_rep_grid(species ~ Chr) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=11),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12,face="italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("Percent repetitive region") + theme(plot.title = element_text(size=12))


#############
#############
#############
#############
#section "Gene diversity is negatively correlated with repeat content in all species but C. inopinata"
#############
#############
#############
#############


#Figure 6


dat <- read.table("gene_density.tsv", sep="\t", header=TRUE)
	# (this is the same as file "all_gene_dens_norm_chr_pos.tsv" line 5488 repeats.sh)


dat$species <- factor(dat$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata","inopinata_no_transposon_cds"))


levels(dat$species)[levels(dat$species)=="briggsae"] <- "C. briggsae"
levels(dat$species)[levels(dat$species)=="nigoni"] <- "C. nigoni"
levels(dat$species)[levels(dat$species)=="elegans"] <- "C. elegans"
levels(dat$species)[levels(dat$species)=="inopinata"] <- "C. inopinata"
levels(dat$species)[levels(dat$species)=="remanei"] <- "C. remanei"
levels(dat$species)[levels(dat$species)=="inopinata_no_transposon_cds"] <- "C. inopinata (no transposons)"


ggplot(dat, aes(x = gene_num, y = perc_N)) + geom_point(alpha=0.125, size=0.8,colour="black") + geom_smooth(method ="lm",colour="black",size=0.33, se=FALSE) + facet_rep_wrap( ~ species,ncol=2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12, face="italic"),strip.text.y = element_text(size=12), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Gene number") + ylab("% repetitive elements") + scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70), limits = c(0,70)) + scale_y_continuous(breaks = c(0,25,50,75,100), limits = c(0,100))

#Figure 7


rep_dat <- read.table("gene_density.tsv", sep="\t", header=TRUE)

rep_dat$MB <- rep_dat$BP/1000000

rep_dat$species <- factor(rep_dat$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata", "inopinata_no_transposon_cds"))

levels(rep_dat$species)[levels(rep_dat$species)=="briggsae"] <- "C. briggsae"
levels(rep_dat$species)[levels(rep_dat$species)=="nigoni"] <- "C. nigoni"
levels(rep_dat$species)[levels(rep_dat$species)=="elegans"] <- "C. elegans"
levels(rep_dat$species)[levels(rep_dat$species)=="inopinata"] <- "C. inopinata"
levels(rep_dat$species)[levels(rep_dat$species)=="remanei"] <- "C. remanei"
levels(rep_dat$species)[levels(rep_dat$species)=="inopinata_no_transposon_cds"] <- "C. inopinata (no\nrepetitive genes)"


library(reshape2)
#get long data to plot gene and repeat genomic landscapes simultaneously
rep_dat_melt <- melt(rep_dat, measure.vars=3:4)

names(rep_dat_melt)[names(rep_dat_melt) == 'variable'] <- 'genome_feature'
names(rep_dat_melt)[names(rep_dat_melt) == 'value'] <- 'feature_density'


levels(rep_dat_melt$genome_feature)[levels(rep_dat_melt$genome_feature)=="gene_num"] <- "genes"
levels(rep_dat_melt$genome_feature)[levels(rep_dat_melt$genome_feature)=="perc_N"] <- "repeats"


ggplot(rep_dat_melt, aes(x = MB, y = feature_density)) + geom_point(aes(colour=genome_feature),alpha=0.12, size=0.25) + stat_smooth(aes(colour=genome_feature),size=0.5) + facet_rep_grid(species ~ Chr) + scale_colour_brewer(palette="Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("Feature density\n(Gene count or percent repeat)") + theme(plot.title = element_text(size=12))+labs(colour = "Genome\nfeature",size=12)


#Supplemental Figure 16

rep_dat <- read.table("gene_density.tsv", sep="\t", header=TRUE)

rep_dat$MB <- rep_dat$BP/1000000

rep_dat$species <- factor(rep_dat$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata", "inopinata_no_transposon_cds"))

levels(rep_dat$species)[levels(rep_dat$species)=="briggsae"] <- "C. briggsae"
levels(rep_dat$species)[levels(rep_dat$species)=="nigoni"] <- "C. nigoni"
levels(rep_dat$species)[levels(rep_dat$species)=="elegans"] <- "C. elegans"
levels(rep_dat$species)[levels(rep_dat$species)=="inopinata"] <- "C. inopinata"
levels(rep_dat$species)[levels(rep_dat$species)=="remanei"] <- "C. remanei"
levels(rep_dat$species)[levels(rep_dat$species)=="inopinata_no_transposon_cds"] <- "C. inopinata (no\nrepetitive genes)"

ggplot(rep_dat, aes(x = MB, y = gene_num)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5, se=FALSE) + facet_rep_grid(species ~ Chr) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("Gene number") + theme(plot.title = element_text(size=12))


#Supplemental Figure 17


rep_dat <- read.table("gene_density.tsv", sep="\t", header=TRUE)

rep_dat$chr_str_type <- ifelse(rep_dat$norm_dist_center >= 0.25,"arms", "centers")

rep_dat$MB <- rep_dat$BP/1000000

rep_dat$species <- factor(rep_dat$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata", "inopinata_no_transposon_cds"))

levels(rep_dat$species)[levels(rep_dat$species)=="briggsae"] <- "C. briggsae"
levels(rep_dat$species)[levels(rep_dat$species)=="nigoni"] <- "C. nigoni"
levels(rep_dat$species)[levels(rep_dat$species)=="elegans"] <- "C. elegans"
levels(rep_dat$species)[levels(rep_dat$species)=="inopinata"] <- "C. inopinata"
levels(rep_dat$species)[levels(rep_dat$species)=="remanei"] <- "C. remanei"
levels(rep_dat$species)[levels(rep_dat$species)=="inopinata_no_transposon_cds"] <- "C. inopinata (no\nrepetitive genes)"


ggplot(rep_dat, aes(x=species,fill=,chr_str_type,y=gene_num)) + geom_sina(aes(colour=chr_str_type),size=0.5,alpha=0.2,) + stat_summary(aes(group=chr_str_type),position = position_dodge(width = 0.9),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black") + scale_colour_brewer(palette="Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12, face = "italic"),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12), strip.background = element_blank(), axis.ticks = element_line(colour = "black"),legend.key=element_blank(),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=13)) + xlab("Species") + ylab("Gene number") +labs(colour = "Chromosome\nregion",size=12) + guides(colour = guide_legend(override.aes = list(size=1,alpha=1)))

#Supplemental Figure 18

dat <- read.table("proteins_hit_transposons.tsv", sep="\t", header=T)
	#see repeats.sh line 4687 for where these data come from

dat$species <- factor(dat$species, levels = c("latens","remanei","briggsae","nigoni","sp26","sp40","sinica","tropicalis","wallacei","doughertyi","brenneri","inopinata","elegans","kamaaina","sp28","sp39","sp29","sp32","afra","japonica","sp31","castelli","angaria","sp38","virilis","plicata","sp21","monodelphis", "D_coronatus"))

latens <- expression(paste(italic("C. latens")))
remanei <- expression(paste(italic("C. remanei")))
briggsae <- expression(paste(italic("C. briggsae")))
nigoni <- expression(paste(italic("C. nigoni")))
sp26 <- expression(paste(italic("C. zanzibari")))
sp40 <- expression(paste(italic("C. tribulationis")))
sinica <- expression(paste(italic("C. sinica")))
tropicalis <- expression(paste(italic("C. tropicalis")))
wallacei <- expression(paste(italic("C. wallacei")))
doughertyi <- expression(paste(italic("C. doughertyi")))
brenneri <- expression(paste(italic("C. brenneri")))
inopinata <- expression(paste(italic("C. inopinata")))
elegans <- expression(paste(italic("C. elegans")))
kamaaina <- expression(paste(italic("C. kamaaina")))
sp28 <- expression(paste(italic("C. panamensis")))
sp39 <- expression(paste(italic("C. waitukubuli")))
sp29 <- expression(paste(italic("C. becei")))
sp32 <- expression(paste(italic("C. sulstoni")))
afra <- expression(paste(italic("C. afra")))
japonica <- expression(paste(italic("C. japonica")))
sp31 <- expression(paste(italic("C. uteleia")))
castelli <- expression(paste(italic("C. castelli")))
angaria <- expression(paste(italic("C. angaria")))
sp38 <- expression(paste(italic("C. quiockensis")))
virilis <- expression(paste(italic("C. virilis")))
plicata <- expression(paste(italic("C. plicata")))
sp21 <- expression(paste(italic("C. parvicauda")))
monodelphis <- expression(paste(italic("C. monodelphis")))
coronatus <- expression(paste(italic("D. coronatus")))


species_labels <- c(latens,remanei,briggsae,nigoni,sp26,sp40,sinica,tropicalis,wallacei,doughertyi,brenneri,inopinata,elegans,kamaaina,sp28,sp39,sp29,sp32,afra,japonica,sp31,castelli,angaria,sp38,virilis,plicata,sp21,monodelphis, coronatus)


#% protein coding genes align transposon

ggplot(dat,aes(x=species, y=percent_prot_genes_transposon)) + geom_bar(stat="identity",aes(fill=rep_mode)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", angle = 45, hjust =1, size=10), axis.text.y = element_text(colour="black", size=11),legend.position="none",axis.title=element_text(size=14), strip.background = element_rect(colour="white", fill="white")) + xlab("Species") + ylab("% protiein-coding genes aligning\nto transposon proteins") + scale_fill_brewer(name = "reproductive mode") + scale_x_discrete(labels = species_labels)

ggsave("/home/gavin/genome/genome/repeats_12-18-18/perc_prot_cod_genes_align_transposon_3-18-19.pdf", width=7.5, height=4.875)
print(a)
dev.off


dat_melt <- melt(dat,measure.vars=c(4:5))


labels_1 <- c(total_protein_coding_genes = "Total Number of Protein-coding Genes", non_transposon_protein_coding_genes = "Number of Protein-coding Genes with Transposon-like Proteins Removed")


ggplot(dat_melt,aes(x=species, y=value,fill=variable)) + geom_bar(stat="identity",position="dodge") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", angle = 45, hjust =1, size=10), axis.text.y = element_text(colour="black", size=11),legend.position="none", strip.text.x = element_text(size=12),axis.title=element_text(size=12), strip.background = element_rect(colour="white", fill="white")) + xlab("Species") + ylab("Number of protein-coding genes") + scale_fill_manual(values=c("#99d8c9","#2ca25f")) + scale_x_discrete(labels = species_labels)
	#legend added in illustrator

#Supplemental Figure 19


eff_data <- read.table("effect_sizes.tsv", sep="\t", header=TRUE)
	# see lines 945-1204 statistics.R for where these data come from
eff_data$species <- factor(eff_data$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata", "inopinata_no_transposon_cds"))

levels(eff_data$species)[levels(eff_data$species)=="briggsae"] <- "C. briggsae"
levels(eff_data$species)[levels(eff_data$species)=="nigoni"] <- "C. nigoni"
levels(eff_data$species)[levels(eff_data$species)=="elegans"] <- "C. elegans"
levels(eff_data$species)[levels(eff_data$species)=="inopinata"] <- "C. inopinata"
levels(eff_data$species)[levels(eff_data$species)=="remanei"] <- "C. remanei"
levels(eff_data$species)[levels(eff_data$species)=="inopinata_no_transposon_cds"] <- "C. inopinata (no\nrepetitive genes)"

ggplot(eff_data,aes(x=species, y=feature_eff_size,fill=genome_feature)) + geom_col(position = "dodge")+ geom_errorbar(aes(ymin=feature_eff_size_lower, ymax=feature_eff_size_upper), colour="black", width=.1,position=position_dodge(0.9)) + scale_fill_brewer(palette="Set1",name = "Genomic\nfeature") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black",size=14,face="italic",angle = 45, hjust = 1), axis.text.y = element_text(colour="black", size=15),legend.text=element_text(colour="black", size=15),legend.title=element_text(colour="black", size=16),axis.title=element_text(size=20), strip.background = element_rect(colour="white", fill="white")) + xlab("Species") + ylab("Arm-center difference\n(Cohen's d effect size)") + scale_y_continuous(breaks = c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2), limits = c(-2,2))


#Supplemental Figure 20


eff_data <- read.table("effect_sizes_b.tsv", sep="\t", header=TRUE)
	# see lines 945-1204 statistics.R for where this data come from
eff_data$species <- factor(eff_data$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata", "inopinata_no_transposon_cds"))

levels(eff_data$species)[levels(eff_data$species)=="briggsae"] <- "C. briggsae"
levels(eff_data$species)[levels(eff_data$species)=="nigoni"] <- "C. nigoni"
levels(eff_data$species)[levels(eff_data$species)=="elegans"] <- "C. elegans"
levels(eff_data$species)[levels(eff_data$species)=="inopinata"] <- "C. inopinata"
levels(eff_data$species)[levels(eff_data$species)=="remanei"] <- "C. remanei"
levels(eff_data$species)[levels(eff_data$species)=="inopinata_no_transposon_cds"] <- "C. inopinata (no\nrepetitive genes)"

library(ggrepel)

ggplot(eff_data,aes(x=gene_eff_size, y=repeat_eff_size)) + geom_point(size=1) + geom_smooth(method = "lm", colour="black",size=0.5,linetype="dotted", se = FALSE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black",size=13), axis.text.y = element_text(colour="black", size=15),legend.text=element_text(colour="black", size=15),legend.title=element_text(colour="black", size=16),axis.title=element_text(size=20), strip.background = element_rect(colour="white", fill="white")) + xlab("Genes effect size") + ylab("Repeats effect size") + scale_colour_brewer() + geom_text_repel(aes(label = species), size = 3.5) + scale_y_continuous(breaks = c(-0.5,0,0.5,1,1.5,2), limits = c(-0.5,2)) + scale_x_continuous(breaks = c(-2,-1.5,-1,-0.5,0,0.5), limits = c(-2,0.5)) +geom_vline(xintercept=0) + geom_hline(yintercept=0) 



#Supplemental Figure 21


rep_dat <- read.table("gene_lengths.tsv", sep="\t", header=TRUE)
	#see line 5496 of repeats.sh to see how this was made
rep_dat$MB <- rep_dat$BP/1000000

rep_dat$species <- factor(rep_dat$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata", "inopinata_no_transposons"))

levels(rep_dat$species)[levels(rep_dat$species)=="briggsae"] <- "C. briggsae"
levels(rep_dat$species)[levels(rep_dat$species)=="nigoni"] <- "C. nigoni"
levels(rep_dat$species)[levels(rep_dat$species)=="elegans"] <- "C. elegans"
levels(rep_dat$species)[levels(rep_dat$species)=="inopinata"] <- "C. inopinata"
levels(rep_dat$species)[levels(rep_dat$species)=="remanei"] <- "C. remanei"
levels(rep_dat$species)[levels(rep_dat$species)=="inopinata_no_transposons"] <- "C. inopinata (no\nrepetitive genes)"

rep_dat$gene_length <- as.integer(rep_dat$gene_length)

rep_dat$gene_length_kb <- rep_dat$gene_length/1000

rep_dat$chr_str_type <- ifelse(rep_dat$norm_dist_center >= 0.25,"arms", "centers")


ggplot(rep_dat, aes(x = MB, y = gene_length_kb)) + geom_point(alpha=0.12, size=0.25) + stat_smooth(size=0.5) + facet_rep_grid(species ~ Chr) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("Gene length (KB)") + theme(plot.title = element_text(size=12)) 



#Supplemental Figure 22


rep_dat <- read.table("gene_lengths.tsv", sep="\t", header=TRUE)
	#see line 5496 of repeats.sh to see how this was made
rep_dat$MB <- rep_dat$BP/1000000

rep_dat$species <- factor(rep_dat$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata", "inopinata_no_transposons"))

levels(rep_dat$species)[levels(rep_dat$species)=="briggsae"] <- "C. briggsae"
levels(rep_dat$species)[levels(rep_dat$species)=="nigoni"] <- "C. nigoni"
levels(rep_dat$species)[levels(rep_dat$species)=="elegans"] <- "C. elegans"
levels(rep_dat$species)[levels(rep_dat$species)=="inopinata"] <- "C. inopinata"
levels(rep_dat$species)[levels(rep_dat$species)=="remanei"] <- "C. remanei"
levels(rep_dat$species)[levels(rep_dat$species)=="inopinata_no_transposons"] <- "C. inopinata (no\nrepetitive genes)"

rep_dat$gene_length <- as.integer(rep_dat$gene_length)

rep_dat$gene_length_kb <- rep_dat$gene_length/1000

rep_dat$chr_str_type <- ifelse(rep_dat$norm_dist_center >= 0.25,"arms", "centers")


library(ggforce)

ggplot(rep_dat, aes(x=species,fill=,chr_str_type,y=gene_length_kb)) + geom_sina(aes(colour=chr_str_type),size=0.5,alpha=0.2,) + stat_summary(aes(group=chr_str_type),position = position_dodge(width = 0.9),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black") + scale_colour_brewer(palette="Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12, face = "italic"),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12), strip.background = element_blank(), axis.ticks = element_line(colour = "black"),legend.key=element_blank()) + xlab("Species") + ylab("Average gene length (kb)") +labs(colour = "Chromosome\nregion",size=12) + theme(legend.key = element_blank()) + guides(colour = guide_legend(override.aes = list(size=1,alpha=1)))

#############
#############
#############
#############
#section "Simulations reveal chromosomal heterogeneity in insertion fitness effects is suficient for promoting non-random genomic repetitive organization"
#############
#############
#############
#############

#Figure 8; each plot is one panel of Figure 8



library(ggplot2)



rep_dat <- read.table("simulations_TE_table_for_plots_50_replicates.txt", sep="\t", header=TRUE)

#three recombination domains and three domains of different s (panel 8f right)


sub_dat <- rep_dat[rep_dat$SCENARIO == "SEL_ARM_WEAK_CENTER_STRONG"  & rep_dat$DOM == "HIGH_LOW_HIGH",]

  #without labels

ggplot(sub_dat, aes(x = POS, y = VALUE)) + geom_rect(aes(xmin=0, xmax=1, ymin= -Inf, ymax=Inf),fill="grey85") + geom_rect(aes(xmin=2, xmax=3, ymin= -Inf, ymax=Inf),fill="grey85") + geom_rect(aes(xmin=1, xmax=2, ymin= -Inf, ymax=Inf),fill="grey65") + geom_vline(xintercept=1,colour="black",linetype="dashed",size=1) + geom_vline(xintercept=2,colour="black",linetype="dashed",size=1) + stat_summary(fun.y = mean,geom = 'line',colour = 'red',size = 1.5) + stat_summary(fun.data = 'mean_sdl',geom = 'ribbon',alpha = 0.0002,colour='red') + labs(title = "", x = "Position (MB)", y = "Number of transposable elements") +  theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black",size=1), axis.ticks = element_line(colour = "black",size=1), axis.text.x = element_blank(), axis.text.y = element_blank(),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title = element_blank(), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + ylim(1500,6000) + annotate("text", x=2.25,y=6000, label="Arm-center difference = 2.8",size=2.5)

# no recombination domains and three domains of different s (panel 8f left)


sub_dat <- rep_dat[rep_dat$SCENARIO == "SEL_ARM_WEAK_CENTER_STRONG"  & rep_dat$DOM == "NO_DOMAINS",]

ggplot(sub_dat, aes(x = POS, y = VALUE)) + geom_rect(aes(xmin=0, xmax=1, ymin= -Inf, ymax=Inf),fill="grey85") + geom_rect(aes(xmin=2, xmax=3, ymin= -Inf, ymax=Inf),fill="grey85") + geom_rect(aes(xmin=1, xmax=2, ymin= -Inf, ymax=Inf),fill="grey65") + stat_summary(fun.y = mean,geom = 'line',colour = 'red',size = 1.5) + stat_summary(fun.data = 'mean_sdl',geom = 'ribbon',alpha = 0.0002,colour='red') + labs(title = "", x = "Position (MB)", y = "Number of transposable elements") +  theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black",size=1), axis.ticks = element_line(colour = "black",size=1), axis.text.x = element_blank(), axis.text.y = element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title = element_blank(), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + ylim(1500,6000) + annotate("text", x=2.25,y=6000, label="Arm-center difference = 4.8",size=2.5) 

#selection domains AND recombination domains, but with CUT AND PASTE (DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS) ; (panel 8g right)


sub_dat <- rep_dat[rep_dat$SCENARIO == "SEL_ARM_WEAK_CENTER_STRONG_LOSS"  & rep_dat$DOM == "HIGH_LOW_HIGH",]

ggplot(sub_dat, aes(x = POS, y = VALUE)) + geom_rect(aes(xmin=0, xmax=1, ymin= -Inf, ymax=Inf),fill="grey85") + geom_rect(aes(xmin=2, xmax=3, ymin= -Inf, ymax=Inf),fill="grey85") + geom_rect(aes(xmin=1, xmax=2, ymin= -Inf, ymax=Inf),fill="grey65") + geom_vline(xintercept=1,colour="black",linetype="dashed",size=1) + geom_vline(xintercept=2,colour="black",linetype="dashed",size=1) + stat_summary(fun.y = mean,geom = 'line',colour = 'red',size = 1.5) + stat_summary(fun.data = 'mean_sdl',geom = 'ribbon',alpha = 0.0002,colour='red') + labs(title = "", x = "Position (MB)", y = "Number of transposable elements") +  theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black",size=1), axis.ticks = element_line(colour = "black",size=1), axis.text.x = element_blank(), axis.text.y = element_blank(),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title = element_blank(), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + ylim(100,300) + annotate("text", x=2.25,y=300, label="Arm-center difference = 2.2",size=2.5)

#selection domains but NO recombination domains, and with CUT AND PASTE (DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS) ; (panel 8g left)


sub_dat <- rep_dat[rep_dat$SCENARIO == "SEL_ARM_WEAK_CENTER_STRONG_LOSS"  & rep_dat$DOM == "NO_DOMAINS",]


ggplot(sub_dat, aes(x = POS, y = VALUE)) + geom_rect(aes(xmin=0, xmax=1, ymin= -Inf, ymax=Inf),fill="grey85") + geom_rect(aes(xmin=2, xmax=3, ymin= -Inf, ymax=Inf),fill="grey85") + geom_rect(aes(xmin=1, xmax=2, ymin= -Inf, ymax=Inf),fill="grey65") + stat_summary(fun.y = mean,geom = 'line',colour = 'red',size = 1.5) + stat_summary(fun.data = 'mean_sdl',geom = 'ribbon',alpha = 0.0002,colour='red') + labs(title = "", x = "Position (MB)", y = "Number of transposable elements") +  theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black",size=1), axis.ticks = element_line(colour = "black",size=1), axis.text.x = element_blank(), axis.text.y =  element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title = element_blank(), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + ylim(100,300) + annotate("text", x=2.25,y=300, label="Arm-center difference = 1.4",size=2.5)

#selection domain only in center, strong, three recombination domains (panel 8e right)


sub_dat <- rep_dat[rep_dat$SCENARIO == "SEL_CENTER"  & rep_dat$DOM == "HIGH_LOW_HIGH",]


ggplot(sub_dat, aes(x = POS, y = VALUE)) + geom_rect(aes(xmin=1, xmax=2, ymin= -Inf, ymax=Inf),fill="grey65") + geom_vline(xintercept=1,colour="black",linetype="dashed",size=1) + geom_vline(xintercept=2,colour="black",linetype="dashed",size=1) + stat_summary(fun.y = mean,geom = 'line',colour = 'red',size = 1.5) + stat_summary(fun.data = 'mean_sdl',geom = 'ribbon',alpha = 0.0002,colour='red') + labs(title = "", x = "Position (MB)", y = "Number of transposable elements") +  theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black",size=1), axis.ticks = element_line(colour = "black",size=1), axis.text.x = element_blank(), axis.text.y = element_blank(),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title = element_blank(), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + ylim(3000,9000) + annotate("text", x=2.25,y=9000, label="Arm-center difference = 5.7",size=2.5)

#selection domain center strong only, no recombination domains (panel 8e left)

sub_dat <- rep_dat[rep_dat$SCENARIO == "SEL_CENTER"  & rep_dat$DOM == "NO_DOMAINS",]

ggplot(sub_dat, aes(x = POS, y = VALUE)) + geom_rect(aes(xmin=1, xmax=2, ymin= -Inf, ymax=Inf),fill="grey65") + stat_summary(fun.y = mean,geom = 'line',colour = 'red',size = 1.5) + stat_summary(fun.data = 'mean_sdl',geom = 'ribbon',alpha = 0.0002,colour='red') + labs(title = "", x = "Position (MB)", y = "Number of transposable elements") +  theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black",size=1), axis.ticks = element_line(colour = "black",size=1), axis.text.x = element_blank(), axis.text.y = element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title = element_blank(), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + ylim(3000,9000) + annotate("text", x=2.25,y=9000, label="Arm-center difference = 9.2",size=2.5)


#selection domain only in center, weak, three recombination domains (panel 8d right)

sub_dat <- rep_dat[rep_dat$SCENARIO == "SEL_CENTER_WEAK"  & rep_dat$DOM == "HIGH_LOW_HIGH",]

ggplot(sub_dat, aes(x = POS, y = VALUE)) + geom_rect(aes(xmin=1, xmax=2, ymin= -Inf, ymax=Inf),fill="grey85") + geom_vline(xintercept=1,colour="black",linetype="dashed",size=1) + geom_vline(xintercept=2,colour="black",linetype="dashed",size=1) + stat_summary(fun.y = mean,geom = 'line',colour = 'red',size = 1.5) + stat_summary(fun.data = 'mean_sdl',geom = 'ribbon',alpha = 0.0002,colour='red') + labs(title = "", x = "Position (MB)", y = "Number of transposable elements") +  theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black",size=1), axis.ticks = element_line(colour = "black",size=1), axis.text.x = element_blank(), axis.text.y = element_blank(),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title = element_blank(), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + ylim(10000,27000) + annotate("text", x=2.25,y=27000, label="Arm-center difference = 2.7",size=2.5)


#sel center weak, no rec domains (panel 8d left)


sub_dat <- rep_dat[rep_dat$SCENARIO == "SEL_CENTER_WEAK"  & rep_dat$DOM == "NO_DOMAINS",]

ggplot(sub_dat, aes(x = POS, y = VALUE)) + geom_rect(aes(xmin=1, xmax=2, ymin= -Inf, ymax=Inf),fill="grey85") + stat_summary(fun.y = mean,geom = 'line',colour = 'red',size = 1.5) + stat_summary(fun.data = 'mean_sdl',geom = 'ribbon',alpha = 0.0002,colour='red') + labs(title = "", x = "Position (MB)", y = "Number of transposable elements") +  theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black",size=1), axis.ticks = element_line(colour = "black",size=1), axis.text.x = element_blank(), axis.text.y = element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title = element_blank(), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + ylim(10000,27000) + annotate("text", x=2.25,y=27000, label="Arm-center difference = 5.5",size=2.5)

#selection strong across whole chromosome, three recombination domains (figure 8c right)


sub_dat <- rep_dat[rep_dat$SCENARIO == "SEL_ALL"  & rep_dat$DOM == "HIGH_LOW_HIGH",]

ggplot(sub_dat, aes(x = POS, y = VALUE)) + geom_rect(aes(xmin=0, xmax=3, ymin= -Inf, ymax=Inf),fill="grey65") + geom_vline(xintercept=1,colour="black",linetype="dashed",size=1) + geom_vline(xintercept=2,colour="black",linetype="dashed",size=1) + stat_summary(fun.y = mean,geom = 'line',colour = 'black',size = 1.5) + stat_summary(fun.data = 'mean_sdl',geom = 'ribbon',alpha = 0.0002,colour='black') + labs(title = "", x = "Position (MB)", y = "Number of transposable elements") +  theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black",size=1), axis.ticks = element_line(colour = "black",size=1), axis.text.x = element_blank(), axis.text.y = element_blank(),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title = element_blank(), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + ylim(1100,2400) + annotate("text", x=2.25,y=2400, label="Arm-center difference = -0.24",size=2.5)

#selection strong across whole chromosome, no recombination domains (panel 8c left)

sub_dat <- rep_dat[rep_dat$SCENARIO == "SEL_ALL"  & rep_dat$DOM == "NO_DOMAINS",]


ggplot(sub_dat, aes(x = POS, y = VALUE)) + geom_rect(aes(xmin=0, xmax=3, ymin= -Inf, ymax=Inf),fill="grey65") + stat_summary(fun.y = mean,geom = 'line',colour = 'black',size = 1.5) + stat_summary(fun.data = 'mean_sdl',geom = 'ribbon',alpha = 0.0002,colour='black') + labs(title = "", x = "Position (MB)", y = "Number of transposable elements") +  theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black",size=1), axis.ticks = element_line(colour = "black",size=1), axis.text.x = element_blank(), axis.text.y = element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title = element_blank(), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + ylim(1100,2400) + annotate("text", x=2.25,y=2400, label="Arm-center difference = -0.040",size=2.5)

#selection weak across whole chromosome, three recombination domains (panel 8b right)


sub_dat <- rep_dat[rep_dat$SCENARIO == "SEL_ALL_WEAK"  & rep_dat$DOM == "HIGH_LOW_HIGH",]

  #without labels

ggplot(sub_dat, aes(x = POS, y = VALUE)) + geom_rect(aes(xmin=0, xmax=3, ymin= -Inf, ymax=Inf),fill="grey85") + geom_vline(xintercept=1,colour="black",linetype="dashed",size=1) + geom_vline(xintercept=2,colour="black",linetype="dashed",size=1) + stat_summary(fun.y = mean,geom = 'line',colour = 'black',size = 1.5) + stat_summary(fun.data = 'mean_sdl',geom = 'ribbon',alpha = 0.0002,colour='black') + labs(title = "", x = "Position (MB)", y = "Number of transposable elements") +  theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black",size=1), axis.ticks = element_line(colour = "black",size=1), axis.text.x = element_blank(), axis.text.y = element_blank(),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title = element_blank(), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + ylim(6000,11000) + annotate("text", x=2.25,y=11000, label="Arm-center difference = 0.028",size=2.5)

#selection weak across whole chromosome,  no recombination domains (panel 8b left)

sub_dat <- rep_dat[rep_dat$SCENARIO == "SEL_ALL_WEAK"  & rep_dat$DOM == "NO_DOMAINS",]

ggplot(sub_dat, aes(x = POS, y = VALUE)) + geom_rect(aes(xmin=0, xmax=3, ymin= -Inf, ymax=Inf),fill="grey85") + stat_summary(fun.y = mean,geom = 'line',colour = 'black',size = 1.5) + stat_summary(fun.data = 'mean_sdl',geom = 'ribbon',alpha = 0.0002,colour='black') + labs(title = "", x = "Position (MB)", y = "Number of transposable elements") +  theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black",size=1), axis.ticks = element_line(colour = "black",size=1), axis.text.x = element_blank(), axis.text.y = element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title = element_blank(), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + ylim(6000,11000) + annotate("text", x=2.25,y=11000, label="Arm-center difference = 0.11",size=2.5)

#no selection whole chromosome, three recombination domains (panel 8a right)


sub_dat <- rep_dat[rep_dat$SCENARIO == "NEUTRAL"  & rep_dat$DOM == "HIGH_LOW_HIGH",]

ggplot(sub_dat, aes(x = POS, y = VALUE))  + geom_vline(xintercept=1,colour="black",linetype="dashed",size=1) + geom_vline(xintercept=2,colour="black",linetype="dashed",size=1) + stat_summary(fun.y = mean,geom = 'line',colour = 'black',size = 1.5) + stat_summary(fun.data = 'mean_sdl',geom = 'ribbon',alpha = 0.0002,colour='black') + labs(title = "", x = "Position (MB)", y = "Number of transposable elements") +  theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black",size=1), axis.ticks = element_line(colour = "black",size=1), axis.text.x = element_blank(), axis.text.y = element_blank(),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title = element_blank(), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + ylim(75000,101000) + annotate("text", x=2.25,y=101000, label="Arm-center difference = -0.13",size=2.5)


#no selection whole chromosome,  no recombination domains (panel 8a left)

sub_dat <- rep_dat[rep_dat$SCENARIO == "NEUTRAL"  & rep_dat$DOM == "NO_DOMAINS",]

ggplot(sub_dat, aes(x = POS, y = VALUE)) + stat_summary(fun.y = mean,geom = 'line',colour = 'black',size = 1.5) + stat_summary(fun.data = 'mean_sdl',geom = 'ribbon',alpha = 0.0002,colour='black') + labs(title = "", x = "Position (MB)", y = "Number of transposable elements") +  theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black",size=1), axis.ticks = element_line(colour = "black",size=1), axis.text.x = element_blank(), axis.text.y = element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title = element_blank(), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + ylim(75000,101000) + annotate("text", x=2.5,y=101000, label="Arm-center difference = 0.0038",size=2.25)



#selection weak across whole chromosome , three domains of recombination beneficial mutations (panel 8h right)


sub_dat <- rep_dat[rep_dat$SCENARIO == "SEL_ALL_WEAK_POS_MUT"  & rep_dat$DOM == "HIGH_LOW_HIGH",]
ggplot(sub_dat, aes(x = POS, y = VALUE)) + geom_rect(aes(xmin=0, xmax=3, ymin= -Inf, ymax=Inf),fill="grey85") + geom_vline(xintercept=1,colour="black",linetype="dashed",size=1) + geom_vline(xintercept=2,colour="black",linetype="dashed",size=1) + stat_summary(fun.y = mean,geom = 'line',colour = 'black',size = 1.5) + stat_summary(fun.data = 'mean_sdl',geom = 'ribbon',alpha = 0.0002,colour='black') + labs(title = "", x = "Position (MB)", y = "Number of transposable elements") +  theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title=element_text(size=15), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + ylim(20000,27000) + annotate("text", x=2.5,y=27000, label="Arm-center difference = -0.00062",size=2.5)



##selection weak across whole chromosome ,  no recombination domains , beneficial mutations (panel 8h left)

sub_dat <- rep_dat[rep_dat$SCENARIO == "SEL_ALL_WEAK_POS_MUT"  & rep_dat$DOM == "NO_DOMAINS",]

ggplot(sub_dat, aes(x = POS, y = VALUE)) + geom_rect(aes(xmin=0, xmax=3, ymin= -Inf, ymax=Inf),fill="grey85") + stat_summary(fun.y = mean,geom = 'line',colour = 'black',size = 1.5) + stat_summary(fun.data = 'mean_sdl',geom = 'ribbon',alpha = 0.0002,colour='black') + labs(title = "", x = "Position (MB)", y = "Number of transposable elements") +  theme(panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=15),axis.title=element_text(size=15), strip.background = element_rect(colour="white", fill="white"),strip.text.x = element_text(size = 12, colour = "black"),legend.key=element_blank()) + ylim(20000,27000) + annotate("text", x=2.5,y=27000, label="Arm-center difference = 0.093",size=2.5)
































































































































































































































