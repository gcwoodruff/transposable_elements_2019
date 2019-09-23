#models and hypothesis tests for Woodruff and Teterina 2019.

#############
#############
#############
#############
#section "Repeat density covaries with chromosomal position in all species but C. inopinata"
#############
#############
#############
#############

#get data in there (see lines 706-1095 of repeats.sh to see how this is generated)
rep_data <- read.table("/home/gavin/genome/genome/repeats_12-18-18/preparing_for_deposition_8-9-19/data/global_repeat_denisty_10kb_win_norm_dist_cent.tsv", sep="\t", header=TRUE)

#get readable unit of bp
rep_data$MB <- rep_data$BP/1000000
#turn number of bp repetitive in 10kb window into percentage
rep_data$perc_N <- (rep_data$num_rep/10000)*100
#get species factors right
rep_data$species <- factor(rep_data$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata"))

#define chromosome arms and centers by middle or outer half of chromosome
rep_data$chr_str_type <- ifelse(rep_data$norm_dist_center >= 0.25,"arms", "centers")
#get the factors right
rep_data$chr_str_type <- factor(rep_data$chr_str_type, levels = c("centers","arms"))


#set aside species

cb <- rep_data[rep_data$species == "briggsae",]
ce <- rep_data[rep_data$species == 'elegans',]
ci <- rep_data[rep_data$species == 'inopinata',]
cn <- rep_data[rep_data$species == 'nigoni',]
cr <- rep_data[rep_data$species == 'remanei',]

#mean repeat content in arms, centers, difference in mean repeat content of arms and centers in briggsae

c(mean(cb[cb$chr_str_type == "arms",]$perc_N),mean(cb[cb$chr_str_type == "centers",]$perc_N), mean(cb[cb$chr_str_type == "arms",]$perc_N)-mean(cb[cb$chr_str_type == "centers",]$perc_N))

# [1] 30.688075 21.165629  9.522445

#mean repeat content in arms, centers, difference in mean repeat content of arms and centers in elegans
c(mean(ce[ce$chr_str_type == "arms",]$perc_N),mean(ce[ce$chr_str_type == "centers",]$perc_N), mean(ce[ce$chr_str_type == "arms",]$perc_N)-mean(ce[ce$chr_str_type == "centers",]$perc_N))

# [1] 26.62572 10.70563 15.92008

#mean repeat content in arms, centers, difference in mean repeat content of arms and centers in nigoni

c(mean(cn[cn$chr_str_type == "arms",]$perc_N),mean(cn[cn$chr_str_type == "centers",]$perc_N), mean(cn[cn$chr_str_type == "arms",]$perc_N)-mean(cn[cn$chr_str_type == "centers",]$perc_N))

# [1] 30.668955 22.204559  8.464396


#mean repeat content in arms, centers, difference in mean repeat content of arms and centers in  remanei

c(mean(cr[cr$chr_str_type == "arms",]$perc_N),mean(cr[cr$chr_str_type == "centers",]$perc_N), mean(cr[cr$chr_str_type == "arms",]$perc_N)-mean(cr[cr$chr_str_type == "centers",]$perc_N))

# [1] 32.34117 14.44945 17.89172


#mean repeat content in arms, centers, difference in mean repeat content of arms and centers in  inopinata

c(mean(ci[ci$chr_str_type == "arms",]$perc_N),mean(ci[ci$chr_str_type == "centers",]$perc_N), mean(ci[ci$chr_str_type == "arms",]$perc_N)-mean(ci[ci$chr_str_type == "centers",]$perc_N))

# [1] 30.7642985 29.7667013  0.9975971


#linear models, repeat density ~ chromosome position
#the coefficients were put into file "global_repeat_denisty_chr_pos_lm_coefficients.tsv" for Supplemental figure 3

summary(lm(perc_N ~ norm_dist_center, data=cb))
	#briggsae

#Call:
#lm(formula = perc_N ~ norm_dist_center, data = cb)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-33.693  -9.379  -1.102   8.339  77.674 
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       17.8506     0.2630   67.88   <2e-16 ***
#norm_dist_center  32.2984     0.9107   35.47   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 13.49 on 10519 degrees of freedom
#Multiple R-squared:  0.1068,	Adjusted R-squared:  0.1067 
#F-statistic:  1258 on 1 and 10519 DF,  p-value: < 2.2e-16


summary(lm(perc_N ~ norm_dist_center, data=ce))
	#elegans
#Call:
#lm(formula = perc_N ~ norm_dist_center, data = ce)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-31.763  -8.541  -2.747   5.869  92.852 
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)        5.1101     0.2783   18.36   <2e-16 ***
#norm_dist_center  54.2140     0.9637   56.26   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 13.93 on 10028 degrees of freedom
#Multiple R-squared:  0.2399,	Adjusted R-squared:  0.2398 
#F-statistic:  3165 on 1 and 10028 DF,  p-value: < 2.2e-16



summary(lm(perc_N ~ norm_dist_center, data=cn))
	#nigoni
#Call:
#lm(formula = perc_N ~ norm_dist_center, data = cn)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-34.336 -10.675  -1.911   8.500  81.457 
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       18.5348     0.2887   64.21   <2e-16 ***
#norm_dist_center  31.6031     0.9998   31.61   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 15.67 on 11784 degrees of freedom
#Multiple R-squared:  0.07817,	Adjusted R-squared:  0.07809 
#F-statistic: 999.2 on 1 and 11784 DF,  p-value: < 2.2e-16



summary(lm(perc_N ~ norm_dist_center, data=cr))
	#remanei
#Call:
#lm(formula = perc_N ~ norm_dist_center, data = cr)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-38.482 -12.251  -4.834   8.345  91.616 
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)        8.3053     0.3398   24.44   <2e-16 ***
#norm_dist_center  60.3526     1.1768   51.28   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 18.99 on 12487 degrees of freedom
#Multiple R-squared:  0.174,	Adjusted R-squared:  0.1739 
#F-statistic:  2630 on 1 and 12487 DF,  p-value: < 2.2e-16
#


summary(lm(perc_N ~ norm_dist_center, data=ci))
	#inopinata
#Call:
#lm(formula = perc_N ~ norm_dist_center, data = ci)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-30.637 -16.978  -3.926  12.610  70.030 
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       30.6403     0.3924  78.093   <2e-16 ***
#norm_dist_center  -1.4987     1.3589  -1.103     0.27    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 21.76 on 12300 degrees of freedom
#Multiple R-squared:  9.889e-05,	Adjusted R-squared:  1.76e-05 
#F-statistic: 1.216 on 1 and 12300 DF,  p-value: 0.2701



#############
#############
#############
#############
#section "Divergent geomic repetitive landscapes are driven by diversity in repeat taxon abundance and chromsomal structure"
#############
#############
#############
#############

#packages
library(dplyr)
library(effsize)
library(Rmisc)


#get data in there (see lines 1434-4447 of repeats.sh to see how this is generated)

rep_dat <- read.table("/home/gavin/genome/genome/repeats_12-18-18/preparing_for_deposition_8-9-19/data/all_repeat_taxa_density.tsv", sep="\t", header=TRUE)

#break up by repeat taxonomic rank

rep_class <- rep_dat[rep_dat$taxonomic_rank == "class",]
rep_order <- rep_dat[rep_dat$taxonomic_rank == "order",]
rep_superfamily <- rep_dat[rep_dat$taxonomic_rank == "superfamily",]
rep_family <- rep_dat[rep_dat$taxonomic_rank == "family",]

#classes first

#remove unused factor levels

rep_class$rep_class <- droplevels(rep_class$rep_class)

#replace NA with "NA"

rep_class$rep_class = factor(rep_class$rep_class, levels=c(levels(rep_class$rep_class), "NA"))


rep_class$rep_class[is.na(rep_class$rep_class)] = "NA"

#subset species

cb <- rep_class[rep_class$species == "briggsae",]
ce <- rep_class[rep_class$species == 'elegans',]
ci <- rep_class[rep_class$species == 'inopinata',]
cn <- rep_class[rep_class$species == 'nigoni',]
cr <- rep_class[rep_class$species == 'remanei',]


#set aside repeat classes and get linear model fit statistics ; first, BASE PAIRS of repeats (not repeat counts)


get_models <- function(arg1, arg2, arg3, arg4){
		#make empty variables to fill later
	stat_df <- NULL
	stat_df_order <- NULL
	cf <- NULL
	adj_r_sq <- NULL
	F_stat <- NULL
	sum_repeats <- NULL
	classes_w_models <- NULL
	repeat_total_bp <- NULL
	percent_genome_repeat <- NULL
	
	arm_lower <- NULL
	arm_mean <- NULL
	arm_upper <- NULL
	
	center_lower <- NULL
	center_mean <- NULL
	center_upper <- NULL
	
	arm_cen_effect_size_lower <- NULL
	arm_cen_effect_size <- NULL
	arm_cen_effect_size_upper <- NULL

	#extract statistics regarding repeat density and chromosome position for each repeat taxon (ie, levels in rep_class)
	for (i in levels(arg1$rep_class)){
		dat <- arg1[arg1$rep_class == i,]
		if (nrow(dat) >14){fit <- summary(lm(num_bp_rep ~ norm_dist_chr_center, data=dat)) #excluding repeat taxa that are in <15 genomic windows across the genome
		dat_center <- dat[dat$norm_dist_chr_center < 0.25,] #define chromosome arms and centers
		dat_arm <- dat[dat$norm_dist_chr_center >= 0.25,]
		cf <- rbind(cf,coefficients(fit)[2,c(1:4)]) #linear model coefficients and statistics
		adj_r_sq <- rbind(adj_r_sq,fit$adj.r.squared)
		F_stat <- rbind(F_stat,fit$fstatistic)
		sum_repeats <- rbind(sum_repeats,sum(dat$num_rep_count))
		classes_w_models <- rbind(classes_w_models,i) #repeat taxon ids
		arm_lower <- rbind(arm_lower, CI(dat_arm$num_bp_rep,ci=0.95)[3]) #95% CI and mean of repeat content in arms
		arm_mean <- rbind(arm_mean, CI(dat_arm$num_bp_rep,ci=0.95)[2])
		arm_upper <- rbind(arm_upper, CI(dat_arm$num_bp_rep,ci=0.95)[1])
		center_lower <- rbind(center_lower, CI(dat_center$num_bp_rep,ci=0.95)[3]) #95% CI and mean of repeat content in arms and centers
		center_mean <- rbind(center_mean, CI(dat_center$num_bp_rep,ci=0.95)[2])
		center_upper <- rbind(center_upper, CI(dat_center$num_bp_rep,ci=0.95)[1])
		arm_cen_effect_size_lower <- rbind(arm_cen_effect_size_lower, cohen.d(dat_arm$num_bp_rep,dat_center$num_bp_rep)$conf.int[1]) #cohen's d effect sizes and 95% CI of repeat content in arms compared to centers
		arm_cen_effect_size <- rbind(arm_cen_effect_size, cohen.d(dat_arm$num_bp_rep,dat_center$num_bp_rep)$estimate)
		arm_cen_effect_size_upper <- rbind(arm_cen_effect_size_upper, cohen.d(dat_arm$num_bp_rep,dat_center$num_bp_rep)$conf.int[2]) 
		percent_genome_repeat <- rbind(percent_genome_repeat,((sum(dat$num_bp_rep)/arg3)*100)) #percentage of genome covered by repeat
		repeat_total_bp <- rbind(repeat_total_bp,sum(dat$num_bp_rep))} #total bp covered by repeat
	}
	#these should all be the same
	nrow(classes_w_models)
	nrow(sum_repeats)
	nrow(adj_r_sq)
	nrow(F_stat)
	nrow(cf)
	nrow(repeat_total_bp)
	nrow(arm_lower)
	nrow(arm_mean)
	nrow(arm_upper)
	nrow(center_lower)
	nrow(center_mean)
	nrow(center_upper)
	nrow(arm_cen_effect_size_lower)
	nrow(arm_cen_effect_size)
	nrow(arm_cen_effect_size_upper)
	#put stats together
	stat_df <- cbind(classes_w_models,sum_repeats,adj_r_sq,F_stat,cf,arm_lower,arm_mean,arm_upper,center_lower,center_mean,center_upper,	arm_cen_effect_size_lower,arm_cen_effect_size,arm_cen_effect_size_upper,percent_genome_repeat,repeat_total_bp)
	#getting data right, column names
	rownames(stat_df) <- NULL
	colnames(stat_df) <- c("rep_class","rep_count","adj_r2","F","numdf","dendf","beta_Estimate","beta_std_error","t_value","p_value", "arm_lower","arm_mean","arm_upper","center_lower","center_mean","center_upper","arm_cen_effect_size_lower","arm_cen_effect_size","	arm_cen_effect_size_upper", "percent_genome_repeat","repeat_total_bp")
	#getting data structure right
	stat_df <- as.data.frame(stat_df)
	
	stat_df$rep_count <-  as.numeric(levels(stat_df$rep_count))[stat_df$rep_count]
	stat_df$adj_r2 <-  as.numeric(levels(stat_df$adj_r2))[stat_df$adj_r2]
	stat_df$F <-  as.numeric(levels(stat_df$F))[stat_df$F]
	stat_df$numdf <-  as.numeric(levels(stat_df$numdf))[stat_df$numdf]
	stat_df$dendf <-  as.numeric(levels(stat_df$dendf))[stat_df$dendf]
	stat_df$beta_Estimate <-  as.numeric(levels(stat_df$beta_Estimate))[stat_df$beta_Estimate]
	stat_df$beta_std_error <-  as.numeric(levels(stat_df$beta_std_error))[stat_df$beta_std_error]
	stat_df$t_value <-  as.numeric(levels(stat_df$t_value))[stat_df$t_value]
	stat_df$p_value <-  as.numeric(levels(stat_df$p_value))[stat_df$p_value]
	stat_df$repeat_total_bp <-  as.numeric(levels(stat_df$repeat_total_bp))[stat_df$repeat_total_bp]
	stat_df$arm_lower <-  as.numeric(levels(stat_df$arm_lower))[stat_df$arm_lower]
	stat_df$arm_mean <-  as.numeric(levels(stat_df$arm_mean))[stat_df$arm_mean]
	stat_df$arm_upper <-  as.numeric(levels(stat_df$arm_upper))[stat_df$arm_upper]
	stat_df$center_lower <-  as.numeric(levels(stat_df$center_lower))[stat_df$center_lower]
	stat_df$center_mean <-  as.numeric(levels(stat_df$center_mean))[stat_df$center_mean]
	stat_df$center_upper <-  as.numeric(levels(stat_df$center_upper))[stat_df$center_upper]
	stat_df$arm_cen_effect_size_lower <-  as.numeric(levels(stat_df$arm_cen_effect_size_lower))[stat_df$arm_cen_effect_size_lower]
	stat_df$arm_cen_effect_size <-  as.numeric(levels(stat_df$arm_cen_effect_size))[stat_df$arm_cen_effect_size]
	stat_df$arm_cen_effect_size_upper <-  as.numeric(levels(stat_df$arm_cen_effect_size_upper))[stat_df$arm_cen_effect_size_upper]
	
	#adjust p values for multiple testing
	stat_df$p_adjust <- p.adjust(stat_df$p_value, method = "BH")
	
	#get significant (TRUE if significant)
	stat_df$p_sig <- stat_df$p_adjust < 0.05
	
	#sort by p
	
	stat_df_order <- stat_df[order(stat_df$p_value),]
	
	#add column of species id
	
	stat_df_order$species <- rep(arg2,nrow(stat_df_order))
	#print table to file
	write.table(stat_df_order, file = paste(arg2,"_",arg4, sep = "", collapse = NULL),row.names = FALSE, quote = FALSE, sep = "\t")
	#return the table
	return(stat_df_order)
}

#get models
	#arg3 is the genome size
briggsae_models_class <- get_models(cb,"briggsae",106936205, "class")
elegans_models_class <- get_models(ce,"elegans",101943821, "class")
inopinata_models_class <- get_models(ci,"inopinata",122994336, "class")
nigoni_models_class <- get_models(cn,"nigoni",119795280, "class")
remanei_models_class <-get_models(cr,"remanei",126105038, "class")

#put all models together

all_models_class <- rbind(briggsae_models_class,elegans_models_class,inopinata_models_class,nigoni_models_class,remanei_models_class, deparse.level = 1, make.row.names = TRUE, stringsAsFactors = default.stringsAsFactors())
	#All reported statistics regarding chromosome position and repeat content for given repeat classes are here! This is also the data Supplemental Figures 4-5.
write.table(all_models_class, file = "all_models_class.tsv",row.names = FALSE, quote = FALSE, sep = "\t")

#do for all other taxa

#orders
#remove unused factor levels

rep_order$rep_class <- droplevels(rep_order$rep_class)

#replace NA with "NA"

rep_order$rep_class = factor(rep_order$rep_class, levels=c(levels(rep_order$rep_class), "NA"))


rep_order$rep_class[is.na(rep_order$rep_class)] = "NA"


#subset species again

cb <- rep_order[rep_order$species == "briggsae",]
ce <- rep_order[rep_order$species == 'elegans',]
ci <- rep_order[rep_order$species == 'inopinata',]
cn <- rep_order[rep_order$species == 'nigoni',]
cr <- rep_order[rep_order$species == 'remanei',]

#get models
briggsae_models_order <- get_models(cb,"briggsae",106936205, "order")
elegans_models_order <- get_models(ce,"elegans",101943821, "order")
inopinata_models_order <- get_models(ci,"inopinata",122994336, "order")
nigoni_models_order <- get_models(cn,"nigoni",119795280, "order")
remanei_models_order <-get_models(cr,"remanei",126105038, "order")


all_models_order <- rbind(briggsae_models_order,elegans_models_order,inopinata_models_order,nigoni_models_order,remanei_models_order, deparse.level = 1, make.row.names = TRUE, stringsAsFactors = default.stringsAsFactors())

write.table(all_models_order, file = "all_models_order.tsv",row.names = FALSE, quote = FALSE, sep = "\t")
	#All reported statistics regarding chromosome position and repeat content for given repeat orders are here! This is also the data Supplemental Figures 6-7.

#superfamilies
#remove unused factor levels

rep_superfamily$rep_class <- droplevels(rep_superfamily$rep_class)

#replace NA with "NA"

rep_superfamily$rep_class = factor(rep_superfamily$rep_class, levels=c(levels(rep_superfamily$rep_class), "NA"))


rep_superfamily$rep_class[is.na(rep_superfamily$rep_class)] = "NA"


#subset species again

cb <- rep_superfamily[rep_superfamily$species == "briggsae",]
ce <- rep_superfamily[rep_superfamily$species == 'elegans',]
ci <- rep_superfamily[rep_superfamily$species == 'inopinata',]
cn <- rep_superfamily[rep_superfamily$species == 'nigoni',]
cr <- rep_superfamily[rep_superfamily$species == 'remanei',]

#get models
briggsae_models_superfamily <- get_models(cb,"briggsae",106936205, "superfamily")
elegans_models_superfamily <- get_models(ce,"elegans",101943821, "superfamily")
inopinata_models_superfamily <- get_models(ci,"inopinata",122994336, "superfamily")
nigoni_models_superfamily <- get_models(cn,"nigoni",119795280, "superfamily")
remanei_models_superfamily <-get_models(cr,"remanei",126105038, "superfamily")


all_models_superfamily <- rbind(briggsae_models_superfamily,elegans_models_superfamily,inopinata_models_superfamily,nigoni_models_superfamily,remanei_models_superfamily, deparse.level = 1, make.row.names = TRUE, stringsAsFactors = default.stringsAsFactors())

write.table(all_models_superfamily, file = "all_models_superfamily.tsv",row.names = FALSE, quote = FALSE, sep = "\t")
	#All reported statistics regarding chromosome position and repeat content for given repeat superfamiles are here! This is also the data for Figures 4c-4d, Supplemental Figures 8-9, Supplemental Figure 11.


#families
#remove unused factor levels

rep_family$rep_class <- droplevels(rep_family$rep_class)

#replace NA with "NA"

rep_family$rep_class = factor(rep_family$rep_class, levels=c(levels(rep_family$rep_class), "NA"))


rep_family$rep_class[is.na(rep_family$rep_class)] = "NA"


#subset species again

cb <- rep_family[rep_family$species == "briggsae",]
ce <- rep_family[rep_family$species == 'elegans',]
ci <- rep_family[rep_family$species == 'inopinata',]
cn <- rep_family[rep_family$species == 'nigoni',]
cr <- rep_family[rep_family$species == 'remanei',]

#get models
briggsae_models_family <- get_models(cb,"briggsae",106936205, "family")
elegans_models_family <- get_models(ce,"elegans",101943821, "family")
inopinata_models_family <- get_models(ci,"inopinata",122994336, "family")
nigoni_models_family <- get_models(cn,"nigoni",119795280, "family")
remanei_models_family <-get_models(cr,"remanei",126105038, "family")


all_models_family <- rbind(briggsae_models_family,elegans_models_family,inopinata_models_family,nigoni_models_family,remanei_models_family, deparse.level = 1, make.row.names = TRUE, stringsAsFactors = default.stringsAsFactors())
	#All reported statistics regarding chromosome position and repeat content for given repeat familes are here! This is also the data for Supplemental Figure 10.

write.table(all_models_family, file = "all_models_family.tsv",row.names = FALSE, quote = FALSE, sep = "\t")


#Relationship between repeat superfamily abundance and chromosomal structure (percent_genome_repeat and arm_cen_effect_size, LM for Figure 4c)


#get data right
all_models_superfamily$percent_genome_repeat <- as.numeric(levels(all_models_superfamily$percent_genome_repeat))[all_models_superfamily$percent_genome_repeat]


#subset species again, this time the data are the statistics for repeat superfamilies
cb <- all_models_superfamily[all_models_superfamily$species == "briggsae",]
ce <- all_models_superfamily[all_models_superfamily$species == 'elegans',]
ci <- all_models_superfamily[all_models_superfamily$species == 'inopinata',]
cn <- all_models_superfamily[all_models_superfamily$species == 'nigoni',]
cr <- all_models_superfamily[all_models_superfamily$species == 'remanei',]


#get linear models
	#log transformed

summary(lm(log(percent_genome_repeat+1) ~ log(arm_cen_effect_size+1), data=cb))
	#briggsae
#Call:
#lm(formula = log(percent_genome_repeat + 1) ~ log(arm_cen_effect_size + 
#    1), data = cb)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-0.4292 -0.2504 -0.1851  0.0463  1.9383 
#
#Coefficients:
#                             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)                    0.2011     0.1393   1.444    0.165  
#log(arm_cen_effect_size + 1)   2.4907     1.1172   2.229    0.038 *
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.5229 on 19 degrees of freedom
#Multiple R-squared:  0.2074,	Adjusted R-squared:  0.1656 
#F-statistic: 4.971 on 1 and 19 DF,  p-value: 0.03805



summary(lm(log(percent_genome_repeat+1) ~ log(arm_cen_effect_size+1), data=ce))
	#elegans
#Call:
#lm(formula = log(percent_genome_repeat + 1) ~ log(arm_cen_effect_size + 
#    1), data = ce)
#
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-0.31634 -0.10353 -0.03131  0.06810  0.65590 
#
#Coefficients:
#                             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                   0.01766    0.06083   0.290    0.774    
#log(arm_cen_effect_size + 1)  2.10221    0.34171   6.152 4.19e-06 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.2213 on 21 degrees of freedom
#Multiple R-squared:  0.6432,	Adjusted R-squared:  0.6262 
#F-statistic: 37.85 on 1 and 21 DF,  p-value: 4.194e-06



summary(lm(log(percent_genome_repeat+1) ~ log(arm_cen_effect_size+1), data=ci))
	#inopinata
#Call:
#lm(formula = log(percent_genome_repeat + 1) ~ log(arm_cen_effect_size + 
#    1), data = ci)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-0.5338 -0.4428 -0.3365  0.2786  1.6483 
#
#Coefficients:
#                             Estimate Std. Error t value Pr(>|t|)   
#(Intercept)                    0.4854     0.1545   3.142  0.00563 **
#log(arm_cen_effect_size + 1)  -1.5691     1.4466  -1.085  0.29237   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.69 on 18 degrees of freedom
#Multiple R-squared:  0.06136,	Adjusted R-squared:  0.009208 
#F-statistic: 1.177 on 1 and 18 DF,  p-value: 0.2924




summary(lm(log(percent_genome_repeat+1) ~ log(arm_cen_effect_size+1), data=cn))
	#nigoni

#Call:
#lm(formula = log(percent_genome_repeat + 1) ~ log(arm_cen_effect_size + 
#    1), data = cn)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-0.4747 -0.2749 -0.1229  0.1361  1.8310 
#
#Coefficients:
#                             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)                    0.2095     0.1268   1.652   0.1141  
#log(arm_cen_effect_size + 1)   4.5134     1.9176   2.354   0.0289 *
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.5091 on 20 degrees of freedom
#Multiple R-squared:  0.2169,	Adjusted R-squared:  0.1778 
#F-statistic:  5.54 on 1 and 20 DF,  p-value: 0.02893
#


summary(lm(log(percent_genome_repeat+1) ~ log(arm_cen_effect_size+1), data=cr))
	#remanei

#Call:
#lm(formula = log(percent_genome_repeat + 1) ~ log(arm_cen_effect_size + 
#    1), data = cr)
#
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-0.53022 -0.15078 -0.03027  0.10910  0.55554 
#
#Coefficients:
#                             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                   0.06868    0.07018   0.979    0.337    
#log(arm_cen_effect_size + 1)  2.67827    0.44651   5.998 3.42e-06 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.2654 on 24 degrees of freedom
#Multiple R-squared:  0.5999,	Adjusted R-squared:  0.5832 
#F-statistic: 35.98 on 1 and 24 DF,  p-value: 3.422e-06


#remove four superfamilies from C. inopinata to see if there is an impact on the global genomic landscape

#get data in there
rep_data <- read.table("/home/gavin/genome/genome/repeats_12-18-18/preparing_for_deposition_8-9-19/data/global_repeat_denisty_remove_four_superfamilies.tsv", sep="\t", header=TRUE)
#readable measure of bp
rep_data$MB <- rep_data$BP/1000000
#percentage of window repetitive
rep_data$perc_N <- (rep_data$num_rep/10000)*100
#order levels
rep_data$species <- factor(rep_data$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata"))

#just set aside inopinata
ci <- rep_data[rep_data$species == "inopinata",]

#get the model (for 10kb window, percentage repetitive ~ normalized distance from chromosome center)

#summary(lm(perc_N ~ norm_dist_center, data=ci))
#
#
#Call:
#lm(formula = perc_N ~ norm_dist_center, data = ci)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-15.478  -7.562  -3.331   4.364  91.019 
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)        6.6291     0.2051   32.32   <2e-16 ***
#norm_dist_center  17.6988     0.7104   24.91   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 11.38 on 12300 degrees of freedom
#Multiple R-squared:  0.04803,	Adjusted R-squared:  0.04796 
#F-statistic: 620.6 on 1 and 12300 DF,  p-value: < 2.2e-16

#Supplemental Figures 12-15



#############
#############
#############
#############
#section "Gene diversity is negatively correlated with repeat content in all species but C. inopinata"
#############
#############
#############
#############


#get gene density data in there

rep_dat <- read.table("/home/gavin/genome/genome/repeats_12-18-18/preparing_for_deposition_8-9-19/data/gene_density.tsv", sep="\t", header=TRUE)
	# (this is the same as file "all_gene_dens_norm_chr_pos.tsv" line 5488 repeats.sh)

#define arms and centers

rep_dat$chr_str_type <- ifelse(rep_dat$norm_dist_center >= 0.25,"arms", "centers")

#split by species

cb <- rep_dat[rep_dat$species == "briggsae",]
ce <- rep_dat[rep_dat$species == "elegans",]
ci <- rep_dat[rep_dat$species == "inopinata",]
cn <- rep_dat[rep_dat$species == "nigoni",]
cr <- rep_dat[rep_dat$species == "remanei",]
ci_no_tsp_cds <- rep_dat[rep_dat$species == "inopinata_no_transposon_cds",]

#lm for gene count ~ normalized chromosome position


summary(lm(gene_num ~ norm_dist_center, data=cb))
	#briggsae
#> summary(lm(gene_num ~ norm_dist_center, data=cb))
#
#Call:
#lm(formula = gene_num ~ norm_dist_center, data = cb)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-16.756  -5.303  -0.550   4.227  32.293 
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       24.8259     0.4407  56.338  < 2e-16 ***
#norm_dist_center -10.3903     1.5233  -6.821 1.53e-11 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 7.153 on 1052 degrees of freedom
#Multiple R-squared:  0.04235,	Adjusted R-squared:  0.04144 
#F-statistic: 46.52 on 1 and 1052 DF,  p-value: 1.525e-11


summary(lm(gene_num ~ norm_dist_center, data=cn))
	#nigoni
#> summary(lm(gene_num ~ norm_dist_center, data=cn))
#
#Call:
#lm(formula = gene_num ~ norm_dist_center, data = cn)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-18.447  -5.149  -0.520   4.509  39.334 
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       25.7438     0.4111  62.627  < 2e-16 ***
#norm_dist_center -10.6259     1.4207  -7.479 1.46e-13 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 7.063 on 1179 degrees of freedom
#Multiple R-squared:  0.0453,	Adjusted R-squared:  0.04449 
#F-statistic: 55.94 on 1 and 1179 DF,  p-value: 1.456e-13


summary(lm(gene_num ~ norm_dist_center, data=cr))
	#remanei
#> summary(lm(gene_num ~ norm_dist_center, data=cr))
#
#Call:
#lm(formula = gene_num ~ norm_dist_center, data = cr)
#
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-20.7976  -5.2958  -0.4358   4.8470  28.0769 
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       26.7030     0.4384   60.91   <2e-16 ***
#norm_dist_center -20.2768     1.5144  -13.39   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 7.755 on 1250 degrees of freedom
#Multiple R-squared:  0.1254,	Adjusted R-squared:  0.1247 
#F-statistic: 179.3 on 1 and 1250 DF,  p-value: < 2.2e-16


summary(lm(gene_num ~ norm_dist_center, data=ce))
	#elegans
#> summary(lm(gene_num ~ norm_dist_center, data=ce))
#
#Call:
#lm(formula = gene_num ~ norm_dist_center, data = ce)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-16.034  -4.610  -0.662   4.293  22.649 
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       24.8415     0.4283   58.00   <2e-16 ***
#norm_dist_center -16.4190     1.4803  -11.09   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 6.789 on 1003 degrees of freedom
#Multiple R-squared:  0.1093,	Adjusted R-squared:  0.1084 
#F-statistic:   123 on 1 and 1003 DF,  p-value: < 2.2e-16


summary(lm(gene_num ~ norm_dist_center, data=ci))
	#inopinata
#> summary(lm(gene_num ~ norm_dist_center, data=ci))
#
#Call:
#lm(formula = gene_num ~ norm_dist_center, data = ci)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-15.461  -3.838  -0.345   3.601  18.329 
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)        19.623      0.311  63.102  < 2e-16 ***
#norm_dist_center   -6.164      1.075  -5.737 1.21e-08 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 5.46 on 1231 degrees of freedom
#Multiple R-squared:  0.02604,	Adjusted R-squared:  0.02525 
#F-statistic: 32.91 on 1 and 1231 DF,  p-value: 1.214e-08


summary(lm(gene_num ~ norm_dist_center, data=ci_no_tsp_cds))
	#inopinata, no transposon-aligning genes
#> summary(lm(gene_num ~ norm_dist_center, data=ci_no_tsp_cds))
#
#Call:
#lm(formula = gene_num ~ norm_dist_center, data = ci_no_tsp_cds)
#
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-13.8017  -3.1808  -0.1751   3.0203  16.1597 
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       16.2504     0.2758  58.918   <2e-16 ***
#norm_dist_center  -0.8989     0.9531  -0.943    0.346    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 4.843 on 1231 degrees of freedom
#Multiple R-squared:  0.0007222,	Adjusted R-squared:  -8.96e-05 
#F-statistic: 0.8896 on 1 and 1231 DF,  p-value: 0.3458



#lm for repeat density ~ gene density


summary(lm(perc_N ~ gene_num, data=cb))
	#briggsae
#Call:
#lm(formula = perc_N ~ gene_num, data = cb)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-31.763  -7.262  -0.974   7.246  27.190 
#
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 36.53040    1.00191  36.461  < 2e-16 ***
#gene_num    -0.35057    0.04283  -8.185 7.82e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 10.15 on 1052 degrees of freedom
#Multiple R-squared:  0.05987,	Adjusted R-squared:  0.05897 
#F-statistic: 66.99 on 1 and 1052 DF,  p-value: 7.822e-16


summary(lm(perc_N ~ gene_num, data=cn))
	#nigoni
#Call:
#lm(formula = perc_N ~ gene_num, data = cn)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-29.937  -7.032  -0.230   6.813  37.164 
#
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 35.53012    0.94227  37.707   <2e-16 ***
#gene_num    -0.38860    0.03896  -9.974   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 9.671 on 1179 degrees of freedom
#Multiple R-squared:  0.07781,	Adjusted R-squared:  0.07703 
#F-statistic: 99.48 on 1 and 1179 DF,  p-value: < 2.2e-16



summary(lm(perc_N ~ gene_num, data=cr))
	#remanei
#Call:
#lm(formula = perc_N ~ gene_num, data = cr)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-39.277 -10.940  -0.805   9.185  69.108 
#
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 39.27686    1.10527   35.54   <2e-16 ***
#gene_num    -0.73195    0.04774  -15.33   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 14 on 1250 degrees of freedom
#Multiple R-squared:  0.1583,	Adjusted R-squared:  0.1576 
#F-statistic: 235.1 on 1 and 1250 DF,  p-value: < 2.2e-16



summary(lm(perc_N ~ gene_num, data=ce))
	#elegans
#Call:
#lm(formula = perc_N ~ gene_num, data = ce)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-19.921  -8.689  -2.457   7.815  47.771 
#
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 29.05977    1.07561   27.02   <2e-16 ***
#gene_num    -0.49833    0.04903  -10.16   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 11.17 on 1003 degrees of freedom
#Multiple R-squared:  0.09337,	Adjusted R-squared:  0.09247 
#F-statistic: 103.3 on 1 and 1003 DF,  p-value: < 2.2e-16


summary(lm(perc_N ~ gene_num, data=ci))
	#inopinata
#
#Call:
#lm(formula = perc_N ~ gene_num, data = ci)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-28.096  -8.457  -1.869   6.853  56.045 
#
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 20.97081    1.25256  16.742  < 2e-16 ***
#gene_num     0.53516    0.06626   8.077 1.57e-15 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 12.86 on 1231 degrees of freedom
#Multiple R-squared:  0.05033,	Adjusted R-squared:  0.04955 
#F-statistic: 65.23 on 1 and 1231 DF,  p-value: 1.573e-15
#


summary(lm(perc_N ~ gene_num, data=ci_no_tsp_cds))
	#inopinata, no transposon-aligning genes

#Call:
#lm(formula = perc_N ~ gene_num, data = ci_no_tsp_cds)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-27.493  -8.634  -1.958   6.321  60.321 
#
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 28.39372    1.29811  21.873   <2e-16 ***
#gene_num     0.14050    0.07754   1.812   0.0702 .  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 13.18 on 1231 degrees of freedom
#Multiple R-squared:  0.00266,	Adjusted R-squared:  0.00185 
#F-statistic: 3.283 on 1 and 1231 DF,  p-value: 0.07025


#get the arm center effect size (arms-centers) for gene density in all species (including inopinata without transposon cds)

#briggsae

summary(cb[cb$chr_str_type =="arms",]$gene_num)
summary(cb[cb$chr_str_type =="centers",]$gene_num)
cohen.d(cb[cb$chr_str_type =="arms",]$gene_num,cb[cb$chr_str_type =="centers",]$gene_num)

#summary(cb[cb$chr_str_type =="arms",]$gene_num)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   4.00   15.00   20.00   20.77   25.00   54.00 
#
#> summary(cb[cb$chr_str_type =="centers",]$gene_num)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   6.00   19.00   23.00   23.68   28.00   53.00 
# Cohen's d
# 
# d estimate: -0.4055431 (small)
# 95 percent confidence interval:
#      lower      upper 
# -0.5276602 -0.2834259 


#elegans

summary(ce[ce$chr_str_type =="arms",]$gene_num)
summary(ce[ce$chr_str_type =="centers",]$gene_num)
cohen.d(ce[ce$chr_str_type =="arms",]$gene_num,ce[ce$chr_str_type =="centers",]$gene_num)

#> summary(ce[ce$chr_str_type =="arms",]$gene_num)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   2.00   14.00   17.00   18.64   23.00   41.00 
#> summary(ce[ce$chr_str_type =="centers",]$gene_num)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   5.00   18.00   22.00   22.83   27.00   41.00 
#> cohen.d(ce[ce$chr_str_type =="arms",]$gene_num,ce[ce$chr_str_type =="centers",]$gene_num)
#
#Cohen's d
#
#d estimate: -0.6088739 (medium)
#95 percent confidence interval:
#     lower      upper 
#-0.7355101 -0.4822377 

#nigoni

summary(cn[cn$chr_str_type =="arms",]$gene_num)
summary(cn[cn$chr_str_type =="centers",]$gene_num)
cohen.d(cn[cn$chr_str_type =="arms",]$gene_num,cn[cn$chr_str_type =="centers",]$gene_num)

#> summary(cn[cn$chr_str_type =="arms",]$gene_num)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   2.00   16.00   21.00   21.61   26.00   62.00 
#> summary(cn[cn$chr_str_type =="centers",]$gene_num)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   8.00   20.00   24.00   24.55   29.00   47.00 
#> cohen.d(cn[cn$chr_str_type =="arms",]$gene_num,cn[cn$chr_str_type =="centers",]$gene_num)
#
#Cohen's d
#
#d estimate: -0.415401 (small)
#95 percent confidence interval:
#     lower      upper 
#-0.5308083 -0.2999936 



#remanei

summary(cr[cr$chr_str_type =="arms",]$gene_num)
summary(cr[cr$chr_str_type =="centers",]$gene_num)
cohen.d(cr[cr$chr_str_type =="arms",]$gene_num,cr[cr$chr_str_type =="centers",]$gene_num)

#> summary(cr[cr$chr_str_type =="arms",]$gene_num)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0.00   13.00   18.00   18.84   24.00   47.00 
#> summary(cr[cr$chr_str_type =="centers",]$gene_num)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   3.00   19.00   24.00   24.41   29.00   51.00 
#> cohen.d(cr[cr$chr_str_type =="arms",]$gene_num,cr[cr$chr_str_type =="centers",]$gene_num)
#
#Cohen's d
#
#d estimate: -0.7123222 (medium)
#95 percent confidence interval:
#     lower      upper 
#-0.8266760 -0.5979684 

#inopinata

summary(ci[ci$chr_str_type =="arms",]$gene_num)
summary(ci[ci$chr_str_type =="centers",]$gene_num)
cohen.d(ci[ci$chr_str_type =="arms",]$gene_num,ci[ci$chr_str_type =="centers",]$gene_num)

#> summary(ci[ci$chr_str_type =="arms",]$gene_num)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   2.00   13.00   17.00   17.18   21.00   35.00 
#> summary(ci[ci$chr_str_type =="centers",]$gene_num)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   4.00   15.00   18.50   18.99   22.00   37.00 
#> cohen.d(ci[ci$chr_str_type =="arms",]$gene_num,ci[ci$chr_str_type =="centers",]$gene_num)
#
#Cohen's d
#
#d estimate: -0.3318094 (small)
#95 percent confidence interval:
#     lower      upper 
#-0.4443205 -0.2192984 


#inopinata no transposon cds

summary(ci_no_tsp_cds[ci_no_tsp_cds$chr_str_type =="arms",]$gene_num)
summary(ci_no_tsp_cds[ci_no_tsp_cds$chr_str_type =="centers",]$gene_num)
cohen.d(ci_no_tsp_cds[ci_no_tsp_cds$chr_str_type =="arms",]$gene_num,ci_no_tsp_cds[ci_no_tsp_cds$chr_str_type =="centers",]$gene_num)

#> summary(ci_no_tsp_cds[ci_no_tsp_cds$chr_str_type =="arms",]$gene_num)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   2.00   12.00   15.00   15.86   19.00   32.00 
#> summary(ci_no_tsp_cds[ci_no_tsp_cds$chr_str_type =="centers",]$gene_num)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   3.00   13.00   16.00   16.19   19.00   31.00 
#> cohen.d(ci_no_tsp_cds[ci_no_tsp_cds$chr_str_type =="arms",]$gene_num,ci_no_tsp_cds[ci_no_tsp_cds$chr_str_type =="centers",]$gene_num)
#
#Cohen's d
#
#d estimate: -0.06872575 (negligible)
#95 percent confidence interval:
#      lower       upper 
#-0.18050353  0.04305204 


#get effect sizes of repeat density as well for all species

summary(cb[cb$chr_str_type =="arms",]$perc_N)
summary(cb[cb$chr_str_type =="centers",]$perc_N)
cohen.d(cb[cb$chr_str_type =="arms",]$perc_N,cb[cb$chr_str_type =="centers",]$perc_N)

#> summary(cb[cb$chr_str_type =="arms",]$perc_N)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  2.664  26.864  35.019  33.924  41.424  59.289 
#> summary(cb[cb$chr_str_type =="centers",]$perc_N)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  6.262  18.390  22.673  23.556  27.552  49.059 
#> cohen.d(cb[cb$chr_str_type =="arms",]$perc_N,cb[cb$chr_str_type =="centers",]$perc_N)
#
#Cohen's d
#
#d estimate: 1.139733 (large)
#95 percent confidence interval:
#   lower    upper 
#1.009407 1.270059 


summary(ce[ce$chr_str_type =="arms",]$perc_N)
summary(ce[ce$chr_str_type =="centers",]$perc_N)
cohen.d(ce[ce$chr_str_type =="arms",]$perc_N,ce[ce$chr_str_type =="centers",]$perc_N)

#> summary(ce[ce$chr_str_type =="arms",]$perc_N)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  6.045  18.390  25.617  26.612  34.586  63.874 
#> summary(ce[ce$chr_str_type =="centers",]$perc_N)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  2.382   7.093   9.683  10.802  12.927  54.704 
#> cohen.d(ce[ce$chr_str_type =="arms",]$perc_N,ce[ce$chr_str_type =="centers",]$perc_N)
#
#Cohen's d
#
#d estimate: 1.825696 (large)
#95 percent confidence interval:
#   lower    upper 
#1.678346 1.973047 



summary(ci[ci$chr_str_type =="arms",]$perc_N)
summary(ci[ci$chr_str_type =="centers",]$perc_N)
cohen.d(ci[ci$chr_str_type =="arms",]$perc_N,ci[ci$chr_str_type =="centers",]$perc_N)

#> summary(ci[ci$chr_str_type =="arms",]$perc_N)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  3.711  22.384  29.901  30.943  38.167  88.254 
#> summary(ci[ci$chr_str_type =="centers",]$perc_N)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  4.386  22.099  27.601  30.345  35.737  91.946 
#> cohen.d(ci[ci$chr_str_type =="arms",]$perc_N,ci[ci$chr_str_type =="centers",]$perc_N)
#
#Cohen's d
#
#d estimate: 0.04538246 (negligible)
#95 percent confidence interval:
#      lower       upper 
#-0.06637672  0.15714165 



summary(cn[cn$chr_str_type =="arms",]$perc_N)
summary(cn[cn$chr_str_type =="centers",]$perc_N)
cohen.d(cn[cn$chr_str_type =="arms",]$perc_N,cn[cn$chr_str_type =="centers",]$perc_N)

#> summary(cn[cn$chr_str_type =="arms",]$perc_N)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  2.743  24.915  30.908  30.729  36.913  64.922 
#> summary(cn[cn$chr_str_type =="centers",]$perc_N)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  5.427  15.649  21.117  22.386  28.113  60.105 
#> cohen.d(cn[cn$chr_str_type =="arms",]$perc_N,cn[cn$chr_str_type =="centers",]$perc_N)
#
#Cohen's d
#
#d estimate: 0.9103546 (large)
#95 percent confidence interval:
#    lower     upper 
#0.7904036 1.0303056 



summary(cr[cr$chr_str_type =="arms",]$perc_N)
summary(cr[cr$chr_str_type =="centers",]$perc_N)
cohen.d(cr[cr$chr_str_type =="arms",]$perc_N,cr[cr$chr_str_type =="centers",]$perc_N)

#> summary(cr[cr$chr_str_type =="arms",]$perc_N)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0.00   24.02   32.34   32.32   39.92   87.16 
#> summary(cr[cr$chr_str_type =="centers",]$perc_N)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.802   5.593   9.651  14.558  20.843  60.956 
#> cohen.d(cr[cr$chr_str_type =="arms",]$perc_N,cr[cr$chr_str_type =="centers",]$perc_N)
#
#Cohen's d
#
#d estimate: 1.432227 (large)
#95 percent confidence interval:
#   lower    upper 
#1.307930 1.556525 



summary(ci_no_tsp_cds[ci_no_tsp_cds$chr_str_type =="arms",]$perc_N)
summary(ci_no_tsp_cds[ci_no_tsp_cds$chr_str_type =="centers",]$perc_N)
cohen.d(ci_no_tsp_cds[ci_no_tsp_cds$chr_str_type =="arms",]$perc_N,ci_no_tsp_cds[ci_no_tsp_cds$chr_str_type =="centers",]$perc_N)


#> summary(ci_no_tsp_cds[ci_no_tsp_cds$chr_str_type =="arms",]$perc_N)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  3.711  22.384  29.901  30.943  38.167  88.254 
#> summary(ci_no_tsp_cds[ci_no_tsp_cds$chr_str_type =="centers",]$perc_N)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  4.386  22.099  27.601  30.345  35.737  91.946 
#> cohen.d(ci_no_tsp_cds[ci_no_tsp_cds$chr_str_type =="arms",]$perc_N,ci_no_tsp_cds[ci_no_tsp_cds$chr_str_type =="centers",]$perc_N)
#
#Cohen's d
#
#d estimate: 0.04538246 (negligible)
#95 percent confidence interval:
#      lower       upper 
#-0.06637672  0.15714165 


#all effect size stats for genes and repeats were put into effect_sizes.tsv and effect_sizes_b.tsv

#lm gene arm-cen effect size and repeats arm-cen effect size


eff_data <- read.table("/home/gavin/genome/genome/repeats_12-18-18/preparing_for_deposition_8-9-19/data/effect_sizes_b.tsv", sep="\t", header=TRUE)


summary(lm(repeat_eff_size ~ gene_eff_size, data=eff_data))
#
#Call:
#lm(formula = repeat_eff_size ~ gene_eff_size, data = eff_data)
#
#Residuals:
#       1        2        3        4        5        6 
# 0.29068  0.41082 -0.59848  0.03387 -0.27052  0.13362 
#
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)  
#(Intercept)    -0.2795     0.3930  -0.711    0.516  
#gene_eff_size  -2.7828     0.8348  -3.334    0.029 *
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.4194 on 4 degrees of freedom
#Multiple R-squared:  0.7353,	Adjusted R-squared:  0.6692 
#F-statistic: 11.11 on 1 and 4 DF,  p-value: 0.02901


#without inopinata (the one that includes transposon-aligning genes)
eff_data_no_inopinata_with_tsp <- eff_data[eff_data$species != "inopinata",]

summary(lm(repeat_eff_size ~ gene_eff_size, data=eff_data_no_inopinata_with_tsp))

#
#Call:
#lm(formula = repeat_eff_size ~ gene_eff_size, data = eff_data_no_inopinata_with_tsp)
#
#Residuals:
#       1        2        4        5        6 
# 0.16100  0.33657 -0.09312 -0.31657 -0.08788 
#
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)  
#(Intercept)   -0.03925    0.29303  -0.134   0.9019  
#gene_eff_size -2.51017    0.59324  -4.231   0.0242 *
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.292 on 3 degrees of freedom
#Multiple R-squared:  0.8565,	Adjusted R-squared:  0.8086 
#F-statistic:  17.9 on 1 and 3 DF,  p-value: 0.02415

#############
#############
#############
#############
#section "Simulations reveal chromosomal heterogeneity in insertion fitness effects is suficient for promoting non-random genomic repetitive organization"
#############
#############
#############
#############

#see folder "SLiM_scripts" for how simulations and data were generated.

#prepare SLiM data for plotting and arm-center effect size statistics

################################################
## Combine the results of TE simulations #######
################################################

setwd("~/Documents/Phillips_lab/drafts/Gavin_Repeats/SLIM_TEs/EXTRA")

#load files
my_files<-list.files(pattern="TE")

#combine them
TE<-c();
for (i in 1:(length(my_files))) {
  print(my_files[i])
  B = read.csv(file = my_files[i],
               header = F,
               sep = " ")
  TE <- rbind(TE, data.frame(B))
}
colnames(TE)<-c("ID","DISARM","DISCENT","ACTARM","ACTCENT")

TE$ARM<-TE$DISARM + TE$ACTARM
TE$CENTER<-TE$DISCENT + TE$ACTCENT

#edit the labels
TE$SIM<-gsub("[0-9]+_TE_DISARM_DISCENT_ACTARM_ACTCENT_","",perl=T,TE$ID)
TE$SIM<-gsub(".txt","",perl=T,TE$SIM)
TE$DOM<-gsub("NO_DOM_.*","NO_DOMAINS", perl=T,TE$SIM)
TE$DOM<-gsub("^DOM_.*","HIGH_LOW_HIGH", perl=T,TE$DOM)
TE$SCENARIO<-gsub("NO_DOM_(.*)","\\1", perl=T,TE$SIM)
TE$SCENARIO<-gsub("DOM_(.*)","\\1", perl=T,TE$SCENARIO)
TE$SCENARIO<-gsub("CENTR", "CENTER",TE$SCENARIO)
TE$SCENARIO<-gsub("-LOSS", "_LOSS",TE$SCENARIO)

#remove some simulations
'%!in%' <- function(x,y)!('%in%'(x,y))
TE<-TE[TE$ID %!in% as.character(TE[TE$SIM=="DOM_SEL_ALL_WEAK_POS_MUT",]$ID[51:59]),]
TE<-TE[TE$ID %!in% as.character(TE[TE$SIM=="NO_DOM_SEL_ALL_WEAK_POS_MUT",]$ID[51:66]),]


#now it's ok
#table(TE[,c("DOM","SCENARIO")])
#SCENARIO
#DOM             NEUTRAL SEL_ALL SEL_ALL_WEAK SEL_ALL_WEAK_POS_MUT
#HIGH_LOW_HIGH      50      50           50                   50
#NO_DOMAINS         50      50           50                   50
#SCENARIO
#DOM             SEL_ARM_WEAK_CENTER_STRONG SEL_ARM_WEAK_CENTER_STRONG_LOSS SEL_CENTER
#HIGH_LOW_HIGH                         50                              50         50
#NO_DOMAINS                            50                              50         50
#SCENARIO
#DOM             SEL_CENTER_WEAK
#HIGH_LOW_HIGH              50
#NO_DOMAINS                 50


#create a data.frame for plots and statistics
NEWTE<-c();
for (i in 1:nrow(TE)) {
  NEWTE <-
    rbind(NEWTE,data.frame(
      POS = 0,
      VALUE = TE$ARM[i],
      DOM = TE$DOM[i],
      SCENARIO = TE$SCENARIO[i],
      SIM = TE$ID[i]
    ));
  NEWTE <-
    rbind(NEWTE,data.frame(
      POS = 0.9,
      VALUE = TE$ARM[i],
      DOM = TE$DOM[i],
      SCENARIO = TE$SCENARIO[i],
      SIM = TE$ID[i]
    ));
  NEWTE <-
    rbind(NEWTE,data.frame(
      POS = 1.1,
      VALUE = TE$CENTER[i],
      DOM = TE$DOM[i],
      SCENARIO = TE$SCENARIO[i],
      SIM = TE$ID[i]
    ));
  NEWTE <-
    rbind(NEWTE,data.frame(
      POS = 1.9,
      VALUE = TE$CENTER[i],
      DOM = TE$DOM[i],
      SCENARIO = TE$SCENARIO[i],
      SIM = TE$ID[i]
    ));
  NEWTE <-
    rbind(NEWTE,data.frame(
      POS = 2.1,
      VALUE = TE$ARM[i],
      DOM = TE$DOM[i],
      SCENARIO = TE$SCENARIO[i],
      SIM = TE$ID[i]
    ));
  NEWTE <-
    rbind(NEWTE,data.frame(
      POS = 3,
      VALUE = TE$ARM[i],
      DOM = TE$DOM[i],
      SCENARIO = TE$SCENARIO[i],
      SIM = TE$ID[i]
    ));
  
}

#reorder scenarious
#NEWTE$SCENARIO<-factor(NEWTE$SCENARIO,levels=c("NEUTRAL","SEL_ALL_WEAK", "SEL_ALL", "SEL_CENTER_WEAK","SEL_CENTER","SEL_ARM_WEAK_CENTER_STRONG","SEL_ARM_WEAK_CENTER_STRONG_LOSS","SEL_ALL_WEAK_POS_MUT" ))


#save the tables
write.table(
  NEWTE,
  file = "simulations_TE_table_for_plots_50_replicates.txt",
  append = FALSE,
  quote = FALSE,
  sep = "\t",
  eol = "\n",
  na = "NA",
  dec = ".",
  row.names = FALSE,
  col.names = TRUE
)
write.table(
  TE,
  file = "simulations_TE_full_table_50_replicates.txt",
  append = FALSE,
  quote = FALSE,
  sep = "\t",
  eol = "\n",
  na = "NA",
  dec = ".",
  row.names = FALSE,
  col.names = TRUE
)


#get arm-center effect sizes (regarding repeat density) for all evolutionary scenarios


#get data in there
rep_dat <- read.table("/home/gavin/genome/genome/repeats_12-18-18/preparing_for_deposition_8-9-19/data/simulations_TE_full_table_50_replicates.txt", sep="\t", header=TRUE)


#get arm-center effect sizes

#make empty things to be filled up later
scen_type <- NULL
arm_cen_effect_size_lower <- NULL
arm_cen_effect_size <- NULL
arm_cen_effect_size_upper <- NULL
wilcox_w <- NULL
wilcox_p <- NULL
mean_tot_te <- NULL
	#get stats for every evolutionary scenario
for (i in levels(rep_dat$SIM)){
	dat <- rep_dat[rep_dat$SIM == i,]
	mean_tot_te <- rbind(mean_tot_te,mean(dat$ARM+dat$CENTER)) #average total number of TE's
	arm_cen_effect_size_lower <- rbind(arm_cen_effect_size_lower, cohen.d(dat$ARM,dat$CENTER)$conf.int[1]) #arm-center effect size stats
	arm_cen_effect_size <- rbind(arm_cen_effect_size, cohen.d(dat$ARM,dat$CENTER)$estimate)
	arm_cen_effect_size_upper <- rbind(arm_cen_effect_size_upper, cohen.d(dat$ARM,dat$CENTER)$conf.int[2])
	scen_type <- rbind(scen_type,i) #scenario id
	wilcox_w <- rbind(wilcox_w, wilcox.test(dat$ARM,dat$CENTER,exact=FALSE)$statistic[[1]]) #wilcox rank sum test statistics
	wilcox_p <- rbind(wilcox_p, wilcox.test(dat$ARM,dat$CENTER,exact=FALSE)$p.value[[1]])
	stat_df <- cbind(scen_type,mean_tot_te,arm_cen_effect_size_lower,arm_cen_effect_size,arm_cen_effect_size_upper,wilcox_w,wilcox_p) #put all together
}

#get data right
rownames(stat_df) <- NULL
colnames(stat_df) <- c("SCENARIO","mean_tot_te","arm_cen_effect_size_lower","arm_cen_effect_size","arm_cen_effect_size_upper","wilcox_w","wilcox_p")

stat_df <- as.data.frame(stat_df)
stat_df$mean_tot_te <-  as.numeric(levels(stat_df$mean_tot_te))[stat_df$mean_tot_te]
stat_df$arm_cen_effect_size_lower <-  as.numeric(levels(stat_df$arm_cen_effect_size_lower))[stat_df$arm_cen_effect_size_lower]
stat_df$arm_cen_effect_size <-  as.numeric(levels(stat_df$arm_cen_effect_size))[stat_df$arm_cen_effect_size]
stat_df$arm_cen_effect_size_upper <-  as.numeric(levels(stat_df$arm_cen_effect_size_upper))[stat_df$arm_cen_effect_size_upper]
stat_df$wilcox_w <-  as.numeric(levels(stat_df$wilcox_w))[stat_df$wilcox_w]
stat_df$wilcox_p <-  as.numeric(levels(stat_df$wilcox_p))[stat_df$wilcox_p]

#write the data to a table
write.table(stat_df, file = "simulations_effect_sizes_wilcox_tests",row.names = FALSE, quote = FALSE, sep = "\t")

#see the stats
stat_df

