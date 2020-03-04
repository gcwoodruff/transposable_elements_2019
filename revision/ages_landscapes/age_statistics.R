#get data in there 
rep_data <- read.table("kimura_distances.tsv", sep="\t", header=TRUE)

#get species factors right
rep_data$species <- factor(rep_data$species, levels = c("briggsae","nigoni","remanei","elegans","inopinata"))


#get readable unit of bp
rep_data$MB <- rep_data$BP/1000000


#define chromosome arms and centers by middle or outer half of chromosome
rep_data$chr_str_type <- ifelse(rep_data$norm_dist_center >= 0.25,"arms", "centers")
#get the factors right
rep_data$chr_str_type <- factor(rep_data$chr_str_type, levels = c("centers","arms"))



library(dplyr)
library(effsize)
library(Rmisc)



rep_class <- rep_data[rep_data$rep_taxononmic_rank == "class",]
rep_order <- rep_data[rep_data$rep_taxononmic_rank == "order",]
rep_superfamily <- rep_data[rep_data$rep_taxononmic_rank == "superfamily",]
rep_family <- rep_data[rep_data$rep_taxononmic_rank == "family",]

#classes first

#remove unused factor levels

rep_class$rep_type <- droplevels(rep_class$rep_type)

#replace NA with "NA"

rep_class$rep_type = factor(rep_class$rep_type, levels=c(levels(rep_class$rep_type), "NA"))


rep_class$rep_type[is.na(rep_class$rep_type)] = "NA"

#subset species

cb <- rep_class[rep_class$species == "briggsae",]
ce <- rep_class[rep_class$species == 'elegans',]
ci <- rep_class[rep_class$species == 'inopinata',]
cn <- rep_class[rep_class$species == 'nigoni',]
cr <- rep_class[rep_class$species == 'remanei',]

# get_models(cb,"briggsae", "class")

get_models <- function(arg1, arg2, arg3){
		#make empty variables to fill later
	stat_df <- NULL
	stat_df_order <- NULL
	cf <- NULL
	adj_r_sq <- NULL
	F_stat <- NULL
	classes_w_models <- NULL
	
	arm_lower <- NULL
	arm_mean <- NULL
	arm_upper <- NULL
	
	center_lower <- NULL
	center_mean <- NULL
	center_upper <- NULL
	
	arm_cen_effect_size_lower <- NULL
	arm_cen_effect_size <- NULL
	arm_cen_effect_size_upper <- NULL

	#extract statistics regarding repeat density and chromosome position for each repeat taxon (ie, levels in rep_type)
	for (i in levels(arg1$rep_type)){
		dat <- arg1[arg1$rep_type == i,]
		if (nrow(dat) >14){fit <- summary(lm(kimura_distance ~ norm_dist_center, data=dat)) #excluding repeat taxa that are in <15 genomic windows across the genome
		dat_center <- dat[dat$norm_dist_center < 0.25,] #define chromosome arms and centers
		dat_arm <- dat[dat$norm_dist_center >= 0.25,]
		cf <- rbind(cf,coefficients(fit)[2,c(1:4)]) #linear model coefficients and statistics
		adj_r_sq <- rbind(adj_r_sq,fit$adj.r.squared)
		F_stat <- rbind(F_stat,fit$fstatistic)
		classes_w_models <- rbind(classes_w_models,i) #repeat taxon ids
		arm_lower <- rbind(arm_lower, CI(dat_arm$kimura_distance,ci=0.95)[3]) #95% CI and mean of repeat content in arms
		arm_mean <- rbind(arm_mean, CI(dat_arm$kimura_distance,ci=0.95)[2])
		arm_upper <- rbind(arm_upper, CI(dat_arm$kimura_distance,ci=0.95)[1])
		center_lower <- rbind(center_lower, CI(dat_center$kimura_distance,ci=0.95)[3]) #95% CI and mean of repeat content in arms and centers
		center_mean <- rbind(center_mean, CI(dat_center$kimura_distance,ci=0.95)[2])
		center_upper <- rbind(center_upper, CI(dat_center$kimura_distance,ci=0.95)[1])
		arm_cen_effect_size_lower <- rbind(arm_cen_effect_size_lower, cohen.d(dat_arm$kimura_distance,dat_center$kimura_distance)$conf.int[1]) #cohen's d effect sizes and 95% CI of repeat content in arms compared to centers
		arm_cen_effect_size <- rbind(arm_cen_effect_size, cohen.d(dat_arm$kimura_distance,dat_center$kimura_distance)$estimate)
		arm_cen_effect_size_upper <- rbind(arm_cen_effect_size_upper, cohen.d(dat_arm$kimura_distance,dat_center$kimura_distance)$conf.int[2])}
	}
	#these should all be the same
	nrow(classes_w_models)
	nrow(adj_r_sq)
	nrow(F_stat)
	nrow(cf)
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
	stat_df <- cbind(classes_w_models,adj_r_sq,F_stat,cf,arm_lower,arm_mean,arm_upper,center_lower,center_mean,center_upper,arm_cen_effect_size_lower,arm_cen_effect_size,arm_cen_effect_size_upper)
	#getting data right, column names
	rownames(stat_df) <- NULL
	colnames(stat_df) <- c("rep_class","adj_r2","F","numdf","dendf","beta_Estimate","beta_std_error","t_value","p_value", "arm_lower","arm_mean","arm_upper","center_lower","center_mean","center_upper","arm_cen_effect_size_lower","arm_cen_effect_size","arm_cen_effect_size_upper")
	#getting data structure right
	stat_df <- as.data.frame(stat_df)
	
	stat_df$adj_r2 <-  as.numeric(levels(stat_df$adj_r2))[stat_df$adj_r2]
	stat_df$F <-  as.numeric(levels(stat_df$F))[stat_df$F]
	stat_df$numdf <-  as.numeric(levels(stat_df$numdf))[stat_df$numdf]
	stat_df$dendf <-  as.numeric(levels(stat_df$dendf))[stat_df$dendf]
	stat_df$beta_Estimate <-  as.numeric(levels(stat_df$beta_Estimate))[stat_df$beta_Estimate]
	stat_df$beta_std_error <-  as.numeric(levels(stat_df$beta_std_error))[stat_df$beta_std_error]
	stat_df$t_value <-  as.numeric(levels(stat_df$t_value))[stat_df$t_value]
	stat_df$p_value <-  as.numeric(levels(stat_df$p_value))[stat_df$p_value]
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
	write.table(stat_df_order, file = paste(arg2,"_",arg3, sep = "", collapse = NULL),row.names = FALSE, quote = FALSE, sep = "\t")
	#return the table
	return(stat_df_order)
}


briggsae_models_class <- get_models(cb,"briggsae","class")
elegans_models_class <- get_models(ce,"elegans", "class")
inopinata_models_class <- get_models(ci,"inopinata","class")
nigoni_models_class <- get_models(cn,"nigoni","class")
remanei_models_class <-get_models(cr,"remanei","class")

#put all models together

all_models_class <- rbind(briggsae_models_class,elegans_models_class,inopinata_models_class,nigoni_models_class,remanei_models_class, deparse.level = 1, make.row.names = TRUE, stringsAsFactors = default.stringsAsFactors())
	#All reported statistics regarding chromosome position and repeat content for given repeat classes are here! This is also the data Supplemental Figures 4-5.
write.table(all_models_class, file = "kimura_distances_models_class.tsv",row.names = FALSE, quote = FALSE, sep = "\t")



#do for all other taxa

#orders
#remove unused factor levels

rep_order$rep_type <- droplevels(rep_order$rep_type)

#replace NA with "NA"

rep_order$rep_type = factor(rep_order$rep_type, levels=c(levels(rep_order$rep_type), "NA"))


rep_order$rep_type[is.na(rep_order$rep_type)] = "NA"


#subset species again

cb <- rep_order[rep_order$species == "briggsae",]
ce <- rep_order[rep_order$species == 'elegans',]
ci <- rep_order[rep_order$species == 'inopinata',]
cn <- rep_order[rep_order$species == 'nigoni',]
cr <- rep_order[rep_order$species == 'remanei',]

#get models
briggsae_models_order <- get_models(cb,"briggsae","order")
elegans_models_order <- get_models(ce,"elegans","order")
inopinata_models_order <- get_models(ci,"inopinata","order")
nigoni_models_order <- get_models(cn,"nigoni","order")
remanei_models_order <-get_models(cr,"remanei","order")


all_models_order <- rbind(briggsae_models_order,elegans_models_order,inopinata_models_order,nigoni_models_order,remanei_models_order, deparse.level = 1, make.row.names = TRUE, stringsAsFactors = default.stringsAsFactors())

write.table(all_models_order, file = "kimura_distances_models_order.tsv",row.names = FALSE, quote = FALSE, sep = "\t")
	#All reported statistics regarding chromosome position and repeat content for given repeat orders are here! This is also the data Supplemental Figures 6-7.




cb <- rep_superfamily[rep_superfamily$species == "briggsae",]
ce <- rep_superfamily[rep_superfamily$species == 'elegans',]
ci <- rep_superfamily[rep_superfamily$species == 'inopinata',]
cn <- rep_superfamily[rep_superfamily$species == 'nigoni',]
cr <- rep_superfamily[rep_superfamily$species == 'remanei',]

#get models
briggsae_models_superfamily <- get_models(cb,"briggsae","superfamily")
elegans_models_superfamily <- get_models(ce,"elegans","superfamily")
inopinata_models_superfamily <- get_models(ci,"inopinata","superfamily")
nigoni_models_superfamily <- get_models(cn,"nigoni","superfamily")
remanei_models_superfamily <-get_models(cr,"remanei","superfamily")


all_models_superfamily <- rbind(briggsae_models_superfamily,elegans_models_superfamily,inopinata_models_superfamily,nigoni_models_superfamily,remanei_models_superfamily, deparse.level = 1, make.row.names = TRUE, stringsAsFactors = default.stringsAsFactors())

write.table(all_models_superfamily, file = "kimura_distances_lm_ac_effect_sizes/kimura_distances_models_superfamily.tsv",row.names = FALSE, quote = FALSE, sep = "\t")
	#All reported statistics regarding chromosome position and repeat content for given repeat superfamilys are here! This is also the data Supplemental Figures 6-7.



cb <- rep_family[rep_family$species == "briggsae",]
ce <- rep_family[rep_family$species == 'elegans',]
ci <- rep_family[rep_family$species == 'inopinata',]
cn <- rep_family[rep_family$species == 'nigoni',]
cr <- rep_family[rep_family$species == 'remanei',]

#get models
briggsae_models_family <- get_models(cb,"briggsae","family")
elegans_models_family <- get_models(ce,"elegans","family")
inopinata_models_family <- get_models(ci,"inopinata","family")
nigoni_models_family <- get_models(cn,"nigoni","family")
remanei_models_family <-get_models(cr,"remanei","family")


all_models_family <- rbind(briggsae_models_family,elegans_models_family,inopinata_models_family,nigoni_models_family,remanei_models_family, deparse.level = 1, make.row.names = TRUE, stringsAsFactors = default.stringsAsFactors())

write.table(all_models_family, file = "kimura_distances_models_family.tsv",row.names = FALSE, quote = FALSE, sep = "\t")
	#All reported statistics regarding chromosome position and repeat content for given repeat familys are here! This is also the data Supplemental Figures 6-7.






#models, sina plot for GLOBAL genomic distribution of kimura distances

rep_data <- read.table("kimura_distances_global_landscape_norm_dist_center.tsv", sep="\t", header=TRUE)


library(effsize)



#set aside species
cb <- rep_data[rep_data$species == "briggsae",]
ce <- rep_data[rep_data$species == 'elegans',]
ci <- rep_data[rep_data$species == 'inopinata',]
cn <- rep_data[rep_data$species == 'nigoni',]
cr <- rep_data[rep_data$species == 'remanei',]

#linear models of kimura distance and normalized chromsome position

summary(lm(kimura_distance ~ norm_dist_center, data=cb))

	#briggsae

#	Call:
#lm(formula = kimura_distance ~ norm_dist_center, data = cb)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-20.243  -2.875   0.373   3.231  36.324 
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)      18.38597    0.09822   187.2   <2e-16 ***
#norm_dist_center  4.44927    0.33956    13.1   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 4.967 on 10363 degrees of freedom
#Multiple R-squared:  0.0163,	Adjusted R-squared:  0.0162 
#F-statistic: 171.7 on 1 and 10363 DF,  p-value: < 2.2e-16



summary(lm(kimura_distance ~ norm_dist_center, data=ce))

	#elegans
#> summary(lm(kimura_distance ~ norm_dist_center, data=ce))
#
#Call:
#lm(formula = kimura_distance ~ norm_dist_center, data = ce)
#
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-18.5202  -3.2351  -0.1678   3.2120  31.3192 
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       18.5492     0.1105 167.800   <2e-16 ***
#norm_dist_center  -3.6912     0.3779  -9.768   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 5.326 on 9653 degrees of freedom
#Multiple R-squared:  0.009788,	Adjusted R-squared:  0.009685 
#F-statistic: 95.42 on 1 and 9653 DF,  p-value: < 2.2e-16




summary(lm(kimura_distance ~ norm_dist_center, data=cn))

	#nigoni

#> summary(lm(kimura_distance ~ norm_dist_center, data=cn))
#
#Call:
#lm(formula = kimura_distance ~ norm_dist_center, data = cn)
#
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-20.8200  -2.7646   0.1681   2.9211  29.2314 
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)      19.98701    0.08731 228.931  < 2e-16 ***
#norm_dist_center  1.72363    0.30210   5.705 1.19e-08 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 4.66 on 11557 degrees of freedom
#Multiple R-squared:  0.002809,	Adjusted R-squared:  0.002723 
#F-statistic: 32.55 on 1 and 11557 DF,  p-value: 1.189e-08



summary(lm(kimura_distance ~ norm_dist_center, data=cr))

	#remanei

#> summary(lm(kimura_distance ~ norm_dist_center, data=cr))
#
#Call:
#lm(formula = kimura_distance ~ norm_dist_center, data = cr)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-18.093  -3.205  -0.109   3.026  37.002 
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       14.0657     0.1117  125.93   <2e-16 ***
#norm_dist_center   8.0627     0.3695   21.82   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 5.447 on 11087 degrees of freedom
#Multiple R-squared:  0.04117,	Adjusted R-squared:  0.04108 
#F-statistic:   476 on 1 and 11087 DF,  p-value: < 2.2e-16


summary(lm(kimura_distance ~ norm_dist_center, data=ci))

	#inopinata


#> summary(lm(kimura_distance ~ norm_dist_center, data=ci))
#
#Call:
#lm(formula = kimura_distance ~ norm_dist_center, data = ci)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-14.058  -2.923  -0.572   2.375  31.666 
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)      13.93122    0.08962 155.445   <2e-16 ***
#norm_dist_center  0.25480    0.30997   0.822    0.411    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 4.872 on 11829 degrees of freedom
#Multiple R-squared:  5.712e-05,	Adjusted R-squared:  -2.741e-05 
#F-statistic: 0.6757 on 1 and 11829 DF,  p-value: 0.4111

cohen.d(cb[cb$chr_str_type == "arms",]$kimura_distance,cb[cb$chr_str_type == "centers",]$kimura_distance)

#Cohen's d
#
#d estimate: 0.2785744 (small)
#95 percent confidence interval:
#    lower     upper 
#0.2398798 0.3172690 
#


cohen.d(ce[ce$chr_str_type == "arms",]$kimura_distance,ce[ce$chr_str_type == "centers",]$kimura_distance)

#
#Cohen's d
#
#d estimate: -0.2060447 (small)
#95 percent confidence interval:
#     lower      upper 
#-0.2460694 -0.1660201 



cohen.d(ci[ci$chr_str_type == "arms",]$kimura_distance,ci[ci$chr_str_type == "centers",]$kimura_distance)
#
#Cohen's d
#
#d estimate: -0.002668666 (negligible)
#95 percent confidence interval:
#      lower       upper 
#-0.03871126  0.03337393 
#
#

cohen.d(cn[cn$chr_str_type == "arms",]$kimura_distance,cn[cn$chr_str_type == "centers",]$kimura_distance)

#Cohen's d
#
#d estimate: 0.1168299 (negligible)
#95 percent confidence interval:
#     lower      upper 
#0.08033428 0.15332548 
#



cohen.d(cr[cr$chr_str_type == "arms",]$kimura_distance,cr[cr$chr_str_type == "centers",]$kimura_distance)


#Cohen's d
#
#d estimate: 0.3469468 (small)
#95 percent confidence interval:
#    lower     upper 
#0.3092000 0.3846936 


#mean, all




#variance
#arms

rep_arms <- rep_data[rep_data$chr_str_type == "arms",]
rep_centers <- rep_data[rep_data$chr_str_type == "centers",]

aggregate(rep_arms$kimura_distance ~ rep_arms$species, FUN = function(x) var(x))

#  rep_arms$species rep_arms$kimura_distance
#1         briggsae                 21.67267
#2           nigoni                 21.01363
#3          remanei                 20.01023
#4          elegans                 20.77421
#5        inopinata                 20.49256


#centers

aggregate(rep_centers$kimura_distance ~ rep_centers$species, FUN = function(x) var(x))


#  rep_centers$species rep_centers$kimura_distance
#1            briggsae                    27.57867
#2              nigoni                    22.40377
#3             remanei                    42.65313
#4             elegans                    36.41553
#5           inopinata                    27.01096
	
	#in all cases, variance in distance is greater in centers than in arms

#hypothesis test for variance



var.test(cb[cb$chr_str_type == "arms",]$kimura_distance,cb[cb$chr_str_type == "centers",]$kimura_distance)


#	F test to compare two variances
#
#data:  cb[cb$chr_str_type == "arms", ]$kimura_distance and cb[cb$chr_str_type == "centers", ]$kimura_distance
#F = 0.78585, num df = 5216, denom df = 5147, p-value < 2.2e-16
#alternative hypothesis: true ratio of variances is not equal to 1
#95 percent confidence interval:
# 0.7441895 0.8298320
#sample estimates:
#ratio of variances 
#         0.7858488 


var.test(ce[ce$chr_str_type == "arms",]$kimura_distance,ce[ce$chr_str_type == "centers",]$kimura_distance)

#	F test to compare two variances
#
#data:  ce[ce$chr_str_type == "arms", ]$kimura_distance and ce[ce$chr_str_type == "centers", ]$kimura_distance
#F = 0.57048, num df = 4981, denom df = 4672, p-value < 2.2e-16
#alternative hypothesis: true ratio of variances is not equal to 1
#95 percent confidence interval:
# 0.5391465 0.6035966
#sample estimates:
#ratio of variances 
#         0.5704768 
#


var.test(ci[ci$chr_str_type == "arms",]$kimura_distance,ci[ci$chr_str_type == "centers",]$kimura_distance)

#	F test to compare two variances
#
#data:  ci[ci$chr_str_type == "arms", ]$kimura_distance and ci[ci$chr_str_type == "centers", ]$kimura_distance
#F = 0.75868, num df = 5940, denom df = 5889, p-value < 2.2e-16
#alternative hypothesis: true ratio of variances is not equal to 1
#95 percent confidence interval:
# 0.7209687 0.7983504
#sample estimates:
#ratio of variances 
#         0.7586758 
#



var.test(cn[cn$chr_str_type == "arms",]$kimura_distance,cn[cn$chr_str_type == "centers",]$kimura_distance)

#	F test to compare two variances
#
#data:  cn[cn$chr_str_type == "arms", ]$kimura_distance and cn[cn$chr_str_type == "centers", ]$kimura_distance
#F = 0.93795, num df = 5811, denom df = 5746, p-value = 0.01491
#alternative hypothesis: true ratio of variances is not equal to 1
#95 percent confidence interval:
# 0.8908006 0.9875885
#sample estimates:
#ratio of variances 
#         0.9379505 



var.test(cr[cr$chr_str_type == "arms",]$kimura_distance,cr[cr$chr_str_type == "centers",]$kimura_distance)




#> var.test(cr[cr$chr_str_type == "arms",]$kimura_distance,cr[cr$chr_str_type == "centers",]$kimura_distance)
#
#	F test to compare two variances
#
#data:  cr[cr$chr_str_type == "arms", ]$kimura_distance and cr[cr$chr_str_type == "centers", ]$kimura_distance
#F = 0.46914, num df = 6171, denom df = 4916, p-value < 2.2e-16
#alternative hypothesis: true ratio of variances is not equal to 1
#95 percent confidence interval:
# 0.4448883 0.4946312
#sample estimates:
#ratio of variances 
#         0.4691387 




#clusters



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


summary(rep_inopinata_no_na_superfamily$insertion_count)


summary(lm(log(insertion_count+1) ~ log(average+1),data=rep_inopinata_no_na_superfamily))

#Call:
#lm(formula = log(insertion_count + 1) ~ log(average + 1), data = rep_inopinata_no_na_superfamily)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-2.9985 -1.1073 -0.3556  0.7873  6.0936 
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       1.04870    0.08498   12.34   <2e-16 ***
#log(average + 1)  0.69514    0.03533   19.68   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 1.556 on 1903 degrees of freedom
#Multiple R-squared:  0.1691,	Adjusted R-squared:  0.1686 
#F-statistic: 387.2 on 1 and 1903 DF,  p-value: < 2.2e-16

#do this again, but exclude clusters with less than 30 insertions


ino_greater_than_30_insertions <- rep_inopinata_no_na_superfamily[rep_inopinata_no_na_superfamily$insertion_count > 29,]


summary(lm(insertion_count ~ average,data=ino_greater_than_30_insertions))


#> summary(lm(insertion_count ~ average,data=ino_greater_than_30_insertions))
#
#Call:
#lm(formula = insertion_count ~ average, data = ino_greater_than_30_insertions)
#
#Residuals:
#   Min     1Q Median     3Q    Max 
#-289.3 -257.7 -218.0  -72.3 8293.7 
#
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  324.945     71.900   4.519 7.63e-06 ***
#average       -1.397      4.825  -0.289    0.772    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 661.1 on 539 degrees of freedom
#Multiple R-squared:  0.0001554,	Adjusted R-squared:  -0.0017 
#F-statistic: 0.08377 on 1 and 539 DF,  p-value: 0.7724
#


summary(lm(log(insertion_count+1) ~ log(average+1),data=ino_greater_than_30_insertions))

#> summary(lm(log(insertion_count+1) ~ log(average+1),data=ino_greater_than_30_insertions))
#
#Call:
#lm(formula = log(insertion_count + 1) ~ log(average + 1), data = ino_greater_than_30_insertions)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-1.4025 -0.9506 -0.3500  0.6701  4.2489 
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       4.68569    0.28164  16.637   <2e-16 ***
#log(average + 1)  0.04525    0.10699   0.423    0.672    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 1.185 on 539 degrees of freedom
#Multiple R-squared:  0.0003318,	Adjusted R-squared:  -0.001523 
#F-statistic: 0.1789 on 1 and 539 DF,  p-value: 0.6725
#

summary(lm(log(insertion_count+1) ~ log(average+1),data=ino_greater_than_30_insertions))

# > summary(lm(log(insertion_count+1) ~ log(average+1),data=ino_greater_than_30_insertions))
# 
# Call:
# lm(formula = log(insertion_count + 1) ~ log(average + 1), data = ino_greater_than_30_insertions)
# 
# Residuals:
#     Min      1Q  Median      3Q     Max 
# -1.4025 -0.9506 -0.3500  0.6701  4.2489 
# 
# Coefficients:
#                  Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       4.68569    0.28164  16.637   <2e-16 ***
# log(average + 1)  0.04525    0.10699   0.423    0.672    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.185 on 539 degrees of freedom
# Multiple R-squared:  0.0003318,	Adjusted R-squared:  -0.001523 
# F-statistic: 0.1789 on 1 and 539 DF,  p-value: 0.6725


#tukey tests


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


#class

cb_class <- cb[cb$rep_taxononmic_rank == 'class',]
ce_class <- ce[ce$rep_taxononmic_rank == 'class',]
ci_class <- ci[ci$rep_taxononmic_rank == 'class',]
cn_class <- cn[cn$rep_taxononmic_rank == 'class',]
cr_class <- cr[cr$rep_taxononmic_rank == 'class',]

#remove unused factor levels

cb_class$rep_type <- droplevels(cb_class$rep_type)
ce_class$rep_type <- droplevels(ce_class$rep_type)
ci_class$rep_type <- droplevels(ci_class$rep_type)
cn_class$rep_type <- droplevels(cn_class$rep_type)
cr_class$rep_type <- droplevels(cr_class$rep_type)



#brigggsae
#exclude classes with <15 observations for tukey test
rep_tbl <- table(cb_class$rep_type)
cb_class_tukey <- droplevels(cb_class[cb_class$rep_type %in% names(rep_tbl)[rep_tbl >= 15],,drop=FALSE])

#cb tukey test class

data.mctp <- mctp(kimura_distance ~ rep_type, data = cb_class_tukey, type = "Tukey", conf.level = 0.95, asy.method = "fisher", info = FALSE)
summary(data.mctp)


#elegans
rep_tbl <- table(ce_class$rep_type)
ce_class_tukey <- droplevels(ce_class[ce_class$rep_type %in% names(rep_tbl)[rep_tbl >= 15],,drop=FALSE])

#cb tukey test class

data.mctp <- mctp(kimura_distance ~ rep_type, data = ce_class_tukey, type = "Tukey", conf.level = 0.95, asy.method = "fisher", info = FALSE)
summary(data.mctp)



#inopinata


rep_tbl <- table(ci_class$rep_type)
ci_class_tukey <- droplevels(ci_class[ci_class$rep_type %in% names(rep_tbl)[rep_tbl >= 15],,drop=FALSE])

#cb tukey test class

data.mctp <- mctp(kimura_distance ~ rep_type, data = ci_class_tukey, type = "Tukey", conf.level = 0.95, asy.method = "fisher", info = FALSE)
summary(data.mctp)


#nigoni

rep_tbl <- table(cn_class$rep_type)
cn_class_tukey <- droplevels(cn_class[cn_class$rep_type %in% names(rep_tbl)[rep_tbl >= 15],,drop=FALSE])

#cb tukey test class

data.mctp <- mctp(kimura_distance ~ rep_type, data = cn_class_tukey, type = "Tukey", conf.level = 0.95, asy.method = "fisher", info = FALSE)
summary(data.mctp)


#remanei

rep_tbl <- table(cr_class$rep_type)
cr_class_tukey <- droplevels(cr_class[cr_class$rep_type %in% names(rep_tbl)[rep_tbl >= 15],,drop=FALSE])

#cb tukey test class

data.mctp <- mctp(kimura_distance ~ rep_type, data = cr_class_tukey, type = "Tukey", conf.level = 0.95, asy.method = "fisher", info = FALSE)
summary(data.mctp)






#superfamily

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


#brigggsae
#exclude classes with <15 observations for tukey test
rep_tbl <- table(cb_superfamily$rep_type)
cb_superfamily_tukey <- droplevels(cb_superfamily[cb_superfamily$rep_type %in% names(rep_tbl)[rep_tbl >= 15],,drop=FALSE])

#cb tukey test superfamily

data.mctp <- mctp(kimura_distance ~ rep_type, data = cb_superfamily_tukey, type = "Tukey", conf.level = 0.95, asy.method = "fisher", info = FALSE)
summary(data.mctp)


#elegans
rep_tbl <- table(ce_superfamily$rep_type)
ce_superfamily_tukey <- droplevels(ce_superfamily[ce_superfamily$rep_type %in% names(rep_tbl)[rep_tbl >= 15],,drop=FALSE])

#cb tukey test superfamily

data.mctp <- mctp(kimura_distance ~ rep_type, data = ce_superfamily_tukey, type = "Tukey", conf.level = 0.95, asy.method = "fisher", info = FALSE)
summary(data.mctp)



#inopinata


rep_tbl <- table(ci_superfamily$rep_type)
ci_superfamily_tukey <- droplevels(ci_superfamily[ci_superfamily$rep_type %in% names(rep_tbl)[rep_tbl >= 15],,drop=FALSE])

#ci tukey test superfamily

data.mctp <- mctp(kimura_distance ~ rep_type, data = ci_superfamily_tukey, type = "Tukey", conf.level = 0.95, asy.method = "fisher", info = FALSE)
summary(data.mctp)


#nigoni

rep_tbl <- table(cn_superfamily$rep_type)
cn_superfamily_tukey <- droplevels(cn_superfamily[cn_superfamily$rep_type %in% names(rep_tbl)[rep_tbl >= 15],,drop=FALSE])

#cn tukey test superfamily

data.mctp <- mctp(kimura_distance ~ rep_type, data = cn_superfamily_tukey, type = "Tukey", conf.level = 0.95, asy.method = "fisher", info = FALSE)
summary(data.mctp)


#remanei

rep_tbl <- table(cr_superfamily$rep_type)
cr_superfamily_tukey <- droplevels(cr_superfamily[cr_superfamily$rep_type %in% names(rep_tbl)[rep_tbl >= 15],,drop=FALSE])

#cr tukey test superfamily

data.mctp <- mctp(kimura_distance ~ rep_type, data = cr_superfamily_tukey, type = "Tukey", conf.level = 0.95, asy.method = "fisher", info = FALSE)
summary(data.mctp)



#tests for multimodality

library(silvermantest)

#classes

#briggsae

cb_class_k_1 <- aggregate( cb_class_tukey$kimura_distance ~ cb_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=1,R=999, adjust = FALSE, digits = 6)@p_value)
cb_class_k_2 <- aggregate( cb_class_tukey$kimura_distance ~ cb_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=2,R=999, adjust = FALSE, digits = 6)@p_value)
cb_class_k_3 <- aggregate( cb_class_tukey$kimura_distance ~ cb_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=3,R=999, adjust = FALSE, digits = 6)@p_value)
cb_class_k_4 <- aggregate( cb_class_tukey$kimura_distance ~ cb_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=4,R=999, adjust = FALSE, digits = 6)@p_value)
cb_class_k_5 <- aggregate( cb_class_tukey$kimura_distance ~ cb_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=5,R=999, adjust = FALSE, digits = 6)@p_value)
cb_class_k_6 <- aggregate( cb_class_tukey$kimura_distance ~ cb_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=6,R=999, adjust = FALSE, digits = 6)@p_value)

colnames(cb_class_k_1)[colnames(cb_class_k_1)=="cb_class_tukey$rep_type"] <- "rep_type"
colnames(cb_class_k_1)[colnames(cb_class_k_1)=="cb_class_tukey$kimura_distance"] <- "p_1"

colnames(cb_class_k_2)[colnames(cb_class_k_2)=="cb_class_tukey$rep_type"] <- "rep_type"
colnames(cb_class_k_2)[colnames(cb_class_k_2)=="cb_class_tukey$kimura_distance"] <- "p_2"

colnames(cb_class_k_3)[colnames(cb_class_k_3)=="cb_class_tukey$rep_type"] <- "rep_type"
colnames(cb_class_k_3)[colnames(cb_class_k_3)=="cb_class_tukey$kimura_distance"] <- "p_3"


colnames(cb_class_k_4)[colnames(cb_class_k_4)=="cb_class_tukey$rep_type"] <- "rep_type"
colnames(cb_class_k_4)[colnames(cb_class_k_4)=="cb_class_tukey$kimura_distance"] <- "p_4"

colnames(cb_class_k_5)[colnames(cb_class_k_5)=="cb_class_tukey$rep_type"] <- "rep_type"
colnames(cb_class_k_5)[colnames(cb_class_k_5)=="cb_class_tukey$kimura_distance"] <- "p_5"

colnames(cb_class_k_6)[colnames(cb_class_k_6)=="cb_class_tukey$rep_type"] <- "rep_type"
colnames(cb_class_k_6)[colnames(cb_class_k_6)=="cb_class_tukey$kimura_distance"] <- "p_6"

cb_class_silverman <- data.frame(rep_type = cb_class_k_1$rep_type, p_1 = cb_class_k_1$p_1, p_2 = cb_class_k_2$p_2, p_3 = cb_class_k_3$p_3, p_4 = cb_class_k_4$p_4, p_5 = cb_class_k_5$p_5, p_6 = cb_class_k_6$p_6)

cb_class_silverman$p1_sig <- cb_class_silverman$p_1 < 0.05
cb_class_silverman$p2_sig <- cb_class_silverman$p_2 < 0.05
cb_class_silverman$p3_sig <- cb_class_silverman$p_3 < 0.05
cb_class_silverman$p4_sig <- cb_class_silverman$p_4 < 0.05
cb_class_silverman$p5_sig <- cb_class_silverman$p_5 < 0.05
cb_class_silverman$p6_sig <- cb_class_silverman$p_6 < 0.05

cb_class_silverman$one_mode <- ifelse(cb_class_silverman$p1_sig == FALSE,"yes","no")
cb_class_silverman$two_modes <- ifelse(cb_class_silverman$p1_sig == TRUE & cb_class_silverman$p2_sig == FALSE,"yes","no")
cb_class_silverman$three_modes <- ifelse(cb_class_silverman$p1_sig == TRUE & cb_class_silverman$p2_sig == TRUE & cb_class_silverman$p3_sig == FALSE,"yes","no")


cb_class_silverman$four_modes <- ifelse(cb_class_silverman$p1_sig == TRUE & cb_class_silverman$p2_sig == TRUE & cb_class_silverman$p3_sig == TRUE & cb_class_silverman$p4_sig == FALSE, "yes","no")

cb_class_silverman$five_modes <- ifelse(cb_class_silverman$p1_sig == TRUE & cb_class_silverman$p2_sig == TRUE & cb_class_silverman$p3_sig == TRUE & cb_class_silverman$p4_sig == TRUE & cb_class_silverman$p5_sig == FALSE, "yes","no")

cb_class_silverman$six_modes <- ifelse(cb_class_silverman$p1_sig == TRUE & cb_class_silverman$p2_sig == TRUE & cb_class_silverman$p3_sig == TRUE & cb_class_silverman$p4_sig == TRUE & cb_class_silverman$p5_sig == TRUE  & cb_class_silverman$p6_sig == FALSE, "yes","no")

cb_class_silverman$more_modes <- ifelse(cb_class_silverman$p1_sig == TRUE & cb_class_silverman$p2_sig == TRUE & cb_class_silverman$p3_sig == TRUE & cb_class_silverman$p4_sig == TRUE & cb_class_silverman$p5_sig == TRUE  & cb_class_silverman$p6_sig == TRUE, "yes","no")


write.table(cb_class_silverman, file = "briggsae_classes",row.names = FALSE, quote = FALSE, sep = "\t")



#elegans class


ce_class_k_1 <- aggregate( ce_class_tukey$kimura_distance ~ ce_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=1,R=999, adjust = FALSE, digits = 6)@p_value)
ce_class_k_2 <- aggregate( ce_class_tukey$kimura_distance ~ ce_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=2,R=999, adjust = FALSE, digits = 6)@p_value)
ce_class_k_3 <- aggregate( ce_class_tukey$kimura_distance ~ ce_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=3,R=999, adjust = FALSE, digits = 6)@p_value)
ce_class_k_4 <- aggregate( ce_class_tukey$kimura_distance ~ ce_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=4,R=999, adjust = FALSE, digits = 6)@p_value)
ce_class_k_5 <- aggregate( ce_class_tukey$kimura_distance ~ ce_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=5,R=999, adjust = FALSE, digits = 6)@p_value)
ce_class_k_6 <- aggregate( ce_class_tukey$kimura_distance ~ ce_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=6,R=999, adjust = FALSE, digits = 6)@p_value)

colnames(ce_class_k_1)[colnames(ce_class_k_1)=="ce_class_tukey$rep_type"] <- "rep_type"
colnames(ce_class_k_1)[colnames(ce_class_k_1)=="ce_class_tukey$kimura_distance"] <- "p_1"

colnames(ce_class_k_2)[colnames(ce_class_k_2)=="ce_class_tukey$rep_type"] <- "rep_type"
colnames(ce_class_k_2)[colnames(ce_class_k_2)=="ce_class_tukey$kimura_distance"] <- "p_2"

colnames(ce_class_k_3)[colnames(ce_class_k_3)=="ce_class_tukey$rep_type"] <- "rep_type"
colnames(ce_class_k_3)[colnames(ce_class_k_3)=="ce_class_tukey$kimura_distance"] <- "p_3"


colnames(ce_class_k_4)[colnames(ce_class_k_4)=="ce_class_tukey$rep_type"] <- "rep_type"
colnames(ce_class_k_4)[colnames(ce_class_k_4)=="ce_class_tukey$kimura_distance"] <- "p_4"

colnames(ce_class_k_5)[colnames(ce_class_k_5)=="ce_class_tukey$rep_type"] <- "rep_type"
colnames(ce_class_k_5)[colnames(ce_class_k_5)=="ce_class_tukey$kimura_distance"] <- "p_5"

colnames(ce_class_k_6)[colnames(ce_class_k_6)=="ce_class_tukey$rep_type"] <- "rep_type"
colnames(ce_class_k_6)[colnames(ce_class_k_6)=="ce_class_tukey$kimura_distance"] <- "p_6"

ce_class_silverman <- data.frame(rep_type = ce_class_k_1$rep_type, p_1 = ce_class_k_1$p_1, p_2 = ce_class_k_2$p_2, p_3 = ce_class_k_3$p_3, p_4 = ce_class_k_4$p_4, p_5 = ce_class_k_5$p_5, p_6 = ce_class_k_6$p_6)

ce_class_silverman$p1_sig <- ce_class_silverman$p_1 < 0.05
ce_class_silverman$p2_sig <- ce_class_silverman$p_2 < 0.05
ce_class_silverman$p3_sig <- ce_class_silverman$p_3 < 0.05
ce_class_silverman$p4_sig <- ce_class_silverman$p_4 < 0.05
ce_class_silverman$p5_sig <- ce_class_silverman$p_5 < 0.05
ce_class_silverman$p6_sig <- ce_class_silverman$p_6 < 0.05

ce_class_silverman$one_mode <- ifelse(ce_class_silverman$p1_sig == FALSE,"yes","no")
ce_class_silverman$two_modes <- ifelse(ce_class_silverman$p1_sig == TRUE & ce_class_silverman$p2_sig == FALSE,"yes","no")
ce_class_silverman$three_modes <- ifelse(ce_class_silverman$p1_sig == TRUE & ce_class_silverman$p2_sig == TRUE & ce_class_silverman$p3_sig == FALSE,"yes","no")


ce_class_silverman$four_modes <- ifelse(ce_class_silverman$p1_sig == TRUE & ce_class_silverman$p2_sig == TRUE & ce_class_silverman$p3_sig == TRUE & ce_class_silverman$p4_sig == FALSE, "yes","no")

ce_class_silverman$five_modes <- ifelse(ce_class_silverman$p1_sig == TRUE & ce_class_silverman$p2_sig == TRUE & ce_class_silverman$p3_sig == TRUE & ce_class_silverman$p4_sig == TRUE & ce_class_silverman$p5_sig == FALSE, "yes","no")

ce_class_silverman$six_modes <- ifelse(ce_class_silverman$p1_sig == TRUE & ce_class_silverman$p2_sig == TRUE & ce_class_silverman$p3_sig == TRUE & ce_class_silverman$p4_sig == TRUE & ce_class_silverman$p5_sig == TRUE  & ce_class_silverman$p6_sig == FALSE, "yes","no")

ce_class_silverman$more_modes <- ifelse(ce_class_silverman$p1_sig == TRUE & ce_class_silverman$p2_sig == TRUE & ce_class_silverman$p3_sig == TRUE & ce_class_silverman$p4_sig == TRUE & ce_class_silverman$p5_sig == TRUE  & ce_class_silverman$p6_sig == TRUE, "yes","no")


write.table(ce_class_silverman, file = "elegans_classes",row.names = FALSE, quote = FALSE, sep = "\t")


#inopinata



ci_class_k_1 <- aggregate( ci_class_tukey$kimura_distance ~ ci_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=1,R=999, adjust = FALSE, digits = 6)@p_value)
ci_class_k_2 <- aggregate( ci_class_tukey$kimura_distance ~ ci_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=2,R=999, adjust = FALSE, digits = 6)@p_value)
ci_class_k_3 <- aggregate( ci_class_tukey$kimura_distance ~ ci_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=3,R=999, adjust = FALSE, digits = 6)@p_value)
ci_class_k_4 <- aggregate( ci_class_tukey$kimura_distance ~ ci_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=4,R=999, adjust = FALSE, digits = 6)@p_value)
ci_class_k_5 <- aggregate( ci_class_tukey$kimura_distance ~ ci_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=5,R=999, adjust = FALSE, digits = 6)@p_value)
ci_class_k_6 <- aggregate( ci_class_tukey$kimura_distance ~ ci_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=6,R=999, adjust = FALSE, digits = 6)@p_value)

colnames(ci_class_k_1)[colnames(ci_class_k_1)=="ci_class_tukey$rep_type"] <- "rep_type"
colnames(ci_class_k_1)[colnames(ci_class_k_1)=="ci_class_tukey$kimura_distance"] <- "p_1"

colnames(ci_class_k_2)[colnames(ci_class_k_2)=="ci_class_tukey$rep_type"] <- "rep_type"
colnames(ci_class_k_2)[colnames(ci_class_k_2)=="ci_class_tukey$kimura_distance"] <- "p_2"

colnames(ci_class_k_3)[colnames(ci_class_k_3)=="ci_class_tukey$rep_type"] <- "rep_type"
colnames(ci_class_k_3)[colnames(ci_class_k_3)=="ci_class_tukey$kimura_distance"] <- "p_3"


colnames(ci_class_k_4)[colnames(ci_class_k_4)=="ci_class_tukey$rep_type"] <- "rep_type"
colnames(ci_class_k_4)[colnames(ci_class_k_4)=="ci_class_tukey$kimura_distance"] <- "p_4"

colnames(ci_class_k_5)[colnames(ci_class_k_5)=="ci_class_tukey$rep_type"] <- "rep_type"
colnames(ci_class_k_5)[colnames(ci_class_k_5)=="ci_class_tukey$kimura_distance"] <- "p_5"

colnames(ci_class_k_6)[colnames(ci_class_k_6)=="ci_class_tukey$rep_type"] <- "rep_type"
colnames(ci_class_k_6)[colnames(ci_class_k_6)=="ci_class_tukey$kimura_distance"] <- "p_6"

ci_class_silverman <- data.frame(rep_type = ci_class_k_1$rep_type, p_1 = ci_class_k_1$p_1, p_2 = ci_class_k_2$p_2, p_3 = ci_class_k_3$p_3, p_4 = ci_class_k_4$p_4, p_5 = ci_class_k_5$p_5, p_6 = ci_class_k_6$p_6)

ci_class_silverman$p1_sig <- ci_class_silverman$p_1 < 0.05
ci_class_silverman$p2_sig <- ci_class_silverman$p_2 < 0.05
ci_class_silverman$p3_sig <- ci_class_silverman$p_3 < 0.05
ci_class_silverman$p4_sig <- ci_class_silverman$p_4 < 0.05
ci_class_silverman$p5_sig <- ci_class_silverman$p_5 < 0.05
ci_class_silverman$p6_sig <- ci_class_silverman$p_6 < 0.05

ci_class_silverman$one_mode <- ifelse(ci_class_silverman$p1_sig == FALSE,"yes","no")
ci_class_silverman$two_modes <- ifelse(ci_class_silverman$p1_sig == TRUE & ci_class_silverman$p2_sig == FALSE,"yes","no")
ci_class_silverman$three_modes <- ifelse(ci_class_silverman$p1_sig == TRUE & ci_class_silverman$p2_sig == TRUE & ci_class_silverman$p3_sig == FALSE,"yes","no")


ci_class_silverman$four_modes <- ifelse(ci_class_silverman$p1_sig == TRUE & ci_class_silverman$p2_sig == TRUE & ci_class_silverman$p3_sig == TRUE & ci_class_silverman$p4_sig == FALSE, "yes","no")

ci_class_silverman$five_modes <- ifelse(ci_class_silverman$p1_sig == TRUE & ci_class_silverman$p2_sig == TRUE & ci_class_silverman$p3_sig == TRUE & ci_class_silverman$p4_sig == TRUE & ci_class_silverman$p5_sig == FALSE, "yes","no")

ci_class_silverman$six_modes <- ifelse(ci_class_silverman$p1_sig == TRUE & ci_class_silverman$p2_sig == TRUE & ci_class_silverman$p3_sig == TRUE & ci_class_silverman$p4_sig == TRUE & ci_class_silverman$p5_sig == TRUE  & ci_class_silverman$p6_sig == FALSE, "yes","no")

ci_class_silverman$more_modes <- ifelse(ci_class_silverman$p1_sig == TRUE & ci_class_silverman$p2_sig == TRUE & ci_class_silverman$p3_sig == TRUE & ci_class_silverman$p4_sig == TRUE & ci_class_silverman$p5_sig == TRUE  & ci_class_silverman$p6_sig == TRUE, "yes","no")


write.table(ci_class_silverman, file = "inopinata_classes",row.names = FALSE, quote = FALSE, sep = "\t")


#nigoni



cn_class_k_1 <- aggregate( cn_class_tukey$kimura_distance ~ cn_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=1,R=999, adjust = FALSE, digits = 6)@p_value)
cn_class_k_2 <- aggregate( cn_class_tukey$kimura_distance ~ cn_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=2,R=999, adjust = FALSE, digits = 6)@p_value)
cn_class_k_3 <- aggregate( cn_class_tukey$kimura_distance ~ cn_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=3,R=999, adjust = FALSE, digits = 6)@p_value)
cn_class_k_4 <- aggregate( cn_class_tukey$kimura_distance ~ cn_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=4,R=999, adjust = FALSE, digits = 6)@p_value)
cn_class_k_5 <- aggregate( cn_class_tukey$kimura_distance ~ cn_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=5,R=999, adjust = FALSE, digits = 6)@p_value)
cn_class_k_6 <- aggregate( cn_class_tukey$kimura_distance ~ cn_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=6,R=999, adjust = FALSE, digits = 6)@p_value)

colnames(cn_class_k_1)[colnames(cn_class_k_1)=="cn_class_tukey$rep_type"] <- "rep_type"
colnames(cn_class_k_1)[colnames(cn_class_k_1)=="cn_class_tukey$kimura_distance"] <- "p_1"

colnames(cn_class_k_2)[colnames(cn_class_k_2)=="cn_class_tukey$rep_type"] <- "rep_type"
colnames(cn_class_k_2)[colnames(cn_class_k_2)=="cn_class_tukey$kimura_distance"] <- "p_2"

colnames(cn_class_k_3)[colnames(cn_class_k_3)=="cn_class_tukey$rep_type"] <- "rep_type"
colnames(cn_class_k_3)[colnames(cn_class_k_3)=="cn_class_tukey$kimura_distance"] <- "p_3"


colnames(cn_class_k_4)[colnames(cn_class_k_4)=="cn_class_tukey$rep_type"] <- "rep_type"
colnames(cn_class_k_4)[colnames(cn_class_k_4)=="cn_class_tukey$kimura_distance"] <- "p_4"

colnames(cn_class_k_5)[colnames(cn_class_k_5)=="cn_class_tukey$rep_type"] <- "rep_type"
colnames(cn_class_k_5)[colnames(cn_class_k_5)=="cn_class_tukey$kimura_distance"] <- "p_5"

colnames(cn_class_k_6)[colnames(cn_class_k_6)=="cn_class_tukey$rep_type"] <- "rep_type"
colnames(cn_class_k_6)[colnames(cn_class_k_6)=="cn_class_tukey$kimura_distance"] <- "p_6"

cn_class_silverman <- data.frame(rep_type = cn_class_k_1$rep_type, p_1 = cn_class_k_1$p_1, p_2 = cn_class_k_2$p_2, p_3 = cn_class_k_3$p_3, p_4 = cn_class_k_4$p_4, p_5 = cn_class_k_5$p_5, p_6 = cn_class_k_6$p_6)

cn_class_silverman$p1_sig <- cn_class_silverman$p_1 < 0.05
cn_class_silverman$p2_sig <- cn_class_silverman$p_2 < 0.05
cn_class_silverman$p3_sig <- cn_class_silverman$p_3 < 0.05
cn_class_silverman$p4_sig <- cn_class_silverman$p_4 < 0.05
cn_class_silverman$p5_sig <- cn_class_silverman$p_5 < 0.05
cn_class_silverman$p6_sig <- cn_class_silverman$p_6 < 0.05

cn_class_silverman$one_mode <- ifelse(cn_class_silverman$p1_sig == FALSE,"yes","no")
cn_class_silverman$two_modes <- ifelse(cn_class_silverman$p1_sig == TRUE & cn_class_silverman$p2_sig == FALSE,"yes","no")
cn_class_silverman$three_modes <- ifelse(cn_class_silverman$p1_sig == TRUE & cn_class_silverman$p2_sig == TRUE & cn_class_silverman$p3_sig == FALSE,"yes","no")


cn_class_silverman$four_modes <- ifelse(cn_class_silverman$p1_sig == TRUE & cn_class_silverman$p2_sig == TRUE & cn_class_silverman$p3_sig == TRUE & cn_class_silverman$p4_sig == FALSE, "yes","no")

cn_class_silverman$five_modes <- ifelse(cn_class_silverman$p1_sig == TRUE & cn_class_silverman$p2_sig == TRUE & cn_class_silverman$p3_sig == TRUE & cn_class_silverman$p4_sig == TRUE & cn_class_silverman$p5_sig == FALSE, "yes","no")

cn_class_silverman$six_modes <- ifelse(cn_class_silverman$p1_sig == TRUE & cn_class_silverman$p2_sig == TRUE & cn_class_silverman$p3_sig == TRUE & cn_class_silverman$p4_sig == TRUE & cn_class_silverman$p5_sig == TRUE  & cn_class_silverman$p6_sig == FALSE, "yes","no")

cn_class_silverman$more_modes <- ifelse(cn_class_silverman$p1_sig == TRUE & cn_class_silverman$p2_sig == TRUE & cn_class_silverman$p3_sig == TRUE & cn_class_silverman$p4_sig == TRUE & cn_class_silverman$p5_sig == TRUE  & cn_class_silverman$p6_sig == TRUE, "yes","no")


write.table(cn_class_silverman, file = "nigoni_classes",row.names = FALSE, quote = FALSE, sep = "\t")


#remanei classes






cr_class_k_1 <- aggregate( cr_class_tukey$kimura_distance ~ cr_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=1,R=999, adjust = FALSE, digits = 6)@p_value)
cr_class_k_2 <- aggregate( cr_class_tukey$kimura_distance ~ cr_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=2,R=999, adjust = FALSE, digits = 6)@p_value)
cr_class_k_3 <- aggregate( cr_class_tukey$kimura_distance ~ cr_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=3,R=999, adjust = FALSE, digits = 6)@p_value)
cr_class_k_4 <- aggregate( cr_class_tukey$kimura_distance ~ cr_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=4,R=999, adjust = FALSE, digits = 6)@p_value)
cr_class_k_5 <- aggregate( cr_class_tukey$kimura_distance ~ cr_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=5,R=999, adjust = FALSE, digits = 6)@p_value)
cr_class_k_6 <- aggregate( cr_class_tukey$kimura_distance ~ cr_class_tukey$rep_type, FUN = function(x) silverman.test(x,k=6,R=999, adjust = FALSE, digits = 6)@p_value)

colnames(cr_class_k_1)[colnames(cr_class_k_1)=="cr_class_tukey$rep_type"] <- "rep_type"
colnames(cr_class_k_1)[colnames(cr_class_k_1)=="cr_class_tukey$kimura_distance"] <- "p_1"

colnames(cr_class_k_2)[colnames(cr_class_k_2)=="cr_class_tukey$rep_type"] <- "rep_type"
colnames(cr_class_k_2)[colnames(cr_class_k_2)=="cr_class_tukey$kimura_distance"] <- "p_2"

colnames(cr_class_k_3)[colnames(cr_class_k_3)=="cr_class_tukey$rep_type"] <- "rep_type"
colnames(cr_class_k_3)[colnames(cr_class_k_3)=="cr_class_tukey$kimura_distance"] <- "p_3"


colnames(cr_class_k_4)[colnames(cr_class_k_4)=="cr_class_tukey$rep_type"] <- "rep_type"
colnames(cr_class_k_4)[colnames(cr_class_k_4)=="cr_class_tukey$kimura_distance"] <- "p_4"

colnames(cr_class_k_5)[colnames(cr_class_k_5)=="cr_class_tukey$rep_type"] <- "rep_type"
colnames(cr_class_k_5)[colnames(cr_class_k_5)=="cr_class_tukey$kimura_distance"] <- "p_5"

colnames(cr_class_k_6)[colnames(cr_class_k_6)=="cr_class_tukey$rep_type"] <- "rep_type"
colnames(cr_class_k_6)[colnames(cr_class_k_6)=="cr_class_tukey$kimura_distance"] <- "p_6"

cr_class_silverman <- data.frame(rep_type = cr_class_k_1$rep_type, p_1 = cr_class_k_1$p_1, p_2 = cr_class_k_2$p_2, p_3 = cr_class_k_3$p_3, p_4 = cr_class_k_4$p_4, p_5 = cr_class_k_5$p_5, p_6 = cr_class_k_6$p_6)

cr_class_silverman$p1_sig <- cr_class_silverman$p_1 < 0.05
cr_class_silverman$p2_sig <- cr_class_silverman$p_2 < 0.05
cr_class_silverman$p3_sig <- cr_class_silverman$p_3 < 0.05
cr_class_silverman$p4_sig <- cr_class_silverman$p_4 < 0.05
cr_class_silverman$p5_sig <- cr_class_silverman$p_5 < 0.05
cr_class_silverman$p6_sig <- cr_class_silverman$p_6 < 0.05

cr_class_silverman$one_mode <- ifelse(cr_class_silverman$p1_sig == FALSE,"yes","no")
cr_class_silverman$two_modes <- ifelse(cr_class_silverman$p1_sig == TRUE & cr_class_silverman$p2_sig == FALSE,"yes","no")
cr_class_silverman$three_modes <- ifelse(cr_class_silverman$p1_sig == TRUE & cr_class_silverman$p2_sig == TRUE & cr_class_silverman$p3_sig == FALSE,"yes","no")


cr_class_silverman$four_modes <- ifelse(cr_class_silverman$p1_sig == TRUE & cr_class_silverman$p2_sig == TRUE & cr_class_silverman$p3_sig == TRUE & cr_class_silverman$p4_sig == FALSE, "yes","no")

cr_class_silverman$five_modes <- ifelse(cr_class_silverman$p1_sig == TRUE & cr_class_silverman$p2_sig == TRUE & cr_class_silverman$p3_sig == TRUE & cr_class_silverman$p4_sig == TRUE & cr_class_silverman$p5_sig == FALSE, "yes","no")

cr_class_silverman$six_modes <- ifelse(cr_class_silverman$p1_sig == TRUE & cr_class_silverman$p2_sig == TRUE & cr_class_silverman$p3_sig == TRUE & cr_class_silverman$p4_sig == TRUE & cr_class_silverman$p5_sig == TRUE  & cr_class_silverman$p6_sig == FALSE, "yes","no")

cr_class_silverman$more_modes <- ifelse(cr_class_silverman$p1_sig == TRUE & cr_class_silverman$p2_sig == TRUE & cr_class_silverman$p3_sig == TRUE & cr_class_silverman$p4_sig == TRUE & cr_class_silverman$p5_sig == TRUE  & cr_class_silverman$p6_sig == TRUE, "yes","no")


write.table(cr_class_silverman, file = "remanei_classes",row.names = FALSE, quote = FALSE, sep = "\t")





#superfamilies age multimodality





cb_superfamily_k_1 <- aggregate( cb_superfamily_tukey$kimura_distance ~ cb_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=1,R=999, adjust = FALSE, digits = 6)@p_value)
cb_superfamily_k_2 <- aggregate( cb_superfamily_tukey$kimura_distance ~ cb_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=2,R=999, adjust = FALSE, digits = 6)@p_value)
cb_superfamily_k_3 <- aggregate( cb_superfamily_tukey$kimura_distance ~ cb_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=3,R=999, adjust = FALSE, digits = 6)@p_value)
cb_superfamily_k_4 <- aggregate( cb_superfamily_tukey$kimura_distance ~ cb_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=4,R=999, adjust = FALSE, digits = 6)@p_value)
cb_superfamily_k_5 <- aggregate( cb_superfamily_tukey$kimura_distance ~ cb_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=5,R=999, adjust = FALSE, digits = 6)@p_value)
cb_superfamily_k_6 <- aggregate( cb_superfamily_tukey$kimura_distance ~ cb_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=6,R=999, adjust = FALSE, digits = 6)@p_value)

colnames(cb_superfamily_k_1)[colnames(cb_superfamily_k_1)=="cb_superfamily_tukey$rep_type"] <- "rep_type"
colnames(cb_superfamily_k_1)[colnames(cb_superfamily_k_1)=="cb_superfamily_tukey$kimura_distance"] <- "p_1"

colnames(cb_superfamily_k_2)[colnames(cb_superfamily_k_2)=="cb_superfamily_tukey$rep_type"] <- "rep_type"
colnames(cb_superfamily_k_2)[colnames(cb_superfamily_k_2)=="cb_superfamily_tukey$kimura_distance"] <- "p_2"

colnames(cb_superfamily_k_3)[colnames(cb_superfamily_k_3)=="cb_superfamily_tukey$rep_type"] <- "rep_type"
colnames(cb_superfamily_k_3)[colnames(cb_superfamily_k_3)=="cb_superfamily_tukey$kimura_distance"] <- "p_3"


colnames(cb_superfamily_k_4)[colnames(cb_superfamily_k_4)=="cb_superfamily_tukey$rep_type"] <- "rep_type"
colnames(cb_superfamily_k_4)[colnames(cb_superfamily_k_4)=="cb_superfamily_tukey$kimura_distance"] <- "p_4"

colnames(cb_superfamily_k_5)[colnames(cb_superfamily_k_5)=="cb_superfamily_tukey$rep_type"] <- "rep_type"
colnames(cb_superfamily_k_5)[colnames(cb_superfamily_k_5)=="cb_superfamily_tukey$kimura_distance"] <- "p_5"

colnames(cb_superfamily_k_6)[colnames(cb_superfamily_k_6)=="cb_superfamily_tukey$rep_type"] <- "rep_type"
colnames(cb_superfamily_k_6)[colnames(cb_superfamily_k_6)=="cb_superfamily_tukey$kimura_distance"] <- "p_6"

cb_superfamily_silverman <- data.frame(rep_type = cb_superfamily_k_1$rep_type, p_1 = cb_superfamily_k_1$p_1, p_2 = cb_superfamily_k_2$p_2, p_3 = cb_superfamily_k_3$p_3, p_4 = cb_superfamily_k_4$p_4, p_5 = cb_superfamily_k_5$p_5, p_6 = cb_superfamily_k_6$p_6)

cb_superfamily_silverman$p1_sig <- cb_superfamily_silverman$p_1 < 0.05
cb_superfamily_silverman$p2_sig <- cb_superfamily_silverman$p_2 < 0.05
cb_superfamily_silverman$p3_sig <- cb_superfamily_silverman$p_3 < 0.05
cb_superfamily_silverman$p4_sig <- cb_superfamily_silverman$p_4 < 0.05
cb_superfamily_silverman$p5_sig <- cb_superfamily_silverman$p_5 < 0.05
cb_superfamily_silverman$p6_sig <- cb_superfamily_silverman$p_6 < 0.05

cb_superfamily_silverman$one_mode <- ifelse(cb_superfamily_silverman$p1_sig == FALSE,"yes","no")
cb_superfamily_silverman$two_modes <- ifelse(cb_superfamily_silverman$p1_sig == TRUE & cb_superfamily_silverman$p2_sig == FALSE,"yes","no")
cb_superfamily_silverman$three_modes <- ifelse(cb_superfamily_silverman$p1_sig == TRUE & cb_superfamily_silverman$p2_sig == TRUE & cb_superfamily_silverman$p3_sig == FALSE,"yes","no")


cb_superfamily_silverman$four_modes <- ifelse(cb_superfamily_silverman$p1_sig == TRUE & cb_superfamily_silverman$p2_sig == TRUE & cb_superfamily_silverman$p3_sig == TRUE & cb_superfamily_silverman$p4_sig == FALSE, "yes","no")

cb_superfamily_silverman$five_modes <- ifelse(cb_superfamily_silverman$p1_sig == TRUE & cb_superfamily_silverman$p2_sig == TRUE & cb_superfamily_silverman$p3_sig == TRUE & cb_superfamily_silverman$p4_sig == TRUE & cb_superfamily_silverman$p5_sig == FALSE, "yes","no")

cb_superfamily_silverman$six_modes <- ifelse(cb_superfamily_silverman$p1_sig == TRUE & cb_superfamily_silverman$p2_sig == TRUE & cb_superfamily_silverman$p3_sig == TRUE & cb_superfamily_silverman$p4_sig == TRUE & cb_superfamily_silverman$p5_sig == TRUE  & cb_superfamily_silverman$p6_sig == FALSE, "yes","no")

cb_superfamily_silverman$more_modes <- ifelse(cb_superfamily_silverman$p1_sig == TRUE & cb_superfamily_silverman$p2_sig == TRUE & cb_superfamily_silverman$p3_sig == TRUE & cb_superfamily_silverman$p4_sig == TRUE & cb_superfamily_silverman$p5_sig == TRUE  & cb_superfamily_silverman$p6_sig == TRUE, "yes","no")


write.table(cb_superfamily_silverman, file = "briggsae_superfamilies",row.names = FALSE, quote = FALSE, sep = "\t")


#elegans superfamily





ce_superfamily_k_1 <- aggregate( ce_superfamily_tukey$kimura_distance ~ ce_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=1,R=999, adjust = FALSE, digits = 6)@p_value)
ce_superfamily_k_2 <- aggregate( ce_superfamily_tukey$kimura_distance ~ ce_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=2,R=999, adjust = FALSE, digits = 6)@p_value)
ce_superfamily_k_3 <- aggregate( ce_superfamily_tukey$kimura_distance ~ ce_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=3,R=999, adjust = FALSE, digits = 6)@p_value)
ce_superfamily_k_4 <- aggregate( ce_superfamily_tukey$kimura_distance ~ ce_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=4,R=999, adjust = FALSE, digits = 6)@p_value)
ce_superfamily_k_5 <- aggregate( ce_superfamily_tukey$kimura_distance ~ ce_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=5,R=999, adjust = FALSE, digits = 6)@p_value)
ce_superfamily_k_6 <- aggregate( ce_superfamily_tukey$kimura_distance ~ ce_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=6,R=999, adjust = FALSE, digits = 6)@p_value)

colnames(ce_superfamily_k_1)[colnames(ce_superfamily_k_1)=="ce_superfamily_tukey$rep_type"] <- "rep_type"
colnames(ce_superfamily_k_1)[colnames(ce_superfamily_k_1)=="ce_superfamily_tukey$kimura_distance"] <- "p_1"

colnames(ce_superfamily_k_2)[colnames(ce_superfamily_k_2)=="ce_superfamily_tukey$rep_type"] <- "rep_type"
colnames(ce_superfamily_k_2)[colnames(ce_superfamily_k_2)=="ce_superfamily_tukey$kimura_distance"] <- "p_2"

colnames(ce_superfamily_k_3)[colnames(ce_superfamily_k_3)=="ce_superfamily_tukey$rep_type"] <- "rep_type"
colnames(ce_superfamily_k_3)[colnames(ce_superfamily_k_3)=="ce_superfamily_tukey$kimura_distance"] <- "p_3"


colnames(ce_superfamily_k_4)[colnames(ce_superfamily_k_4)=="ce_superfamily_tukey$rep_type"] <- "rep_type"
colnames(ce_superfamily_k_4)[colnames(ce_superfamily_k_4)=="ce_superfamily_tukey$kimura_distance"] <- "p_4"

colnames(ce_superfamily_k_5)[colnames(ce_superfamily_k_5)=="ce_superfamily_tukey$rep_type"] <- "rep_type"
colnames(ce_superfamily_k_5)[colnames(ce_superfamily_k_5)=="ce_superfamily_tukey$kimura_distance"] <- "p_5"

colnames(ce_superfamily_k_6)[colnames(ce_superfamily_k_6)=="ce_superfamily_tukey$rep_type"] <- "rep_type"
colnames(ce_superfamily_k_6)[colnames(ce_superfamily_k_6)=="ce_superfamily_tukey$kimura_distance"] <- "p_6"

ce_superfamily_silverman <- data.frame(rep_type = ce_superfamily_k_1$rep_type, p_1 = ce_superfamily_k_1$p_1, p_2 = ce_superfamily_k_2$p_2, p_3 = ce_superfamily_k_3$p_3, p_4 = ce_superfamily_k_4$p_4, p_5 = ce_superfamily_k_5$p_5, p_6 = ce_superfamily_k_6$p_6)

ce_superfamily_silverman$p1_sig <- ce_superfamily_silverman$p_1 < 0.05
ce_superfamily_silverman$p2_sig <- ce_superfamily_silverman$p_2 < 0.05
ce_superfamily_silverman$p3_sig <- ce_superfamily_silverman$p_3 < 0.05
ce_superfamily_silverman$p4_sig <- ce_superfamily_silverman$p_4 < 0.05
ce_superfamily_silverman$p5_sig <- ce_superfamily_silverman$p_5 < 0.05
ce_superfamily_silverman$p6_sig <- ce_superfamily_silverman$p_6 < 0.05

ce_superfamily_silverman$one_mode <- ifelse(ce_superfamily_silverman$p1_sig == FALSE,"yes","no")
ce_superfamily_silverman$two_modes <- ifelse(ce_superfamily_silverman$p1_sig == TRUE & ce_superfamily_silverman$p2_sig == FALSE,"yes","no")
ce_superfamily_silverman$three_modes <- ifelse(ce_superfamily_silverman$p1_sig == TRUE & ce_superfamily_silverman$p2_sig == TRUE & ce_superfamily_silverman$p3_sig == FALSE,"yes","no")


ce_superfamily_silverman$four_modes <- ifelse(ce_superfamily_silverman$p1_sig == TRUE & ce_superfamily_silverman$p2_sig == TRUE & ce_superfamily_silverman$p3_sig == TRUE & ce_superfamily_silverman$p4_sig == FALSE, "yes","no")

ce_superfamily_silverman$five_modes <- ifelse(ce_superfamily_silverman$p1_sig == TRUE & ce_superfamily_silverman$p2_sig == TRUE & ce_superfamily_silverman$p3_sig == TRUE & ce_superfamily_silverman$p4_sig == TRUE & ce_superfamily_silverman$p5_sig == FALSE, "yes","no")

ce_superfamily_silverman$six_modes <- ifelse(ce_superfamily_silverman$p1_sig == TRUE & ce_superfamily_silverman$p2_sig == TRUE & ce_superfamily_silverman$p3_sig == TRUE & ce_superfamily_silverman$p4_sig == TRUE & ce_superfamily_silverman$p5_sig == TRUE  & ce_superfamily_silverman$p6_sig == FALSE, "yes","no")

ce_superfamily_silverman$more_modes <- ifelse(ce_superfamily_silverman$p1_sig == TRUE & ce_superfamily_silverman$p2_sig == TRUE & ce_superfamily_silverman$p3_sig == TRUE & ce_superfamily_silverman$p4_sig == TRUE & ce_superfamily_silverman$p5_sig == TRUE  & ce_superfamily_silverman$p6_sig == TRUE, "yes","no")


write.table(ce_superfamily_silverman, file = "elegans_superfamilies",row.names = FALSE, quote = FALSE, sep = "\t")


#inopinata





ci_superfamily_k_1 <- aggregate( ci_superfamily_tukey$kimura_distance ~ ci_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=1,R=999, adjust = FALSE, digits = 6)@p_value)
ci_superfamily_k_2 <- aggregate( ci_superfamily_tukey$kimura_distance ~ ci_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=2,R=999, adjust = FALSE, digits = 6)@p_value)
ci_superfamily_k_3 <- aggregate( ci_superfamily_tukey$kimura_distance ~ ci_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=3,R=999, adjust = FALSE, digits = 6)@p_value)
ci_superfamily_k_4 <- aggregate( ci_superfamily_tukey$kimura_distance ~ ci_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=4,R=999, adjust = FALSE, digits = 6)@p_value)
ci_superfamily_k_5 <- aggregate( ci_superfamily_tukey$kimura_distance ~ ci_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=5,R=999, adjust = FALSE, digits = 6)@p_value)
ci_superfamily_k_6 <- aggregate( ci_superfamily_tukey$kimura_distance ~ ci_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=6,R=999, adjust = FALSE, digits = 6)@p_value)

colnames(ci_superfamily_k_1)[colnames(ci_superfamily_k_1)=="ci_superfamily_tukey$rep_type"] <- "rep_type"
colnames(ci_superfamily_k_1)[colnames(ci_superfamily_k_1)=="ci_superfamily_tukey$kimura_distance"] <- "p_1"

colnames(ci_superfamily_k_2)[colnames(ci_superfamily_k_2)=="ci_superfamily_tukey$rep_type"] <- "rep_type"
colnames(ci_superfamily_k_2)[colnames(ci_superfamily_k_2)=="ci_superfamily_tukey$kimura_distance"] <- "p_2"

colnames(ci_superfamily_k_3)[colnames(ci_superfamily_k_3)=="ci_superfamily_tukey$rep_type"] <- "rep_type"
colnames(ci_superfamily_k_3)[colnames(ci_superfamily_k_3)=="ci_superfamily_tukey$kimura_distance"] <- "p_3"


colnames(ci_superfamily_k_4)[colnames(ci_superfamily_k_4)=="ci_superfamily_tukey$rep_type"] <- "rep_type"
colnames(ci_superfamily_k_4)[colnames(ci_superfamily_k_4)=="ci_superfamily_tukey$kimura_distance"] <- "p_4"

colnames(ci_superfamily_k_5)[colnames(ci_superfamily_k_5)=="ci_superfamily_tukey$rep_type"] <- "rep_type"
colnames(ci_superfamily_k_5)[colnames(ci_superfamily_k_5)=="ci_superfamily_tukey$kimura_distance"] <- "p_5"

colnames(ci_superfamily_k_6)[colnames(ci_superfamily_k_6)=="ci_superfamily_tukey$rep_type"] <- "rep_type"
colnames(ci_superfamily_k_6)[colnames(ci_superfamily_k_6)=="ci_superfamily_tukey$kimura_distance"] <- "p_6"

ci_superfamily_silverman <- data.frame(rep_type = ci_superfamily_k_1$rep_type, p_1 = ci_superfamily_k_1$p_1, p_2 = ci_superfamily_k_2$p_2, p_3 = ci_superfamily_k_3$p_3, p_4 = ci_superfamily_k_4$p_4, p_5 = ci_superfamily_k_5$p_5, p_6 = ci_superfamily_k_6$p_6)

ci_superfamily_silverman$p1_sig <- ci_superfamily_silverman$p_1 < 0.05
ci_superfamily_silverman$p2_sig <- ci_superfamily_silverman$p_2 < 0.05
ci_superfamily_silverman$p3_sig <- ci_superfamily_silverman$p_3 < 0.05
ci_superfamily_silverman$p4_sig <- ci_superfamily_silverman$p_4 < 0.05
ci_superfamily_silverman$p5_sig <- ci_superfamily_silverman$p_5 < 0.05
ci_superfamily_silverman$p6_sig <- ci_superfamily_silverman$p_6 < 0.05

ci_superfamily_silverman$one_mode <- ifelse(ci_superfamily_silverman$p1_sig == FALSE,"yes","no")
ci_superfamily_silverman$two_modes <- ifelse(ci_superfamily_silverman$p1_sig == TRUE & ci_superfamily_silverman$p2_sig == FALSE,"yes","no")
ci_superfamily_silverman$three_modes <- ifelse(ci_superfamily_silverman$p1_sig == TRUE & ci_superfamily_silverman$p2_sig == TRUE & ci_superfamily_silverman$p3_sig == FALSE,"yes","no")


ci_superfamily_silverman$four_modes <- ifelse(ci_superfamily_silverman$p1_sig == TRUE & ci_superfamily_silverman$p2_sig == TRUE & ci_superfamily_silverman$p3_sig == TRUE & ci_superfamily_silverman$p4_sig == FALSE, "yes","no")

ci_superfamily_silverman$five_modes <- ifelse(ci_superfamily_silverman$p1_sig == TRUE & ci_superfamily_silverman$p2_sig == TRUE & ci_superfamily_silverman$p3_sig == TRUE & ci_superfamily_silverman$p4_sig == TRUE & ci_superfamily_silverman$p5_sig == FALSE, "yes","no")

ci_superfamily_silverman$six_modes <- ifelse(ci_superfamily_silverman$p1_sig == TRUE & ci_superfamily_silverman$p2_sig == TRUE & ci_superfamily_silverman$p3_sig == TRUE & ci_superfamily_silverman$p4_sig == TRUE & ci_superfamily_silverman$p5_sig == TRUE  & ci_superfamily_silverman$p6_sig == FALSE, "yes","no")

ci_superfamily_silverman$more_modes <- ifelse(ci_superfamily_silverman$p1_sig == TRUE & ci_superfamily_silverman$p2_sig == TRUE & ci_superfamily_silverman$p3_sig == TRUE & ci_superfamily_silverman$p4_sig == TRUE & ci_superfamily_silverman$p5_sig == TRUE  & ci_superfamily_silverman$p6_sig == TRUE, "yes","no")


write.table(ci_superfamily_silverman, file = "inopinata_superfamilies",row.names = FALSE, quote = FALSE, sep = "\t")


#nigoni





cn_superfamily_k_1 <- aggregate( cn_superfamily_tukey$kimura_distance ~ cn_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=1,R=999, adjust = FALSE, digits = 6)@p_value)
cn_superfamily_k_2 <- aggregate( cn_superfamily_tukey$kimura_distance ~ cn_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=2,R=999, adjust = FALSE, digits = 6)@p_value)
cn_superfamily_k_3 <- aggregate( cn_superfamily_tukey$kimura_distance ~ cn_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=3,R=999, adjust = FALSE, digits = 6)@p_value)
cn_superfamily_k_4 <- aggregate( cn_superfamily_tukey$kimura_distance ~ cn_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=4,R=999, adjust = FALSE, digits = 6)@p_value)
cn_superfamily_k_5 <- aggregate( cn_superfamily_tukey$kimura_distance ~ cn_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=5,R=999, adjust = FALSE, digits = 6)@p_value)
cn_superfamily_k_6 <- aggregate( cn_superfamily_tukey$kimura_distance ~ cn_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=6,R=999, adjust = FALSE, digits = 6)@p_value)

colnames(cn_superfamily_k_1)[colnames(cn_superfamily_k_1)=="cn_superfamily_tukey$rep_type"] <- "rep_type"
colnames(cn_superfamily_k_1)[colnames(cn_superfamily_k_1)=="cn_superfamily_tukey$kimura_distance"] <- "p_1"

colnames(cn_superfamily_k_2)[colnames(cn_superfamily_k_2)=="cn_superfamily_tukey$rep_type"] <- "rep_type"
colnames(cn_superfamily_k_2)[colnames(cn_superfamily_k_2)=="cn_superfamily_tukey$kimura_distance"] <- "p_2"

colnames(cn_superfamily_k_3)[colnames(cn_superfamily_k_3)=="cn_superfamily_tukey$rep_type"] <- "rep_type"
colnames(cn_superfamily_k_3)[colnames(cn_superfamily_k_3)=="cn_superfamily_tukey$kimura_distance"] <- "p_3"


colnames(cn_superfamily_k_4)[colnames(cn_superfamily_k_4)=="cn_superfamily_tukey$rep_type"] <- "rep_type"
colnames(cn_superfamily_k_4)[colnames(cn_superfamily_k_4)=="cn_superfamily_tukey$kimura_distance"] <- "p_4"

colnames(cn_superfamily_k_5)[colnames(cn_superfamily_k_5)=="cn_superfamily_tukey$rep_type"] <- "rep_type"
colnames(cn_superfamily_k_5)[colnames(cn_superfamily_k_5)=="cn_superfamily_tukey$kimura_distance"] <- "p_5"

colnames(cn_superfamily_k_6)[colnames(cn_superfamily_k_6)=="cn_superfamily_tukey$rep_type"] <- "rep_type"
colnames(cn_superfamily_k_6)[colnames(cn_superfamily_k_6)=="cn_superfamily_tukey$kimura_distance"] <- "p_6"

cn_superfamily_silverman <- data.frame(rep_type = cn_superfamily_k_1$rep_type, p_1 = cn_superfamily_k_1$p_1, p_2 = cn_superfamily_k_2$p_2, p_3 = cn_superfamily_k_3$p_3, p_4 = cn_superfamily_k_4$p_4, p_5 = cn_superfamily_k_5$p_5, p_6 = cn_superfamily_k_6$p_6)

cn_superfamily_silverman$p1_sig <- cn_superfamily_silverman$p_1 < 0.05
cn_superfamily_silverman$p2_sig <- cn_superfamily_silverman$p_2 < 0.05
cn_superfamily_silverman$p3_sig <- cn_superfamily_silverman$p_3 < 0.05
cn_superfamily_silverman$p4_sig <- cn_superfamily_silverman$p_4 < 0.05
cn_superfamily_silverman$p5_sig <- cn_superfamily_silverman$p_5 < 0.05
cn_superfamily_silverman$p6_sig <- cn_superfamily_silverman$p_6 < 0.05

cn_superfamily_silverman$one_mode <- ifelse(cn_superfamily_silverman$p1_sig == FALSE,"yes","no")
cn_superfamily_silverman$two_modes <- ifelse(cn_superfamily_silverman$p1_sig == TRUE & cn_superfamily_silverman$p2_sig == FALSE,"yes","no")
cn_superfamily_silverman$three_modes <- ifelse(cn_superfamily_silverman$p1_sig == TRUE & cn_superfamily_silverman$p2_sig == TRUE & cn_superfamily_silverman$p3_sig == FALSE,"yes","no")


cn_superfamily_silverman$four_modes <- ifelse(cn_superfamily_silverman$p1_sig == TRUE & cn_superfamily_silverman$p2_sig == TRUE & cn_superfamily_silverman$p3_sig == TRUE & cn_superfamily_silverman$p4_sig == FALSE, "yes","no")

cn_superfamily_silverman$five_modes <- ifelse(cn_superfamily_silverman$p1_sig == TRUE & cn_superfamily_silverman$p2_sig == TRUE & cn_superfamily_silverman$p3_sig == TRUE & cn_superfamily_silverman$p4_sig == TRUE & cn_superfamily_silverman$p5_sig == FALSE, "yes","no")

cn_superfamily_silverman$six_modes <- ifelse(cn_superfamily_silverman$p1_sig == TRUE & cn_superfamily_silverman$p2_sig == TRUE & cn_superfamily_silverman$p3_sig == TRUE & cn_superfamily_silverman$p4_sig == TRUE & cn_superfamily_silverman$p5_sig == TRUE  & cn_superfamily_silverman$p6_sig == FALSE, "yes","no")

cn_superfamily_silverman$more_modes <- ifelse(cn_superfamily_silverman$p1_sig == TRUE & cn_superfamily_silverman$p2_sig == TRUE & cn_superfamily_silverman$p3_sig == TRUE & cn_superfamily_silverman$p4_sig == TRUE & cn_superfamily_silverman$p5_sig == TRUE  & cn_superfamily_silverman$p6_sig == TRUE, "yes","no")


write.table(cn_superfamily_silverman, file = "nigoni_superfamilies",row.names = FALSE, quote = FALSE, sep = "\t")


#remanei superfamily






cr_superfamily_k_1 <- aggregate( cr_superfamily_tukey$kimura_distance ~ cr_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=1,R=999, adjust = FALSE, digits = 6)@p_value)
cr_superfamily_k_2 <- aggregate( cr_superfamily_tukey$kimura_distance ~ cr_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=2,R=999, adjust = FALSE, digits = 6)@p_value)
cr_superfamily_k_3 <- aggregate( cr_superfamily_tukey$kimura_distance ~ cr_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=3,R=999, adjust = FALSE, digits = 6)@p_value)
cr_superfamily_k_4 <- aggregate( cr_superfamily_tukey$kimura_distance ~ cr_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=4,R=999, adjust = FALSE, digits = 6)@p_value)
cr_superfamily_k_5 <- aggregate( cr_superfamily_tukey$kimura_distance ~ cr_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=5,R=999, adjust = FALSE, digits = 6)@p_value)
cr_superfamily_k_6 <- aggregate( cr_superfamily_tukey$kimura_distance ~ cr_superfamily_tukey$rep_type, FUN = function(x) silverman.test(x,k=6,R=999, adjust = FALSE, digits = 6)@p_value)

colnames(cr_superfamily_k_1)[colnames(cr_superfamily_k_1)=="cr_superfamily_tukey$rep_type"] <- "rep_type"
colnames(cr_superfamily_k_1)[colnames(cr_superfamily_k_1)=="cr_superfamily_tukey$kimura_distance"] <- "p_1"

colnames(cr_superfamily_k_2)[colnames(cr_superfamily_k_2)=="cr_superfamily_tukey$rep_type"] <- "rep_type"
colnames(cr_superfamily_k_2)[colnames(cr_superfamily_k_2)=="cr_superfamily_tukey$kimura_distance"] <- "p_2"

colnames(cr_superfamily_k_3)[colnames(cr_superfamily_k_3)=="cr_superfamily_tukey$rep_type"] <- "rep_type"
colnames(cr_superfamily_k_3)[colnames(cr_superfamily_k_3)=="cr_superfamily_tukey$kimura_distance"] <- "p_3"


colnames(cr_superfamily_k_4)[colnames(cr_superfamily_k_4)=="cr_superfamily_tukey$rep_type"] <- "rep_type"
colnames(cr_superfamily_k_4)[colnames(cr_superfamily_k_4)=="cr_superfamily_tukey$kimura_distance"] <- "p_4"

colnames(cr_superfamily_k_5)[colnames(cr_superfamily_k_5)=="cr_superfamily_tukey$rep_type"] <- "rep_type"
colnames(cr_superfamily_k_5)[colnames(cr_superfamily_k_5)=="cr_superfamily_tukey$kimura_distance"] <- "p_5"

colnames(cr_superfamily_k_6)[colnames(cr_superfamily_k_6)=="cr_superfamily_tukey$rep_type"] <- "rep_type"
colnames(cr_superfamily_k_6)[colnames(cr_superfamily_k_6)=="cr_superfamily_tukey$kimura_distance"] <- "p_6"

cr_superfamily_silverman <- data.frame(rep_type = cr_superfamily_k_1$rep_type, p_1 = cr_superfamily_k_1$p_1, p_2 = cr_superfamily_k_2$p_2, p_3 = cr_superfamily_k_3$p_3, p_4 = cr_superfamily_k_4$p_4, p_5 = cr_superfamily_k_5$p_5, p_6 = cr_superfamily_k_6$p_6)

cr_superfamily_silverman$p1_sig <- cr_superfamily_silverman$p_1 < 0.05
cr_superfamily_silverman$p2_sig <- cr_superfamily_silverman$p_2 < 0.05
cr_superfamily_silverman$p3_sig <- cr_superfamily_silverman$p_3 < 0.05
cr_superfamily_silverman$p4_sig <- cr_superfamily_silverman$p_4 < 0.05
cr_superfamily_silverman$p5_sig <- cr_superfamily_silverman$p_5 < 0.05
cr_superfamily_silverman$p6_sig <- cr_superfamily_silverman$p_6 < 0.05

cr_superfamily_silverman$one_mode <- ifelse(cr_superfamily_silverman$p1_sig == FALSE,"yes","no")
cr_superfamily_silverman$two_modes <- ifelse(cr_superfamily_silverman$p1_sig == TRUE & cr_superfamily_silverman$p2_sig == FALSE,"yes","no")
cr_superfamily_silverman$three_modes <- ifelse(cr_superfamily_silverman$p1_sig == TRUE & cr_superfamily_silverman$p2_sig == TRUE & cr_superfamily_silverman$p3_sig == FALSE,"yes","no")


cr_superfamily_silverman$four_modes <- ifelse(cr_superfamily_silverman$p1_sig == TRUE & cr_superfamily_silverman$p2_sig == TRUE & cr_superfamily_silverman$p3_sig == TRUE & cr_superfamily_silverman$p4_sig == FALSE, "yes","no")

cr_superfamily_silverman$five_modes <- ifelse(cr_superfamily_silverman$p1_sig == TRUE & cr_superfamily_silverman$p2_sig == TRUE & cr_superfamily_silverman$p3_sig == TRUE & cr_superfamily_silverman$p4_sig == TRUE & cr_superfamily_silverman$p5_sig == FALSE, "yes","no")

cr_superfamily_silverman$six_modes <- ifelse(cr_superfamily_silverman$p1_sig == TRUE & cr_superfamily_silverman$p2_sig == TRUE & cr_superfamily_silverman$p3_sig == TRUE & cr_superfamily_silverman$p4_sig == TRUE & cr_superfamily_silverman$p5_sig == TRUE  & cr_superfamily_silverman$p6_sig == FALSE, "yes","no")

cr_superfamily_silverman$more_modes <- ifelse(cr_superfamily_silverman$p1_sig == TRUE & cr_superfamily_silverman$p2_sig == TRUE & cr_superfamily_silverman$p3_sig == TRUE & cr_superfamily_silverman$p4_sig == TRUE & cr_superfamily_silverman$p5_sig == TRUE  & cr_superfamily_silverman$p6_sig == TRUE, "yes","no")


write.table(cr_superfamily_silverman, file = "remanei_superfamilies",row.names = FALSE, quote = FALSE, sep = "\t")


