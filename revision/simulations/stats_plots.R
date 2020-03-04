

library(effsize)

#get data in there
rep_dat <- read.table("all_dec_27_simulations_TE_sites.tsv", sep="\t", header=TRUE)


#define chromosome arms and centers by middle or outer half of chromosome
rep_dat$chr_str_type <- ifelse(rep_dat$norm_dist_center >= 0.25,"arms", "centers")
#get the factors right
rep_dat$chr_str_type <- factor(rep_dat$chr_str_type, levels = c("centers","arms"))


#get arm-center effect sizes

#make empty things to be filled up later
scen_type <- NULL
mean_arms <- NULL
mean_centers <-NULL
arm_cen_effect_size_lower <- NULL
arm_cen_effect_size <- NULL
arm_cen_effect_size_upper <- NULL
wilcox_w <- NULL
wilcox_p <- NULL
mean_tot_te <- NULL
	#get stats for every evolutionary scenario
for (i in levels(rep_dat$scenario)){
	dat <- rep_dat[rep_dat$scenario == i,]
	dat_arms <- dat[dat$chr_str_type == "arms",]
	dat_centers <- dat[dat$chr_str_type == "centers",]
	mean_arms <- rbind(mean_arms,mean(aggregate(dat_arms$tot_TE ~ dat_arms$sim_id, FUN = function(x) sum(x))[,2]))
	mean_centers <- rbind(mean_centers,mean(aggregate(dat_centers$tot_TE ~ dat_centers$sim_id, FUN = function(x) sum(x))[,2]))
	mean_tot_te <- rbind(mean_tot_te,mean(aggregate(dat$tot_TE ~ dat$sim_id, FUN = function(x) sum(x))[,2])) #average total number of TE's among simulations
	arm_cen_effect_size_lower <- rbind(arm_cen_effect_size_lower, cohen.d(dat_arms$tot_TE,dat_centers$tot_TE)$conf.int[1]) #arm-center effect size stats
	arm_cen_effect_size <- rbind(arm_cen_effect_size, cohen.d(dat_arms$tot_TE,dat_centers$tot_TE)$estimate)
	arm_cen_effect_size_upper <- rbind(arm_cen_effect_size_upper, cohen.d(dat_arms$tot_TE,dat_centers$tot_TE)$conf.int[2])
	scen_type <- rbind(scen_type,i) #scenario id
	wilcox_w <- rbind(wilcox_w, wilcox.test(dat_arms$tot_TE,dat_centers$tot_TE,exact=FALSE)$statistic[[1]]) #wilcox rank sum test statistics
	wilcox_p <- rbind(wilcox_p, wilcox.test(dat_arms$tot_TE,dat_centers$tot_TE,exact=FALSE)$p.value[[1]])
	stat_df <- cbind(scen_type,mean_tot_te,mean_arms,mean_centers,arm_cen_effect_size_lower,arm_cen_effect_size,arm_cen_effect_size_upper,wilcox_w,wilcox_p) #put all together
}

#get data right
rownames(stat_df) <- NULL
colnames(stat_df) <- c("SCENARIO","mean_tot_te","mean_arms","mean_centers","arm_cen_effect_size_lower","arm_cen_effect_size","arm_cen_effect_size_upper","wilcox_w","wilcox_p")

stat_df <- as.data.frame(stat_df)
stat_df$mean_tot_te <-  as.numeric(levels(stat_df$mean_tot_te))[stat_df$mean_tot_te]
stat_df$mean_arms <-  as.numeric(levels(stat_df$mean_arms))[stat_df$mean_arms]
stat_df$mean_centers <-  as.numeric(levels(stat_df$mean_centers))[stat_df$mean_centers]
stat_df$arm_cen_effect_size_lower <-  as.numeric(levels(stat_df$arm_cen_effect_size_lower))[stat_df$arm_cen_effect_size_lower]
stat_df$arm_cen_effect_size <-  as.numeric(levels(stat_df$arm_cen_effect_size))[stat_df$arm_cen_effect_size]
stat_df$arm_cen_effect_size_upper <-  as.numeric(levels(stat_df$arm_cen_effect_size_upper))[stat_df$arm_cen_effect_size_upper]
stat_df$wilcox_w <-  as.numeric(levels(stat_df$wilcox_w))[stat_df$wilcox_w]
stat_df$wilcox_p <-  as.numeric(levels(stat_df$wilcox_p))[stat_df$wilcox_p]

#write the data to a table
write.table(stat_df, file = "/home/gavin/genome/genome/repeats_12-18-18/revisions/simulations_12-20-19/work_on_macbook_12-22-19/06_SIMULATIONS_DEC_27_2019/local_00_awk_TE_sites/simulations_effect_sizes_wilcox_tests_12-27-19",row.names = FALSE, quote = FALSE, sep = "\t")

#see the stats
stat_df




#do it again with older simulations
#get data in there
rep_dat <- read.table("early_dec_simulations_TE_sites.tsv", sep="\t", header=TRUE)


#define chromosome arms and centers by middle or outer half of chromosome
rep_dat$chr_str_type <- ifelse(rep_dat$norm_dist_center >= 0.25,"arms", "centers")
#get the factors right
rep_dat$chr_str_type <- factor(rep_dat$chr_str_type, levels = c("centers","arms"))


#get arm-center effect sizes

#make empty things to be filled up later
scen_type <- NULL
mean_arms <- NULL
mean_centers <-NULL
arm_cen_effect_size_lower <- NULL
arm_cen_effect_size <- NULL
arm_cen_effect_size_upper <- NULL
wilcox_w <- NULL
wilcox_p <- NULL
mean_tot_te <- NULL
	#get stats for every evolutionary scenario
for (i in levels(rep_dat$scenario)){
	dat <- rep_dat[rep_dat$scenario == i,]
	dat_arms <- dat[dat$chr_str_type == "arms",]
	dat_centers <- dat[dat$chr_str_type == "centers",]
	mean_arms <- rbind(mean_arms,mean(aggregate(dat_arms$tot_TE ~ dat_arms$sim_id, FUN = function(x) sum(x))[,2]))
	mean_centers <- rbind(mean_centers,mean(aggregate(dat_centers$tot_TE ~ dat_centers$sim_id, FUN = function(x) sum(x))[,2]))
	mean_tot_te <- rbind(mean_tot_te,mean(aggregate(dat$tot_TE ~ dat$sim_id, FUN = function(x) sum(x))[,2])) #average total number of TE's among simulations
	arm_cen_effect_size_lower <- rbind(arm_cen_effect_size_lower, cohen.d(dat_arms$tot_TE,dat_centers$tot_TE)$conf.int[1]) #arm-center effect size stats
	arm_cen_effect_size <- rbind(arm_cen_effect_size, cohen.d(dat_arms$tot_TE,dat_centers$tot_TE)$estimate)
	arm_cen_effect_size_upper <- rbind(arm_cen_effect_size_upper, cohen.d(dat_arms$tot_TE,dat_centers$tot_TE)$conf.int[2])
	scen_type <- rbind(scen_type,i) #scenario id
	wilcox_w <- rbind(wilcox_w, wilcox.test(dat_arms$tot_TE,dat_centers$tot_TE,exact=FALSE)$statistic[[1]]) #wilcox rank sum test statistics
	wilcox_p <- rbind(wilcox_p, wilcox.test(dat_arms$tot_TE,dat_centers$tot_TE,exact=FALSE)$p.value[[1]])
	stat_df <- cbind(scen_type,mean_tot_te,mean_arms,mean_centers,arm_cen_effect_size_lower,arm_cen_effect_size,arm_cen_effect_size_upper,wilcox_w,wilcox_p) #put all together
}

#get data right
rownames(stat_df) <- NULL
colnames(stat_df) <- c("SCENARIO","mean_tot_te","mean_arms","mean_centers","arm_cen_effect_size_lower","arm_cen_effect_size","arm_cen_effect_size_upper","wilcox_w","wilcox_p")

stat_df <- as.data.frame(stat_df)
stat_df$mean_tot_te <-  as.numeric(levels(stat_df$mean_tot_te))[stat_df$mean_tot_te]
stat_df$mean_arms <-  as.numeric(levels(stat_df$mean_arms))[stat_df$mean_arms]
stat_df$mean_centers <-  as.numeric(levels(stat_df$mean_centers))[stat_df$mean_centers]
stat_df$arm_cen_effect_size_lower <-  as.numeric(levels(stat_df$arm_cen_effect_size_lower))[stat_df$arm_cen_effect_size_lower]
stat_df$arm_cen_effect_size <-  as.numeric(levels(stat_df$arm_cen_effect_size))[stat_df$arm_cen_effect_size]
stat_df$arm_cen_effect_size_upper <-  as.numeric(levels(stat_df$arm_cen_effect_size_upper))[stat_df$arm_cen_effect_size_upper]
stat_df$wilcox_w <-  as.numeric(levels(stat_df$wilcox_w))[stat_df$wilcox_w]
stat_df$wilcox_p <-  as.numeric(levels(stat_df$wilcox_p))[stat_df$wilcox_p]

#write the data to a table
write.table(stat_df, file = "early_dec_simulations_effect_sizes_wilcox_tests_12-27-19",row.names = FALSE, quote = FALSE, sep = "\t")

#see the stats
stat_df

#supplemental figure 29
rep_dat <- read.table("ac_effect_size_data.tsv", sep="\t", header=TRUE)


notable_scen <- rep_dat[rep_dat$scenario == "DOM_NEUTRAL" | rep_dat$scenario == "DOM_sel_all_WEAK_SEL_-0.0002" | rep_dat$scenario == "DOM_sel_all_WEAK_SEL_-0.0005" | rep_dat$scenario == "DOM_sel_all_WEAK_SEL_-0.001" | rep_dat$scenario == "DOM_sel_all_WEAK_SEL_-0.0015" | rep_dat$scenario == "DOM_sel_all_WEAK_SEL_-0.002",]

notable_scen$scenario <- droplevels(notable_scen$scenario)


#define chromosome arms and centers by middle or outer half of chromosome
notable_scen$chr_str_type <- ifelse(notable_scen$norm_dist_center >= 0.25,"arms", "centers")
#get the factors right
notable_scen$chr_str_type <- factor(notable_scen$chr_str_type, levels = c("centers","arms"))

#Ns = 0
Ns_0 <- notable_scen[notable_scen$scenario == "DOM_NEUTRAL",]

scen_type <- NULL
sim_ID <- NULL
arm_cen_effect_size <- NULL
Ns_0_df <- NULL

Ns_0$sim_id <- droplevels(Ns_0$sim_id)
Ns_0$scenario <- droplevels(Ns_0$scenario)


for (i in levels(Ns_0$sim_id)){
	dat <- Ns_0[Ns_0$sim_id == i,]
	dat_arms <- dat[dat$chr_str_type == "arms",]
	dat_centers <- dat[dat$chr_str_type == "centers",]
	arm_cen_effect_size <- rbind(arm_cen_effect_size, cohen.d(dat_arms$tot_TE,dat_centers$tot_TE)$estimate)
	scen_type <- rbind(scen_type,unique(levels(Ns_0$scenario)))
	sim_ID <- rbind(sim_ID,i)
	Ns_0_df <- cbind(scen_type,sim_ID,arm_cen_effect_size)
}

rownames(Ns_0_df) <- NULL
colnames(Ns_0_df) <- c("SCENARIO","sim_id","arm_cen_effect_size")
Ns_0_df <- as.data.frame(Ns_0_df)
Ns_0_df$arm_cen_effect_size <-  as.numeric(levels(Ns_0_df$arm_cen_effect_size))[Ns_0_df$arm_cen_effect_size]
Ns_0_df


#Ns = 1

Ns_1 <- notable_scen[notable_scen$scenario == "DOM_sel_all_WEAK_SEL_-0.0002",]

scen_type <- NULL
sim_ID <- NULL
arm_cen_effect_size <- NULL
Ns_1_df <- NULL

Ns_1$sim_id <- droplevels(Ns_1$sim_id)
Ns_1$scenario <- droplevels(Ns_1$scenario)


for (i in levels(Ns_1$sim_id)){
	dat <- Ns_1[Ns_1$sim_id == i,]
	dat_arms <- dat[dat$chr_str_type == "arms",]
	dat_centers <- dat[dat$chr_str_type == "centers",]
	arm_cen_effect_size <- rbind(arm_cen_effect_size, cohen.d(dat_arms$tot_TE,dat_centers$tot_TE)$estimate)
	scen_type <- rbind(scen_type,unique(levels(Ns_1$scenario)))
	sim_ID <- rbind(sim_ID,i)
	Ns_1_df <- cbind(scen_type,sim_ID,arm_cen_effect_size)
}

rownames(Ns_1_df) <- NULL
colnames(Ns_1_df) <- c("SCENARIO","sim_id","arm_cen_effect_size")
Ns_1_df <- as.data.frame(Ns_1_df)
Ns_1_df$arm_cen_effect_size <-  as.numeric(levels(Ns_1_df$arm_cen_effect_size))[Ns_1_df$arm_cen_effect_size]
Ns_1_df




#Ns = 2.5

Ns_2.5 <- notable_scen[notable_scen$scenario == "DOM_sel_all_WEAK_SEL_-0.0005",]

scen_type <- NULL
sim_ID <- NULL
arm_cen_effect_size <- NULL
Ns_2.5_df <- NULL

Ns_2.5$sim_id <- droplevels(Ns_2.5$sim_id)
Ns_2.5$scenario <- droplevels(Ns_2.5$scenario)


for (i in levels(Ns_2.5$sim_id)){
	dat <- Ns_2.5[Ns_2.5$sim_id == i,]
	dat_arms <- dat[dat$chr_str_type == "arms",]
	dat_centers <- dat[dat$chr_str_type == "centers",]
	arm_cen_effect_size <- rbind(arm_cen_effect_size, cohen.d(dat_arms$tot_TE,dat_centers$tot_TE)$estimate)
	scen_type <- rbind(scen_type,unique(levels(Ns_2.5$scenario)))
	sim_ID <- rbind(sim_ID,i)
	Ns_2.5_df <- cbind(scen_type,sim_ID,arm_cen_effect_size)
}

rownames(Ns_2.5_df) <- NULL
colnames(Ns_2.5_df) <- c("SCENARIO","sim_id","arm_cen_effect_size")
Ns_2.5_df <- as.data.frame(Ns_2.5_df)
Ns_2.5_df$arm_cen_effect_size <-  as.numeric(levels(Ns_2.5_df$arm_cen_effect_size))[Ns_2.5_df$arm_cen_effect_size]
Ns_2.5_df



#Ns = 5

Ns_5 <- notable_scen[notable_scen$scenario == "DOM_sel_all_WEAK_SEL_-0.001",]

scen_type <- NULL
sim_ID <- NULL
arm_cen_effect_size <- NULL
Ns_5_df <- NULL

Ns_5$sim_id <- droplevels(Ns_5$sim_id)
Ns_5$scenario <- droplevels(Ns_5$scenario)


for (i in levels(Ns_5$sim_id)){
	dat <- Ns_5[Ns_5$sim_id == i,]
	dat_arms <- dat[dat$chr_str_type == "arms",]
	dat_centers <- dat[dat$chr_str_type == "centers",]
	arm_cen_effect_size <- rbind(arm_cen_effect_size, cohen.d(dat_arms$tot_TE,dat_centers$tot_TE)$estimate)
	scen_type <- rbind(scen_type,unique(levels(Ns_5$scenario)))
	sim_ID <- rbind(sim_ID,i)
	Ns_5_df <- cbind(scen_type,sim_ID,arm_cen_effect_size)
}

rownames(Ns_5_df) <- NULL
colnames(Ns_5_df) <- c("SCENARIO","sim_id","arm_cen_effect_size")
Ns_5_df <- as.data.frame(Ns_5_df)
Ns_5_df$arm_cen_effect_size <-  as.numeric(levels(Ns_5_df$arm_cen_effect_size))[Ns_5_df$arm_cen_effect_size]
Ns_5_df



#Ns = 7.5

Ns_7.5 <- notable_scen[notable_scen$scenario == "DOM_sel_all_WEAK_SEL_-0.0015",]

scen_type <- NULL
sim_ID <- NULL
arm_cen_effect_size <- NULL
Ns_7.5_df <- NULL

Ns_7.5$sim_id <- droplevels(Ns_7.5$sim_id)
Ns_7.5$scenario <- droplevels(Ns_7.5$scenario)


for (i in levels(Ns_7.5$sim_id)){
	dat <- Ns_7.5[Ns_7.5$sim_id == i,]
	dat_arms <- dat[dat$chr_str_type == "arms",]
	dat_centers <- dat[dat$chr_str_type == "centers",]
	arm_cen_effect_size <- rbind(arm_cen_effect_size, cohen.d(dat_arms$tot_TE,dat_centers$tot_TE)$estimate)
	scen_type <- rbind(scen_type,unique(levels(Ns_7.5$scenario)))
	sim_ID <- rbind(sim_ID,i)
	Ns_7.5_df <- cbind(scen_type,sim_ID,arm_cen_effect_size)
}

rownames(Ns_7.5_df) <- NULL
colnames(Ns_7.5_df) <- c("SCENARIO","sim_id","arm_cen_effect_size")
Ns_7.5_df <- as.data.frame(Ns_7.5_df)
Ns_7.5_df$arm_cen_effect_size <-  as.numeric(levels(Ns_7.5_df$arm_cen_effect_size))[Ns_7.5_df$arm_cen_effect_size]
Ns_7.5_df




#Ns = 10

Ns_10 <- notable_scen[notable_scen$scenario == "DOM_sel_all_WEAK_SEL_-0.002",]

scen_type <- NULL
sim_ID <- NULL
arm_cen_effect_size <- NULL
Ns_10_df <- NULL

Ns_10$sim_id <- droplevels(Ns_10$sim_id)
Ns_10$scenario <- droplevels(Ns_10$scenario)


for (i in levels(Ns_10$sim_id)){
	dat <- Ns_10[Ns_10$sim_id == i,]
	dat_arms <- dat[dat$chr_str_type == "arms",]
	dat_centers <- dat[dat$chr_str_type == "centers",]
	arm_cen_effect_size <- rbind(arm_cen_effect_size, cohen.d(dat_arms$tot_TE,dat_centers$tot_TE)$estimate)
	scen_type <- rbind(scen_type,unique(levels(Ns_10$scenario)))
	sim_ID <- rbind(sim_ID,i)
	Ns_10_df <- cbind(scen_type,sim_ID,arm_cen_effect_size)
}

rownames(Ns_10_df) <- NULL
colnames(Ns_10_df) <- c("SCENARIO","sim_id","arm_cen_effect_size")
Ns_10_df <- as.data.frame(Ns_10_df)
Ns_10_df$arm_cen_effect_size <-  as.numeric(levels(Ns_10_df$arm_cen_effect_size))[Ns_10_df$arm_cen_effect_size]
Ns_10_df


all_ac_ef_df <- rbind(Ns_0_df,Ns_1_df,Ns_2.5_df,Ns_5_df,Ns_7.5_df,Ns_10_df)

all_ac_ef_df$Ns <- ifelse(all_ac_ef_df$SCENARIO == "DOM_NEUTRAL",0, 1)

all_ac_ef_df[all_ac_ef_df$SCENARIO == "DOM_sel_all_WEAK_SEL_-0.0002",]$Ns <- 1
all_ac_ef_df[all_ac_ef_df$SCENARIO == "DOM_sel_all_WEAK_SEL_-0.0005",]$Ns <- 2.5
all_ac_ef_df[all_ac_ef_df$SCENARIO == "DOM_sel_all_WEAK_SEL_-0.001",]$Ns <- 5
all_ac_ef_df[all_ac_ef_df$SCENARIO == "DOM_sel_all_WEAK_SEL_-0.0015",]$Ns <- 7.5
all_ac_ef_df[all_ac_ef_df$SCENARIO == "DOM_sel_all_WEAK_SEL_-0.002",]$Ns <- 10


summary(all_ac_ef_df$arm_cen_effect_size)

#this is supplemental figure 29

ggplot(all_ac_ef_df, aes(x=Ns, y=arm_cen_effect_size)) + geom_hline(yintercept=0, linetype="dashed") + geom_hline(yintercept=0.04538246,colour="red",size=0.5) + geom_hline(yintercept=1.139733,colour="blue",size=0.5) + geom_smooth(size=0.5, colour="darkgray", alpha=0.125) + geom_point(size=1,alpha=0.25) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),axis.title=element_text(size=14), axis.ticks = element_line(size = 0.5)) + xlab("Ns")+ ylab("Arm-center difference (effect size)") + annotate("text", x = 9, y = 0.1, label = "C. inopinata", colour="red") + annotate("text", x = 9, y = 1.2, label = "C. briggsae", colour="blue") + scale_y_continuous(limits = c(-1.5,1.5),breaks = c(-1.5,-1.2,-0.9,-0.6,-0.3,0,0.3,0.6,0.9,1.2,1.5)) + scale_x_continuous(limits=c(-0,10),breaks=c(0:10)) 




#figures


library(ggplot2)

#get data in there (see lines 706-1095 of repeats.sh to see how this is generated)

options(scipen=999)

rep_data <- read.table("DOM_NEUTRAL", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = tot_TE)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("Number of TE's") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_NEUTRAL")






rep_data <- read.table("NO_DOM_NEUTRAL", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = tot_TE)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("Number of TE's") + theme(plot.title = element_text(size=12)) + ggtitle("NO_DOM_NEUTRAL")









rep_data <- read.table("DOM_SEL_ALL", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = tot_TE)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("Number of TE's") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_SEL_ALL")









rep_data <- read.table("NO_DOM_SEL_ALL", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = tot_TE)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("Number of TE's") + theme(plot.title = element_text(size=12)) + ggtitle("NO_DOM_SEL_ALL")









rep_data <- read.table("DOM_SEL_ALL_WEAK", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = tot_TE)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("Number of TE's") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_SEL_ALL_WEAK")








rep_data <- read.table("NO_DOM_SEL_ALL_WEAK", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = tot_TE)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("Number of TE's") + theme(plot.title = element_text(size=12)) + ggtitle("NO_DOM_SEL_ALL_WEAK")







rep_data <- read.table("DOM_SEL_ALL_WEAK_POS_MUT", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = tot_TE)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("Number of TE's") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_SEL_ALL_WEAK_POS_MUT")









rep_data <- read.table("NO_DOM_SEL_ALL_WEAK_POS_MUT", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = tot_TE)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("Number of TE's") + theme(plot.title = element_text(size=12)) + ggtitle("NO_DOM_SEL_ALL_WEAK_POS_MUT")







rep_data <- read.table("DOM_SEL_ARM_WEAK_CENTER_STRONG", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = tot_TE)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("Number of TE's") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_SEL_ARM_WEAK_CENTER_STRONG")










rep_data <- read.table("NO_DOM_SEL_ARM_WEAK_CENTR_STRONG", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = tot_TE)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("Number of TE's") + theme(plot.title = element_text(size=12)) + ggtitle("NO_DOM_SEL_ARM_WEAK_CENTR_STRONG")






rep_data <- read.table("DOM_SEL_CENTR_WEAK", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = tot_TE)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("Number of TE's") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_SEL_CENTR_WEAK")







rep_data <- read.table("NO_DOM_SEL_CENTR_WEAK", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = tot_TE)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("Number of TE's") + theme(plot.title = element_text(size=12)) + ggtitle("NO_DOM_SEL_CENTR_WEAK")





library(ggplot2)

#get data in there (see lines 706-1095 of repeats.sh to see how this is generated)

options(scipen=999)

rep_data <- read.table("DOM_sel_all_WEAK_SEL_-0.0002", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = tot_TE)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_sel_all_WEAK_SEL_-0.0002")









rep_data <- read.table("DOM_sel_all_WEAK_SEL_-0.0005", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = tot_TE)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_sel_all_WEAK_SEL_-0.0005")







rep_data <- read.table("DOM_sel_all_WEAK_SEL_-0.001", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = tot_TE)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_sel_all_WEAK_SEL_-0.001")







rep_data <- read.table("DOM_sel_all_WEAK_SEL_-0.0015", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = tot_TE)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_sel_all_WEAK_SEL_-0.0015")











rep_data <- read.table("DOM_sel_all_WEAK_SEL_-0.002", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = tot_TE)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_sel_all_WEAK_SEL_-0.002")











rep_data <- read.table("no_dom_sel_all_weak_SELFING", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = tot_TE)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("no_dom_sel_all_weak_SELFING")









rep_data <- read.table("DOM_sel_all_WEAK_SELFING", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = tot_TE)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_sel_all_WEAK_SELFING")










rep_data <- read.table("no_dom_sel_arm_weak_Center_Strong_SELFING", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = tot_TE)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("no_dom_sel_arm_weak_Center_Strong_SELFING")











rep_data <- read.table("DOM_sel_arm_WEAK_center_Strong_SELFING", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = tot_TE)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_sel_arm_WEAK_center_Strong_SELFING")








rep_data <- read.table("no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = tot_TE)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN")










rep_data <- read.table("DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = tot_TE)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN")






#ages plots


library(ggplot2)

#get data in there (see lines 706-1095 of repeats.sh to see how this is generated)

options(scipen=999)

rep_data <- read.table("DOM_sel_all_WEAK_SEL_-0.0002", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_sel_all_WEAK_SEL_-0.0002")









rep_data <- read.table("DOM_sel_all_WEAK_SEL_-0.0005", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_sel_all_WEAK_SEL_-0.0005")







rep_data <- read.table("DOM_sel_all_WEAK_SEL_-0.001", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_sel_all_WEAK_SEL_-0.001")







rep_data <- read.table("DOM_sel_all_WEAK_SEL_-0.0015", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_sel_all_WEAK_SEL_-0.0015")











rep_data <- read.table("DOM_sel_all_WEAK_SEL_-0.002", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_sel_all_WEAK_SEL_-0.002")











rep_data <- read.table("no_dom_sel_all_weak_SELFING", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("no_dom_sel_all_weak_SELFING")









rep_data <- read.table("DOM_sel_all_WEAK_SELFING", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_sel_all_WEAK_SELFING")










rep_data <- read.table("no_dom_sel_arm_weak_Center_Strong_SELFING", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("no_dom_sel_arm_weak_Center_Strong_SELFING")











rep_data <- read.table("DOM_sel_arm_WEAK_center_Strong_SELFING", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_sel_arm_WEAK_center_Strong_SELFING")








rep_data <- read.table("no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("no_dom_sel_arm_weak_Center_Strong-LOSS_0.5BROKEN")










rep_data <- read.table("DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_sel_arm_WEAK_center_Strong_LOSS_0.5BROKEN")







library(ggplot2)

#get data in there (see lines 706-1095 of repeats.sh to see how this is generated)

options(scipen=999)

rep_data <- read.table("DOM_NEUTRAL", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_NEUTRAL")






rep_data <- read.table("NO_DOM_NEUTRAL", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("NO_DOM_NEUTRAL")









rep_data <- read.table("DOM_SEL_ALL", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_SEL_ALL")









rep_data <- read.table("NO_DOM_SEL_ALL", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("NO_DOM_SEL_ALL")









rep_data <- read.table("DOM_SEL_ALL_WEAK", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_SEL_ALL_WEAK")








rep_data <- read.table("NO_DOM_SEL_ALL_WEAK", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("NO_DOM_SEL_ALL_WEAK")







rep_data <- read.table("DOM_SEL_ALL_WEAK_POS_MUT", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_SEL_ALL_WEAK_POS_MUT")









rep_data <- read.table("NO_DOM_SEL_ALL_WEAK_POS_MUT", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("NO_DOM_SEL_ALL_WEAK_POS_MUT")







rep_data <- read.table("DOM_SEL_ARM_WEAK_CENTER_STRONG", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_SEL_ARM_WEAK_CENTER_STRONG")










rep_data <- read.table("NO_DOM_SEL_ARM_WEAK_CENTR_STRONG", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("NO_DOM_SEL_ARM_WEAK_CENTR_STRONG")










rep_data <- read.table("DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS")







rep_data <- read.table("NO_DOM_SEL_ARM_WEAK_CENTR_STRONG-LOSS", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("NO_DOM_SEL_ARM_WEAK_CENTER_STRONG_LOSS")







rep_data <- read.table("DOM_SEL_CENTR_WEAK", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("DOM_SEL_CENTR_WEAK")







rep_data <- read.table("NO_DOM_SEL_CENTR_WEAK", sep="\t", header=TRUE)

#bp to MB
rep_data$MB <- rep_data$BP/1000000

ggplot(rep_data, aes(x = MB, y = age)) + geom_point(alpha=0.12, size=0.25,colour="black") + stat_smooth(size=0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face = "italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("TE age (generations old)") + theme(plot.title = element_text(size=12)) + ggtitle("NO_DOM_SEL_CENTR_WEAK")




















#supplemental figure 37

library(ggplot2)

rep_dat <- read.table("early_dec_neutral_simulations_over_time.tsv", sep="\t", header=TRUE)


library(lemon)


levels(rep_dat$scenario)[levels(rep_dat$scenario)=="NO_DOM_NEUTRAL"] <- "Uniform recombination (β=7.7)"
levels(rep_dat$scenario)[levels(rep_dat$scenario)=="DOM_NEUTRAL"] <- "Three recombination domains (β=13)"

rep_dat$scenario <- factor(rep_dat$scenario, levels = c("Uniform recombination (β=7.7)","Three recombination domains (β=13)"))


ggplot(rep_dat, aes(x=generation, y=tot_TE, group=sim_id)) + geom_line(alpha=0.5) + facet_rep_wrap(~scenario,ncol=1) + stat_smooth(aes(group=scenario),method=lm, linetype="dashed") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),axis.title=element_text(size=14), axis.ticks = element_line(size = 0.5), strip.background = element_rect(colour="white", fill="white"), strip.text = element_text(size=14, colour="black")) + xlab("Generation")+ ylab("Total number of TE's") + ylim(150000,200000)





geom_smooth(data=rep_dat, method = lm, formula = lm(tot_TE ~ generation, data=rep_dat)) 

summary(lm(tot_TE ~ generation, data=rep_dat))

#> summary(lm(tot_TE ~ generation, data=rep_dat))
#
#Call:
#lm(formula = tot_TE ~ generation, data = rep_dat)
#
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-12342.5  -4728.5    -77.3   3527.0  16173.1 
#
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept) -3.467e+05  1.043e+05  -3.325 0.000889 ***
#generation   1.040e+01  2.083e+00   4.992 6.09e-07 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 5969 on 9824 degrees of freedom
#Multiple R-squared:  0.00253,	Adjusted R-squared:  0.002428 
#F-statistic: 24.92 on 1 and 9824 DF,  p-value: 6.089e-07
