######################################
###            METADATA            ###
######################################
#
# Author: Alexander Skeels
#
# Date: 11.06.2021
#
# Projects: BIGEST: Evolutionary Speed Hypothesis (ESH)
#
# Description: script to see how the distribution of summary statsitics differs between simulated and empirical data
#
# Contact: alexander.skeels@gmail.com
#
######################################


############################################################################
############################################################################
###                                                                      ###
###                              SECTION 1:                              ###
###                          LIBRARIES AND DATA                          ###
###                                                                      ###
############################################################################
############################################################################

# directories
setwd("PATH/TO/DATA")

# libraries
library(ggplot2)
library(reshape)

# data
m0_ss <- read.table("simulation_m0_summary_stats.txt", header=T)
m1_ss <- read.table("simulation_m1_summary_stats.txt", header=T)
m2_ss <- read.table("simulation_m2_summary_stats.txt", header=T)
m3_ss <- read.table("simulation_m3_summary_stats.txt", header=T)

# give the model name
m0_ss$m <- "m0"
m1_ss$m <- "m1"
m2_ss$m <- "m2"
m3_ss$m <- "m3"

# combine four models
m_ss <-  rbind(m0_ss,m1_ss,m2_ss,m3_ss)
# subset the variables extracted from both empiricfal and simulated datasets
m_ss <- m_ss[,c(142, 30:31, 40:41, 48:49, 52:74, 93:104, 114:125)]

# now load empirical data
e_ss <- read.csv("order_empirical_summary_statistics.csv")

# just want to look at diverse clades
e_ss <- na.omit(e_ss[which(e_ss$n_species >= 20),])
e_ss$taxon <- tolower(e_ss$taxon)
colnames(e_ss)[which(colnames(e_ss)=="rs_kutosis")] <- "rs_kurtosis"
colnames(e_ss)[which(colnames(e_ss)=="n_species")] <- "n_extant_diversity"
colnames(e_ss) <- gsub("_p_cor", "_cor",colnames(e_ss)) # change _p_cor for posterior samplescould also change _m_cor to use MCC samples
colnames(e_ss) <- gsub("DivRate", "DR",colnames(e_ss))
colnames(e_ss)[which(colnames(e_ss) == "taxon")] <- "m"
colnames(e_ss)[which(colnames(e_ss) == "collessI_post")] <- "collessI"
colnames(e_ss)[which(colnames(e_ss) == "sackinI_mcc")] <- "sackinI"
colnames(e_ss)[which(colnames(e_ss) == "gamma_mcc")] <- "gamma"


# match empirical and simulated data frames
e_ss <- e_ss[, which(colnames(e_ss) %in% colnames(m_ss))]
m_ss <- m_ss[, which(colnames(m_ss) %in% colnames(e_ss))]
e_ss <- e_ss[,match(colnames(m_ss), colnames(e_ss))]
# check
colnames(e_ss) == colnames(m_ss) 



############################################################################
############################################################################
###                                                                      ###
###     SECTION 2: TRENDS IN SIMULATED SUMMARY STATITICS (DR ~ TEMP)     ###
###                                                                      ###
############################################################################
############################################################################

# Fisher transofrom Spearman correlation coefficient
# Find confidence intervals
# Repeat for each submodel

#M1
m1_ss_redux$temp_DR_cor_fisher   <- FisherZ(m1_ss_redux$temp_DR_cor)
temp_DR_fisher_se    <- CorCI(m1_ss_redux$temp_DR_cor,   m1_ss_redux$n_extant_diversity, conf.level = 0.95, alternative = c("two.sided"))
m1_ss_redux$temp_DR_fisher_se        <-      temp_DR_fisher_se[grepl("cor", names(temp_DR_fisher_se      ))]          

#M2
m2_ss_redux$temp_DR_cor_fisher   <- FisherZ(m2_ss_redux$temp_DR_cor)
temp_DR_fisher_se    <- CorCI(m2_ss_redux$temp_DR_cor,   m2_ss_redux$n_extant_diversity, conf.level = 0.95, alternative = c("two.sided"))
m2_ss_redux$temp_DR_fisher_se        <-      temp_DR_fisher_se[grepl("cor", names(temp_DR_fisher_se      ))]          

#M3
m3_ss_redux$temp_DR_cor_fisher   <- FisherZ(m3_ss_redux$temp_DR_cor)
temp_DR_fisher_se    <- CorCI(m3_ss_redux$temp_DR_cor,   m3_ss_redux$n_extant_diversity, conf.level = 0.95, alternative = c("two.sided"))
m3_ss_redux$temp_DR_fisher_se        <-      temp_DR_fisher_se[grepl("cor", names(temp_DR_fisher_se      ))]          

#M0
m0_ss_redux$temp_DR_cor_fisher   <- FisherZ(m0_ss_redux$temp_DR_cor)
temp_DR_fisher_se    <- CorCI(m0_ss_redux$temp_DR_cor,   m0_ss_redux$n_extant_diversity, conf.level = 0.95, alternative = c("two.sided"))
m0_ss_redux$temp_DR_fisher_se        <-      temp_DR_fisher_se[grepl("cor", names(temp_DR_fisher_se      ))]          

# Fit linear mixed effect models to perform meta-analyses
meta_temp_dr0 <- rma.uni(yi=temp_DR_cor_fisher,   sei=temp_DR_fisher_se,    data=m0_ss_redux)
meta_temp_dr1 <- rma.uni(yi=temp_DR_cor_fisher,   sei=temp_DR_fisher_se,    data=m1_ss_redux)
meta_temp_dr2 <- rma.uni(yi=temp_DR_cor_fisher,   sei=temp_DR_fisher_se,    data=m2_ss_redux)
meta_temp_dr3 <- rma.uni(yi=temp_DR_cor_fisher,   sei=temp_DR_fisher_se,    data=m3_ss_redux)

# summarise
summary(meta_temp_dr0)
summary(meta_temp_dr1)
summary(meta_temp_dr2)
summary(meta_temp_dr3)

#############################################################################
#############################################################################
###                                                                       ###
###  SECTION 3: COMPARISON OF EMPIRICAL AND SIMULATED SUMMARY STATISTICS  ###
###                                                                       ###
#############################################################################
#############################################################################

# add additional class variable
e_ss$class <- "empirical"
m_ss$class <- "simulated"
# combine and 'melt' the data for plotting
data <- rbind(e_ss, m_ss)
data[, 2:(ncol(data)-1)] <- scale(data[, 2:(ncol(data)-1)])
melt_data <- melt(data)

# summary function for the mean and range
data_summary <- function(x) {
  m <- mean(x)
  ymin <- min(x)
  ymax <- max(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# plot the mean and range of each summary statistic
# Figure S9
S9 <- ggplot(data=melt_data,aes(y=value, x=variable,fill=class, colour=class))+ 
  stat_summary(fun.data=data_summary,position=position_dodge(0.8))+
  coord_flip()+ 
  theme_light()+
  scale_color_manual(values=c( "#E69F00", "#56B4E9"))

# now lets quantify the overlap for each summary statsitic
# get the proportion of the empirical data that falls within the range of the simulated data
withinD <- function(x, y){
  in_between <- which(x[which(y %in% "empirical")] < max(x[which(y %in% "simulated")], na.rm=T) &
        x[which(y %in% "empirical")] > min(x[which(y %in% "simulated")], na.rm=T) )
  return(length(in_between)/ length(which(!is.na(x[which(y %in% "empirical")]))))
}

# find proportion of overlap
prop_inbetween <- sapply(data, FUN=function(x){withinD(x, y=data$class)})

# find number of variables that have > 95% overlap
length(which(prop_inbetween ==1))
length(which(prop_inbetween > 0.95))
length(which(prop_inbetween < 0.95))
prop_inbetween[which(prop_inbetween < 0.95)] # which variables?


