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
# Description: Script to perform meta-analysis of key correlations in tetrapods
#
# Contact: alexander.skeels@gmail.com
#
######################################

###########################################################################
###########################################################################
###                                                                     ###
###                    SECTION 1: LIBRARIES AND DATA                    ###
###                                                                     ###
###########################################################################
###########################################################################

# libraries
library(metafor)
library(DescTools)

# read in empirical summary statistics
empirical_data <- read.csv("empirical_summary_statistics.csv")

# just want to look at diverse clades
empirical_data <- na.omit(empirical_data[which(empirical_data$n_species > 20),])
empirical_data$taxon <- tolower(empirical_data$taxon)


############################################################################
############################################################################
###                                                                      ###
###               SECTION 2: META-ANALYSIS OF CORRELATIONS               ###
###                                                                      ###
############################################################################
############################################################################
meta_df <- empirical_data

# Fisher transfromation of rho values
meta_df$bodysize_DivRate_fisher <- FisherZ(meta_df$bodysize_DR_cor)
meta_df$temp_DivRate_fisher <- FisherZ(meta_df$temp_DR_cor)
meta_df$richness_lat_cor_fisher <- FisherZ(meta_df$richness_lat_cor)
meta_df$richness_temp_cor_fisher <- FisherZ(meta_df$richness_temp0mya_cor)

# get CIs from cor and n
bodysize_DivRate_fisher_se <- CorCI(meta_df$bodysize_DR_cor,   meta_df$n_species, conf.level = 0.95, alternative = c("two.sided"))
temp_DivRate_fisher_se <- CorCI(meta_df$temp_DR_cor,       meta_df$n_species, conf.level = 0.95, alternative = c("two.sided"))
richness_lat_cor_fisher_se <- CorCI(meta_df$richness_lat_cor,      meta_df$n_species, conf.level = 0.95, alternative = c("two.sided"))
richness_temp_cor_fisher_se <- CorCI(meta_df$richness_temp0mya_cor, meta_df$n_species, conf.level = 0.95, alternative = c("two.sided"))

meta_df$bodysize_DivRate_fisher_se <- bodysize_DivRate_fisher_se[grepl("cor", names(bodysize_DivRate_fisher_se  ))]   
meta_df$temp_DivRate_fisher_se <- temp_DivRate_fisher_se[grepl("cor", names(temp_DivRate_fisher_se      ))]          
meta_df$richness_lat_cor_fisher_se <- richness_lat_cor_fisher_se[grepl("cor", names(richness_lat_cor_fisher_se  ))]
meta_df$richness_temp_cor_fisher_se <- richness_temp_cor_fisher_se[grepl("cor", names(richness_temp_cor_fisher_se ))]

# fit 
meta_bs_dr <- rma.uni(yi=bodysize_DivRate_fisher, sei=bodysize_DivRate_fisher_se, data=meta_df)
meta_t_dr <- rma.uni(yi=temp_DivRate_fisher, sei=temp_DivRate_fisher_se,data=meta_df)
meta_r_l  <- rma.uni(yi=richness_lat_cor_fisher,  sei = richness_lat_cor_fisher_se ,  data=meta_df)
meta_r_t  <- rma.uni(yi=richness_temp_cor_fisher, sei = richness_temp_cor_fisher_se,  data=meta_df)

summary(meta_bs_dr)
summary(meta_t_dr )
summary(meta_r_l )
summary(meta_r_t )
