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
# Description: Script to perform sensitivity analysis following Prowse et al (2016, Ecosphere, 7:3, e01238, doi.org/10.1002/ecs2.1238)
#
# Contact: alexander.skeels@gmail.com
#
######################################

###########################################################################
###########################################################################
###                                                                     ###
###                    SECTION 1: DATA AND LIBRARIES                    ###
###                                                                     ###
###########################################################################
###########################################################################

# directories
setwd("PATH/TO/DATA")

#libraries
library(caret)
library(compiler)
library(dismo)
library(ggplot2)
library(reshape)
library(scico)
library(gridExtra)
library(MDM)
library(corrplot)

# data
m0_ss <- read.table("simulation_m0_summary_stats.txt", header=T)
m1_ss <- read.table("simulation_m1_summary_stats.txt", header=T)
m2_ss <- read.table("simulation_m2_summary_stats.txt", header=T)
m3_ss <- read.table("simulation_m3_summary_stats.txt", header=T)

# give models an identifier
m0_ss$m <- "m0"
m1_ss$m <- "m1"
m2_ss$m <- "m2"
m3_ss$m <- "m3"

# combine
dataset <- rbind(m0_ss, m1_ss, m2_ss, m3_ss)
dataset$m <- as.factor(dataset$m)

# subset only the variables we can extract for both empirical and simulated datasets as well as the model parameters
# also remove simulations with missing observations as these were those which did not run successfully  (na.omit)
dataset <- na.omit(dataset[,c(142,1:6,21, 30:31, 40:41, 48:49, 52:74, 93:104, 114:125)])

# create a datset to look at correlations
corrplot_dataset <- dataset[8:61][sort(names(dataset)[8:61])]
# Figure S13.
corrplot(cor(corrplot_dataset), method="square", type="lower")

############################################################################
############################################################################
###                                                                      ###
###    SECTION 2: SENSITIVITY ANALYSIS USING BOOSTED REGRESSION TREES    ###
###                                                                      ###
############################################################################
############################################################################

# using boosted regression trees (BRT) in sensitivity analysis
# lets loop over a number of sampling levels to see the effect of N
# for each sample size (1: j_sample) we will randomly subset the simulated dataset and fit BRT for each summary statistic ~ 7 model parameters

sensitivity_list <- vector("list", 20)
sample_list <- vector("list", 20)

sample_number <- round(seq(from=50, to=nrow(dataset), length.out = 20))

for(j_sample in 1:20){
  
  # print to keep track
  print(paste("sample  ::  ", j_sample))
  
  # model parameters are predictor variables
  predictor_vars <- c("divergence_threshold", "lambda", "omega", "sigma_bs", "sigma_t", "dispersal", "m")
  var_list <- list()
  
  # subsample the dataset
  dataset_tmp <- dataset[runif(sample_number[j_sample], 0, nrow(dataset)),]
  
  # for each variable we will fit a BRT with a learning rate of 0.01 and a tree complexity of 1
  for(i_var in 8:61) {
    
    learn_rate <- 0.01
    # response var is summary statistic i
    response_var <- names(dataset_tmp)[i_var]
    brt_fit <- NULL
    x_col <- which(names(dataset_tmp)%in% predictor_vars)
    y_col <- which(names(dataset_tmp)==response_var)
    brt_dataset <- dataset_tmp[!is.na(dataset_tmp[,y_col]),]
    # fit the BRT model
    brt_fit <- try(gbm.step( learning.rate = learn_rate, data=brt_dataset, gbm.x=x_col, gbm.y=y_col,family="gaussian",tree.complexity=1,n.folds=5,tolerance.method='auto',max.trees=200000))
    
    # increase learning rate for unsuccessful fits
    while(is.null(brt_fit)){
      learn_rate <- learn_rate-0.001
      brt_fit <- try(gbm.step( learning.rate = learn_rate, data=brt_dataset, gbm.x=x_col, gbm.y=y_col,family="gaussian",tree.complexity=1,n.folds=5,tolerance.method='auto',max.trees=200000))
    }
    var_list[[i_var]] <- brt_fit
  }
  var_list <-  var_list[8:61]
  
  sample_list[[j_sample]] <- var_list
  # for each sample size, we will also extract model sensitivity information including validation statistics and R2 
  # first create a data frame to store the infromation
  sensitivity_df <- data.frame(var= names(dataset_tmp)[8:60], omega=NA, sigma_squared_t=NA, sigma_squared_bs=NA,
                                dispersal=NA, lambda=NA, divergence_threshold=NA, m=NA, R2=NA, cvDev=NA)
  # then loop over each summary statitic (i_var) and calultae the sensitivity infromation
  for(i_var in 1:54){
    
    if(class(var_list[[i_var]]) == 'try-error'){next}
    cont <- var_list[[i_var]]$contributions # variable contributions
    if(is.null(cont)){next}
    sensitivity_df$omega[i_var]  <- cont[which(cont$var == "omega"), 2]
    sensitivity_df$sigma_squared_t[i_var]  <- cont[which(cont$var == "sigma_t"), 2]
    sensitivity_df$sigma_squared_bs[i_var]  <- cont[which(cont$var == "sigma_bs"), 2]
    sensitivity_df$dispersal[i_var]  <- cont[which(cont$var == "dispersal"), 2]
    sensitivity_df$lambda[i_var]  <- cont[which(cont$var == "lambda"), 2]
    sensitivity_df$divergence_threshold[i_var]  <- cont[which(cont$var == "divergence_threshold"), 2]
    sensitivity_df$m[i_var]  <- cont[which(cont$var == "m"), 2]
    sensitivity_df$R2[i_var]  <- (cor(var_list[[i_var]]$fit, var_list[[i_var]]$data$y))^2
    sensitivity_df$cvDev[i_var] <- var_list[[i_var]]$cv.statistics$deviance.mean
  }
  
  sensitivity_list[[j_sample]] <- sensitivity_df
}

###########################################################################
###########################################################################
###                                                                     ###
###                SECTION 3: STABILITY WITH SAMPLE SIZE                ###
###                                                                     ###
###########################################################################
###########################################################################

# here we want to see how our estimates of relative variable contributions to the BRT models change with increasing sample size
# this is measured using beta diversity metrics of relative importance factors

# variable names
varss <- names(dataset_tmp)[8:ncol(dataset_tmp)]
beta_div_list <- list() 

# loop over summary statistics
for (i_var in 1:54) {
  
  var <-varss[i_var]
  tmp_df <- data.frame(var=colnames(sensitivity_df)[2:8])
  
  for(j_sample in 1:20){
    
    sensitivity_df <- sensitivity_list[[j_sample]]
    if(all(is.na(as.numeric(sensitivity_df[i_var,2:8])))){
      tmp_df <- cbind(tmp_df, as.numeric(tmp_df[, ncol(tmp_df)]))
    } else {
      tmp_df <- cbind(tmp_df, data.frame(rel.inf= as.numeric(sensitivity_df[i_var,2:8])))
    }
  }
  
  betaDiv.dummy <- data.frame(subsample=sample_number,resp=NA)
  for (j in 2:length(sample_number)) {
    if(all(is.na(tmp_df[,j]))){betaDiv.dummy[j,2] <- NA; next}
    beta.div <- as.numeric(ed(t(as.matrix(tmp_df[,j:(j+1)])),q=1,retq=T)['beta'])
    betaDiv.dummy[j,2] <- beta.div    
  }
  beta_div_list[[i_var]] <- betaDiv.dummy
}
names(beta_div_list) <- varss


###########################################################################
###########################################################################
###                                                                     ###
###                     SECTION 4: PLOTTING RESULTS                     ###
###                                                                     ###
###########################################################################
###########################################################################

# firstly look at the contribution of each model parameter to variance in summary statsitics
sensitivity_df <- sensitivity_df[which(!sensitivity_df$var %in% c("lat_bssd_cor","bs_skewness","bs_kurtosis",
                                                                  "lat_bsm_cor", "pd_bssd_cor",  "mpd_bssd_cor", "mntd_bssd_cor")),]
# melt data for plotting
cont_melt <- melt(sensitivity_df[, 1:8])

# plot the variable contributions for each statistic
var_cont_p <- ggplot(cont_melt, aes(x=var, y=value))+
  geom_bar(aes(fill=variable),position="fill", stat="identity")+
  scale_fill_scico_d(palette = 'hawaii')+
  theme(axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# plot the R2 of each summary statistic as  bar plot
r2_p <- ggplot(sensitivity_df, aes(x=var, y=R2))+
  geom_bar(stat="identity")+
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# plot together - Figure S10
grid.arrange(var_cont_p, r2_p)


# look at the validation scores (R2) to see how much variation in the biodiversity patterns we can attribute to the parameters
# see how this varies with sample size

rsquared_df <- data.frame(sensitivity_list[[1]]$var, lapply(sensitivity_list, FUN=function(x){x$R2}))
names(rsquared_df) <- c("var",paste0(sample_number))
melt_rsquared <- melt(rsquared_df)

# plot them all together
s11_p1 <- ggplot(melt_rsquared, aes(y=value, x=variable, fill=var, group=var))+
  geom_point(alpha=0.5)+
  geom_line(alpha=0.7)+
  theme_light()+
  scale_color_scico_d(palette = 'hawaii')+
  xlab("Sample size")+
  ylab("Validation metric (pseudo R-squared)")+
  theme(legend.title = element_blank(), legend.position = "none")

s11_p2 <- ggplot(melt_rsquared, aes(y=value, x=variable, colour=variable))+
  geom_boxplot()+
  scale_colour_viridis_d(alpha=0.7)+
  geom_jitter(width=0.15)+
  xlab("Sample size")+
  ylab("Validation metric (pseudo R-squared)")+
  theme_light()+
  theme(legend.title = element_blank(), legend.position = "none", axis.text.x =element_text(angle = 90))

grid.arrange(s11_p1, s11_p2)


# See how the stability of the contribution of model parameters to summary statsitics varies with sample size

for(i in 1:54){
  beta_div_list[[i]]$var <- varss[[i]]
}

relative_importance_df <-  na.omit(do.call(rbind, beta_div_list))

s12_p1 <- ggplot(relative_importance_df , aes(y=resp, x=as.factor(subsample), group=var))+
  scale_colour_viridis_d(alpha=0.5)+
  geom_jitter(width=0.15, alpha=0.2)+
  geom_line(alpha=0.5)+
  xlab("Sample size")+
  ylab("Stability")+
  theme_light()+
  theme(legend.title = element_blank(), legend.position = "none", axis.text.x =element_text(angle = 90))

s12_p2 <- ggplot(relative_importance_df , aes(y=resp, x=as.factor(subsample), colour=as.factor(subsample)))+
  scale_colour_viridis_d(alpha=0.5)+
  geom_boxplot()+
  geom_jitter(width=0.15)+
  xlab("Sample size")+
  ylab("Stability")+
  theme_light()+
  theme(legend.title = element_blank(), legend.position = "none", axis.text.x =element_text(angle = 90))

grid.arrange(s12_p1, s12_p2)




