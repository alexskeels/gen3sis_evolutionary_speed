######################################
###            METADATA            ###
######################################
#
# Author: Alexander Skeels
#
# Date: 10.01.2022
#
# Projects: BIGEST: Evolutionary Speed Hypothesis (ESH)
#
# Description: Script to perform simulation based inference model selection
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

# directories
setwd("PATH/TO/DATA")

# libraries
library(caret)
library(corrplot)
library(ggplot2)
library(reshape2)
library(caret)
library(MASS)
library(klaR)
library(gplots)
library(randomForest)
require(caTools)
library(fields)
library(gridExtra)
library(caret)
library(gbm)
library(metafor)
library(DescTools)

# data
m0_ss <- read.table("simulation_m0_summary_stats.txt", header=T)
m1_ss <- read.table("simulation_m1_summary_stats.txt", header=T)
m2_ss <- read.table("simulation_m2_summary_stats.txt", header=T)
m3_ss <- read.table("simulation_m3_summary_stats.txt", header=T)

# read in empirical summary statistics
empirical_data <- read.csv("order_empirical_summary_statistics.csv")

############################################################################
############################################################################
###                                                                      ###
###                      SECTION 2: DATA PROCESSING                      ###
###                                                                      ###
############################################################################
############################################################################

# following CARET Package and advice from author Max Kuhn 
# https://topepo.github.io/caret/
# https://www.r-project.org/conferences/useR-2013/Tutorials/kuhn/user_caret_2up.pdf

# give models an identifier
m0_ss$m <- "m0"
m1_ss$m <- "m1"
m2_ss$m <- "m2"
m3_ss$m <- "m3"

#first subset to only use parameters that were successful for each model
m0_ss_redux <- m0_ss[which(m0_ss$running_time_step==0 & m1_ss$running_time_step==0 & m2_ss$running_time_step==0 & m3_ss$running_time_step==0),]
m1_ss_redux <- m1_ss[which(m0_ss$running_time_step==0 & m1_ss$running_time_step==0 & m2_ss$running_time_step==0 & m3_ss$running_time_step==0),]
m2_ss_redux <- m2_ss[which(m0_ss$running_time_step==0 & m1_ss$running_time_step==0 & m2_ss$running_time_step==0 & m3_ss$running_time_step==0),]
m3_ss_redux <- m3_ss[which(m0_ss$running_time_step==0 & m1_ss$running_time_step==0 & m2_ss$running_time_step==0 & m3_ss$running_time_step==0),]

# combine model data frames
m_ss <- rbind(m0_ss_redux, m1_ss_redux, m2_ss_redux, m3_ss_redux)

# subset to the 54 variables which we can also estimate for empirical data
m_ss_subset <- na.omit(m_ss[,c(142,21, 30:31, 40:41, 48:49, 52:74, 93:104, 114:125)])

# remove the 7 summary statsitics that showed poor match to the empirical data based on validation procedure (ESH_simulation_validation.R script for details)
m_ss_subset <- m_ss_subset[, which(!colnames(m_ss_subset) %in% c("lat_bssd_cor","bs_skewness","bs_kurtosis",
                                                                 "lat_bsm_cor", "pd_bssd_cor",  "mpd_bssd_cor", "mntd_bssd_cor"))]

ncol(m_ss_subset) # comparing 47 summary statistics
nrow(m_ss_subset) # 1384 simulations

# remove highly colinear variables 
correlations <- cor(m_ss_subset[, 2:ncol(m_ss_subset)])
highly_correlated <- sort(findCorrelation(correlations, cutoff=0.90, names=T)) 
m_ss_subset <- m_ss_subset[,which(!colnames(m_ss_subset) %in% highly_correlated)]
ncol(m_ss_subset) # 26 summarys statistics

### apply classification models ###

# formula - using all non-correlated variables
f1 <- formula(paste("m ~ ", paste(names(m_ss_subset)[2:c(length(names(m_ss_subset)))], collapse=" + ")))

# index train and test data
set.seed(666)
m_ss_subset$m <- as.factor(m_ss_subset$m)
train_index <- createDataPartition(m_ss_subset$m, p = .66, list = FALSE,  times = 1)
train_data <- m_ss_subset[ train_index ,]
test_data  <- m_ss_subset[-train_index ,]

# preprocess values - we will scale and center values
preprocessed_values <- preProcess(train_data, method = c("center", "scale"))
train_transformed   <- predict(preprocessed_values, train_data)
test_transformed    <- predict(preprocessed_values, test_data )


###########################################################################
###########################################################################
###                                                                     ###
###               SECTION 3: CLASSIFICATION MODEL FITTING               ###
###                                                                     ###
###########################################################################
###########################################################################
        
# configure the cross-validation paramaters
train_control <- trainControl( method = "repeatedcv", number = 10, repeats = 10, classProbs = TRUE, savePredictions = TRUE)

# FIT MODELS

## LINEAR DISCRIMINANT ANALYSIS
lda_train        <- train(f1, data=train_transformed, method = "lda", trControl = train_control, verbose = T)

## AVERAGED NEURAL NET ALGORITHM
avNNet_train     <- train(f1, data=train_transformed,  method = "avNNet", trControl = train_control, verbose = T)

## RECURSIVE PARITIONING AND REGRESSION TREES
rpart_train      <- train(f1, data=train_transformed,  method = "rpart1SE", trControl = train_control)

## GRADIENT BOOSTED MODEL
gbm_train        <- train(f1, data=train_transformed,  method = "gbm", trControl = train_control, verbose = T)

## RANDOM FOREST
rf_train         <- train(f1, data=train_transformed,  method = "rf", trControl = train_control, verbose = T)

## SUPPORT-VECTOR MACHINES
svmLinear_train  <- train(f1, data=train_transformed,  method = "svmLinear", trControl = train_control, verbose = T)

## NAIVE BAYES
naivebayes_train <- train(f1, data=train_transformed,  method = "naive_bayes", trControl = train_control, verbose = T)

#predict on test data
lda_test        <- predict(lda_train, test_transformed)
rf_test         <- predict(rf_train, test_transformed)
rpart_test      <- predict(rpart_train, test_transformed)
avNNet_test     <- predict(avNNet_train, test_transformed)
gbm_test        <- predict(gbm_train, test_transformed)
svmLinear_test  <- predict(svmLinear_train, test_transformed)
naivebayes_test <- predict(naivebayes_train, test_transformed)

# classification accuracy
lda_cm        <- confusionMatrix(data = lda_test, reference = test_transformed$m, mode = "prec_recall")
rf_cm         <- confusionMatrix(data = rf_test, reference = test_transformed$m, mode = "prec_recall")
rpart_cm      <- confusionMatrix(data = rpart_test, reference = test_transformed$m, mode = "prec_recall")
avNNet_cm     <- confusionMatrix(data = avNNet_test, reference = test_transformed$m, mode = "prec_recall")
gbm_cm        <- confusionMatrix(data = gbm_test, reference = test_transformed$m, mode = "prec_recall")
svmLinear_cm  <- confusionMatrix(data = svmLinear_test , reference = test_transformed$m, mode = "prec_recall")
naivebayes_cm <- confusionMatrix(data = naivebayes_test, reference = test_transformed$m, mode = "prec_recall")


# how well did each model perform
cvValues <- resamples(list(LDA = lda_train, 
                           RF  = rf_train,    
                           rpart = rpart_train, 
                           NeuralNet = avNNet_train,  
                           GBM = gbm_train,
                           SVM = svmLinear_train,
                           NBayes = naivebayes_train))
summary(cvValues)

# Figure S14
dotplot(cvValues, metric="Kappa") 
dotplot(cvValues, metric="Accuracy") 

# differences in kappa from training data
kappaDiffs <- diff(cvValues, metric = "Kappa")
summary(kappaDiffs)

# Variable importance
lda_varI <- varImp(lda_train, scale = T)
rf_varI <- varImp(rf_train, scale = T)
rpart_varI <- varImp(rpart_train, scale = T)
avNNet_varI <- varImp(avNNet_train, scale = T)
gbm_varI <- varImp(gbm_train, scale = T)
svmLinear_varI <- varImp(svmLinear_train, scale = T)
naivebayes_varI <- varImp(naivebayes_train, scale =T)

# model weights
# weight the models by their kappa values
classification_df <- data.frame(model = c("lda", "rf", "rpart","avNNet", "gbm", "SVMLinear","naiave_bayes"), 
                                accuracy = c(lda_cm$overall[1], rf_cm$overall[1], rpart_cm$overall[1], avNNet_cm$overall[1], gbm_cm$overall[1], svmLinear_cm$overall[1], naivebayes_cm$overall[1]), 
                                kappa = c(lda_cm$overall[2], rf_cm$overall[2], rpart_cm$overall[2], avNNet_cm$overall[2], gbm_cm$overall[2], svmLinear_cm$overall[2], naivebayes_cm$overall[2]))

weights <- classification_df$kappa / sum(classification_df$kappa)
 
# weight the variable importance by Kappa
lda_varI_weights <- rowMeans(lda_varI$importance * weights[1])
rf_varI_weights <- rowMeans(rf_varI$importance * weights[2])
rpart_varI_weights <- rowMeans(rpart_varI$importance * weights[3])
avNNet_varI_weights <- rowMeans(avNNet_varI$importance * weights[4])
gbm_varI_weights <- rowMeans(gbm_varI$importance * weights[5])
svmLinear_varI_weights <- rowMeans(svmLinear_varI$importance * weights[6])
naivebayes_varI_weights <- rowMeans(naivebayes_varI$importance* weights[7])

# order the vectors
lda_varI_weights <- lda_varI_weights[order(names(lda_varI_weights))]
rf_varI_weights <- rf_varI_weights[order(names(rf_varI_weights))]
rpart_varI_weights <- rpart_varI_weights[order(names(rpart_varI_weights))]
avNNet_varI_weights <- avNNet_varI_weights[order(names(avNNet_varI_weights))]
gbm_varI_weights <- gbm_varI_weights[order(names(gbm_varI_weights))]
svmLinear_varI_weights <- svmLinear_varI_weights [order(names(svmLinear_varI_weights ))]
naivebayes_varI_weights <- naivebayes_varI_weights[order(names(naivebayes_varI_weights))]

# combine and place in data frame
weightedVIF <- c(lda_varI_weights, rf_varI_weights, rpart_varI_weights, avNNet_varI_weights,
                 gbm_varI_weights, svmLinear_varI_weights, naivebayes_varI_weights)
weightedVIF <- data.frame(VIF=weightedVIF[order(weightedVIF, decreasing=T)])


###########################################################################
###########################################################################
###                                                                     ###
###     SECTION 4: CLASSIFICATION MODEL SELECTION ON EMPIRICAL DATA     ###
###                                                                     ###
###########################################################################
###########################################################################

# first subset and clean the colnames of the empirical data
# just want to look at diverse clades
empirical_data <- na.omit(empirical_data[which(empirical_data$n_species >= 20),])
empirical_data$taxon <- tolower(empirical_data$taxon)
colnames(empirical_data)[which(colnames(empirical_data) == "taxon")] <- "m"
colnames(empirical_data)[which(colnames(empirical_data)=="rs_kutosis")] <- "rs_kurtosis"
colnames(empirical_data)[which(colnames(empirical_data)=="n_species")] <- "n_extant_diversity"
colnames(empirical_data) <- gsub("_p_cor", "_cor",colnames(empirical_data)) # change _p_cor for posterior samplescould also change _m_cor to use MCC samples
colnames(empirical_data) <- gsub("DivRate", "DR",colnames(empirical_data))
colnames(empirical_data)[which(colnames(empirical_data) == "taxon")] <- "m"
colnames(empirical_data)[which(colnames(empirical_data) == "collessI_post")] <- "collessI"
colnames(empirical_data)[which(colnames(empirical_data) == "sackinI_mcc")] <- "sackinI"
colnames(empirical_data)[which(colnames(empirical_data) == "gamma_mcc")] <- "gamma"

empirical_data <- empirical_data[, which(colnames(empirical_data) %in% colnames(m_ss_subset))]
empirical_subset <- empirical_data[, match(colnames(m_ss_subset), colnames(empirical_data))]

# check names match
colnames(m_ss_subset)[which(!colnames(m_ss_subset) %in% colnames(empirical_subset))]
colnames(empirical_subset)[which(!colnames(empirical_subset) %in% colnames(m_ss_subset))]
colnames(empirical_subset) == colnames(m_ss_subset)


# Process empirical data in the same way as for the simulated data
simulated_transformed <- predict(preprocessed_values, m_ss_subset )
empirical_transformed <- predict(preprocessed_values, empirical_subset )

model_set <- list(lda_train, rf_train, rpart_train, avNNet_train,
                  gbm_train, svmLinear_train, naivebayes_train)

# predict on empirical
class_predictions   <- predict(model_set, newdata = empirical_transformed, type = "raw", na.action = na.omit)
class_probabilities <- predict(model_set, newdata = empirical_transformed, type = "prob", na.action = na.omit)

# weight by kappa
class_prediction_table <- lapply(class_predictions, table)

# get the model probabilities based on their weights
for(i in 1:7){
  class_prediction_table[[i]] <- class_prediction_table[[i]]* weights[i]
}

# sum the classes
colSums(do.call(rbind, class_prediction_table))

# get the model probabilities based on their weights
for(i in 1:7){
  class_probabilities[[i]] <- class_probabilities[[i]]* weights[i]
}

weighted_probabilities <- Reduce('+', class_probabilities)


