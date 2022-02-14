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

# set up data wd
data_home <- "PATH/TO/DATA"

# libraries
library(metafor)
library(DescTools)
library(nlme)
library(ape)
library(phylobase)
library(phylosignal)
library(raster)

# read in empirical summary statistics
empirical_data <- read.csv("empirical_summary_statistics.csv")

# just want to look at diverse clades
empirical_data <- na.omit(empirical_data[which(empirical_data$n_species >= 20),])
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


###########################################################################
###########################################################################
###                                                                     ###
###       SECTION 3: AUTOCORRELATION SENSITIVITY AND VERIFICATION       ###
###                                                                     ###
###########################################################################
###########################################################################

# get clade file names
setwd(file.path(data_home))
taxon_files <- list.files()

# seperate the DR stats and the spatial stats extracted at 220km
taxon_files_220 <- taxon_files[grepl("220", taxon_files)] 
taxon_files_DR  <- taxon_files[grepl("_DR", taxon_files)] 

# load phylo data
setwd(data_home)
mam  <- readRDS("spatial_phylo_cleaned_mammal_data.rds")
bird <- readRDS("spatial_phylo_cleaned_bird_data.rds")
amph <- readRDS("spatial_phylo_cleaned_amphibian_data.rds")
squa <- readRDS("spatial_phylo_cleaned_squamate_data.rds")
croc <- readRDS("spatial_phylo_cleaned_crocturt_data.rds")

# set up df
#taxon_data <- data.frame(taxon=NA, n_species=NA)
setwd(data_store)
taxon_data <- read.csv("logBS_pgls_data.csv")

for(i in 1:length(taxon_files_DR)){
  
  ##~~~~~~~ LOAD AND CLEAN DATA ~~~~~~~~#
  print(i)
  setwd(data_store)
  data_set_220 <- read.csv(taxon_files_220[i])
  data_set_DR  <- read.csv(taxon_files_DR[i])
  
  clade <- NA
  if("Mass.g" %in% colnames(data_set_DR)){clade <- "mammal"}
  if("Body_length_mm" %in% colnames(data_set_DR)){clade <- "amphibian"}
  if("Mass" %in% colnames(data_set_DR)){clade <- "croc"}
  if("maximum_SVL" %in% colnames(data_set_DR)){clade <- "lizard"}
  if("BodyMass.Value" %in% colnames(data_set_DR)){clade <- "bird"}
  if("body_mass_g" %in% colnames(data_set_DR)){clade <- "snake"}
  
  # change the body mass / body size column name
  colnames(data_set_DR)[which(colnames(data_set_DR) %in% c("Mass.g", "Body_length_mm", "BodyMass.Value", "maximum_SVL", "Mass", "body_mass_g"))] <- "body_size"
  
  # scale SVL to body mass in squamates
  if(any(colnames(data_set_DR) %in% c("intercept"))){
    data_set_DR$body_size[which(data_set_DR$body_size %in% "no measurements")] <- NA
    data_set_DR$body_size[which(data_set_DR$body_size %in% "juveniles_only")] <- NA
    data_set_DR$body_size[which(data_set_DR$body_size %in% "known from a single specimen, posterior part of the body missing (Ribeiro et al. 2009)")] <- NA
    data_set_DR$body_size <- data_set_DR$intercept + data_set_DR$slope * as.numeric(data_set_DR$body_size)
  }
  
  data_set_DR <- data_set_DR[, c(2,9:ncol(data_set_DR))]
  data_set <- merge(data_set_DR, data_set_220, by="species")
  data_set <- data_set[!duplicated(data_set$species),]
  rownames(data_set) <- data_set$species
  
  # taxon name and number of species
  taxon_data[i, "taxon"] <- strsplit(taxon_files_DR[i], "_")[[1]][1]
  taxon_data[i, "n_species"] <- nrow(data_set)
  taxon_data$clade <- clade
  
  # skip baby clades
  if(nrow(data_set) < 20){next}
  
  # standardise
  data_set_tmp <- data_set
  test1 <- try(scale(data_set_tmp[,c("temp_mean", "DR", "body_size")]))
  
  if(class( test1)=='try-error' | any(is.na(data_set_tmp$body_size))){
    data_set_tmp <- data_set_tmp[which(!is.na(data_set_tmp$body_size)),]
  } 
  data_set_tmp$body_size <- log(data_set_tmp$body_size)
  data_set_tmp[,c("temp_mean", "DR", "body_size")] <- scale(data_set_tmp[,c("temp_mean", "DR", "body_size")])
  
  ## load phylos
  if(clade == "mammal"){    phy_set <- mam$phy}
  if(clade == "bird"){      phy_set <- bird$phy}
  if(clade == "amphibian"){ phy_set <- amph$phy}
  if(clade == "lizard"){    phy_set <- squa$phy}
  if(clade == "snake"){     phy_set <- squa$phy}
  if(clade == "croc"){      phy_set <- croc$phy}
  
  # set seed to sub sample tree
  set.seed(666)
  phy_index <- sample(1:100, 1)
  phy <- phy_set[[phy_index]]
  phy <- drop.tip(phy, phy$tip.label[which(!phy$tip.label %in% data_set_tmp$species)])
  phy <- multi2di(phy)
  
  # make data set match tips in the tree
  data_set_tmp <- data_set_tmp[which(data_set_tmp$species %in% phy$tip.label),]
  
  # Fit models
  pglsModel_1 <- gls(DR ~ temp_mean, correlation = corBrownian(phy = phy, form=~species), data = data_set_tmp, method = "ML")
  pglsModel_2 <- gls(DR ~ body_size, correlation = corBrownian(1, phy = phy, form=~species), data = data_set_tmp, method = "ML")
  
  # also get correlations
  cor_1 <- cor(data_set_tmp$DR, data_set_tmp$temp_mean, method="spearman")
  cor_2 <- cor(data_set_tmp$DR, data_set_tmp$body_size, method="spearman")

  # also calculate phylogenetic signal using Blomberg's K
  phy4d <- phylo4d(phy, tip.data=data_set_tmp[, c("DR", "body_size", "temp_mean")])
  signal_test <- phyloSignal(phy4d, methods = c("K"), reps = 999, W = NULL)
  
  # store model parameters
  taxon_data$DR_TE_pgls_int[i] <-   pglsModel_1$coefficients[1]
  taxon_data$DR_BS_pgls_int[i] <-   pglsModel_2$coefficients[1]
  
  taxon_data$DR_TE_pgls_slope[i] <-   pglsModel_1$coefficients[2]
  taxon_data$DR_BS_pgls_slope[i] <-   pglsModel_2$coefficients[2]
  
  taxon_data$DR_TE_pgls_se_int[i] <-   summary(pglsModel_1)$tTable[1,2]
  taxon_data$DR_BS_pgls_se_int[i] <-   summary(pglsModel_2)$tTable[1,2]
  
  taxon_data$DR_TE_pgls_se_slope[i] <-   summary(pglsModel_1)$tTable[2,2]
  taxon_data$DR_BS_pgls_se_slope[i] <-   summary(pglsModel_2)$tTable[2,2]
  
  taxon_data$DR_TE_pgls_t_int[i] <-   summary(pglsModel_1)$tTable[1,3]
  taxon_data$DR_BS_pgls_t_int[i] <-   summary(pglsModel_2)$tTable[1,3]
  
  taxon_data$DR_TE_pgls_t_slope[i] <-   summary(pglsModel_1)$tTable[2,3]
  taxon_data$DR_BS_pgls_t_slope[i] <-   summary(pglsModel_2)$tTable[2,3]
  
  taxon_data$DR_TE_pgls_P_int[i] <-   summary(pglsModel_1)$tTable[1,4]
  taxon_data$DR_BS_pgls_P_int[i] <-   summary(pglsModel_2)$tTable[1,4]
  
  taxon_data$DR_TE_pgls_P_slope[i] <-   summary(pglsModel_1)$tTable[2,4]
  taxon_data$DR_BS_pgls_P_slope[i] <-   summary(pglsModel_2)$tTable[2,4]
  
  taxon_data$DR_TE_pgls_sigma[i] <-   pglsModel_1$sigma
  taxon_data$DR_BS_pgls_sigma[i] <-   pglsModel_2$sigma
  
  taxon_data$DR_K[i] <- signal_test$stat[1,1] 
  taxon_data$BS_K[i] <- signal_test$stat[2,1] 
  taxon_data$TE_K[i] <- signal_test$stat[3,1] 
  
  taxon_data$DR_Kp[i] <- signal_test$pvalue[1,1] 
  taxon_data$BS_Kp[i] <- signal_test$pvalue[2,1] 
  taxon_data$TE_Kp[i] <- signal_test$pvalue[3,1] 
  
  write.csv(taxon_data, "logBS_pgls_data.csv")
}



# do same meta-analysis using the PGLS slopes

library(metafor)
meta_dr_te <- rma.uni(yi=taxon_data$DR_TE_pgls_slope, sei=taxon_data$DR_TE_pgls_se_slope)
meta_dr_bs <- rma.uni(yi=taxon_data$DR_BS_pgls_slope, sei=taxon_data$DR_BS_pgls_se_slope)

summary(meta_dr_bs)
summary(meta_dr_te)

plot(taxon_data$DR_TE_cor, taxon_data$DR_TE_pgls_slope)
summary(lm(taxon_data$DR_TE_cor ~ taxon_data$DR_TE_pgls_slope))
abline(lm(taxon_data$DR_TE_cor ~ taxon_data$DR_TE_pgls_slope))

plot(taxon_data$DR_TE_lm_slope, taxon_data$DR_TE_pgls_slope)
summary(lm(taxon_data$DR_TE_lm_slope ~ taxon_data$DR_TE_pgls_slope))
abline(lm(taxon_data$DR_TE_lm_slope ~ taxon_data$DR_TE_pgls_slope))

plot(taxon_data$DR_BS_cor, taxon_data$DR_BS_pgls_slope)
summary(lm(taxon_data$DR_BS_cor ~ taxon_data$DR_BS_pgls_slope))
abline(lm(taxon_data$DR_BS_cor ~ taxon_data$DR_BS_pgls_slope))

plot(taxon_data$DR_BS_lm_slope, taxon_data$DR_BS_pgls_slope)
summary(lm(taxon_data$DR_BS_lm_slope ~ taxon_data$DR_BS_pgls_slope))
abline(lm(taxon_data$DR_BS_lm_slope ~ taxon_data$DR_BS_pgls_slope))


#### now perform spatial autocorrelation analyses

# read WorldClim Bio01
bio01 <-  raster("bio01.tif")

# make it equal area at 220km x 220km  resolution
bio01_equal_area  <- projectRaster(bio01, 
                                   crs=crs('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +no_defs +towgs84=0,0,0'), 
                                   res=220000, 
                                   over=T)

for(i in 1:length(taxon_files_DR)){
  
  ##~~~~~~~ LOAD AND CLEAN DATA ~~~~~~~~#
  print(i)
  setwd(data_store)
  data_set_220 <- read.csv(taxon_files_220[i])
  data_set_DR  <- read.csv(taxon_files_DR[i])
  
  clade <- NA
  if("Mass.g" %in% colnames(data_set_DR)){clade <- "mammal"}
  if("Body_length_mm" %in% colnames(data_set_DR)){clade <- "amphibian"}
  if("Mass" %in% colnames(data_set_DR)){clade <- "croc"}
  if("maximum_SVL" %in% colnames(data_set_DR)){clade <- "lizard"}
  if("BodyMass.Value" %in% colnames(data_set_DR)){clade <- "bird"}
  if("body_mass_g" %in% colnames(data_set_DR)){clade <- "snake"}
  
  # change the body mass / body size column name
  colnames(data_set_DR)[which(colnames(data_set_DR) %in% c("Mass.g", "Body_length_mm", "BodyMass.Value", "maximum_SVL", "Mass", "body_mass_g"))] <- "body_size"
  # scale SVL to body mass in squamates
  if(any(colnames(data_set_DR) %in% c("intercept"))){
    data_set_DR$body_size[which(data_set_DR$body_size %in% "no measurements")] <- NA
    data_set_DR$body_size[which(data_set_DR$body_size %in% "juveniles_only")] <- NA
    data_set_DR$body_size[which(data_set_DR$body_size %in% "known from a single specimen, posterior part of the body missing (Ribeiro et al. 2009)")] <- NA
    data_set_DR$body_size <- data_set_DR$intercept + data_set_DR$slope * as.numeric(data_set_DR$body_size)
  }
  data_set_DR <- data_set_DR[, c(2,9:ncol(data_set_DR))]
  data_set <- merge(data_set_DR, data_set_220, by="species")
  data_set <- data_set[!duplicated(data_set$species),]
  rownames(data_set) <- data_set$species
  
  # taxon name and number of species
  taxon_data[i, "taxon"] <- strsplit(taxon_files_DR[i], "_")[[1]][1]
  taxon_data[i, "n_species"] <- nrow(data_set)
  taxon_data$clade <- clade
  
  # skip baby clades
  if(nrow(data_set) < 20){next}
  
  # standardise
  data_set_tmp <- data_set
  test1 <- try(scale(data_set_tmp[,c("temp_mean", "DR", "body_size")]))
  
  if(class( test1)=='try-error' | any(is.na(data_set_tmp$body_size))){
    data_set_tmp <- data_set_tmp[which(!is.na(data_set_tmp$body_size)),]
  } 
  data_set_tmp$body_size <- log(data_set_tmp$body_size)
  data_set_tmp[,c("temp_mean", "DR", "body_size")] <- scale(data_set_tmp[,c("temp_mean", "DR", "body_size")])
  
  ## load phylos
  if(clade == "mammal"){  sp <- mam$sp}
  if(clade == "bird"){   sp <- bird$sp}
  if(clade == "amphibian"){ sp <- amph$sp}
  if(clade == "lizard"){  sp <- as.data.frame(squa$sp)}
  if(clade == "snake"){   sp <-  as.data.frame(squa$sp)}
  if(clade == "croc"){   sp <- croc$sp}
  
  sp <- sp[c(1:2, which(colnames(sp) %in% data_set_tmp$species))]
  sp <- sp[which(rowSums(sp[, 3:ncol(sp)]) > 0),]
  
  if(clade %in% c("mammal", "bird", "amphibian")){
    x <- sp$x
    y <- sp$y
    spatial_coordinates <- SpatialPoints(cbind(x, y))
    crs(spatial_coordinates) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
    spatial_coordinates_equal_area <- spTransform(spatial_coordinates, '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +no_defs +towgs84=0,0,0')
  }
  
  if(clade %in% c("croc")){
    x <- sp$Longitude_x_
    y <- sp$Latitude_y_
    spatial_coordinates_equal_area <- SpatialPoints(cbind(x, y))
    crs(spatial_coordinates_equal_area) <- '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +no_defs +towgs84=0,0,0'
  }
  if(clade %in% c("lizard", "snake")){
    x <- sp$'Longitude(x)'
    y <- sp$'Latitude(y)'
    spatial_coordinates_equal_area <- SpatialPoints(cbind(x, y))
    crs(spatial_coordinates_equal_area) <- '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +no_defs +towgs84=0,0,0'
    
  }
  
  richness <- rowSums(sp[, 3:ncol(sp)])
  bio01_extraction <- extract(bio01_equal_area, spatial_coordinates_equal_area, df=T)
  df <- cbind(x, y, richness, bio01_extraction)
  df <- na.omit(df)
  # estimate autocorrelation here
  #.....
  df.dists <- as.matrix(dist(cbind(df$x, df$y)))
  df.dists.inv <- 1/df.dists
  diag(df.dists.inv) <- 0
  
  moran_richness <- Moran.I(log(df$richness), df.dists.inv)
  moran_temperature <- Moran.I(df$bio01_equal_area, df.dists.inv)
  
  taxon_data$moran_richness[i] <- moran_richness$observed
  taxon_data$moran_richnessP[i] <-  moran_richness$p.value
  
  taxon_data$moran_temp[i] <- moran_temperature$observed
  taxon_data$moran_tempP[i] <- moran_temperature$p.value
  
  # then subsample and do gls
  n <- 1000
  dfss <- try(df[sample(1:nrow(df), n, replace=F), ])
  while(class(dfss) =="try-error" ){
    n <- n-100
    dfss <- try(df[sample(1:nrow(df), n, replace=F), ])
  }
  
  df[, 3:ncol(df)] <- scale(df[, 3:ncol(df)])
  m1 <- try(gls(log(richness) ~ y, correlation = corGaus(form = ~x + y, nugget = TRUE), data = dfss))
  m2 <- try(gls(log(richness) ~ bio01_equal_area, correlation = corGaus(form = ~x + y, nugget = TRUE), data = dfss))
  
  if(class(m1) != 'try-error'){
    taxon_data$SR_LA_sgls_int[i] <-   m1$coefficients[1]
    taxon_data$SR_LA_sgls_slope[i] <-   m1$coefficients[2]
    taxon_data$SR_LA_sgls_se_int[i] <-   summary(m1)$tTable[1,2]
    taxon_data$SR_LA_sgls_se_slope[i] <-   summary(m1)$tTable[2,2]
    taxon_data$SR_LA_sgls_t_int[i] <-   summary(m1)$tTable[1,3]
    taxon_data$SR_LA_sgls_t_slope[i] <-   summary(m1)$tTable[2,3]
    taxon_data$SR_LA_sgls_P_int[i] <-   summary(m1)$tTable[1,4]
    taxon_data$SR_LA_sgls_P_slope[i] <-   summary(m1)$tTable[2,4]
    taxon_data$SR_LA_sgls_sigma[i] <-   m1$sigma
  } else {
    taxon_data$SR_LA_sgls_int[i] <-   NA
    taxon_data$SR_LA_sgls_slope[i] <-  NA
    taxon_data$SR_LA_sgls_se_int[i] <- NA
    taxon_data$SR_LA_sgls_se_slope[i] <NA
    taxon_data$SR_LA_sgls_t_int[i] <-  NA
    taxon_data$SR_LA_sgls_t_slope[i] <-NA
    taxon_data$SR_LA_sgls_P_int[i] <-  NA
    taxon_data$SR_LA_sgls_P_slope[i] <-NA
    taxon_data$SR_LA_sgls_sigma[i] <-  NA
  }
  
  if(class(m2) != 'try-error'){
    taxon_data$SR_TE_sgls_int[i] <-   m2$coefficients[1]
    taxon_data$SR_TE_sgls_slope[i] <-   m2$coefficients[2]
    taxon_data$SR_TE_sgls_se_int[i] <-   summary(m2)$tTable[1,2]
    taxon_data$SR_TE_sgls_se_slope[i] <-   summary(m2)$tTable[2,2]
    taxon_data$SR_TE_sgls_t_int[i] <-   summary(m2)$tTable[1,3]
    taxon_data$SR_TE_sgls_t_slope[i] <-   summary(m2)$tTable[2,3]
    taxon_data$SR_TE_sgls_P_int[i] <-   summary(m2)$tTable[1,4]
    taxon_data$SR_TE_sgls_P_slope[i] <-   summary(m2)$tTable[2,4]
    taxon_data$SR_TE_sgls_sigma[i] <-   m2$sigma
  } else{
    taxon_data$SR_TE_sgls_int[i] <-   NA
    taxon_data$SR_TE_sgls_slope[i] <-   NA
    taxon_data$SR_TE_sgls_se_int[i] <-  NA
    taxon_data$SR_TE_sgls_se_slope[i] <-NA
    taxon_data$SR_TE_sgls_t_int[i] <-   NA
    taxon_data$SR_TE_sgls_t_slope[i] <- NA
    taxon_data$SR_TE_sgls_P_int[i] <-   NA
    taxon_data$SR_TE_sgls_P_slope[i] <- NA
    taxon_data$SR_TE_sgls_sigma[i] <-   NA
  }

  write.csv(taxon_data, file="spatial_auto_data.csv")
}


meta_sr_la <- rma.uni(yi=taxon_data$SR_LA_sgls_slope, sei=taxon_data$SR_LA_sgls_se_slope)
meta_sr_te <- rma.uni(yi=taxon_data$SR_TE_sgls_slope, sei=taxon_data$SR_TE_sgls_se_slope)

summary(meta_sr_la)
summary(meta_sr_te)

tdf <- taxon_data[which(taxon_data$n_species > 20),]


length(which(tdf$DR_TE_pgls_slope > 0 & tdf$DR_TE_pgls_P_slope < 0.05)) # 4
length(which(tdf$DR_TE_pgls_slope < 0 & tdf$DR_TE_pgls_P_slope < 0.05)) # 13

length(which(tdf$DR_BS_pgls_slope > 0 & tdf$DR_BS_pgls_P_slope < 0.05)) # 8
length(which(tdf$DR_BS_pgls_slope < 0 & tdf$DR_BS_pgls_P_slope < 0.05)) # 9

length(which(tdf$SR_LA_sgls_slope > 0 & tdf$SR_LA_sgls_P_slope < 0.05)) # 3
length(which(tdf$SR_LA_sgls_slope < 0 & tdf$SR_LA_sgls_P_slope < 0.05)) # 34

length(which(tdf$SR_TE_sgls_slope > 0 & tdf$SR_TE_sgls_P_slope < 0.05)) # 29
length(which(tdf$SR_TE_sgls_slope < 0 & tdf$SR_TE_sgls_P_slope < 0.05)) # 4
