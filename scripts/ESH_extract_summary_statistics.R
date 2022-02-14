######################################
###            METADATA            ###
######################################
# Version: 1.0
#
# Author: Alexander Skeels
#
# Date: 11.06.2021
#
# Projects: BIGEST: Evolutionary Speed Hypothesis (ESH)
#
# Description: Scripts to extract summary statistics from simulation experiment
#
# Contact: alexander.skeels@gmail.com
#
######################################

# load libraries
library(ape)
library(raster)
library(rgeos)
library(gen3sis)
library(picante)
library(apTreeshape)
library(moments)
library(adephylo)
library(phytools)
library(dplyr)
library(PhyloMeasures)

# rescale tree function
rescaleTree<-function(tree,scale){
  tree$edge.length<-
    tree$edge.length/max(nodeHeights(tree)[,2])*scale
  return(tree)
}

# summarise traits in sites function
summarisePresenceAbsence <- function(p, y, z){
  
  y[y == 0] <- 1e-38
  
  for (i in 1:ncol(p)) {
    pos <- z[i]
    p[, i] <- p[, i] * y[i,]
    pos2 <- p[, i] == 0
    p[pos2, i] <- NA
  }
  func2 <- function(x) {
    pos <- is.na(x)
    resu <- mean(x[!pos])
  }
  func3 <- function(x) {
    pos <- is.na(x)
    resu <- sd(x[!pos])
  }
  mean_bodysize <- apply(p, 1, func2)
  sd_bodysize <- apply(p, 1, func3)
  
  return(data.frame(site=rownames(p),mean_bs= mean_bodysize, sd_bs=sd_bodysize))
}

# set directory to where simulation output is stored
# each submodel (M0-3) should be stored in a seperate folder names either m1, m2, m3, or m0
setwd("PATH/To/OUTPUT/FOLDER")

# start by extracting from m0, change this line to extract from other models
model="m0"

# read in parameter data generated in the script "ESH_config_generator.R" and used to run the simulations
params <- read.table(paste0(model,"_parameters.txt"), header=T)

# add extra columns for the summary statistics
# time step
params$running_time_step <- NA

# diversity summary
params$n_total_diversity   <- NA
params$n_extant_diversity  <- NA
params$n_extinct_diversity <- NA
params$turnover <- NA

# grid cell richnes summaries
params$richness_mean <- NA
params$richness_max <-NA
params$richness_median <- NA
params$richness_sd <- NA
params$richness_skewness <- NA
params$richness_kurtosis <- NA

# grid cell richness ~ latitude/temperature/paleo-temperature correlations
params$richness_lat_cor <-    NA    
params$richness_temp0mya_cor <- NA
params$richness_temp10mya_cor <-NA
params$richness_temp20mya_cor <-NA
params$richness_temp30mya_cor <-NA
params$richness_temp40mya_cor <-NA
params$richness_temp50mya_cor <-NA
params$richness_temp60mya_cor <-NA

# trait summaries
params$temp_mean <- NA
params$temp_sd <- NA
params$temp_skewness <- NA
params$temp_kurtosis <- NA

params$temp_trait_mean <- NA
params$temp_trait_sd <- NA
params$temp_trait_skewness <- NA
params$temp_trait_kurtosis <- NA

params$bs_mean <- NA
params$bs_sd <- NA
params$bs_skewness <- NA
params$bs_kurtosis <- NA

#range summaries
params$rs_mean <- NA
params$rs_sd <- NA
params$rs_skewness <- NA
params$rs_kurtosis <- NA

# trait ~ trait correlations
params$bodysize_temp_cor<- NA
params$bodysize_rangesize_cor<- NA
params$rangesize_temp_cor<- NA

# trait ~ tip metrics correlations
params$bodysize_MRDa_cor<- NA
params$bodysize_MRDn_cor<- NA
params$bodysize_MRDs_cor<- NA
params$bodysize_ED_cor<- NA
params$bodysize_ES_cor<- NA
params$bodysize_DR_cor<- NA
params$temp_MRDa_cor<- NA
params$temp_MRDn_cor<- NA
params$temp_MRDs_cor<- NA
params$temp_ED_cor <- NA
params$temp_ES_cor<- NA
params$temp_DR_cor<- NA
params$rangesize_MRDa_cor <-   NA
params$rangesize_MRDn_cor <-  NA 
params$rangesize_MRDs_cor <-NA   
params$rangesize_ED_cor <-NA      
params$rangesize_ES_cor <-  NA    
params$rangesize_DR_cor <- NA

params$bodysize_MRDa_cor_p<- NA
params$bodysize_MRDn_cor_p<- NA
params$bodysize_MRDs_cor_p<- NA
params$bodysize_ED_cor_p<- NA
params$bodysize_ES_cor_p<- NA
params$bodysize_DR_cor_p<- NA
params$temp_MRDa_cor_p<- NA
params$temp_MRDn_cor_p<- NA
params$temp_MRDs_cor_p<- NA
params$temp_ED_cor_p <- NA
params$temp_ES_cor_p<- NA
params$temp_DR_cor_p<- NA
params$rangesize_MRDa_cor_p <-   NA
params$rangesize_MRDn_cor_p <-  NA 
params$rangesize_MRDs_cor_p <-NA   
params$rangesize_ED_cor_p <-NA      
params$rangesize_ES_cor_p <-  NA    
params$rangesize_DR_cor_p <- NA

# tree shape metrics
params$collessI<- NA
params$sackinI<- NA
params$gamma<- NA

# spatial phylogenetic metrics
params$lat_pd_cor   <-  NA
params$lat_mpd_cor  <-  NA
params$lat_mntd_cor <-  NA
params$temp_pd_cor  <- NA
params$temp_mpd_cor  <- NA
params$temp_mntd_cor <- NA
params$richness_pd_cor  <- NA
params$richness_mpd_cor <- NA
params$richness_mntd_cor <- NA

params$lat_pd_cor_p   <-  NA
params$lat_mpd_cor_p  <-  NA
params$lat_mntd_cor_p <-  NA
params$temp_pd_cor_p  <- NA
params$temp_mpd_cor_p  <- NA
params$temp_mntd_cor_p <- NA
params$richness_pd_cor_p  <- NA
params$richness_mpd_cor_p <- NA
params$richness_mntd_cor_p <- NA

# spatial functional metrics
params$lat_bsm_cor  <-      NA
params$lat_bssd_cor <-      NA
params$temp_bsm_cor <-      NA
params$temp_bssd_cor <-     NA
params$richness_bsm_cor <-  NA
params$richness_bssd_cor <- NA
params$pd_bsm_cor <- NA
params$mpd_bsm_cor <- NA
params$mntd_bsm_cor <- NA
params$pd_bssd_cor <- NA
params$mpd_bssd_cor <- NA
params$mntd_bssd_cor <- NA

params$lat_bsm_cor_p <-     NA
params$lat_bssd_cor_p <-    NA
params$temp_bsm_cor_p <-    NA
params$temp_bssd_cor_p <-   NA
params$richness_bsm_cor_p <-  NA
params$richness_bssd_cor_p <- NA
params$pd_bsm_cor_p <- NA
params$mpd_bsm_cor_p <- NA
params$mntd_bsm_cor_p <- NA
params$pd_bssd_cor_p <- NA
params$mpd_bssd_cor_p <- NA
params$mntd_bssd_cor_p <- NA


# load the landscape object
landscape <- readRDS("PATH/TO/landscapes.rds")
landscape_t <- landscape$temp
landscape_t$site <- rownames(landscape_t)
  
# go into the model output folder (either m0, m1, m2 , m3)
setwd(paste0("output/",model))

# list the files
m_files <- list.files(pattern=model)
m_files <- m_files[order(as.numeric(sapply(m_files, FUN=function(m)strsplit(m, "_")[[1]][6])))]


#mega loop that collects all the suammry stats and saves the output!
for(i in c(1:length(m_files))){
  
  setwd(file.path(m_files[i]))

  setwd(file.path("richness"))
  print(m_files[i])
  files <- list.files()
  
  most_recent <- min(as.numeric(sapply(files, FUN=function(x){strsplit(strsplit(x, "_")[[1]][3], ".rds")[[1]]})))
  if(most_recent != 0){ most_recent<- most_recent+1}
  params$running_time_step[i] <- most_recent
 
  # RICHNESS SUMMARY STATS
  richness <- readRDS(paste0("richness_t_", most_recent,".rds"))
  params$richness_mean[i] <- mean(richness, na.rm=T)
  params$richness_max[i] <- max(richness, na.rm=T)
  params$richness_median[i] <- median(richness, na.rm=T)
  params$richness_sd[i] <- sd(richness, na.rm=T)
  params$richness_skewness[i] <- skewness(richness, na.rm=T)
  params$richness_kurtosis[i] <- kurtosis(richness, na.rm=T)
  
  richness <- data.frame(richness)
  richness$site <- rownames(richness)
  
  richness_l <- left_join(richness, landscape_t, by='site', match='all')
  
  # if total extinction move on to next simulation
  if(sum(richness$richness)==0){
    params$n_extant_diversity[i] <- 0
    setwd("../../")
    next
  }
  
  params$richness_lat_cor[i] <-         cor(richness_l$richness, abs(richness_l$y), use='complete.obs', method = "spearman")
  params$richness_temp0mya_cor[i] <-  cor(richness_l$richness, richness_l$'0', use='complete.obs',   method = "spearman")
  params$richness_temp10mya_cor[i] <- cor(richness_l$richness, richness_l$'60', use='complete.obs',  method = "spearman")
  params$richness_temp20mya_cor[i] <- cor(richness_l$richness, richness_l$'120', use='complete.obs', method = "spearman")
  params$richness_temp30mya_cor[i] <- cor(richness_l$richness, richness_l$'180', use='complete.obs', method = "spearman")
  params$richness_temp40mya_cor[i] <- cor(richness_l$richness, richness_l$'240', use='complete.obs', method = "spearman")
  params$richness_temp50mya_cor[i] <- cor(richness_l$richness, richness_l$'300', use='complete.obs', method = "spearman")
  params$richness_temp60mya_cor[i] <- cor(richness_l$richness, richness_l$'360', use='complete.obs', method = "spearman")
  
  #SPECIES-LEVEL SUMMARY STATS
  setwd("../species")
  species <- readRDS(paste0("species_t_", most_recent,".rds"))
  
  params$n_extant_diversity[i] <- sum(sapply(species, FUN=function(x){if(length(x$abundance)>0){return(1)}else{0}}))
  
  # don't really want to waste time on calulating further metrics on runs with too many species (will take forever for the phylogenetic ones)
  if(most_recent > 0){
    setwd("../../")
    next
    }
  
  setwd("../")
  sim <- readRDS("sgen3sis.rds")
  params$occupancy[i] <- tail(sim$summary$occupancy, 1)
  
  #TRAIT SUMMARY STATS
  setwd("traits")
  traits <- readRDS(paste0("traits_t_", most_recent,".rds"))
  
  # realised temperature niche
  
  species_realised_temp_niche <-  do.call(rbind, lapply(species, FUN=function(x){
    realised_temp <- landscape_t[which(landscape_t$site %in% names(x$abundance)),"0"]
    return(data.frame(Species=paste0("species", x$id), niche_trait=mean(realised_temp, na.rm=T)))
  }))
  
  params$temp_mean[i]   <-  mean(species_realised_temp_niche$niche_trait, na.rm=T) # mean temp niche
  params$temp_sd[i]     <- sd(species_realised_temp_niche$niche_trait, na.rm=T)# sd temp niche
  params$temp_skewness[i] <- skewness(species_realised_temp_niche$niche_trait, na.rm=T)
  params$temp_kurtosis[i] <- kurtosis(species_realised_temp_niche$niche_trait, na.rm=T)  
  # trait values
  temp_trait_means <- unlist(lapply(traits, FUN=function(x){y=mean(x[,"temp"], na.rm=T);return(y)}))
  bs_trait_means <- unlist(lapply(traits, FUN=function(x){y=mean(x[,"body_size"], na.rm=T);return(y)}))
  
  params$temp_trait_mean[i]<- mean(temp_trait_means, na.rm=T) # mean temp niche
  params$temp_trait_sd[i] <- sd(temp_trait_means, na.rm=T) # sd temp niche
  params$temp_trait_skewness[i] <- skewness(temp_trait_means, na.rm=T)
  params$temp_trait_kurtosis[i] <- kurtosis(temp_trait_means, na.rm=T)
  
  params$bs_mean[i]<- mean(bs_trait_means, na.rm=T) # mean temp niche
  params$bs_sd[i] <- sd(bs_trait_means, na.rm=T) # sd temp niche
  params$bs_skewness[i] <- skewness(bs_trait_means, na.rm=T)
  params$bs_kurtosis[i] <- kurtosis(bs_trait_means, na.rm=T)
  
  # RANGE SIZE FREQUENCY DISTRIBUTIONS
  range_size_l <- lapply(traits, nrow)
  range_size <- data.frame(Species = paste0("species",names(range_size_l)), range_size =unlist(range_size_l)) 
  
  params$rs_mean[i] <- mean(unlist(range_size_l), na.rm=T)
  params$rs_sd[i] <- sd(unlist(range_size_l), na.rm=T)
  params$rs_skewness[i] <- skewness(unlist(range_size_l))
  params$rs_kurtosis[i] <- kurtosis(unlist(range_size_l))
  
  #PHYLO SUMMARY STATS
  setwd("../phylogeny")
  phy <- try(read.nexus(paste0("phylogeny_t_", most_recent,".nex")))
  if(class(phy)=="try-error"){setwd("../../") ; next}  
  
  # phy's will differ based on age of clade - so try to make comparable 
  phy_og <- phy
  
  phy <- try(drop.fossil(phy))
  # tree manually dropping extincts
  if(class(phy) == "try-error"){
    extinct_species <- paste0("species",which(sapply(species, FUN=function(x){length(x$abundance)==0})))
    phy <- phy_og
    phy <- try(drop.tip(phy, extinct_species))
  }
  # if that fails too move on
  if(class(phy) == "try-error"){setwd("../../") ; next}
  phy <- rescaleTree(phy, 1)
  
  # Tips and Diversity
  params$n_total_diversity[i] <- length(phy_og$tip.label)
  params$n_extinct_diversity[i] <- length(phy_og$tip.label) - length(phy$tip.label)
  
  # can't do uch if no species so move on to next sim
  if(length(phy$tip) < 3){setwd("../../") ; next} else { 

    if(any(phy$edge.length < 0 )){ phy$edge.length <- phy$edge.length + abs(min(phy$edge.length))}
    
    # tree shape parameters
    phy_treeshape <- as.treeshape(phy)
    params$collessI[i] <- colless(phy_treeshape, "yule")
    #params$beta[i] <- maxlik.betasplit(phy_treeshape, up = 10, remove.outgroup = FALSE, confidence.interval = "none", conf.level = 0.95, size.bootstrap = 100)$max_lik
    params$sackinI[i] <- sackin(phy_treeshape, "yule")
    params$gamma[i] <- gammaStat(phy)
    
    
    # Evolutionary Distinctiveness
    ED <- data.frame(phyloregion::evol_distinct(phy, type = c("fair.proportion"),scale = FALSE, use.branch.lengths = TRUE))
    ED$Species <- rownames(ED)
    colnames(ED)[1] <- "ED"
    # Equal Splits
    ES <- data.frame(phyloregion::evol_distinct(phy, type = c("equal.splits"),scale = FALSE, use.branch.lengths = TRUE))
    ES$Species <- rownames(ES)
    colnames(ES)[1] <- "ES"
    
    # Div rates
    tree_metrics <- left_join(ES,ED, by ="Species", match = "all")
    tree_metrics$DivRate <- 1/tree_metrics$ES
    
    #MRD three ways
    MRD_n <- data.frame(MRD_n = distRoot(phy, tips = "all", method = c("nNodes")))
    MRD_n$Species <- row.names(MRD_n)
    MRD_A <- data.frame(MRD_a = distRoot(phy, tips = "all", method = c("Abouheif")))
    MRD_A$Species <- row.names(MRD_A)
    MRD_s <- data.frame(MRD_s = distRoot(phy, tips = "all", method = c("sumDD")))
    MRD_s$Species <- row.names(MRD_s)
    MRD <- left_join(MRD_n, MRD_A, by = "Species", match = "all")
    MRD <- left_join(MRD, MRD_s,  by = "Species", match = "all")
    names(MRD) <- c("MRD_nodes", "Species",  "MRD_Abouheif", "MRD_sumDD")
    
    # join em all together
    tree_metrics <- left_join(tree_metrics, MRD, by = "Species", match = "all")
    
    
    # TIP RATE CORRELATIONS
    
    # now link em with trait data
    niche_trait <- species_realised_temp_niche
    body_size_l <- lapply(traits, FUN=function(x){y=mean(x[,"body_size"], na.rm=T);return(y)})
    body_size <- data.frame(Species = paste0("species",names(body_size_l)), body_size =unlist(body_size_l)) 
    
    traits_means <- left_join(body_size, niche_trait, by = "Species", match = "all")
    
    trait_trees <- left_join(traits_means, tree_metrics, by = "Species", match = "all")
    
    trait_trees <- left_join(traits_means, tree_metrics, by = "Species", match = "all")
    trait_trees_ranges <- left_join(trait_trees, range_size, by = "Species", match = "all")
    
    
    correlations <- try(cor.table(na.omit(trait_trees_ranges[, 2:ncol(trait_trees_ranges)]), cor.method = "spearman"))
    
    if(class(correlations) == "try-error"){setwd("../../"); next}
    
    params$bodysize_MRDa_cor[i] <-   correlations$r["body_size","MRD_Abouheif"]
    params$bodysize_MRDn_cor[i] <-   correlations$r["body_size","MRD_nodes"]
    params$bodysize_MRDs_cor[i] <-   correlations$r["body_size","MRD_sumDD"]
    params$bodysize_ED_cor[i] <-      correlations$r["body_size","ED"]
    params$bodysize_ES_cor[i] <-      correlations$r["body_size","ES"]
    params$bodysize_DR_cor[i] <- correlations$r["body_size","DivRate"]
    
    params$temp_MRDa_cor[i] <-   correlations$r["niche_trait","MRD_Abouheif"]
    params$temp_MRDn_cor[i] <-   correlations$r["niche_trait","MRD_nodes"]
    params$temp_MRDs_cor[i] <-   correlations$r["niche_trait","MRD_sumDD"]
    params$temp_ED_cor[i] <-      correlations$r["niche_trait","ED"]
    params$temp_ES_cor[i] <-      correlations$r["niche_trait","ES"]
    params$temp_DR_cor[i] <- correlations$r["niche_trait","DivRate"]
    
    params$rangesize_MRDa_cor[i] <-   correlations$r["range_size","MRD_Abouheif"]
    params$rangesize_MRDn_cor[i] <-   correlations$r["range_size","MRD_nodes"]
    params$rangesize_MRDs_cor[i] <-   correlations$r["range_size","MRD_sumDD"]
    params$rangesize_ED_cor[i] <-      correlations$r["range_size","ED"]
    params$rangesize_ES_cor[i] <-      correlations$r["range_size","ES"]
    params$rangesize_DR_cor[i] <- correlations$r["range_size","DivRate"]
    
    params$bodysize_temp_cor[i]       <-   correlations$r["body_size","niche_trait"]
    params$bodysize_rangesize_cor[i]  <-   correlations$r["body_size","range_size"]
    params$rangesize_temp_cor[i]      <-   correlations$r["range_size","niche_trait"]
    
    ############################################# P
    params$bodysize_MRDa_cor_p[i] <-   correlations$P["body_size","MRD_Abouheif"]
    params$bodysize_MRDn_cor_p[i] <-   correlations$P["body_size","MRD_nodes"]
    params$bodysize_MRDs_cor_p[i] <-   correlations$P["body_size","MRD_sumDD"]
    params$bodysize_ED_cor_p[i]   <-   correlations$P["body_size","ED"]
    params$bodysize_ES_cor_p[i]   <-   correlations$P["body_size","ES"]
    params$bodysize_DR_cor_p[i]   <-   correlations$P["body_size","DivRate"]
    
    params$temp_MRDa_cor_p[i] <-  correlations$P["niche_trait","MRD_Abouheif"]
    params$temp_MRDn_cor_p[i] <-  correlations$P["niche_trait","MRD_nodes"]
    params$temp_MRDs_cor_p[i] <-  correlations$P["niche_trait","MRD_sumDD"]
    params$temp_ED_cor_p[i]   <-  correlations$P["niche_trait","ED"]
    params$temp_ES_cor_p[i]   <-  correlations$P["niche_trait","ES"]
    params$temp_DR_cor_p[i]   <-  correlations$P["niche_trait","DivRate"]
    
    params$rangesize_MRDa_cor_p[i] <-  correlations$P["range_size","MRD_Abouheif"]
    params$rangesize_MRDn_cor_p[i] <-  correlations$P["range_size","MRD_nodes"]
    params$rangesize_MRDs_cor_p[i] <-  correlations$P["range_size","MRD_sumDD"]
    params$rangesize_ED_cor_p[i] <-    correlations$P["range_size","ED"]
    params$rangesize_ES_cor_p[i] <-    correlations$P["range_size","ES"]
    params$rangesize_DR_cor_p[i] <-    correlations$P["range_size","DivRate"]
    
    params$bodysize_temp_cor_p[i]       <-   correlations$P["body_size","niche_trait"]
    params$bodysize_rangesize_cor_p[i]  <-   correlations$P["body_size","range_size"]
    params$rangesize_temp_cor_p[i]      <-   correlations$P["range_size","niche_trait"]
    
    # SPATIAL PHYLOGENETIC CALCULATIONS AND CORRELATIONS
    
    setwd("../pa_matrices")
    
    #fix p/a matrix to match tree
    species_pa <- read.table(paste0("species_pa_", most_recent,".txt"), header=T)
    species_pa_tmp <- species_pa[, 3:ncol(species_pa)]
    colnames( species_pa_tmp) <- gsub("X", "species", colnames(species_pa_tmp))
    species_pa_tmp <- species_pa_tmp[, which(colSums(species_pa_tmp)>0)]
    
  
    #estimate PD
    pd_estimate <- pd.query(phy, species_pa_tmp, standardize = T)
    # estimate MPD
    mpd_estimate <- mpd.query(phy, species_pa_tmp, standardize = T)
    # estimate MPNTD
    mntd_estimate <- mntd.query(phy, species_pa_tmp, standardize = T)
    
    spatial_phylogenetic_metrics <- data.frame(site=rownames(species_pa_tmp), x=species_pa$x, y=species_pa$y,
                                               pd=pd_estimate, mpd=mpd_estimate, 
                                               mntd=mntd_estimate)
    
    spatial_phylogenetic_metrics$x <- round(spatial_phylogenetic_metrics$x)
    richness_l$x <- round(richness_l$x)
    spatial_phylogenetic_metrics$y <- round(spatial_phylogenetic_metrics$y)
    richness_l$y <- round(richness_l$y)
    richness_l$site <- rownames(richness_l)
    
    spatial_phylogenetic_metrics <- left_join(spatial_phylogenetic_metrics, richness_l[, c("site", "y","richness","0")], by=c("site", "y"))
    spatial_phylogenetic_metrics$y <- abs(spatial_phylogenetic_metrics$y)
    pd_correlations <- try(cor.table(na.omit(spatial_phylogenetic_metrics[, 2:ncol(spatial_phylogenetic_metrics)]), cor.method = "spearman"))
    
    if(class( pd_correlations)=="try-error"){
      
    } else {
      
      params$lat_pd_cor[i]   <-   pd_correlations$r["y","pd"]
      params$lat_mpd_cor[i]  <-   pd_correlations$r["y","mpd"]
      params$lat_mntd_cor[i] <-   pd_correlations$r["y","mntd"]
      
      params$temp_pd_cor[i]   <-   pd_correlations$r["0","pd"]
      params$temp_mpd_cor[i]  <-   pd_correlations$r["0","mpd"]
      params$temp_mntd_cor[i] <-   pd_correlations$r["0","mntd"]
      
      params$richness_pd_cor[i]   <-   pd_correlations$r["richness","pd"]
      params$richness_mpd_cor[i]  <-   pd_correlations$r["richness","mpd"]
      params$richness_mntd_cor[i] <-   pd_correlations$r["richness","mntd"]
      
      
      params$lat_pd_cor_p[i]   <-   pd_correlations$P["y","pd"]
      params$lat_mpd_cor_p[i]  <-   pd_correlations$P["y","mpd"]
      params$lat_mntd_cor_p[i] <-   pd_correlations$P["y","mntd"]
      
      params$temp_pd_cor_p[i]   <-   pd_correlations$P["0","pd"]
      params$temp_mpd_cor_p[i]  <-   pd_correlations$P["0","mpd"]
      params$temp_mntd_cor_p[i] <-   pd_correlations$P["0","mntd"]
      
      params$richness_pd_cor_p[i]   <-   pd_correlations$P["richness","pd"]
      params$richness_mpd_cor_p[i]  <-   pd_correlations$P["richness","mpd"]
      params$richness_mntd_cor_p[i] <-   pd_correlations$P["richness","mntd"]
    }
    
    # SPATIAL FUNCTIONAL CALCULATIONS AND CORRELATIONS

    traits_FD <- na.omit(traits_means)
    traits_FD <- traits_FD[order(traits_FD$Species),]
    species_tmp <- traits_FD$Species
    traits_FD <- data.frame(body_size = traits_FD$body_size)
    rownames(traits_FD) <- species_tmp
    # also scale FD so comparable between clades
    traits_FD$body_size <- scale(traits_FD$body_size, center = T)
    
    species_pa_tmp <- species_pa_tmp[, order(colnames(species_pa_tmp))]
    species_pa_tmp <- species_pa_tmp[, which(colnames(species_pa_tmp) %in% rownames(traits_FD))]
    species_pa_tmp <-  species_pa_tmp[which(rowSums(species_pa_tmp)>0), which(colSums(species_pa_tmp)>0)]
    
    # modified from maplizer in LetsR
    body_size_community_msd <- summarisePresenceAbsence(species_pa_tmp, traits_FD, rownames(traits_FD))
    spatial_functional_phylogenetic_metrics <- left_join(body_size_community_msd, spatial_phylogenetic_metrics[, c("site", "y","richness","0", "pd", "mpd", "mntd")], by=c("site"))
    spatial_functional_phylogenetic_metrics$y <- abs(spatial_functional_phylogenetic_metrics$y)
    
    
    fd_correlations <- try(cor.table(na.omit(spatial_functional_phylogenetic_metrics[, 2:ncol(spatial_functional_phylogenetic_metrics)]), cor.method = "spearman"))
    
    if(class( fd_correlations)=="try-error"){
      
    } else {
      
      params$lat_bsm_cor[i]   <-   fd_correlations$r["y","mean_bs"]
      params$lat_bssd_cor[i]  <-   fd_correlations$r["y","sd_bs"]

      params$temp_bsm_cor[i]   <-   fd_correlations$r["0","mean_bs"]
      params$temp_bssd_cor[i]  <-   fd_correlations$r["0","sd_bs"]

      params$richness_bsm_cor[i]   <-   fd_correlations$r["richness","mean_bs"]
      params$richness_bssd_cor[i]  <-   fd_correlations$r["richness","sd_bs"]

      params$pd_bsm_cor[i]     <-   fd_correlations$r["pd","mean_bs"]
      params$mpd_bsm_cor[i]    <-   fd_correlations$r["mpd","mean_bs"]
      params$mntd_bsm_cor[i]   <-   fd_correlations$r["mntd","mean_bs"]
      
      params$pd_bssd_cor[i]    <-   fd_correlations$r["pd","sd_bs"]
      params$mpd_bssd_cor[i]   <-   fd_correlations$r["mpd","sd_bs"]
      params$mntd_bssd_cor[i]  <-   fd_correlations$r["mntd","sd_bs"]
      
      
      params$lat_bsm_cor_p[i]   <-   fd_correlations$P["y","mean_bs"]
      params$lat_bssd_cor_p[i]  <-   fd_correlations$P["y","sd_bs"]

      params$temp_bsm_cor_p[i]   <-   fd_correlations$P["0","mean_bs"]
      params$temp_bssd_cor_p[i]  <-   fd_correlations$P["0","sd_bs"]

      params$richness_bsm_cor_p[i]   <-   fd_correlations$P["richness","mean_bs"]
      params$richness_bssd_cor_p[i]  <-   fd_correlations$P["richness","sd_bs"]
      
      params$pd_bsm_cor_p[i]     <-   fd_correlations$P["pd","mean_bs"]
      params$mpd_bsm_cor_p[i]    <-   fd_correlations$P["mpd","mean_bs"]
      params$mntd_bsm_cor_p[i]   <-   fd_correlations$P["mntd","mean_bs"]
      
      params$pd_bssd_cor_p[i]    <-   fd_correlations$P["pd","sd_bs"]
      params$mpd_bssd_cor_p[i]   <-   fd_correlations$P["mpd","sd_bs"]
      params$mntd_bssd_cor_p[i]  <-   fd_correlations$P["mntd","sd_bs"]
    }
    
    params$turnover[i]  <-   params$n_extinct_diversity[i] / params$n_extant_diversity[i] 
    
    # save everything so don't have to calulate again
    setwd("../")
    results_list <- list("richness"=richness,
                         "traits"=traits,
                         "phylogeny"=phy,
                         "species" = species,
                         "trait_range_tree_metrics"=trait_trees_ranges,
                         "spatial_functional_phylogenetic_metrics" = spatial_functional_phylogenetic_metrics)
    saveRDS(results_list, file="phylogeny_traits_richness_spatial_and_tip_metrics.rds")
    setwd("../")
    write.table(params, file=paste0(model,"_summary_stats.txt"), row.names = F, col.names = T)
  }
}

setwd("../")

write.table(params, file=paste0(model,"_summary_stats.txt"), row.names = F, col.names = T)