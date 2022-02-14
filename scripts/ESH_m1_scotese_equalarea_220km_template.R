######################################
###            METADATA            ###
######################################
# Version: 1.0
#
# Author: Alexander Skeels
#
# Date: 11.06.2021
#
# Landscape: Scotese_Behrmann_world_220x220km
#
# Projects: BIGEST: Evolutionary Speed Hypothesis (ESH)
#
# Description: M1 - temperature dependant population divergence model of the ESH
#
# Contact: alexander.skeels@gmail.com
#
######################################


########################
### General settings ###
########################


random_seed = params$seed
start_time = params$start_time
end_time = 0
max_number_of_species = params$max_species
max_number_of_coexisting_species = 1e6
initial_abundance =  params$initial_abundance

# a list of traits to include with each species
trait_names = c("temp","prec", "dispersal", "body_size")

# ranges to scale the input environemts with:
# not listed variable:         no scaling takes place
# listed, set to NA:           the environmental variable will be scaled from [min, max] to [0, 1]
# lsited with a given range r: the environmental variable will be scaled from [r1, r2] to [0, 1]
# Scotese 220km range: -37.50240, 26.67235
# Straume 220km range: -29.56649  33.37760

environmental_ranges = list("temp" = c(-38,  34), "prec"=NA)


#########################
### Observer Function ###
#########################


end_of_timestep_observer = function(data, vars, config){
  #save stuff
  plot_richness(data$all_species, data$landscape)
  save_richness()
  save_phylogeny()
  save_species()
  save_traits()
  # make p/a matrices
  if(!file.exists(file.path(config$directories$output, "pa_matrices"))){dir.create(file.path(config$directories$output, "pa_matrices"))}
  pa_matrix <- data.frame(matrix(0,nrow=length(data$landscape$coordinates[,1]), ncol=(length(data$all_species)+2)))
  pa_matrix[,1:2]<-data$landscape$coordinates
  rownames(pa_matrix) <- rownames(data$landscape$coordinates)
  names(pa_matrix)[1:2] <- c("x", "y")
  names(pa_matrix)[3:length(pa_matrix)] <-unlist(lapply(data$all_species, FUN=function(x){x$id}))
  for(i in 3:(length(pa_matrix[1,]))){
    pa_matrix[names(data$all_species[[i-2]]$abundance),i] <- 1
  }
  write.table(pa_matrix, file.path(config$directories$output,"pa_matrices",  paste0("species_pa_",vars$ti, ".txt")), row.names = F, col.names = T)
}


######################
### Initialization ###
######################


create_ancestor_species <- function(landscape, config) {
  #browser()
  # set up a raster of global temperatures 
  t_world <- raster::rasterFromXYZ(cbind(landscape$coordinates,landscape$environment[, 1, drop = F]))
  t_world <- extend(t_world, landscape$extent)
  
  # Will the simulation sample the ancestor from a partcular continent or randomly from anywhere?
  start_cells <- as.character(Which(t_world, cells=T))
  
  # now fill the list of new species
  all_species <- list()
  
  new_species <- create_species(as.character(start_cells), config)
  new_species$traits[ , "dispersal"] <- 1
  new_species$traits[ , "temp"] <- landscape$environment[start_cells, "temp"]
  new_species$traits[ , "prec"] <- landscape$environment[start_cells, "prec"]
  new_species$traits[ , "body_size"] <- 0.2
  
  all_species <- append(all_species, list(new_species))
  
  
  return(all_species)
  
}


#################
### Dispersal ###
#################


grid_cell_distance <- params$grid_cell_distance  # resolution of raster (distance between grid cells)
dispersal_multiplier <- params$dispersal   # 


get_dispersal_values <- function(n, species, landscape, config) {
  #values <- rep(c( params$grid_cell_distance * params$dispersal ),n)
  values <- rweibull(n, shape=2, scale=(grid_cell_distance *  dispersal_multiplier) )
  return(values)
}


##################
### Divergence ###
##################


# threshold for genetic distance after which a speciation event takes place
divergence_threshold =  params$divergence_threshold
lambda =  params$lambda 

get_divergence_factor <- function(species, cluster_indices, landscape, config) {
  divergences <- vector("numeric", length(unique(cluster_indices)))
  
  divergence_matrix <- matrix(0, ncol=length(unique(cluster_indices)),  nrow=length(unique(cluster_indices)))
  cluster_indices_order <- unique(cluster_indices)[order(unique(cluster_indices))]
  # order them, find pairwise temp dist
  for(i_x in 1:length(cluster_indices_order)){
    mean_temp_i <- mean(landscape$environment[names(species$abundance[which(cluster_indices %in% cluster_indices_order[i_x])]),"temp"], na.rm=T)
    for(j_x in 1:length(cluster_indices_order)){
      mean_temp_j <- mean(landscape$environment[names(species$abundance[which(cluster_indices %in% cluster_indices_order[j_x])]),"temp"], na.rm=T)
      divergence_matrix[i_x, j_x] <- (sum(mean_temp_i, mean_temp_j)/(2))^lambda
    }
  }
  rownames(divergence_matrix) <- cluster_indices_order
  colnames(divergence_matrix) <- cluster_indices_order
  
  if(length(divergence_matrix) == 1){ divergence_matrix <- as.numeric(divergence_matrix)}else{
    
    diag(divergence_matrix) <- 0
  }
  
  return(divergence_matrix)
}



#######################
### Trait Evolution ###
#######################


sigma_bs <-  params$sigma_bs
sigma_t <-  params$sigma_t
psi <-  params$psi

apply_evolution <- function(species, cluster_indices, landscape, config) {
  
  # cell names
  traits <- species[["traits"]]
  cells <- rownames( traits )
  
  # evolve traits for each cluster
  for(cluster_index in unique(cluster_indices)){
    cells_cluster <- cells[which(cluster_indices == cluster_index)]
    t_theta_cluster <- mean(landscape$environment[cells_cluster,"temp"], na.rm=T)
    # evolve temperature
    traits[cells_cluster, "temp"] <- traits[cells_cluster, "temp"] + ( psi * (t_theta_cluster - traits[cells_cluster, "temp"]) ) + rnorm(1, mean = 0, sd = sigma_t)
    # evolve body size
    traits[cells_cluster, "body_size"] <- traits[cells_cluster, "body_size"] + rnorm(1, mean = 0, sd = sigma_bs)
  }
  
  # set bounds between 0 and 1 so the species can;t evolve a niche beyond that present in the data (all temp is scaled between 0 and 1)
  if(any(traits[, "temp"] > 1)){traits[which(traits[,"temp"]>1), "temp"] <- 1}
  if(any(traits[, "temp"] < 0)){traits[which(traits[,"temp"]<0), "temp"] <- 0}
  if(any(traits[, "body_size"] > 1)){traits[which(traits[,"body_size"]>1), "body_size"] <- 1}
  if(any(traits[, "body_size"] < 0)){traits[which(traits[,"body_size"]< 0), "body_size"] <- 0}
  
  return(traits)
  
}

#################################################
### Environmental and Ecological Interactions ###
#################################################


K_opt_max <-  params$K    # environmental filtering parameter (maximum abundance of all individuals in a mesic cell)
omega <- params$omega     # environmental filtering parameter (higher values equal less restrictive niche)
aridity_cost <-  params$aridity_cost # environmental filtering parameter (higher values equal less abundance in arid cells)

x0 <- params$inflexion # extinction parameter: inflexion point of extinction probability curve
decay <- params$decay  # extinction parameter: inflexion point of extinction probability curve

apply_ecology <- function(abundance, traits, landscape, config) {
  
  # Realised carrying capacity (reduced in arid cells)
  K <- K_opt_max * exp(1-aridity_cost*landscape[, "prec"])
  
  # difference between landscape and species optimum temperature niche
  diff_t <- abs(traits[, "temp"] - landscape[, "temp"])

  # potential population size based on resource use efficieny 
  # omega is strength of environmental filtering
  # will equal K when tthe species is perfectly adapted
  # modified from McPeek 2008/2007
  Nij <- K * exp(-(diff_t/omega)^2)
  
  # These lines make carrying capcaity a zero-sum game - species abundance is scaled based on their resourse use efficiency (Nij)
  Nij_sum <- sum(Nij)
  Nij_hat <- Nij * (sort(c(Nij_sum, K),partial=1)[1] / Nij_sum)
  Nij_hat[which(is.na(Nij_hat))] <- 0
  # now do an extinction filter based on population size
  prob_extinction <- (1/(1+exp(-decay*(x0 -  Nij_hat))))
  Nij_hat[which(sapply(prob_extinction, FUN=function(x){sample(c(0,1),1, prob= c(x, 1-(x)))})==0)] <- 0 
  
  # absolute precipitation threshold 
  #if(landscape[, "prec"] > 0.5){Nij <- rep(0, length(Nij))}
  
  return(Nij_hat)
}



