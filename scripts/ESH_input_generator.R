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
# Description: Function to generate gen3sis inputs 
#
# Contact: alexander.skeels@gmail.com
#
######################################

library(gen3sis)
library(raster)


# read in the data
# tmperature and aridity reconstructions at 1Degree resolution
landscapes <- readRDS("data/scotese_1D_landscape.rds")

# get number of time-steps
timesteps <- colnames(landscapes$temp)[3:length(colnames(landscapes$temp))]

# create rasters of temperature and aridity at each time step by converting 1Degree to 220km x220km resolution
# do ardity index
prec_list<-vector("list", (length(landscapes$prec[1,])-2))
for(i in 3:length(landscapes$prec[1,])){
  r_tmp <- rasterFromXYZ(landscapes$prec[,c(1,2,i)], crs ="+proj=longlat +datum=WGS84") 
  r_tmp_proj <- projectRaster(r_tmp,res=220000, crs = crs('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +no_defs'))
  prec_list[[i-2]] <- r_tmp_proj
}

# do temprature
temp_list <- vector("list", (length(landscapes$temp[1,])-2))
for(i in 3:length(landscapes$temp[1,])){
  r_tmp <- rasterFromXYZ(landscapes$temp[,c(1,2,i)], crs ="+proj=longlat +datum=WGS84") 
  r_tmp_proj <- projectRaster(r_tmp,res=220000, crs = crs('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +no_defs'))
  temp_list[[i-2]] <- r_tmp_proj
}

landscape_list <- list(temp=temp_list, prec=prec_list)

cost_function <- function(source, habitable_src, dest, habitable_dest) {
  if(!all(habitable_src, habitable_dest)) {
    return(2/1000)
  } else {
    return(1/1000)
  }
}

overwrite_output <- T
directions <- 8
output_directory <- "Scotese_Behrmann_world_220x220km"

create_input_landscape(landscapes = landscape_list,
                       timesteps = timesteps,
                       cost_function = cost_function,
                       directions = directions,
                       output_directory = output_directory,
                       overwrite_output = overwrite_output,
                       crs = '+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +no_defs',
                       calculate_full_distance_matrices = T)

