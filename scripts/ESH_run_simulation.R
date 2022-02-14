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
# Description: Function to run gen3sis inputs 
#
# Contact: alexander.skeels@gmail.com
#
######################################


library(gen3sis)
setwd(file.path(HOME, "scripts", "gen3sis", "config"))

# config files
setwd("m0")
config_M0 <- list.files()
setwd("../m1")
config_M1 <- list.files()
setwd("../m2")
config_M2 <- list.files()
setwd("../m3")
config_M3 <- list.files()

# input file
input_dir <- file.path(HOME, "data", "Scotese_Behrmann_world_220x220km")

# output
setwd(file.path(HOME, "scripts", "gen3sis"))
if(!dir.exists("output")){dir.create("output"); setwd("output"); dir.create("m0");dir.create("m1");dir.create("m2");dir.create("m3")}
output_dir <- file.path(HOME, "scripts", "gen3sis", "output")

# run single simulation # NOTE: to run large batchs of these parameters will need HPC
setwd(file.path(HOME))

sim <- run_simulation(config = config_M0[1],
                      landscape = input_dir,
                      output_directory = output_dir,
                      save_state = NA,
                      call_observer = "all",
                      verbose = 1)

