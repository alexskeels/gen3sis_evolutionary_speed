# Gen3sis: Evolutionary Speed Hypothesis

Skeels, A., Bach, W., Hagen, O., Jetz, W., Pellissier, L. in prep. 

R code and data from a study of temperature and body size dependant rates of population divergence in tetrapods using gen3sis.

All scripts written by Alexander Skeels: alexander.skeels@gmail.com

Feel free to ask me any questions!

#~~~~~~ data

data used and generated in this study

      m0_summary_statistics.txt         # summary statistics for m0
  
      m1_summary_statistics.txt         # summary statistics for m1
  
      m2_summary_statistics.txt         # summary statistics for m2
  
      m3_summary_statistics.txt         # summary statistics for m3
      
      empirical_summary_statistics.csv  # summary statistics for empirical tetrapod data
  
      scotese_1D_landscapes.rds         # paleo-temperature and ariditylandscape data 
      

  
#~~~~~~ scripts/gen3sis

scripts to set up and run simulation experiment using gen3sis
  
    ESH_input_generator.R                       # generates landscape in gen3sis format from paleoenvironmental reconstructions
  
    ESH_config_generator.R                      # generates 500 configs under each model from sobol seuqences - stored in scripts/gen3sis/config
    
    ESH_m0_scotese_equalarea_220km_template.R   # genesis config templates for m0 submodels

    ESH_m1_scotese_equalarea_220km_template.R   # genesis config templates for m1 submodels
  
    ESH_m2_scotese_equalarea_220km_template.R   # genesis config templates for m2 submodels
  
    ESH_m3_scotese_equalarea_220km_template.R   # genesis config templates for m3 submodels
  
    ESH_run_simulation.R                        # script to run a single gen3sis simulation
    
 scripts for data analysis
    
    ESH_empirical_summary_statistic_metaanalysis.R   # script to run metanalysis on empirical summary statistics
    
    ESH_simulation_validation.R                      # script to validate siumulated patterns against empirical patterns
    
    ESH_simulation_sensitity_analysis.R              # script to explore simulated patterns and model parameters
    
    ESH_simulation_based_inference.R                 # script to perform model selection using simulation based inference methods
    
