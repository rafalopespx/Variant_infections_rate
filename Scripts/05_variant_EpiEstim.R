## Cleaning the ambient
rm(list = ls())
gc()

## Loading Libraries
packs = c("tidyverse", "vroom", "jsonlite", "httr", "geofacet")
lapply(packs,require, character.only = TRUE)

## Loading functions
source("Scripts/Functions/functions.R")
# source("Scripts/estimate_rt_ro_fun.R")

infections_variants_weekly<-vroom("Data/infections_estimates_variants_weekly.csv.xz")

infections_variants_daily<-vroom("Data/infections_estimates_variants_daily.csv.xz")

results<-c()

estimates_df<-infections_variants_daily

# For each state-variant combination, calculate the Rt
for (i in unique(estimates_df$name_states)) {
  for (j in unique(estimates_df$variant)) {
    
    # Subset the data for the current state-variant combination
    subset_data <- estimates_df[estimates_df$name_states == i & estimates_df$variant == j, ]
    
    # Estimate the Rt using the EpiEstim package
    epi_result <- tryCatch(estimate_R(subset_data, 
                             method = "parametric_si", 
                             config = make_config(list(mean_si = 3.5, 
                                                       std_si = 1))), 
                           error = function(err) NA)
    
    epi_df<-epi_result$R |> 
      as.data.frame() |> 
      mutate(name_states = i, 
             variant = j)
    
    # Store the results in a new dataframe
    results[[i]][[j]] <- epi_df
  }
}


