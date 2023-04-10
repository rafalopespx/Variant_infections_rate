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

#run for everything -Other
rt_fun <- function(df){
  suppressPackageStartupMessages(require(EpiEstim))
  
  #1st day with infections of variant to start the R estimate otherwise R estimate artificially high
  non0 <- min(which(df$I > 0)) 
  #last day with infections of variant to start the R estimate otherwise R estimate artificially high
  end0 <- max(which(df$I > 0))
  #dataframe filtered where there is the first case of variant
  df2 <- df[non0:end0,] 
  
  #input of interval for R estimate
  #estimate burn in period for R estimate: t_end - t_start
  t_start<-seq(2, nrow(df2)-15) 
  t_end<- t_start + 15
  
  # configuration of input data for R estimate
  # essentially we are estimating the serial intervals of SARS-CoV-2 (Omicron variant specific)
  # by drawing from two ( truncated normal ) distributions 
  # for the mean and standard deviation of the serial interval
  config <- make_config(list(mean_si = 3.5, 
                             std_mean_si = 1, 
                             min_mean_si = 1, 
                             max_mean_si = 6, # estimates for SARS-CoV-2 serial interval
                             std_si = 1, 
                             std_std_si = 0.5, 
                             min_std_si = 0.5, 
                             max_std_si = 1.5,
                             n1= 80, 
                             n2=20, 
                             t_start=t_start, 
                             t_end=t_end)
  )
  
  # main R estimate function:
  # will search for column named I which was created in the ci_fun but explicitly named here
  mean_Rt <- estimate_R(df2$I, 
                        method="uncertain_si",
                        config = config)
  
  #adds back in days that were filtered out to match the days in the main dataframe
  mean_Rt$R$t_start <- mean_Rt$R$t_start + non0 
  mean_Rt$R$t_end <- mean_Rt$R$t_end + non0
  
  # #binds them into a dataframe
  rt_df<-cbind.data.frame(day = mean_Rt$R$t_start,
                          Rt = mean_Rt$R$`Mean(R)`,
                          lower = mean_Rt$R$`Quantile.0.05(R)`,
                          upper = mean_Rt$R$`Quantile.0.95(R)`)
  
  #merges the Rt value with the other variant data and renames Rt to have variant suffix
  merge <- df %>%
    # mutate(date = as.Date(date)) |> 
    arrange(days) %>% #keep in week so that the day variable lines up with the first week
    mutate(day = 1:nrow(df)) %>% #used to merge with the estimate_R variable output for the day
    left_join(rt_df, by = c("day")) |> 
    dplyr::mutate(date_start = seq.Date(from = min(days) + non0, 
                                        length.out = nrow(df),
                                        by = "days"),
                  date_end = seq.Date(from = min(days) + non0 + 15, 
                                      length.out = nrow(df),
                                      by = "days")) |> 
    select(days, date_start, date_end, name_states, variant, I, Rt, lower, upper)
  #renames the smooth_spline output to have variant prefix
  #rename_with(.fn = ~paste0(name,"_",.), .cols = c("Rt", "rtlowci", "rtupci")) 
  
  return(merge)
  
}

## rt_fun tryCatch
rt_safe<-function(x){
  result <- tryCatch(rt_fun(df = x), error = function(err) NA)
  return(result)
}

## Creating estimates for Rt
## Try do it in parallel
rt_list<-c()

## Estimates to use
estimates_df<-infections_variants_daily

## Vectors to for loops
variants<-unique(estimates_df$variant)

states<-unique(estimates_df$name_states)

for (i in variants) {
  for (j in states) {
    
    tmp <- estimates_df |> 
      filter(variant == i, 
             name_states == j)
    
    ## Try to handle when any Rt fails and continue it 
    rt_list[[i]][[j]]<-rt_safe(x = tmp)
    rm(tmp)
    gc()
    cat("Finished state: ", j, "\n")
  }
  # rt_list[[i]]<-rt_list[[i]] |> 
  #   reduce(rbind)
  
  ## Prompting messages, to monitor progress
  cat("Finished variant: ", i, "over all states \n")
}

