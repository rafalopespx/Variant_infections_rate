## Cleaning the ambient
rm(list = ls())
gc()

## Loading Libraries
packs = c("tidyverse", "vroom", "jsonlite", "httr", "geofacet", "EpiEstim")
lapply(packs,require, character.only = TRUE)

## Loading functions
source("Scripts/Functions/functions.R")
# source("Scripts/estimate_rt_ro_fun.R")

infections_variants_weekly<-vroom("Data/infections_estimates_variants_weekly.csv.xz")

infections_variants_weekly<-infections_variants_weekly |> 
  group_by(name_states, days, variant) |> 
  summarise(I = sum(I, na.rm = T))

#run for everything -Other
rt_fun_week <- function(df, wallinga_teunis = FALSE, ...){
  suppressPackageStartupMessages(require(EpiEstim))
  
  dots<-list(...)
  
  ## Making the time series have all weeks in the period
  df2<-df |> 
    ungroup() |> 
    complete(days = full_seq(df$days, period = 7)) |> 
    fill(c("name_states", "variant"), .direction = "down") |> 
    mutate(across(is.numeric, ~round(zoo::na.approx(.x), 0)))
  
  # main R estimate function:
  # will search for column named I which was created in the ci_fun but explicitly named here
  if(wallinga_teunis){
    config<-list(t_start = t_start, 
                 t_end = t_end, 
                 mean_si = 3.5, 
                 std_si = 1,
                 n_sim = 10)
    
    mean_Rt<-wallinga_teunis(incid = df2$I, 
                             method = "parametric_si",
                             config = config)
    
    #adds back in days that were filtered out to match the days in the main dataframe
    mean_Rt$R$t_start <- mean_Rt$R$t_start + non0 
    mean_Rt$R$t_end <- mean_Rt$R$t_end + non0
    
    # #binds them into a dataframe
    rt_df<-cbind.data.frame(day = mean_Rt$R$t_start,
                            Rt = mean_Rt$R$`Mean(R)`,
                            lower = mean_Rt$R$`Quantile.0.025(R)`,
                            upper = mean_Rt$R$`Quantile.0.975(R)`)
    
  }else{
    # configuration of input data for R estimate
    # essentially we are estimating the serial intervals of SARS-CoV-2 (Omicron variant specific)
    # by drawing from two ( truncated normal ) distributions 
    # for the mean and standard deviation of the serial interval
    config <- make_config(list(mean_si = 3.5, 
                               std_si = 1))
    
    mean_Rt <- suppressMessages(EpiEstim::estimate_R(incid = df2$I,
                                    dt = 7L, 
                                    dt_out = 7L,
                                    iter = 50L,
                                    # grid = list(precision = 0.001, min = -1, max = 1),
                                    method="parametric_si",
                                    config = config))
    
    # #binds them into a dataframe
    rt_df<-cbind.data.frame(day = mean_Rt$R$t_start,
                            Rt = round(mean_Rt$R$`Mean(R)`, 3),
                            lower = round(mean_Rt$R$`Quantile.0.025(R)`, 3),
                            upper = round(mean_Rt$R$`Quantile.0.975(R)`, 3))
  }
  
  #merges the Rt value with the other variant data and renames Rt to have variant suffix
  merge <- df2 %>%
    arrange(days) %>% #keep in week so that the day variable lines up with the first week
    rowid_to_column(var = "week") %>% #used to merge with the estimate_R variable output for the day
    ungroup() |> 
    complete(days = full_seq(df2$days, period = 1)) |> 
    fill(c("name_states", "variant", "week"), .direction = "down") |> 
    mutate(across(is.numeric, ~round(zoo::na.spline(.x), 0))) |> 
    rowid_to_column(var = "day") |> 
    left_join(rt_df, by = c("day"))
  
  return(merge)
  
}

## rt_fun tryCatch
rt_safe<-function(x, walling_teunis = FALSE, ...){
  result <- tryCatch(rt_fun_week(df = x, wallinga_teunis = walling_teunis, ...), error = function(err) x)
  return(result)
}

## Creating estimates for Rt
## Try do it in parallel
rt_list<-rt_walling_teunis<-c()

## Estimates to use
estimates_df<-infections_variants_weekly

## Vectors to for loops
variants<-unique(estimates_df$variant)

states<-unique(estimates_df$name_states)

## Cori et al. Method
for (i in variants) {
  for (j in states) {
    
    tmp <- estimates_df |> 
      filter(variant == i, 
             name_states == j)
    
    if(nrow(tmp) == 0) next
    
    ## Try to handle when any Rt fails and continue it 
    rt_list[[i]][[j]]<-rt_safe(x = tmp)
    rm(tmp)
    gc()
    cat("Finished state: ", j, "\n")
  }
  rt_list[[i]]<-bind_rows(rt_list[[i]])
  
  ## Prompting messages, to monitor progress
  cat("Finished variant: ", i, "over all states \n")
}

rt_estimates<-bind_rows(rt_list)

vroom_write(x = rt_estimates, 
            file = "Output/Tables/rt_estimates_cori_method_weekly.tsv.xz")

# ## Walling-Teunis et al. Method
# for (i in variants_reduced) {
#   for (j in states) {
#     
#     tmp <- estimates_df |> 
#       filter(variant_reduced == i, 
#              name_states == j)
#     
#     if(nrow(tmp) == 0) next
#     
#     ## Try to handle when any Rt fails and continue it 
#     rt_walling_teunis[[i]][[j]]<-rt_safe(x = tmp, walling_teunis = TRUE)
#     rm(tmp)
#     gc()
#     cat("Finished state: ", j, "\n")
#   }
#   rt_walling_teunis[[i]]<-bind_rows(rt_walling_teunis[[i]])
#   
#   ## Prompting messages, to monitor progress
#   cat("Finished variant: ", i, "over all states \n")
# }
# 
# rt_estimates_walling_teunis<-bind_rows(rt_walling_teunis)
# 
# vroom_write(x = rt_walling_teunis, 
#             file = "Output/Tables/rt_estimates_walling_teunis_method.csv.xz")
