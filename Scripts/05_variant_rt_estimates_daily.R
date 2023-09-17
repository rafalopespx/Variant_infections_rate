## Cleaning the ambient
rm(list = ls())
gc()

## Loading Libraries
packs = c("tidyverse", "vroom", "jsonlite", "httr", "geofacet", "EpiEstim", "doParallel")
lapply(packs,require, character.only = TRUE)

## Loading functions
source("Scripts/Functions/functions.R")
# source("Scripts/estimate_rt_ro_fun.R")

infections_variants_daily<-vroom("Data/infections_estimates_variants_daily.csv.xz")

#run for everything -Other
rt_fun <- function(df, wallinga_teunis = FALSE){
  suppressPackageStartupMessages(require(EpiEstim))
  
  ## Completing the dates, if there is some missing
  df<-df |> 
    tidyr::complete(days = full_seq(days, period = 1)) |> 
    tidyr::fill(c(name_states, variant), .direction = "down") |> 
    dplyr::mutate(I = round(zoo::na.approx(I), 0))
  
  #1st day with infections of variant to start the R estimate otherwise R estimate artificially high
  non0 <- min(which(df$I > 0)) 
  #last day with infections of variant to start the R estimate otherwise R estimate artificially high
  end0 <- max(which(df$I > 0))
  #dataframe filtered where there is the first case of variant
  df2 <- df[non0:end0,]
  
  
  ## Uncomment to switch to a different way to filter the infections
  # #1st day with infections of variant to start the R estimate otherwise R estimate artificially high
  # non0 <- which.min(df$I > 0)
  # #last day with more than 10 infections, after this is filtered out
  # end0 <- which.max(df$I < 10)
  # #dataframe filtered where there is the first case of variant
  # df2 <- df[non0:end0,]
  
  #input of interval for R estimate
  #estimate burn in period for R estimate: t_end - t_start
  t_start<-seq(2, nrow(df2)-15) 
  t_end<- t_start + 15
  
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
    config <- list(mean_si = 3.5, 
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
    
    mean_Rt <- EpiEstim::estimate_R(df2$I, 
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
  }
  
  #merges the Rt value with the other variant data and renames Rt to have variant suffix
  merge <- df %>%
    dplyr::arrange(days) %>% #keep in week so that the day variable lines up with the first week
    tibble::rowid_to_column(var = "day") %>% #used to merge with the estimate_R variable output for the day
    dplyr::left_join(rt_df, by = c("day"))
  
  return(merge)
  
}

## rt_fun tryCatch
rt_safe<-function(x, walling_teunis = FALSE){
  result <- tryCatch(rt_fun(df = x, wallinga_teunis = walling_teunis), error = function(err) x)
  return(result)
}

## Creating estimates for Rt
## Try do it in parallel
rt_list<-rt_walling_teunis<-c()

## Estimates to use
estimates_df<-infections_variants_daily

## Vectors to for loops
variants<-unique(estimates_df$variant)

variants_reduced<-unique(estimates_df$variant)

states<-unique(estimates_df$name_states)

## Cori et al. Method
# for (i in variants) {
cores<-detectCores() - 1
cl<-makeCluster(cores)
clusterExport(cl, c("rt_safe", "estimates_df"))
registerDoParallel(cl)

rt_list<-
  foreach(i = variants, 
          .combine = "bind_rows") %:%
  foreach(j = states, 
          .packages = c("dplyr", "tidyr"), 
          .combine = "bind_rows") %dopar% {
            
            tmp <- estimates_df |> 
              dplyr::filter(variant == i, 
                            name_states == j)
            
            # if(nrow(tmp) == 0) next
            
            ## Try to handle when any Rt fails and continue it 
            return(rt_safe(x = tmp))
            # rm(tmp)
            # gc()
            # cat("Finished state: ", j, "\n")
          }
stopCluster(cl)

## Prompting messages, to monitor progress
# cat("Finished variant: ", i, "over all states \n")
# }

vroom_write(x = rt_list, 
            file = "Output/Tables/rt_estimates_cori_method_daily.tsv.xz")
#
