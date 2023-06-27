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

rt_splines<-function(data,
                     N_weeks, 
                     N_weeks_before = 28/7,
                     n_spl_rt_knotwidth = 2){
  
  if(missing(N_weeks)){
    N_weeks<-nrow(data)
  }
  
  if(missing(N_weeks_before)){
    N_weeks_before<-1
  }
  
  n_spl_par_rt <- max(4,ceiling((N_weeks + N_weeks_before)/n_spl_rt_knotwidth))
  des_mat_rt <- splines::bs(data$I, 
                            1:(N_weeks + N_weeks_before), 
                            df=n_spl_par_rt, 
                            degree=3, 
                            intercept=T)
  
  n_spl_par_dx <- max(4,ceiling((N_weeks + N_weeks_before)/3)) 
  des_mat_dx <- splines::bs(data$I, 
                            1:(N_weeks + N_weeks_before), 
                            df=n_spl_par_dx, 
                            degree=3, 
                            intercept=T) 
  
  spl_basis_rt = as.matrix(as.data.frame(des_mat_rt))
  spl_basis_dx = as.matrix(as.data.frame(des_mat_dx))
  
  # logRt <- spl_basis_rt*spl_par_rt
  
}

states<-unique(infections_variants$name_states)

variants<-unique(infections_variants$variant)

rt_list<-c()

for (i in states) {
  for (j in variants) {
    
    tmp<-infections_variants |> 
      filter(name_states == i, 
             variant == j)
    
    rt_list[[i]][[j]]<-rt_fun(tmp)
  }
}








