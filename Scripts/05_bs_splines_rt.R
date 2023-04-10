## Cleaning the ambient
rm(list = ls())
gc()

## Loading Libraries
packs = c("tidyverse", "vroom", "jsonlite", "httr", "geofacet")
lapply(packs,require, character.only = TRUE)

## Loading functions
source("Scripts/Functions/functions.R")
# source("Scripts/estimate_rt_ro_fun.R")

## Reading the database
variants_count<-vroom("Data/variant_counts_us.csv.xz")

variants_count_wide<-variants_count |> 
  pivot_wider(id_cols = c("epiweek", "name_states"), 
              names_from = "voc_cdc", 
              values_from = c("freq", "n")) |> 
  {\(.) {replace(.,is.na(.),0)}}() ## Trick to use replace(is.na(.), 0)

## CovidEstim State-level Estimates
url<-GET(paste('https://api2.covidestim.org/latest_runs?geo_type=eq.state&select=*%2Ctimeseries(*)'))

covidestim_state<-fromJSON(rawToChar(url$content))
name_states<-covidestim_state$geo_name
covidestim_state<-covidestim_state[[8]]
names(covidestim_state)<-name_states
covidestim_state<-covidestim_state |> 
  bind_rows(.id = "name_states") |> 
  select(name_states, date, infections, infections_p2_5, infections_p97_5) |> 
  mutate(epiweek = end.of.epiweek(as.Date(date)))

## Removing heavy objects
rm(url)
gc()

## Restriging the dates of analysis on both data-set
estimates_variant<-variants_count_wide |>
  left_join(covidestim_state,
            by = c("epiweek","name_states")) |> #merges CovidEstim data with our data and keeps only 
  filter(!is.na(infections))|>
  select(-c(n_Other,freq_Other))

infections_estimates<-function(data, infections.col, week.col, daily = FALSE){
  
  ## Creating states var to loop over it
  states<-unique(data$name_states)
  
  data_list<-lapply(states, function(x){
    
    data<-data |> 
      filter(name_states == x)
    
    data <- data |> 
      # arrange by week for creating the days
      dplyr::arrange(data[[week.col]]) |> 
      # create variable for days
      dplyr::mutate(days = epiweek) |>
      # cutting off the frequencies less than 0.02, and setting as 0
      dplyr::mutate(across(starts_with("freq"), ~ifelse(.x <= 0.02, 0, .x))) |>
      #calculate variant specific number of infections
      dplyr::mutate(across(starts_with("freq"), list("infections" = ~.x*infections))) |>
      #remove unnecessary columns
      select(days, ends_with("infections")) |>
      #removes variants with no infections (after freq <= -.0.02) and thus no frequency
      select_if(function(col) max(col) != 0) |>
      #round infections to whole individuals 
      dplyr::mutate(across(ends_with("infections"), ~round(.x, 0)))
    
    ## Infections data.frame
    data_infections<- data |> 
      ## Selecting columns of interest
      select(days, ends_with("_infections")) |> 
      #pivot longer to split by variant for estimate_R
      pivot_longer(cols = c(ends_with("_infections")),
                   values_to = "I",
                   names_to = "variant") |>
      #removed because it is annoying
      dplyr::mutate(variant = str_remove(variant, "freq_")) |>
      dplyr::mutate(variant = str_remove(variant, "_infections")) |>
      dplyr::arrange(variant) |> 
      filter(I != 0)
    
  })
  
  names(data_list)<-states
  
  return(data_list)
}

infections_variants<-estimates_variant |> 
  infections_estimates(week.col = "epiweek", 
                       infections.col = "infections") |> 
  bind_rows(.id = "name_states")

vroom_write(x = infections_variants, 
            file = "Data/infections_estimates_variants_weekly.csv.xz")

infections_variants |> 
  filter(name_states == "Alabama") |>
  ggplot(aes(x = days, y = I, col = variant, fill = variant))+
  geom_line()+
  scale_x_date(date_breaks = "4 months", date_labels = "%b-%Y")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom")+
  scale_fill_brewer(palette = "Paired", 
                    name = "VOCs", 
                    aesthetics = c("color", "fill"))+
  labs(x = "Date", y = "Rt \n per variant", subtitle = "Connecticut")


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








