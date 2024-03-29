## Cleaning the ambient
rm(list = ls())
gc()

## Loading Libraries
packs = c("tidyverse", "vroom", "jsonlite", "httr", "geofacet")
lapply(packs,require, character.only = TRUE)

## Loading functions
source("Scripts/Functions/functions.R")

## Reading the database
variant_count<-vroom("Data/variant_counts_us.csv.xz")

## Resetting to use just Omicron descendants subvariants
variant_count<-variant_count |> 
  filter(!voc_cdc %in% c("Alpha*", "Beta*", "Gamma*", "Delta*")) |> 
  mutate(voc_cdc = droplevels(factor(voc_cdc))) |> 
  mutate(voc_cdc = factor(voc_cdc,
                          levels = c("Omicron BA.1*", "Omicron BA.2*", 
                                     "Omicron BA.3*", "Omicron BA.4*", "Omicron BA.5*",
                                     "Omicron XBB*","Recombinant", "Other")))

variants_count_wide<-variant_count |> 
  pivot_wider(id_cols = c("epiweek", "name_states"), 
              names_from = "voc_cdc", 
              values_from = c("freq", "n")) |> 
  {\(.) {replace(.,is.na(.),0)}}() ## Trick to use replace(is.na(.), 0)

# ## CovidEstim State-level Estimates
# url<-GET(paste('https://api2.covidestim.org/latest_runs?geo_type=eq.state&select=*%2Ctimeseries(*)'))
# 
# covidestim_state<-fromJSON(rawToChar(url$content))
# name_states<-covidestim_state$geo_name
# covidestim_state<-covidestim_state[[8]]
# names(covidestim_state)<-name_states
# covidestim_state<-covidestim_state |>
#   bind_rows(.id = "name_states") |>
#   select(name_states, date, infections, infections_p2_5, infections_p97_5) |>
#   mutate(epiweek = end.of.epiweek(as.Date(date)))

# vroom_write(x = covidestim_state, file = "Output/Tables/covidestim_estimates_states.csv.xz")

# ## Removing heavy objects
# rm(url)
# gc()

covidestim_state<-vroom("Output/Tables/covidestim_estimates_states.csv.xz")

## Restriging the dates of analysis on both data-set
estimates_variant<-covidestim_state |>
  left_join(variants_count_wide,
            by = c("epiweek","name_states")) |> #merges CovidEstim data with our data and keeps only 
  filter(!is.na(infections))

infections_estimates<-function(data, daily = FALSE, infections_col){
  
  ## Creating states var to loop over it
  states<-unique(data$name_states)
  
  data_list<-lapply(states, function(x){
    
    data<-data |> 
      filter(name_states == x)
    
    if(daily){
      data2 <- data |> 
        # arrange by week for creating the days
        dplyr::arrange(epiweek) |> 
        # create variable for days
        dplyr::mutate(days = epiweek) |>
        # dividing by 7 to make it an average of the week
        dplyr::mutate(across(starts_with("infections"), ~.x/7)) |>
        # Complete the dates
        tidyr::complete(days = seq.Date(min(epiweek), 
                                 max(epiweek), 
                                 by = "days")) |> 
        tidyr::fill(c(name_states, epiweek), .direction = "down") |> 
        dplyr::mutate(across(where(is.numeric), ~zoo::na.approx(.x, maxgap = 6)))|> 
        # cutting off the frequencies less than 0.02, and setting as 0
        dplyr::mutate(across(starts_with("freq"), ~ifelse(.x <= 0.02, 0, .x))) |>
        #calculate variant specific number of infections, frequency is 0 to 100, not 0 to 1, has to divide by 100
        dplyr::mutate(across(starts_with("freq"), .fns = ~.x*{{ infections_col }}/100, .names = "{.col}_infections")) |>
        #remove unnecessary columns
        dplyr::select(days, name_states, ends_with("_infections")) |>
        #removes variants with no infections (after freq <= -.0.02) and thus no frequency
        dplyr::select_if(function(col) max(col, na.rm = T) != 0) |>
        #round infections to whole individuals
        dplyr::mutate(across(ends_with("infections"), ~round(.x, 0)))
      
    }else{
      
      data2 <- data |> 
        # arrange by week for creating the days
        dplyr::arrange(epiweek) |> 
        # create variable for days
        dplyr::mutate(days = epiweek) |>
        # cutting off the frequencies less than 0.02, and setting as 0
        dplyr::mutate(across(starts_with("freq"), ~ifelse(.x <= 0.02, 0, .x))) |>
        #calculate variant specific number of infections
        dplyr::mutate(across(starts_with("freq"), 
                             list("infections" = ~.x*{{ infections_col }}/100))) |>
        #remove unnecessary columns
        dplyr::select(days, ends_with("infections")) |>
        #removes variants with no infections (after freq <= -.0.02) and thus no frequency
        dplyr::select_if(function(col) max(col, na.rm = T) != 0) |>
        #round infections to whole individuals 
        dplyr::mutate(across(ends_with("infections"), ~round(.x, 0)))
      
    }
    
    ## Infections data.frame
    data_infections<- data2 |> 
      ## Selecting columns of interest
      dplyr::select(days, name_states, ends_with("_infections")) |> 
      #pivot longer to split by variant for estimate_R
      tidyr::pivot_longer(cols = ends_with("_infections"),
                   values_to = "I",
                   names_to = "variant") |>
      #removed because it is annoying
      dplyr::mutate(variant = str_remove(variant, "freq_")) |>
      dplyr::mutate(variant = str_remove(variant, "_infections")) |>
      dplyr::arrange(variant) |> 
      dplyr::filter(I != 0)
    
  })
  
  # names(data_list)<-states
  
  return(data_list)
}


## Creating estimates on weekly and daily basis

# ## Weekly
# infections_variants_weekly<-estimates_variant |> 
#   infections_estimates() |> 
#   bind_rows(.id = "name_states")
# 
# vroom_write(x = infections_variants_weekly, 
#             file = "Data/infections_estimates_variants_weekly.csv.xz")

## Daily
### Central estimate
infections_variants_daily<-estimates_variant |> 
  infections_estimates(daily = TRUE, infections_col = infections) |> 
  bind_rows()

### Upper estimate
infections_variants_daily_upper<-estimates_variant |> 
  infections_estimates(daily = TRUE, infections_col = infections_p97_5) |> 
  bind_rows() |> 
  rename(upper = I)

### Lower estimate
infections_variants_daily_lower<-estimates_variant |> 
  infections_estimates(daily = TRUE, infections_col = infections_p2_5) |> 
  bind_rows() |> 
  rename(lower = I)

### Altogether
infections_variants_daily <- infections_variants_daily |> 
  left_join(infections_variants_daily_upper) |> 
  left_join(infections_variants_daily_lower)

vroom_write(x = infections_variants_daily, 
            file = "Data/infections_estimates_variants_daily.csv.xz")

#

