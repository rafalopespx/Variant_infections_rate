## Cleaning the ambient
rm(list = ls())
gc()

## Loading Libraries
packs = c("tidyverse", "vroom", "jsonlite", "httr")
lapply(packs,require, character.only = TRUE)

## Loading functions
source("Scripts/functions.R")
# source("Scripts/estimate_rt_ro_fun.R")

## Reading the database
variants_count<-vroom("Data/variant_counts_us.csv.xz")

# variants_count_wide<-variants_count |> 
#   pivot_wider(id_cols = c("epiweek", "name_states"), 
#               names_from = "voc_cdc", 
#               values_from = c("freq", "n")) |> 
#   {\(.) {replace(.,is.na(.),0)}}() ## Trick to use replace(is.na(.), 0)

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
estimates_variant<-variants_count |>
  dplyr::rename(variant = voc_cdc) |> 
  left_join(covidestim_state,
            by = c("epiweek","name_states")) |> #merges CovidEstim data with our data and keeps only 
  filter(!is.na(infections))|>
  # select(-c(n.Other,freq.Other)) |> 
  select_if(function(col)max(col) != 0) |> 
  dplyr::mutate(across(starts_with("infections"), ~.x/7)) |> 
  group_by(name_states, variant) |>
  complete(epiweek = seq.Date(min(epiweek), 
                              max(epiweek), 
                              by = "day")) |>
  fill(c(n, freq, date, starts_with("infections")), 
       .direction = "down") |> 
  mutate(I = round(infections, 0))

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
    dplyr::mutate(date = as.Date(date)) |> 
    arrange(date) %>% #keep in week so that the day variable lines up with the first week
    mutate(day = 1:nrow(df)) %>% #used to merge with the estimate_R variable output for the day
    left_join(rt_df, by = "day") |> 
    filter(!is.na(Rt)) |> 
    select(name_states, epiweek, variant, I, Rt, lower, upper)
  
  return(merge)
  
}

rt_safely<-possibly(.f = rt_fun, quiet = T)

## Creating estimates for Rt
## Try do it in parallel
rt_list<-c()

## Vectors to for loops
variants<-unique(estimates_variant$variant)

states<-unique(estimates_variant$name_states)

for (i in variants) {
  for (j in states) {
    
    tmp <- estimates_variant |> 
      filter(variant == i, 
             name_states == j)
    
    ## Try to handle when any Rt fails and continue it 
    rt_list[[i]][[j]]<- rt_safely(tmp)
    
    ## Prompting messages, to monitor progress
    rm(tmp)
    gc()
    cat("Finished state: ", j, "\n")
  }
  
  ## Bind rows to a single data.frame for the Variant over all states
  rt_list[[i]]<-rt_list[[i]] |>
    bind_rows(.id = "name_state")
  
  ## Prompting messages, to monitor progress
  cat("Finished variant: ", i, "over all states \n")
}


test<-rt_list |> 
  bind_rows(.id = "variant")


test |> 
  filter(name_state == "Connecticut") |> 
  ggplot(aes(x = epiweek, y = Rt, 
             ymin = lower, ymax = upper, 
             col = variant, fill = variant))+
  geom_line()+
  geom_ribbon(alpha = .5)+
  scale_x_date(date_breaks = "2 months", 
               date_labels = "%b-%y")+
  theme_minimal()+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90))+
  facet_wrap(variant~., scales = "free_y")
  # facet_geo(name_states~., grid = "us_state_grid1")
  


