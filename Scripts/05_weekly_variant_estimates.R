## Cleaning the ambient
rm(list = ls())
gc()

## Loading Libraries
packs = c("tidyverse", "vroom", "jsonlite", "httr", "EpiEstim")
lapply(packs,require, character.only = TRUE)

## Loading functions
source("Scripts/Functions/functions.R")
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
  mutate(infections = infections*freq) |> 
  mutate(I = round(infections, 0))

# configuration of input data for R estimate
# essentially we are estimating the serial intervals of SARS-CoV-2 (Omicron variant specific)
# by drawing from two ( truncated normal ) distributions 
# for the mean and standard deviation of the serial interval
config <- make_config(list(mean_si = 3.5, 
                           # std_mean_si = 1, 
                           # min_mean_si = 1, 
                           # max_mean_si = 6, # estimates for SARS-CoV-2 serial interval
                           std_si = 1
                           # std_std_si = 0.5, 
                           # min_std_si = 0.5, 
                           # max_std_si = 1.5
                           )
)

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
  
  # main R estimate function:
  # will search for column named I which was created in the ci_fun but explicitly named here
  mean_Rt <- estimate_R(tmp$I, 
                        dt = 7L,
                        dt_out = 7L,
                        iter = 10L,
                        grid = list(precision = 0.001, min = -1, max = 1),
                        method="parametric_si",
                        config = config)
  
  # #binds them into a dataframe
  rt_df<-cbind.data.frame(day = mean_Rt$R$t_start,
                          Rt = mean_Rt$R$`Mean(R)`,
                          lower = mean_Rt$R$`Quantile.0.05(R)`,
                          upper = mean_Rt$R$`Quantile.0.95(R)`)
  
  #merges the Rt value with the other variant data and renames Rt to have variant suffix
  rt_list[[i]][[j]]<- tmp %>%
    dplyr::mutate(date = as.Date(date)) |> 
    arrange(date) %>% #keep in week so that the day variable lines up with the first week
    mutate(day = 1:nrow(tmp)) %>% #used to merge with the estimate_R variable output for the day
    left_join(rt_df, by = "day") |> 
    filter(!is.na(Rt)) |> 
    select(name_states, epiweek, variant, I, Rt, lower, upper)
  
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


rt_safely<-possibly(.f = rt_fun, quiet = T)

rt_list[[1]] |> 
  filter(name_states == "Connecticut") |>
  ggplot(aes(x = epiweek, y = Rt, 
             ymin = lower, ymax = upper, 
             col = variant, fill = variant))+
  geom_line()+
  geom_ribbon(alpha = .5)+
  scale_x_date(date_breaks = "2 months", 
               date_labels = "%b-%y")


