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
  select(-c(n_Other,freq_Other)) |>
  select_if(function(col)max(col) != 0)

## function to convert into daily counts, respecting the frequencies generated
daily_spread<-function(data, infections.col, week.col, freq.col){
  
  ## Creating states var to loop over it
  states<-unique(data$name_states)
  
  data_list<-lapply(states, function(x){
    
    data<-data |> 
      filter(name_states == x)
    
    data<-data[rep(seq_len(nrow(data)), 7),]
    
    data <- data |> 
      #arrange by week for creating the days
      dplyr::arrange(data[[week.col]]) |> 
      #create variable for days
      dplyr::mutate(days = seq.Date(from = min(data[[week.col]]),
                                    # Remember to change this for .data pronouns or something similar
                                    to = min(data[[week.col]])+nrow(data)-1, 
                                    by = "days")) |> 
      # Divide by 7 to convert weekly total to avg daily
      dplyr::mutate(across(all_of(starts_with("infections")), ~.x/7)) |> 
      # Revert freq back to original because freq is for the entire week
      dplyr::mutate(across(starts_with("freq"), ~.x*7)) |> 
      # daily rolling mean
      dplyr::mutate(across(where(is.numeric), ~zoo::rollmean(.x, 
                                                             k = 7, 
                                                             fill = 0))) |> 
      # cutting off the frequencies less than 0.02, and setting as 0
      dplyr::mutate(across(starts_with("freq"), ~ifelse(.x <= 0.02, 0, .x))) |>  
      #calculate variant specific number of infections
      dplyr::mutate(across(starts_with("freq"), list("infections" = ~.x*infections))) |> 
      #round infections to whole individuals 
      dplyr::mutate(across(ends_with("infections"), ~round(.x, 0)))
    
    ## Cleaning and putting on longer format
    data <- data |>
      #remove unnecessary columns
      select(days, ends_with("infections")) |>
      #removes variants with no infections (after freq <= -.0.02) and thus no frequency
      select_if(function(col) max(col) != 0) |>
      #pivot longer to split by variant for estimate_R
      pivot_longer(cols = c(starts_with("freq")),
                   values_to = "I",
                   names_to = "variant") |>
      #removed because it is annoying
      dplyr::mutate(variant = str_remove(variant, "freq_")) |>
      dplyr::mutate(variant = str_remove(variant, "_infections")) |>
      dplyr::arrange(variant)
    
    # ## Create a list by variant ready to calculate estimate_R
    data<-data |>
      group_by(variant) |>
      arrange(variant) |>
      group_split()
    
  })
  
  names(data_list)<-states
  
  return(data_list)
}

estimates_variant_main<-estimates_variant |> 
  daily_spread(week.col = "epiweek", 
               infections.col = "infections")

# estimates_df<-estimates_variant_main |>
# bind_rows(.id = "name_states")

# estimates_df |>
#   # filter(name_states == "Connecticut") |>
#   # filter(days >= "2022-01-01", days <= "2022-01-31") |>
#   ggplot(aes(x = days, y = I, col = variant, fill = variant))+
#   geom_line()+
#   facet_geo(name_states~., grid = "us_state_grid1", scales = "free_y")+
#   scale_x_date(date_breaks = "4 months", date_labels = "%b-%Y")+
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 90),
#         legend.position = "bottom")+
#   scale_fill_brewer(palette = "Spectral", name = "VOCs", aesthetics = c("color", "fill"))+
#   labs(x = "Date", y = "Infections counts \n per variant")

## Creating estimates for Rt
## Try do it in parallel
# Function to estimate Rt for a given state-variant combination
estimate_Rt <- function(data, state, variant) {
  # Subset the data for the current state-variant combination
  subset_data <- data %>%
    filter(name_states == state,
           variant == variant)
  
  # Convert the data to a time series
  ts_data <- zoo::zoo(subset_data$infections, 
                      subset_data$dates)
  
  # Estimate the Rt using the EpiEstim package
  epi_result <- EpiEstim::estimate_R(zoo::na.omit(ts_data), 
                                     method = "parametric_si",
                                     config = EpiEstim::make_config(list(mean_si = 3.5, std_si = 1)))
  
  # Return the estimated Rt as a data frame
  data.frame(state = state, 
             variant = variant, 
             dates = index(epi_result$R),
             Rt = coredata(epi_result$R))
}

## Vectors to for loops
estimates_df<-estimates_variant_main |>
  bind_rows(.id = "name_states") |> 
  group_by(name_states, variant) %>%
  filter(sum(infections) >= 10) %>%
  ungroup()

estimates_Rt<-estimates_variant_main %>% 
  rowwise() %>%
  do.call(rbind, map2(.$name_states, .$variant, estimate_Rt))

results<-c()

# For each state-variant combination, calculate the Rt
for (i in unique(estimates_df$name_states)) {
  for (j in unique(estimates_df$variant)) {
    # Subset the data for the current state-variant combination
    subset_data <- estimates_df[estimates_df$name_states == i & estimates_df$variant == j, ]
    
    # Convert the data to an EpiEstim object
    epi_data <- subset_data
    
    # Estimate the Rt using the EpiEstim package
    epi_result <- estimate_R(epi_data, 
                             method = "parametric_si", 
                             config = make_config(list(mean_si = 3.5, 
                                                       std_si = 1)))
    
    epi_df<-epi_result$R |> 
      as.data.frame() |> 
      mutate(name_states = i, 
             variant = j)
    
    
    # Store the results in a new dataframe
    results <- rbind(results, epi_df)
  }
}


for (i in variants) {
  for (j in states) {
    
    tmp <- estimates_df |> 
      filter(variant == i, 
             name_states == j)
    
    ## Try to handle when any Rt fails and continue it 
    rt_list[[i]][[j]]<-rt_fun(df = tmp)
    rm(tmp)
    gc()
    cat("Finished state: ", j, "\n")
  }
  rt_list[[i]]<-rt_list[[i]] |>
    bind_rows(.id = "name_state")
  
  ## Prompting messages, to monitor progress
  cat("Finished variant: ", i, "over all states \n")
}

test<-rt_list |> 
  bind_rows()

rt_list$`Omicron BA.1*` |> 
  # filter(Rt < 10) |>
  ggplot(aes(x = date_start, y = Rt, 
             ymin = lower, ymax = upper, 
             col = variant, fill = variant))+
  geom_line()+
  geom_ribbon(alpha = .5)+
  scale_x_date(date_breaks = "2 months", 
               date_labels = "%b-%y")+
  facet_geo(name_states~., grid = "us_state_grid1", scales = "free_y")+
  theme(axis.text.x = element_text(angle = 90))

rt_list2<-lapply(estimates_list, function(x){
  x<-x |> 
    filter(variant == x, 
           name_states == y)
  
  z<-rt_fun(x)
  
  return(z)
}, name_states)

