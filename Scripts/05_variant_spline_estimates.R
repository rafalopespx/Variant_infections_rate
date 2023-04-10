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

## function to convert into daily counts, respecting the frequencies generated
daily_spread<-function(data, infections.col, week.col){
  
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
      # dividing by 7 to make it an average of the week
      dplyr::mutate(across(starts_with("infections"), ~.x/7)) |> 
      # Complete the dates
      complete(days = seq.Date(min(data[[week.col]]), 
                               max(data[[week.col]]), 
                               by = "days")) |> 
      fill(name_states, .direction = "down") |> 
      dplyr::mutate(across(where(is.numeric), ~zoo::na.spline(.x))) |> 
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
    # |> 
    #   # Creating Rt
    #   dplyr::mutate(across(ends_with("infections"), list("Rt" = ~round((.x/lag(.x))^(5.8), 2))))
    
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
    
    # ## Rt data.frame
    # data_rt<-data |> 
    #   ## Selecting columns of interest
    #   select(days, ends_with("_Rt")) |> 
    #   #pivot longer to split by variant for estimate_R
    #   pivot_longer(cols = c(ends_with("_Rt")),
    #                values_to = "Rt",
    #                names_to = "variant") |>
    #   #removed because it is annoying
    #   dplyr::mutate(variant = str_remove(variant, "freq_")) |>
    #   dplyr::mutate(variant = str_remove(variant, "_infections_Rt")) |>
    #   dplyr::arrange(variant) |> 
    #   filter(Rt != 0)
    
    ## Joining
    data_join<-data_infections 
    # |> 
    #   left_join(data_rt, by = c("days", "variant"))
    
  })
  
  names(data_list)<-states
  
  return(data_list)
}

estimates_variant_main<-estimates_variant |> 
  daily_spread(week.col = "epiweek", 
               infections.col = "infections")

estimates_df<-estimates_variant_main |> 
  bind_rows(.id = "name_states")

estimates_df |> 
  filter(name_states == "Connecticut") |>
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

## Vectors to for loops
variants<-unique(estimates_df$variant)

states<-unique(estimates_df$name_states)

for (i in variants) {
  for (j in states) {
    
    tmp <- estimates_df |> 
      filter(variant == i, 
             name_states == j)
    
    ## Try to handle when any Rt fails and continue it 
    rt_list[[i]][[j]]<-rt_safe(x = tmp) |> 
      mutate(name_states = j)
    rm(tmp)
    gc()
    cat("Finished state: ", j, "\n")
  }
  rt_list[[i]]<-rt_list[[i]] |> 
    reduce(rbind)
  
  ## Prompting messages, to monitor progress
  cat("Finished variant: ", i, "over all states \n")
}

test<-rt_list

test2<-map(test, reduce(.x = ~.x, .f = rbind))

test2 |> 
  # filter(Rt < 10) |>
  ggplot(aes(x = date_start, y = Rt, 
             ymin = lower, ymax = upper, 
             col = variant, fill = variant))+
  geom_line()+
  geom_ribbon(alpha = .5)+
  scale_x_date(date_breaks = "2 months", 
               date_labels = "%b-%y")+
  facet_geo(name_states~., grid = "us_state_grid1", scales = "free_y")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90))

# rt_list2<-lapply(estimates_list, function(x){
  x<-x |> 
    filter(variant == x, 
           name_states == y)
  
  z<-rt_fun(x)
  
  return(z)
}, name_states)

