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

## Plot daily
plot_infections_daily<-infections_variants_daily |> 
  # filter(variant != "Other") |>
  mutate(variant = droplevels(factor(variant))) |> 
  ggplot(aes(x = days, y = I, 
             col = variant, fill = variant))+
  geom_line()+
  theme_minimal()+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90))+
  facet_geo(name_states~., scales = "free_y")+
  scale_x_date(date_breaks = "4 months", 
               date_labels = "%b %y")+
  labs(x = "Date", 
       y = "Infections number per variant", 
       subtitle = "Daily basis")+
  scale_fill_brewer(aesthetics = c("color", "fill"), 
                    palette = "Spectral", 
                    name = "VOCs")
plot_infections_daily

## Saving infections plot
ggsave(filename = "Output/Plots/plt_infections_daily_states.png", 
       plot = plot_infections_daily, 
       width = 15,
       height = 9, 
       dpi = 100)

## Plot weekly
plot_infections_weekly<-infections_variants_weekly |> 
  # filter(variant != "Other") |>
  mutate(variant = droplevels(factor(variant))) |> 
  ggplot(aes(x = days, y = I, 
             col = variant, fill = variant))+
  geom_line()+
  theme_minimal()+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90))+
  facet_geo(name_states~., scales = "free_y")+
  scale_x_date(date_breaks = "4 months", 
               date_labels = "%b %y")+
  labs(x = "Date", 
       y = "Infections number per variant")+
  scale_fill_brewer(aesthetics = c("color", "fill"), 
                    palette = "Spectral", 
                    name = "VOCs")
plot_infections_weekly

## Saving infections plot
ggsave(filename = "Output/Plots/plt_infections_weekly_states.png", 
       plot = plot_infections_daily, 
       width = 15,
       height = 9, 
       dpi = 100)

## Plots for each state

plt_fun<-function(x, title){
  
  x<-x|>
    mutate(variant = droplevels(factor(variant))) |> 
    ggplot(aes(x = days, y = I, 
               col = variant, fill = variant))+
    geom_line()+
    theme_minimal()+
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90))+
    scale_x_date(date_breaks = "4 months", 
                 date_labels = "%b %y")+
    labs(x = "Date", 
         y = "Infections counts (N) \n  per variant",
         title = title)+
    colorspace::scale_fill_discrete_divergingx(name = "VOCs", 
                                               palette = "Spectral", 
                                               aesthetics = c("color", "fill"))
  
}

## Per states with facet per variant, daily
states<-unique(infections_variants_daily$name_states)

plot_infections_list_daily<-lapply(states, function(x){
  infections_variants_daily |> 
    filter(name_states == x) |>
    plt_fun(title = x)
  
  ggsave(filename = paste0("Output/Plots/States_infections/plt_infections_estimates_", x, "_daily.png"), 
         width = 15, 
         height = 9, 
         dpi = 100)
})

## Per states with facet per variant, weekly
states<-unique(infections_variants_weekly$name_states)

plot_infections_list_weekly<-lapply(states, function(x){
  infections_variants_weekly|> 
    filter(name_states == x) |>
    plt_fun(title = x)
  
  ggsave(filename = paste0("Output/Plots/States_infections/plt_infections_estimates_", x, "_weekly.png"), 
         width = 15, 
         height = 9, 
         dpi = 100)
})


## Split over states and variants
rt_split<-rt_estimates %>% 
  split(list(.$name_states, .$variant))

## Function to estimate the average Rts
mean_rt<-function(x, first_something_days){
  x<-x |> 
    filter(!is.na(Rt))
  
  if(!missing(first_something_days)){
    first_something_days<-as.numeric(first_something_days)
    x<-x[1:first_something_days,]
  }
  
  x<-x |> 
    group_by(name_states, variant) |> 
    summarise(mean_rt = mean(Rt, na.rm = T), 
              mean_lower = mean(lower, na.rm = T),
              mean_upper = mean(upper, na.rm = T))
  
  return(x)
}