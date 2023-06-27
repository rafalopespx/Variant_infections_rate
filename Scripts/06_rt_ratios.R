## Cleaning the ambient
rm(list = ls())
gc()

## Loading Libraries
packs = c("tidyverse", "vroom", "jsonlite", "httr", "geofacet", "EpiEstim")
lapply(packs,require, character.only = TRUE)

## Loading functions
source("Scripts/Functions/functions.R")

## Loading the estimates
rt_estimates<-vroom("Output/Tables/rt_estimates_cori_method_daily.tsv.xz")|>
  filter(!variant %in% c("Alpha*", "Beta*", "Gamma*", "Delta*", "Omicron BA.3*", "Recombinant", "Other")) |>
  mutate(variant = droplevels(factor(variant))) |>
  mutate(variant = factor(variant,
                          levels = c("Omicron BA.1*", "Omicron BA.2*", "Omicron BA.4*", "Omicron BA.5*",
                                     "Omicron XBB*")))

## Population by states
pop_states<-vroom("https://raw.githubusercontent.com/covidestim/covidestim-sources/master/data-sources/statepop.csv")

## Loading the frequencies of variants
variant_count<-vroom("Data/variant_counts_us.csv.xz") |>
  rename(days = epiweek,
         variant = voc_cdc) |>
  filter(!variant %in% c("Alpha*", "Beta*", "Gamma*", "Delta*", "Omicron BA.3*", "Recombinant", "Other")) |>
  mutate(variant = droplevels(factor(variant))) |>
  mutate(variant = factor(variant,
                          levels = c("Omicron BA.1*", "Omicron BA.2*", "Omicron BA.4*", "Omicron BA.5*",
                                     "Omicron XBB*")))

## Joining Rt with variant counts
estimates_rt_incidence<-rt_estimates |>
  left_join(variant_count) |>
  mutate(variant = droplevels(factor(variant))) |> 
  left_join(pop_states, by = c("name_states" = "state")) |> 
  ## Renaming infections
  rename(infections = I) |> 
  ## Creating incidence per 100k
  mutate(incidence = (infections/pop)*1e5)
# |>
#   ## Filtering to greater incidence than 100/100k
#   filter(incidence >= 5)

## Rt ratios
rt_ratios<-estimates_rt_incidence |> 
  rename(median = Rt) |> 
  pivot_wider(id_cols = c("days", "name_states"),
              names_from = variant,
              values_from = c(median, upper, lower))

## Ratios calculation
median_ratios <- rt_ratios |>
  select(days, name_states, starts_with('median')) |> 
  ## BA.2 ratios
  mutate(across(.cols = c(ends_with('BA.4*'):ends_with('XBB*')),
                .fns = ~.x/`median_Omicron BA.2*`, 
                .names = "{.col}/BA.2*_ratio")) |> 
  ## BA.2/BA.1 ratio
  mutate(across(.cols = c(ends_with('BA.2*')), 
                .fns = ~.x/`median_Omicron BA.1*`,
                .names = "{.col}/BA.1*_ratio")) |> 
  ## BA.5 ratios
  mutate(across(.cols = c(ends_with('BA.4*'), ends_with("XBB*")), 
                .fns = ~.x/`median_Omicron BA.5*`,
                .names = "{.col}/BA.5*_ratio")) |> 
  ## BA.5*/BA.4*
  mutate(`median_Omicron BA.5*/BA.4*_ratio` = `median_Omicron BA.5*`/`median_Omicron BA.4*`) |> 
  select(days, name_states, ends_with("_ratio")) |> 
  pivot_longer(cols = ends_with("_ratio"),
               names_to = "ratios",
               values_to = "median") |> 
  mutate(ratios = str_replace(ratios, "median_Omicron ","")) |> 
  mutate(ratios = str_replace(ratios, "_ratio",""))

upper_ratios <- rt_ratios |>
  select(days, name_states, starts_with('upper')) |> 
  ## BA.2 ratios
  mutate(across(.cols = c(ends_with('BA.4*'):ends_with('XBB*')),
                .fns = ~.x/`upper_Omicron BA.2*`, 
                .names = "{.col}/BA.2*_ratio")) |> 
  ## BA.2/BA.1 ratio
  mutate(across(.cols = c(ends_with('BA.2*')), 
                .fns = ~.x/`upper_Omicron BA.1*`,
                .names = "{.col}/BA.1*_ratio")) |> 
  ## BA.4/BA.5 ratio
  mutate(across(.cols = c(ends_with('BA.4*'), ends_with("XBB*")), 
                .fns = ~.x/`upper_Omicron BA.5*`,
                .names = "{.col}/BA.5*_ratio")) |> 
  ## BA.5*/BA.4*
  mutate(`upper_Omicron BA.5*/BA.4*_ratio` = `upper_Omicron BA.5*`/`upper_Omicron BA.4*`) |> 
  select(days, name_states, ends_with("_ratio")) |> 
  pivot_longer(cols = ends_with("_ratio"),
               names_to = "ratios",
               values_to = "upper") |> 
  mutate(ratios = str_replace(ratios, "upper_Omicron ","")) |> 
  mutate(ratios = str_replace(ratios, "_ratio",""))

lower_ratios <- rt_ratios |>
  select(days, name_states, starts_with('lower')) |> 
  ## BA.2 ratios
  mutate(across(.cols = c(ends_with('BA.4*'):ends_with('XBB*')),
                .fns = ~.x/`lower_Omicron BA.2*`, 
                .names = "{.col}/BA.2*_ratio")) |> 
  ## BA.2/BA.1 ratio
  mutate(across(.cols = c(ends_with('BA.2*')), 
                .fns = ~.x/`lower_Omicron BA.1*`,
                .names = "{.col}/BA.1*_ratio")) |> 
  ## BA.4/BA.5 ratio
  mutate(across(.cols = c(ends_with('BA.4*'), ends_with("XBB*")), 
                .fns = ~.x/`lower_Omicron BA.5*`,
                .names = "{.col}/BA.5*_ratio")) |> 
  ## BA.5*/BA.4*
  mutate(`lower_Omicron BA.5*/BA.4*_ratio` = `lower_Omicron BA.5*`/`lower_Omicron BA.4*`) |> 
  select(days, name_states, ends_with("_ratio")) |> 
  pivot_longer(cols = ends_with("_ratio"),
               names_to = "ratios",
               values_to = "lower") |> 
  mutate(ratios = str_replace(ratios, "lower_Omicron ","")) |> 
  mutate(ratios = str_replace(ratios, "_ratio",""))

## Rt ratios
rt_ratios<-median_ratios |> 
  left_join(upper_ratios) |> 
  left_join(lower_ratios) 
# |> 
#   filter(ratios != "XBB*/BA.2*", ratios != "BA.4*/BA.5*") |> 
#   mutate(ratios = factor(droplevels(factor(ratios)), 
#                          levels = c("BA.2*/BA.1*", "BA.4*/BA.2*", "BA.5*/BA.2*", 
#                                     "BA.5*/BA.4*", "XBB*/BA.5*")))

vroom_write(x = rt_ratios, file = "Output/Tables/rt_ratios.csv.xz")

## Joining the ratios
plt_rt_ratios<-function(x){
  plt<-rt_ratios |>
    filter(name_states == x) |> 
    ggplot(aes(x = days, y = median, 
               ymin = lower, ymax = upper,
               col = ratios, fill = ratios))+
    geom_hline(yintercept = 1, show.legend = F, col = "grey9")+
    geom_line()+
    geom_ribbon(alpha = .15)+
    theme_minimal()+
    labs(x = "Date", y = "Rt Ratios value", title = x)+
    theme(legend.title = element_blank(), 
          legend.position = "bottom")+
    facet_wrap(ratios~., ncol = 1, strip.position = "right")+
    scale_x_date(date_breaks = "2 months",
                 date_labels = "%b %y"
    )+
    MetBrewer::scale_color_met_d(palette_name = "Austria", name = "", direction = -1)+
    MetBrewer::scale_fill_met_d(palette_name = "Austria", name = "", direction = -1)
}

states<-unique(rt_ratios$name_states)

plt_rt_ratios_list<-lapply(states, plt_rt_ratios)

names(plt_rt_ratios_list)<-states

plt_rt_ratios_list$California
plt_rt_ratios_list$Connecticut
plt_rt_ratios_list$`New York`
plt_rt_ratios_list$Texas

## Function to estimate the average median, upper, lower
avg_over_something_days<-function(x, first_something_days){
  
  if(!missing(first_something_days)){
    first_something_days<-as.numeric(first_something_days)
  }
  
  
  x<-x |>
    filter(!is.na(median)) |> 
    slice(first_something_days, 
          .by = c("name_states","ratios")) |> 
    summarise(avg_median = ave(median, FUN = mean),
              avg_upper = ave(upper, FUN = mean),
              avg_lower = ave(lower, FUN = mean),
              .by = c("name_states","ratios"))
  
  return(x)
}

## Average for the X first days
## 30 first days
avg_ratios_30_first_days <- rt_ratios |> 
  avg_over_something_days(first_something_days = 30)

## 60 first days
avg_ratios_60_first_days <- rt_ratios |> 
  avg_over_something_days(first_something_days = 60)

## 90 first days
avg_ratios_90_first_days <- rt_ratios |> 
  avg_over_something_days(first_something_days = 90)

## 120 first days
avg_ratios_120_first_days <- rt_ratios |> 
  avg_over_something_days(first_something_days = 120)

avg_ratios_30_first_days |> 
  mutate(name_states2 = tidytext::reorder_within(name_states, avg_median, within = ratios)) |> 
  ggplot(aes(y = name_states2, 
             x = avg_median, 
             xmin = avg_lower, xmax = avg_upper,
             col = ratios))+
  geom_vline(xintercept = 1, show.legend = F, col = "grey9")+
  geom_pointrange()+
  theme_minimal()+
  labs(y = "States", 
       x = "Average Rt ratio", 
       subtitle = "First 30 days of ratios time series",
       caption = "States reorder acording to average ratio")+
  theme(legend.title = element_blank(), 
        legend.position = "bottom")+
  facet_wrap(ratios~., nrow = 1, scales = "free", strip.position = "right")+
  tidytext::scale_y_reordered()+
  MetBrewer::scale_color_met_d(palette_name = "Austria", name = "", direction = -1)+
  MetBrewer::scale_fill_met_d(palette_name = "Austria", name = "", direction = -1)

avg_ratios_60_first_days |> 
  mutate(name_states2 = tidytext::reorder_within(name_states, avg_median, within = ratios)) |> 
  ggplot(aes(y = name_states2, 
             x = avg_median, 
             xmin = avg_lower, xmax = avg_upper,
             col = ratios))+
  geom_vline(xintercept = 1, show.legend = F, col = "grey9")+
  geom_pointrange()+
  theme_minimal()+
  labs(y = "States", 
       x = "Average Rt ratio", 
       subtitle = "First 60 days of ratios time series",
       caption = "States reorder acording to average ratio")+
  theme(legend.title = element_blank(), 
        legend.position = "bottom")+
  facet_wrap(ratios~., nrow = 1, scales = "free", strip.position = "right")+
  tidytext::scale_y_reordered()+
  MetBrewer::scale_color_met_d(palette_name = "Austria", name = "", direction = -1)+
  MetBrewer::scale_fill_met_d(palette_name = "Austria", name = "", direction = -1)

#
