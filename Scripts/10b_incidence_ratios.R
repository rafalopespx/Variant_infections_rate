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
  mutate(incidence = (infections/pop)*1e5, 
         percentual_incidence = round(incidence/pop, 2)) 
# |>
#   ## Filtering to greater incidence than 100/100k
#   filter(incidence >= 100)

## Incidence ratios
incidence_ratios<-estimates_rt_incidence |> 
  reframe(total_incidence = sum(incidence, na.rm = T),
          .by = c("name_states", "variant"))

incidence_ratios |> 
  filter(name_states == "New York") |>
  ggplot(aes(x = reorder(variant, total_incidence, decreasing = T), 
             y = total_incidence, 
             fill = variant))+
  geom_col()+
  geom_text(aes(label = round(total_incidence/1e5, 2),
                size = 1.5),
            vjust = 1.5)+
  # scale_y_continuous(labels = scales::percent)+
  # facet_geo(name_states~., scales = "free_y")+
  labs(y = "Total Incidence \n per 100k", )+
  theme_minimal()+
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = "bottom")

# incidence_ratios<-incidence_ratios |>
#   select(days, name_states, starts_with('incidence')) |> 
#   ## BA.2 ratios
#   mutate(across(.cols = c(ends_with('BA.4*'):ends_with('XBB*')),
#                 .fns = ~.x/`incidence_Omicron BA.2*`, 
#                 .names = "{.col}/BA.2*_ratio")) |> 
#   ## BA.2/BA.1 ratios
#   mutate(across(.cols = c(ends_with('BA.2*')), 
#                 .fns = ~.x/`incidence_Omicron BA.1*`,
#                 .names = "{.col}/BA.1*_ratio")) |> 
#   ## BA.4/BA.5 ratio
#   mutate(across(.cols = c(ends_with('BA.4*')), 
#                 .fns = ~.x/`incidence_Omicron BA.5*`,
#                 .names = "{.col}/BA.5*_ratio")) |> 
#   select(days, name_states, ends_with("_ratio")) |> 
#   pivot_longer(cols = ends_with("_ratio"),
#                names_to = "ratios",
#                values_to = "value") |> 
#   mutate(ratios = str_replace(ratios, "incidence_Omicron ","")) |> 
#   mutate(ratios = str_replace(ratios, "_ratio",""))

## Joining the ratios
plt_incidence_ratios<-function(x){
  plt<-incidence_ratios |>
    filter(name_states == x) |> 
    ggplot(aes(x = days, y = incidence,
               col = variant, fill = variant))+
    geom_hline(yintercept = 1, show.legend = F, col = "grey9")+
    geom_col(position = position_dodge())+
    theme_minimal()+
    labs(x = "Ratios", y = "Total Incidence", title = x)+
    theme(legend.title = element_blank(), 
          legend.position = "bottom")+
    # facet_wrap(ratios~., ncol = 1, strip.position = "right", scales = "free_y")+
    scale_x_date(date_breaks = "2 months",
                 date_labels = "%b %y"
    )
}

states<-unique(incidence_ratios$name_states)

plt_incidence_ratios_list<-lapply(states, plt_incidence_ratios)

names(plt_incidence_ratios_list)<-states

plt_incidence_ratios_list$California
plt_incidence_ratios_list$Connecticut
plt_incidence_ratios_list$`New York`

## Function to estimate the average median, upper, lower
avg_over_something_days<-function(x, first_something_days){
  
  if(!missing(first_something_days)){
    first_something_days<-as.numeric(first_something_days)
  }
  
  
  x<-x |>
    filter(!is.na(value)) |> 
    # slice(first_something_days, 
    #       .by = c("name_states","ratios")) |> 
    summarise(avg_ratio = mean(value, na.rm = T),
              .by = c("name_states","ratios"))
  
  return(x)
}

## Average for the X first days
## 30 first days
avg_ratios_30_first_days <- incidence_ratios |> 
  summarise(avg_ratio = mean(value, na.rm = T),
            .by = c("name_states","ratios"))

# ## 60 first days
# avg_ratios_60_first_days <- incidence_ratios |> 
#   avg_over_something_days(first_something_days = 60)
# 
# ## 90 first days
# avg_ratios_90_first_days <- incidence_ratios |> 
#   avg_over_something_days(first_something_days = 90)
# 
# ## 120 first days
# avg_ratios_120_first_days <- incidence_ratios |> 
#   avg_over_something_days(first_something_days = 120)

avg_ratios_30_first_days |> 
  mutate(name_states2 = tidytext::reorder_within(name_states, avg_ratio, within = ratios)) |> 
  ggplot(aes(y = name_states2, 
             x = avg_ratio, 
             col = ratios))+
  geom_vline(xintercept = 1, show.legend = F, col = "grey9")+
  geom_point()+
  theme_minimal()+
  labs(y = "States", 
       x = "Average Incidence ratio", 
       subtitle = "",
       caption = "States reorder acording to average ratio")+
  theme(legend.title = element_blank(), 
        legend.position = "bottom")+
  facet_wrap(ratios~., nrow = 1, scales = "free", strip.position = "right")+
  tidytext::scale_y_reordered()

# avg_ratios_60_first_days |> 
#   mutate(name_states2 = tidytext::reorder_within(name_states, avg_ratio, within = ratios)) |> 
#   ggplot(aes(y = name_states2, 
#              x = avg_ratio, 
#              col = ratios))+
#   geom_vline(xintercept = 1, show.legend = F, col = "grey9")+
#   geom_point()+
#   theme_minimal()+
#   labs(y = "States", 
#        x = "Average Incidence ratio", 
#        subtitle = "First 60 days of ratios time series",
#        caption = "States reorder acording to average ratio")+
#   theme(legend.title = element_blank(), 
#         legend.position = "bottom")+
#   facet_wrap(ratios~., nrow = 1, scales = "free", strip.position = "right")+
#   tidytext::scale_y_reordered()

