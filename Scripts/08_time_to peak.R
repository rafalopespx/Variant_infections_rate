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

## Days to peak
# ### By Incidence
# days_to_peak_by_incidence<-estimates_rt_incidence |> 
#   slice_max(incidence, 
#             n = 1, 
#             by = c(name_states, variant), 
#             na_rm = T, 
#             with_ties = F) |> 
#   select(day, name_states, variant)

### By Frequency
days_to_peak_by_frequency<-estimates_rt_incidence |> 
  slice_max(freq, 
            n = 1, 
            by = c(name_states, variant), 
            na_rm = T, 
            with_ties = F) |> 
  select(day, name_states, variant)

days_to_peak_by_frequency_ratios<-days_to_peak_by_frequency |> 
  pivot_wider(id_cols = "name_states",
              names_from = variant,
              values_from = day) |> 
  ## BA.2 ratios
  mutate(across(.cols = c(ends_with('BA.4*'):ends_with('XBB*')),
                .fns = ~.x/`Omicron BA.2*`, 
                .names = "{.col}/BA.2*_ratio")) |> 
  ## BA.2/BA.1 ratio
  mutate(across(.cols = c(ends_with('BA.2*')), 
                .fns = ~.x/`Omicron BA.1*`,
                .names = "{.col}/BA.1*_ratio")) |> 
  ## BA.4/BA.5 ratio
  mutate(across(.cols = c(ends_with('BA.4*'), ends_with("XBB*")), 
                .fns = ~.x/`Omicron BA.5*`,
                .names = "{.col}/BA.5*_ratio")) |> 
  select(name_states, ends_with("_ratio"))|> 
  pivot_longer(cols = ends_with("_ratio"),
               names_to = "ratios",
               values_to = "value") |> 
  mutate(ratios = str_replace(ratios, "Omicron ","")) |> 
  mutate(ratios = str_replace(ratios, "_ratio",""))

days_to_peak_by_frequency_ratios <- days_to_peak_by_frequency_ratios |>
  mutate(value = round(value, 2)) |> 
  mutate(ratios2 = tidytext::reorder_within(ratios, -value, within = name_states))
  
days_to_peak_by_frequency_ratios |> 
  # filter(name_states == "New York") |>
  ggplot(aes(y = ratios2,
             x = value,
             fill = ratios))+
  geom_vline(xintercept = 1, col = "gray5")+
  geom_col()+
  geom_text(data = subset(days_to_peak_by_frequency_ratios, 
                          value < 1),
            aes(label = value),
            color = "gray50",
            size = 2.5,
            hjust = -0.5)+
  lims(x = c(NA,10))+
  scale_fill_discrete(name = "VOCs")+
  facet_geo(name_states~., scales = "free_y")+
  labs(x = "Ratio Value",
       title = "Ratio between n. of days \n from first emergence to peak")+
  theme_minimal()+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90),
        legend.position = "bottom")+
  tidytext::scale_y_reordered()

days_to_peak_by_frequency |> 
  # filter(name_states == "New York") |>
  mutate(variant2 = tidytext::reorder_within(variant, -day, within = name_states)) |> 
  ggplot(aes(y = variant2, 
             x = day, 
             fill = variant))+
  geom_col()+
  geom_text(aes(label = day),
            color = "gray90",
            size = 2.5,
            hjust = 1)+
  # lims(x = c(NA,400))+
  scale_fill_discrete(name = "VOCs")+
  facet_geo(name_states~., scales = "free_y")+
  labs(x = "Days to reach max. incidence", 
       title = "Days from first emergence to peak")+
  theme_minimal()+
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = "bottom")+
  tidytext::scale_y_reordered()

