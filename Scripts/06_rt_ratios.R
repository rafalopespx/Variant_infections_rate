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
  ## Creating incidence per 100k
  mutate(incidence = (infections/pop)*1e5)

## Rt ratios
rt_ratios<-estimates_rt_incidence |> 
  rename(median = Rt,
         upper = Rt_upper,
         lower = Rt_lower) |> 
  pivot_wider(id_cols = c("days", "name_states"),
              names_from = variant,
              values_from = c(median, upper, lower))

## Ratios calculation
median_ratios <- rt_ratios |>
  dplyr::select(days, name_states, starts_with('median')) |> 
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
  dplyr::select(days, name_states, ends_with("_ratio")) |> 
  pivot_longer(cols = ends_with("_ratio"),
               names_to = "ratios",
               values_to = "median") |> 
  mutate(ratios = str_replace(ratios, "median_Omicron ","")) |> 
  mutate(ratios = str_replace(ratios, "_ratio",""))

upper_ratios <- rt_ratios |>
  dplyr::select(days, name_states, starts_with('upper')) |> 
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
  dplyr::select(days, name_states, ends_with("_ratio")) |> 
  pivot_longer(cols = ends_with("_ratio"),
               names_to = "ratios",
               values_to = "upper") |> 
  mutate(ratios = str_replace(ratios, "upper_Omicron ","")) |> 
  mutate(ratios = str_replace(ratios, "_ratio",""))

lower_ratios <- rt_ratios |>
  dplyr::select(days, name_states, starts_with('lower')) |> 
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
  dplyr::select(days, name_states, ends_with("_ratio")) |> 
  pivot_longer(cols = ends_with("_ratio"),
               names_to = "ratios",
               values_to = "lower") |> 
  mutate(ratios = str_replace(ratios, "lower_Omicron ","")) |> 
  mutate(ratios = str_replace(ratios, "_ratio",""))

## Rt ratios
rt_ratios<-median_ratios |> 
  left_join(upper_ratios) |> 
  left_join(lower_ratios)

vroom_write(x = rt_ratios, file = "Output/Tables/rt_ratios.csv.xz")

#