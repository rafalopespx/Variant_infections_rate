## Cleaning the ambient
rm(list = ls())
gc()

## Loading Libraries
packs = c("tidyverse", "vroom", "patchwork", "geofacet", "tigris", "ggforce", "ggthemes", "MetBrewer", "ggrepel", "gghighlight", "gtsummary")
lapply(packs,require, character.only = TRUE)

# Loading functions
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

## States Abbreviation
states_abb <- vroom("Data/state_abbreviation.tsv")

## Loading the frequencies of variants
variant_count<-vroom("Data/variant_counts_us.csv.xz") |>
  rename(days = epiweek,
         variant = voc_cdc) |>
  filter(!variant %in% c("Alpha*", "Beta*", "Gamma*", "Delta*", "Omicron BA.3*", "Recombinant", "Other")) |>
  mutate(variant = droplevels(factor(variant))) |>
  mutate(variant = factor(variant,
                          levels = c("Omicron BA.1*", "Omicron BA.2*", "Omicron BA.4*", "Omicron BA.5*",
                                     "Omicron XBB*")))

# ## Joining Rt with variant counts
# estimates_rt_incidence<-rt_estimates |>
#   left_join(variant_count) |>
#   mutate(variant = droplevels(factor(variant))) |> 
#   left_join(pop_states, by = c("name_states" = "state")) |> 
#   ## Renaming infections
#   rename(infections = I) |> 
#   ## Creating incidence per 100k
#   mutate(incidence = (infections/pop)*1e5, 
#          percentual_incidence = round(incidence/pop, 2))
# 
# vroom_write(x = estimates_rt_incidence, 
#             file = "Data/state_full_data.csv.xz")

## Loading estimates_rt_incidence dataset
estimates_rt_incidence <- vroom("Data/state_full_data.csv.xz")

## Peak of infections, by state and variants
peak_infections_by_variant_table <- estimates_rt_incidence |> 
  reframe(sum_infections = sum(infections, na.rm = T),
          .by = c(days, variant)) %>%
  group_by(variant) |> 
  summarise(peak_infections = max(sum_infections, na.rm = T)) %>%
  # add_row(variant = "Peak", summarise(., across(where(is.numeric), max))) |>
  tbl_summary(by = variant, 
              # include = peak_infections,
              type = list(everything() ~ "continuous"),
              statistic = all_continuous() ~ "{mean}",
              label = list(peak_infections ~ "")) |> 
  modify_header(list(label ~ "", all_stat_cols() ~ "**{level}**")) |> 
  modify_footnote(all_stat_cols() ~ NA) |> 
  modify_spanning_header(c("stat_1", "stat_2", "stat_3", "stat_4", "stat_5") ~ "**Variant Categories**") |> 
  modify_caption("**Supplementary Table 2. Peaks of infections a day \n by variant categories**") |> 
  bold_labels()
peak_infections_by_variant_table

peak_infections_by_state_table <- estimates_rt_incidence |> 
  reframe(sum_infections = sum(infections, na.rm = T),
          .by = c(days, name_states, variant)) %>%
  group_by(name_states, variant) |> 
  summarise(peak_infections = max(sum_infections, na.rm = T)) |> 
  # filter(name_states %in% c("California", "New York", "Texas", "Florida")) |> 
  pivot_wider(names_from = name_states,
              values_from = peak_infections) %>%
  # add_row(variant = "Peak", summarise(., across(where(is.numeric), max))) |>
  tbl_summary(by = variant, 
              type = list(everything() ~ "continuous"),
              statistic = list(all_continuous() ~ "{mean}")) |> 
  modify_header(label ~ "**State**", all_stat_cols() ~ "**{level}**") |> 
  modify_footnote(all_stat_cols() ~ NA) |> 
  modify_spanning_header(c("stat_1", "stat_2", "stat_3", "stat_4", "stat_5") ~ "**Variant Categories**") |> 
  modify_caption("**Supplementary Table 2. Peaks of infections a day \n by variant categories**") |> 
  bold_labels()
peak_infections_by_state_table

peak_infections_table <- tbl_stack(list(peak_infections_by_variant_table, 
                                        peak_infections_by_state_table), 
                                   group_header = c("USA", "State")) |> 
  bold_labels()
peak_infections_table

## Saving the table
peak_infections_table |> 
  as_gt() |> 
  gt::gtsave(filename = "Output/Tables/peak_infections_table.png")

peak_infections_table |> 
  as_gt() |> 
  gt::gtsave(filename = "Output/Tables/peak_infections_table.docx")

## Total of infections, by state and variant
total_infections_by_variant_table <- estimates_rt_incidence |> 
  group_by(variant) |> 
  summarise(total_infections = sum(infections, na.rm = T)) %>%
  add_row(variant = "Overall", summarise(., across(where(is.numeric), sum))) |> 
  tbl_summary(by = variant, 
              type = list(everything() ~ "continuous"),
              statistic = all_continuous() ~ "{mean}",
              label = list(total_infections ~ "")) |> 
  modify_header(list(label ~ "", all_stat_cols() ~ "**{level}**")) |> 
  modify_footnote(all_stat_cols() ~ NA) |> 
  modify_spanning_header(c("stat_1", "stat_2", "stat_3", "stat_4", "stat_5") ~ "**Variant Categories**") |> 
  modify_caption("**Supplementary Table 3. Totals of infections \n by variant categories**") |> 
  bold_labels()
total_infections_by_variant_table

total_infections_by_state_table <- estimates_rt_incidence |> 
  group_by(name_states, variant) |> 
  summarise(total_infections = sum(infections, na.rm = T)) |> 
  pivot_wider(names_from = name_states,
              values_from = total_infections) %>%
  add_row(variant = "Total", summarise(., across(where(is.numeric), sum))) |>
  tbl_summary(by = variant, 
              type = list(everything() ~ "continuous"),
              statistic = list(all_continuous() ~ "{mean}")) |> 
  modify_header(label ~ "**State**", all_stat_cols() ~ "**{level}**") |> 
  modify_footnote(all_stat_cols() ~ NA) |> 
  modify_spanning_header(c("stat_1", "stat_2", "stat_3", "stat_4", "stat_5") ~ "**Variant Categories**") |> 
  modify_caption("**Supplementary Table 3. Totals of infections \n by variant categories**") |> 
  bold_labels()
total_infections_by_state_table

total_infections_table <- tbl_stack(list(total_infections_by_variant_table, 
                                         total_infections_by_state_table), 
                                    group_header = c("USA", "State")) |> 
  bold_labels()
total_infections_table

## Saving the table
total_infections_table |> 
  as_gt() |> 
  gt::gtsave(filename = "Output/Tables/total_infections_table.png")

total_infections_table |> 
  as_gt() |> 
  gt::gtsave(filename = "Output/Tables/total_infections_table.docx")

## Attack Rate tables
total_attack_rate_us_table <- estimates_rt_incidence |> 
  group_by(name_states, variant) |>
  summarise(total_infections = (sum(infections, na.rm = T)),
            pop = first(pop)) |> 
  mutate(attack_rate = (total_infections/pop)*1e2) |> 
  dplyr::select(variant, attack_rate) |> 
  group_by(variant) |> 
  summarise(mean_ar = mean(attack_rate, na.rm = T),
            max_ar = max(attack_rate, na.rm = T),
            min_ar = min(attack_rate, na.rm = T)) |> 
  pivot_longer(cols = c(mean_ar, max_ar, min_ar),
               names_to = "interval", 
               values_to = "values") |> 
  dplyr::select(variant, values) |> 
  tbl_summary(by = variant, 
              type = list(everything() ~ "continuous"), 
              statistic = list(all_continuous() ~ "{mean}% ({max}%, {min}%)"),
              label = list(values ~ ""),
              digits = list(all_continuous() ~ 1))   |> 
  add_overall() |> 
  modify_header(label ~ "", 
                all_stat_cols() ~ "**{level}**") |> 
  modify_footnote(all_stat_cols() ~ NA) |> 
  modify_spanning_header(all_stat_cols() ~ "**Variant Categories**") |> 
  modify_caption("**Supplementary Table 5. Attack Rates per state \n by variant categories**") |> 
  bold_labels()

total_attack_rate_us_table

total_attack_rate_by_state_table <- estimates_rt_incidence |> 
  group_by(name_states, variant) |>
  summarise(total_infections = (sum(infections, na.rm = T)),
            pop = first(pop)) |> 
  mutate(attack_rate = (total_infections/pop)*1e2) |> 
  select(name_states, variant, attack_rate) |> 
  pivot_wider(names_from = name_states,
              values_from = attack_rate) %>%
  tbl_summary(by = variant, 
              type = list(everything() ~ "continuous"), 
              statistic = list(all_continuous() ~ "{mean}%"), 
              digits = list(all_continuous() ~ 1)) |> 
  # add_overall() |> 
  modify_header(label ~ "**State**", 
                all_stat_cols() ~ "**{level}**") |> 
  modify_footnote(all_stat_cols() ~ NA) |> 
  modify_spanning_header(all_stat_cols() ~ "**Variants**") |> 
  modify_caption("**Supplementary Table 5. Attack Rates per state \n by variant categories**") |> 
  bold_labels()

total_attack_rate_by_state_table

total_attack_rate_table <- tbl_stack(list(total_attack_rate_us_table, 
                                         total_attack_rate_by_state_table), 
                                    group_header = c("USA \n mean (max, min)", "States")) |> 
  bold_labels()

total_attack_rate_table

## Saving the table
total_attack_rate_table |> 
  as_gt() |> 
  gt::gtsave(filename = "Output/Tables/attack_rates_table.png")

total_attack_rate_table |> 
  as_gt() |> 
  gt::gtsave(filename = "Output/Tables/attack_rates_table.docx")

## median, max and min for Rt estimates
rt_by_variant_table <- estimates_rt_incidence |> 
  rename(median = Rt) |> 
  select(variant, median, upper, lower) |> 
  drop_na() |> 
  tbl_summary(by = variant, 
              type = list(everything() ~ "continuous"), 
              statistic = list(all_continuous() ~ "{max} ({p25}, {p75})"), 
              digits = list(all_continuous() ~ 2)) |> 
  modify_header(label ~ "max (p25, p75)", 
                all_stat_cols() ~ "**{level}**") |> 
  modify_footnote(all_stat_cols() ~ NA) |> 
  modify_spanning_header(all_stat_cols() ~ "**Variant Categories**") |> 
  modify_caption("**Supplementary Table 6. Median, Upper, and Lower interval of Rt by variant categories**") |> 
  bold_labels()
rt_by_variant_table

## Saving the table
rt_by_variant_table |>
  as_gt() |> 
  gt::gtsave(filename = "Output/Tables/rt_by_variant.png")

rt_by_variant_table |>
  as_gt() |> 
  gt::gtsave(filename = "Output/Tables/rt_by_variant.docx")

##US Totals table

us_total_table <- tbl_stack(list(total_attack_rate_us_table,
                                 peak_infections_by_variant_table,
                                 total_infections_by_variant_table),
                            group_header = list("Attack Rate \n mean (max, min)", 
                                             "Peak of Infections", 
                                             "Total of Infections")) |> 
  modify_caption("**Table 1. Variant-specific Attack Rate, Peak and Total of infections**")
us_total_table

us_total_table |> 
  as_gt() |> 
  gt::gtsave(filename = "Output/Tables/usa_ar_peak_total.png")

us_total_table |> 
  as_gt() |> 
  gt::gtsave(filename = "Output/Tables/usa_ar_peak_total.docx")

#
  