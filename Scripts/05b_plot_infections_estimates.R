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

infections_variants_daily<-infections_variants_daily |> 
  mutate(variant_reduced = case_when(variant == "Omicron BA.2.75*" ~ "Omicron BA.2*", 
                                     variant %in% c("Omicron BQ.1*", "Omicron BJ.1*") ~ "Omicron BA.5*",
                                     variant %in% c("XBB.1*" ,"XBB.1.5*") ~ "XBB*", 
                                     variant == "Recombinant" ~ "Other", 
                                     TRUE ~ variant))

infections_variants_daily_reduced<-infections_variants_daily |> 
  group_by(name_states, days, variant_reduced) |> 
  summarise(I = sum(I, na.rm = T))

## Plot
plot_infections<-infections_variants_daily_reduced |> 
  filter(variant_reduced != "Other") |>
  ggplot(aes(x = days, y = I, 
             col = variant_reduced, fill = variant_reduced))+
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
plot_infections

plot_ct<-infections_variants_daily_reduced |> 
  filter(variant_reduced != "Other", 
         name_states == "Connecticut") |>
  ggplot(aes(x = days, y = I, 
             col = variant_reduced, fill = variant_reduced))+
  geom_line()+
  theme_minimal()+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90))+
  # facet_geo(name_states~., scales = "free_y")+
  scale_x_date(date_breaks = "4 months", 
               date_labels = "%b %y")+
  labs(x = "Date", 
       y = "Infections number per variant",
       title = "Connecticut")+
  scale_fill_brewer(aesthetics = c("color", "fill"), 
                    palette = "Spectral", 
                    name = "VOCs")
plot_ct
