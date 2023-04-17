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

infections_variants_daily_reduced<-infections_variants_daily |> 
  group_by(name_states, days, variant) |> 
  summarise(I = sum(I, na.rm = T))

## Plot
plot_infections<-infections_variants_daily_reduced |> 
  filter(variant != "Other") |>
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
plot_infections

## Saving infections plot
ggsave(filename = "Output/Plots/plt_infections_daily_states.png", 
       plot = plot_infections, 
       width = 11,
       height = 9, 
       dpi = 100)

plot_ct<-infections_variants_daily_reduced |> 
  filter(variant != "Other", 
         name_states == "Connecticut") |>
  mutate(variant = droplevels(factor(variant))) |> 
  ggplot(aes(x = days, y = I, 
             col = variant, fill = variant))+
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

## Saving CT plot as individual example
ggsave(filename = "Output/Plots/plt_infections_daily_CT.png", 
       plot = plot_ct, 
       width = 11, 
       height = 9, 
       dpi = 100)
