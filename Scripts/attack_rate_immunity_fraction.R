## Cleaning the ambient
rm(list = ls())
gc()

## Loading Libraries
packs = c("tidyverse", "vroom", "MASS", "ggeffects", "marginaleffects", "broom.helpers", "sf", "tidycensus", "tidyverse", "ggstats", "geofacet", "biscale", "cowplot", "pals", "MetBrewer")
lapply(packs,require, character.only = TRUE)

# Loading functions
source("Scripts/Functions/functions.R")
source("Scripts/get_svi.R")

## Loading data sources
states_fulldata<-vroom("Data/state_full_data.csv.xz")

data_national <- states_fulldata |>
  filter(days <= "2023-04-01",
         infections > 1000) |>
  reframe(infections = sum(infections, na.rm = T),
          .by = c(days))

data_national_variants <- states_fulldata |>
  filter(days <= "2023-04-01",
         infections > 1000) |>
  reframe(infections = sum(infections, na.rm = T),
          .by = c(days, variant))

days_vertical <- as.Date(c("2022-05-01",
                           "2022-09-01",
                           "2022-12-01",
                           "2023-04-01",
                           "2023-04-01"))

days_labels <- c("Omicron BA.1*", 
                 "Omicron BA.2*", 
                 "Omicron BA.4*",
                 "Omicron BA.5*",
                 "Omicron XBB*")

variants <- unique(data_national_variants$variant)

fill_colors <- c(met.brewer(palette_name = "Archambault",
                            type = "discrete",
                            n = 5))

plt_list <- list()

for (i in 1:5) {
  plt_list[[i]] <- data_national_variants |>
    ggplot(aes(x = days, y = infections, group = variant))+
    geom_line(aes(color = variant))+
    geom_area(data = data_national_variants |>
                filter(variant == variants[i]),
              aes(x = days, y = infections),
              fill = fill_colors[i],
              alpha = 0.50)+
    theme_minimal()+
    labs(x = "Date", y = "Infections per day")+
    scale_x_date(date_breaks = "2 months", date_labels = "%b %y")+
    scale_y_continuous(labels = scales::label_scientific())+
    scale_color_manual(values = c(met.brewer(palette_name = "Archambault",
                                             type = "discrete",
                                             n = 5)),
                       name = "",
                       aesthetics = "color")+
    annotate(geom = "vline",
             x = days_vertical[i],
             # y= 0,
             ymax = 2e6,
             xintercept = days_vertical[i],
             linetype = "dashed")+
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank())
}

plt_list[[1]] <- plt_list[[1]]+
  annotate(geom = "text",
           x = days_vertical[1]+75,
           y = 3e6, 
           label = paste0("End of ",days_labels[1], "\n wave"),
           angle = 0)
plt_list[[1]]

##Saving plots
ggsave(plot = plt_list[[1]],
       filename = paste0("Output/Plots/area_under_the_curve_",
                         str_remove(variants[1],
                                    pattern = "Omicron "),
                         ".png"),
       width = 11,
       height = 9,
       dpi = 100)

plt_list[[2]] <- plt_list[[2]]+
  annotate(geom = "text",
           x = days_vertical[2]+75,
           y = 3e6, 
           label = paste0("End of ",days_labels[2], "\n wave"),
           angle = 0)
plt_list[[2]]

##Saving plots
ggsave(plot = plt_list[[2]],
       filename = paste0("Output/Plots/area_under_the_curve_",
                         str_remove(variants[2],
                                    pattern = "Omicron "),
                         ".png"),
       width = 11,
       height = 9,
       dpi = 100)

plt_list[[3]] <- plt_list[[3]]+
  annotate(geom = "text",
           x = days_vertical[3]-75,
           y = 3e6, 
           label = paste0("End of ",days_labels[3], "\n wave"),
           angle = 0)
plt_list[[3]]

##Saving plots
ggsave(plot = plt_list[[3]],
       filename = paste0("Output/Plots/area_under_the_curve_",
                         str_remove(variants[3],
                                    pattern = "Omicron "),
                         ".png"),
       width = 11,
       height = 9,
       dpi = 100)

plt_list[[4]] <- plt_list[[4]]+
  annotate(geom = "text",
           x = days_vertical[4]-75,
           y = 3e6, 
           label = paste0("End of ",days_labels[4], "\n wave"),
           angle = 0)
plt_list[[4]]

##Saving plots
ggsave(plot = plt_list[[4]],
       filename = paste0("Output/Plots/area_under_the_curve_",
                         str_remove(variants[4],
                                    pattern = "Omicron "),
                         ".png"),
       width = 11,
       height = 9,
       dpi = 100)

plt_list[[5]] <- plt_list[[5]]+
  annotate(geom = "text",
           x = days_vertical[5]-75,
           y = 3e6, 
           label = paste0("End of ",days_labels[5], "\n wave"),
           angle = 0)
plt_list[[5]]

##Saving plots
ggsave(plot = plt_list[[5]],
       filename = paste0("Output/Plots/area_under_the_curve_",
                         str_remove(variants[5],
                                    pattern = "Omicron "),
                         ".png"),
       width = 11,
       height = 9,
       dpi = 100)


states_attack_rates <- states_fulldata |> 
  # ## filtering to days before the peak of infections of each variant, all states
  # filter(days <= days[which.max(infections)],
  #        .by = c("name_states", "variant")) |>
  reframe(total_infections = sum(infections, na.rm = T),
          total_incidence = sum(incidence, na.rm = T),
          pop = first(pop),
          .by = c("name_states", "variant")) |> 
  mutate(attack_rate = round((total_infections/pop)*100, 2))

## Loading Immunity estimates

## Pre Omicron model
immunity_fraction_pos_omicron<-vroom("Data/immunity-weekly-state.csv.xz") |> 
  rename(name_states = location,
         days = date) 

## Pos Omicron model
immunity_fraction_pre_omicron<-vroom("Data/immunity-daily-state.csv.xz") |> 
  rename(name_states = location,
         days = date)

immunity_fraction_longer_pos_omicron <- immunity_fraction_pos_omicron|> 
  filter(interval == "median") |> 
  pivot_longer(cols = exposed_infection:protectedSev,
               names_to = "type_protected",
               values_to = "values") |> 
  filter(startsWith(type_protected, "protectedInf")) |> 
  mutate(epiweek = end.of.epiweek(days)) 
# |> 
#   filter(values > 0 ) |> 
#   pivot_wider(names_from = "interval", 
#               values_from = "values")

immunity_fraction_longer_pre_omicron <- immunity_fraction_pre_omicron|> 
  filter(interval == "median") |> 
  pivot_longer(cols = exposed_infection:protectedSev,
               names_to = "type_protected",
               values_to = "values") |> 
  filter(startsWith(type_protected, "protectedInf")) |> 
  mutate(epiweek = end.of.epiweek(days))
# |> 
#   filter(values > 0 ) |> 
#   pivot_wider(names_from = "interval", 
#               values_from = "values")

# ## Reframing pre Omicron estimates into weekly basis
# immunity_fraction_longer_pre_omicron <- immunity_fraction_longer_pre_omicron |> 
#   reframe(low = mean(low, na.rm = T),
#           median = mean(median, na.rm = T),
#           high = mean(high, na.rm = T), 
#           .by = c("name_states", "epiweek", "type_protected"))
# 
# immunity_fraction_longer_pos_omicron <- immunity_fraction_longer_pos_omicron |> 
#   reframe(low = mean(low, na.rm = T),
#           median = mean(median, na.rm = T),
#           high = mean(high, na.rm = T), 
#           .by = c("name_states", "epiweek", "type_protected"))

immunity_fraction_longer_pos_omicron |> 
  ggplot(aes(x = epiweek, y = values,
             fill = type_protected))+
  geom_line()+
  # geom_ribbon(aes(ymin = low, ymax = high), alpha = .5)+
  geofacet::facet_geo(name_states~.)+
  theme(axis.text.x = element_text(angle = 90))+
  scale_x_date(name = "Date", date_breaks = "2 months", date_labels = "%y - %b")

immunity_fraction_longer_pre_omicron |> 
  ggplot(aes(x = epiweek, y = values,
             fill = type_protected))+
  geom_line()+
  # geom_ribbon(aes(ymin = low, ymax = high), alpha = .5)+
  geofacet::facet_geo(name_states~.)+
  theme(axis.text.x = element_text(angle = 90))+
  scale_x_date(name = "Date", date_breaks = "2 months", date_labels = "%y - %b")

immunity_fraction_longer <- full_join(immunity_fraction_longer_pre_omicron, 
                                      immunity_fraction_longer_pos_omicron) 

immunity_fraction_longer |> 
  filter(name_states != "Puerto Rico", 
         # name_states != "District of Columbia", 
         type_protected != "protectedInf_omicron") |> 
  ggplot(aes(x = epiweek, y = values,
             col = type_protected))+
  geom_line()+
  geofacet::facet_geo(name_states~.)+
  theme(axis.text.x = element_text(angle = 90))+
  scale_x_date(name = "Date", date_breaks = "4 months", date_labels = "%y - %b")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90), 
        legend.position = "bottom")+
  scale_color_met_d(name = "", 
                    labels = c("All-kind protection", 
                               "Hybrid protection", 
                               "Protection by infection", 
                               "Protection against Omicron", 
                               "Protection by vaccination"),
                    palette_name = "Juarez")

## Joining data
states_attack_rate_immune<-left_join(states_fulldata, 
                                     immunity_fraction_longer) |> 
  mutate(epiweek = end.of.epiweek(days)) 

states_sum <- states_attack_rate_immune |> 
  # ## filtering to days before the peak of infections of each variant, all states
  # filter(days <= days[which.max(infections)],
  #        .by = c("name_states", "variant")) |>
  # rename(values = median) |> 
  reframe(values = round(sum(values, na.rm = T), 2), 
          .by = c("name_states", "variant", "type_protected")) |>
  filter(!is.na(type_protected)) |> 
  right_join(states_attack_rates) |> 
  mutate(attack_rate = round((total_infections/pop)*100, 2))

# vroom_write(x = states_sum,
#             file = "Data/states_immune_attack_rate.csv.xz")

states_abb<-vroom("Data/states_abbrev_region_division.csv") |> 
  rename(name_states = `State`)

states_map<-tigris::states(cb = TRUE) |> 
  tigris::shift_geometry() |> 
  rename(name_states = NAME) |> 
  filter(name_states %in% states_abb$name_states)

states_immune_map<-left_join(states_map, states_sum)

type_protected<-"protectedInf"

states_immune_map <- states_immune_map |> 
  # filter(type_protected == type_protected) |> 
  filter(name_states != "District of Columbia") |> 
  bi_class(x = values, 
           y = attack_rate, 
           style = "quantile",
           dim = 4, 
           keep_factors = T)

unique(states_immune_map$bi_class)

source('Scripts/custom_pal.R')

protected_types<-unique(states_immune_map |> 
                          filter(!is.na(type_protected)) |> 
                          pull(var = type_protected))
protected_labels<-c("Infection protected estimates",
                    "Vaccination protected estimates",
                    "Hybrid protected estimates",
                    "All-kind protected estimates")

plot_list<- map <- list()

for(i in 1:length(protected_types)){
  for (j in 1:length(variants)) {
    ## Data filtering
    data<-states_immune_map |> 
      filter(type_protected == protected_types[i],
             variant == variants[j])
    
    plt_map<- data |> 
      ggplot()+
      geom_sf(aes(fill = bi_class),
              color = "grey45",
              size = 0.01, 
              show.legend = FALSE)+
      bi_scale_fill(pal = custom_pal4_1, 
                    dim = 4) +
      labs(subtitle = variants[j]
      #   title = "Attack Rate vs. Spike-exposed during",
      #   subtitle = paste0(variants[j], 
      #                     " wave via ", 
      #                     str_remove(protected_labels[i],
      #                                " protected estimates"))
      ) +
      theme_void()
    
    map[[j]] <- (plt_map / plt_list[[j]])+
      plot_layout(nrow = 2, 
                  widths = c(1,1), 
                  heights = c(4,1))
    
  }

## Labels
labels1 <- bi_class_breaks(states_immune_map,
                           x = values,
                           y = attack_rate,
                           style = "quantile",
                           dim = 4,
                           dig_lab = 2,
                           split = FALSE)
## Legend
legend <- bi_legend(pal = custom_pal4_1,
                    breaks = labels1,
                    dim = 4,
                    xlab = "Higher (%) Effectively Protected",
                    ylab = "Higher (%) Attack Rate",
                    size = 12)+
  theme_minimal()+
  # bi_theme()+
  theme(axis.text.x = element_text(angle = 90))

# combine map with legend
plot_list[[i]] <- (map[[1]] | map[[2]] | map[[3]])/(map[[4]] | map[[5]] | legend)

# ## Saving the plots
ggsave(filename = paste0("Output/Plots/bivariate_plots/plot_",protected_types[i], ".png"),
       plot = plot_list[[i]],
       width = 11,
       height = 9,
       dpi = 100)
# finalPlot
}

plot_list[[4]]
plot_list[[3]]
plot_list[[2]]
plot_list[[1]]


## Attack Rates scatterplots

states_omicron_ba1_ar <- states_attack_rates|> 
  filter(variant == "Omicron BA.1*") |> 
  dplyr::select(name_states, attack_rate) |> 
  rename("Omicron BA.1* Attack Rate" = attack_rate)

states_omicron_ba2_ar <- states_attack_rates|> 
  filter(variant == "Omicron BA.2*") |> 
  dplyr::select(name_states, attack_rate) |> 
  rename("Omicron BA.2* Attack Rate" = attack_rate)

states_omicron_ba4_ar <- states_attack_rates|> 
  filter(variant == "Omicron BA.4*") |> 
  dplyr::select(name_states, attack_rate) |> 
  rename("Omicron BA.4* Attack Rate" = attack_rate)

states_omicron_ba5_ar <- states_attack_rates|> 
  filter(variant == "Omicron BA.5*") |> 
  dplyr::select(name_states, attack_rate) |> 
  rename("Omicron BA.5* Attack Rate" = attack_rate)

states_omicron_xbb_ar <- states_attack_rates|> 
  filter(variant == "Omicron XBB*") |> 
  dplyr::select(name_states, attack_rate) |> 
  rename("Omicron XBB* Attack Rate" = attack_rate)

states_joined <- full_join(states_omicron_ba1_ar, 
                           states_omicron_ba2_ar) |>
  left_join(states_omicron_ba4_ar) |> 
  left_join(states_omicron_ba5_ar) |> 
  left_join(states_omicron_xbb_ar) |> 
  left_join(states_abb)

ba1_ba2 <- states_joined |> 
  ggplot(aes(x = `Omicron BA.1* Attack Rate`, 
             y = `Omicron BA.2* Attack Rate`,
             label = `State Code`, 
             col = Region))+
  geom_point()+
  geom_text_repel(show.legend = F)+
  theme_bw()+
  theme(legend.position = "bottom")

ba2_ba4 <- states_joined |> 
  ggplot(aes(x = `Omicron BA.2* Attack Rate`, 
             y = `Omicron BA.4* Attack Rate`,
             label = `State Code`, 
             col = Region))+
  geom_point()+
  geom_text_repel(show.legend = F)+
  theme_bw()+
  theme(legend.position = "bottom")

ba2_ba5 <- states_joined |> 
  ggplot(aes(x = `Omicron BA.2* Attack Rate`, 
             y = `Omicron BA.5* Attack Rate`,
             label = `State Code`, 
             col = Region))+
  geom_point()+
  geom_text_repel(show.legend = F)+
  theme_bw()+
  theme(legend.position = "bottom")

ba4_ba5 <- states_joined |> 
  ggplot(aes(x = `Omicron BA.4* Attack Rate`, 
             y = `Omicron BA.5* Attack Rate`,
             label = `State Code`, 
             col = Region))+
  geom_point()+
  geom_text_repel(show.legend = F)+
  theme_bw()+
  theme(legend.position = "bottom")

ba5_xbb <- states_joined |> 
  ggplot(aes(x = `Omicron BA.5* Attack Rate`, 
             y = `Omicron XBB* Attack Rate`,
             label = `State Code`, 
             col = Region))+
  geom_point()+
  geom_text_repel(show.legend = F)+
  theme_bw()+
  theme(legend.position = "bottom")

patchwork_ar <- (ba1_ba2 | ba2_ba4 | ba2_ba5 | ba4_ba5 | ba5_xbb)+
  plot_layout(guides = 'collect')&
  theme(legend.position = "bottom")
patchwork_ar
