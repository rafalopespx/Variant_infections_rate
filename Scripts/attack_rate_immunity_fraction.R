## Cleaning the ambient
rm(list = ls())
gc()

## Loading Libraries
packs = c("tidyverse", "vroom", "MASS", "ggeffects", "marginaleffects", "broom.helpers", "sf", "tidycensus", "tidyverse", "ggstats")
lapply(packs,require, character.only = TRUE)

# Loading functions
source("Scripts/Functions/functions.R")
source("Scripts/get_svi.R")

## Loading data sources
states_fulldata<-vroom("Data/state_full_data.csv.xz")

states_attack_rates<-vroom("Data/state_attack_rate_variants.csv.xz") 

immunity_fraction<-vroom("Data/immunity-weekly-state.csv.xz") |> 
  rename(name_states = location,
         days = date) |> 
  filter(interval == "median") 

immunity_fraction_longer <- immunity_fraction|> 
  pivot_longer(cols = exposed_infection:protectedSev,
               names_to = "type_protected",
               values_to = "values") |> 
  filter(startsWith(type_protected, "protectedInf"))

immunity_fraction_longer|> 
  ggplot(aes(x = days, y = values,
             col = type_protected))+
  geom_line()+
  geofacet::facet_geo(name_states~.)

## Joining data
states_attack_rate_immune<-left_join(states_fulldata, 
                                     immunity_fraction_longer) |> 
  mutate(epiweek = end.of.epiweek(days)) 

states_sum <- states_attack_rate_immune|> 
  reframe(values = round(sum(values, na.rm = T), 2), 
          .by = c("name_states", "variant", "type_protected")) |>
  filter(!is.na(type_protected)) |> 
  right_join(states_attack_rates) |> 
  mutate(attack_rate = round((total_infections/pop)*100, 2))

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
           dim = 6, 
           keep_factors = T)
unique(states_immune_map$bi_class)

source('Scripts/custom_pal.R')

protected_types<-unique(states_immune_map$type_protected)

plot_list<-list()

for(i in 1:length(protected_types)){
  ## Data filtering
  data<-states_immune_map |> 
    filter(type_protected == protected_types[i])
    
  map<- data |> 
    ggplot()+
    geom_sf(aes(fill = bi_class),
            color = "grey50",
            size = 0.1, 
            show.legend = FALSE)+
    bi_scale_fill(pal = custom_pal6, 
                  dim = 6) +
    # bi_scale_color(pal = custom_pal6, 
    #                dim = 6) +
    labs(
      title = "Attack Rate vs. Effectively Protected",
      subtitle = protected_types[i]
    ) +
    theme_minimal()+
    facet_wrap(variant~.)
  map
  
  ## Labels
  labels1 <- bi_class_breaks(data, 
                             x = values, 
                             y = attack_rate, 
                             style = "quantile", 
                             dim = 6, 
                             dig_lab = 2, 
                             split = FALSE)
  ## Legend
  legend <- bi_legend(pal = custom_pal6,
                      breaks = labels1,
                      dim = 6,
                      xlab = "Higher (%) Effectively Protected",
                      ylab = "Higher (%) Attack Rate",
                      size = 12)+
    theme(axis.text.x = element_text(angle = 90))
  
  # combine map with legend
  plot_list[[i]] <- ggdraw() +
    draw_plot(map, 0, 0, 1, 1) +
    draw_plot(legend, .70, .10, 0.35, 0.35)
  # finalPlot
}

plot_list[[1]]
plot_list[[2]]
plot_list[[3]]
plot_list[[4]]

# variants<-unique(states_immune_map$variant)
# 
# plot_list2<-list()
# 
# for(i in 1:length(variants)){
#   ## Data filtering
#   data<-states_immune_map |> 
#     filter(variant == variants[i])
#   
#   map<- data |> 
#     ggplot()+
#     geom_sf(aes(fill = bi_class),
#             color = "grey50",
#             size = 0.1, 
#             show.legend = FALSE)+
#     bi_scale_fill(pal = custom_pal6, 
#                   dim = 6) +
#     # bi_scale_color(pal = custom_pal6, 
#     #                dim = 6) +
#     labs(
#       title = "Attack Rate vs. Effectively Protected",
#       subtitle = variants[i]
#     ) +
#     theme_minimal()+
#     facet_wrap(type_protected~., nrow = 1)
#   map
#   
#   ## Labels
#   labels1 <- bi_class_breaks(data, 
#                              x = values, 
#                              y = attack_rate, 
#                              style = "quantile", 
#                              dim = 6, 
#                              dig_lab = 2, 
#                              split = FALSE)
#   ## Legend
#   legend <- bi_legend(pal = custom_pal6,
#                       breaks = labels1,
#                       dim = 6,
#                       xlab = "Higher (%) Effectively Protected",
#                       ylab = "Higher (%) Attack Rate",
#                       size = 12)+
#     theme(axis.text.x = element_text(angle = 90))
#   
#   # combine map with legend
#   plot_list2[[i]] <- ggdraw() +
#     draw_plot(map, 0, 0, 1, 1) +
#     draw_plot(legend, .1, .1, .35, 35)
#   # finalPlot
# }
# 
# plot_list2[[1]]
# plot_list2[[2]]
# plot_list2[[3]]
# plot_list2[[4]]
# 
# 
# 
# 
