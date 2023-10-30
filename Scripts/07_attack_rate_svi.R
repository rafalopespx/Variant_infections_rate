## Cleaning the ambient
rm(list = ls())
gc()

## Loading Libraries
packs = c("tidyverse", "vroom", "MASS", "ggeffects", "marginaleffects", "broom.helpers", "sf", "tidycensus", "tidyverse", "ggstats", "geofacet", "biscale", "cowplot", "pals", "MetBrewer", "ggrepel", "patchwork")
lapply(packs,require, character.only = TRUE)

# Loading functions
source("Scripts/Functions/functions.R")
source("Scripts/get_svi.R")

## Loading data sources
states_fulldata<-vroom("Data/state_full_data.csv.xz")

states_abb<-vroom("Data/state_abbreviation.tsv")

states_attack_rates <- states_fulldata |> 
  reframe(total_infections = sum(infections, na.rm = T),
          total_incidence = sum(incidence, na.rm = T),
          pop = first(pop),
          .by = c("name_states", "variant")) |> 
  mutate(attack_rate = round((total_infections/pop)*100, 2))

## Attack Rates scatterplots
states_attack_rates_wider <- states_attack_rates |> 
  dplyr::select(name_states, variant, attack_rate, pop) |> 
  pivot_wider(names_from = "variant", 
              values_from = "attack_rate") |> 
  left_join(states_abb)

# Loading functions
source("Scripts/Functions/functions.R")
source("Scripts/Functions/get_svi.R")

# ## SVI variable vector of all states
# us <- unique(fips_codes$state)[1:51]
# acs_year <- 2020
# geo_unit <- "state"
# svi_df_raw <- map_df(us, function(x) {
#   get_svi(geo_unit, acs_year, x)
# })
# 
# ## States SVI
# states_svi <- svi_df_raw |> 
#   # group_by(ST) %>%
#   mutate(
#     EPL_POV150 = percent_rank(EP_POV150),
#     EPL_UNEMP = percent_rank(EP_UNEMP),
#     EPL_HBURD = percent_rank(EP_HBURD),
#     EPL_NOHSDP = percent_rank(EP_NOHSDP),
#     EPL_UNINSUR = percent_rank(EP_UNINSUR),
#     SPL_THEME1 = EPL_POV150 + EPL_UNEMP +  EPL_HBURD + EPL_NOHSDP + EPL_UNINSUR,
#     RPL_THEME1 = percent_rank(SPL_THEME1),
#     
#     EPL_AGE65 = percent_rank(EP_AGE65),
#     EPL_AGE17 = percent_rank(EP_AGE17),
#     EPL_DISABL = percent_rank(EP_DISABL),
#     EPL_SNGPNT = percent_rank(EP_SNGPNT),
#     EPL_LIMENG = percent_rank(EP_LIMENG),
#     SPL_THEME2 = EPL_AGE65 + EPL_AGE17 + EPL_DISABL + EPL_SNGPNT + EPL_LIMENG,
#     RPL_THEME2 = percent_rank(SPL_THEME2),
#     
#     EPL_MINRTY = percent_rank(EP_MINRTY),
#     SPL_THEME3 = EPL_MINRTY,
#     RPL_THEME3 = percent_rank(SPL_THEME3),
#     
#     EPL_MUNIT = percent_rank(EP_MUNIT),
#     EPL_MOBILE = percent_rank(EP_MOBILE),
#     EPL_CROWD = percent_rank(EP_CROWD),
#     EPL_NOVEH = percent_rank(EP_NOVEH),
#     EPL_GROUPQ = percent_rank(EP_GROUPQ),
#     SPL_THEME4 = EPL_MUNIT + EPL_MOBILE + EPL_CROWD + EPL_NOVEH + EPL_GROUPQ,
#     RPL_THEME4 = percent_rank(SPL_THEME4),
#     
#     SPL_THEMES = SPL_THEME1 + SPL_THEME2 + SPL_THEME3 + SPL_THEME4,
#     RPL_THEMES = percent_rank(SPL_THEMES)) |> 
#   rename(SVI = SPL_THEMES, SVI_rank = RPL_THEMES) |> 
#   dplyr::select(GEOID, name_states, geometry, SVI, SVI_rank)
# 
# ## Saving the states_svi as a .shp file
# st_write(states_svi, "Data/states_svi.shp")

## Read directly from the saved .shp file
states_svi <- st_read("Data/states_svi.shp") |> 
  rename(name_states = nm_stts,
         SVI_rank = SVI_rnk) |> 
  tigris::shift_geometry() |> 
  left_join(states_abb)

SVI_states_map <- states_svi |> 
  ggplot()+
  geom_sf(aes(fill = SVI), 
          show.legend = F, 
          lwd = 0.1)+
  theme_void()+
  scale_fill_met_c(palette = "Tam", 
                   # direction = -1,
                   # show.limits = T,
                   guide = guide_colorbar(title = "SVI",
                                          title.position = "top",
                                          title.hjust = 0.5,
                                          nbin = 10)
  )+
  theme(legend.position = "bottom", 
        legend.key.width = grid::unit(2.5, "cm"))
SVI_states_map

ggsave(filename = "Output/Plots/states_svi_map.png",
       plot = SVI_states_map,
       width = 16, 
       height = 9, 
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures//states_svi_map.pdf",
       plot = SVI_states_map,
       width = 16, 
       height = 9, 
       dpi = 100)

states_ar_svi <- states_attack_rates_wider |> 
  left_join(states_svi) |> 
  mutate(across((starts_with("Omicron ")), 
                .fns = ~.x/100)) |> 
  mutate(Region = factor(Region, levels = c("West", "Midwest", "Northeast", "South")))

states_ar_svi_longer <- states_ar_svi |>
  pivot_longer(cols = c(`Omicron BA.1*`:`Omicron XBB*`),
               names_to = "variant",
               values_to = "attack_rate")

vroom_write(x = states_ar_svi_longer, 
            file = "Data/states_ar_svi_longer.csv.xz")

label_ar_svi <- states_ar_svi_longer |> 
  filter(SVI == max(SVI))

ar_svi <-
  ggplot(data = states_ar_svi_longer,
         aes(x = SVI, 
             y = attack_rate,
             col = variant, size = pop, 
             label = `State Code`))+
  geom_point(alpha = 0.5)+
  geom_smooth(aes(fill = variant, weight = pop),
              method = "lm",
              show.legend = F,
              alpha = 0.25)+
  ggrepel::geom_text_repel(data= label_ar_svi, 
                           aes(label = variant, 
                               x = SVI + 0.10),
                           angle = -90,
                           color = "black",
                           size = 3)+
  ggpubr::stat_cor(method = "pearson",
                   geom = "label",
                   cor.coef.name = "R",
                   aes(label = paste(after_stat(..r.label..),
                                     after_stat(..p.label..),
                                     sep = "~`,`~")), 
                   label.x.npc = 0,
                   label.y.npc = 1,
                   r.digits = 2,
                   r.accuracy = 0.01, 
                   show.legend = F)+
  theme_minimal()+
  scale_color_met_d(palette = "Archambault")+
  scale_fill_met_d(palette = "Archambault")+
  scale_y_continuous(labels = scales::percent,
                     limits = c(0,0.80))+
  scale_x_continuous(breaks = seq(4.5, 11.5, 1))+
  scale_size_binned(labels = scales::comma, 
                    range = c(1,14),
                    # nice.breaks = T,
                    # n.breaks = 7,
                    guide = guide_legend(title = "Population size", 
                                         title.position = "top",
                                         title.hjust = 0.5, 
                                         nrow = 1))+
  labs(y = "Attack Rate \n Percent of state-population ever infected", 
       x = "Social Vulnerability Index \n SVI")+
  theme(legend.position = "bottom", 
        legend.key.width = unit(2.5, "cm"),
        # axis.text.x = element_blank(),
        strip.placement = "inside", 
        strip.text.x.top = element_text(hjust = 0))+
  guides(color = "none")
ar_svi

ar_svi <- ar_svi + 
  inset_element(p = SVI_states_map+
                  theme(legend.position = "none"), 
                left = 0.65, 
                bottom = 0.75, 
                right = 1, top = 1)
ar_svi

ba1_ba2 <- states_ar_svi |> 
  ggplot(aes(x = `Omicron BA.1*`, 
             y = `Omicron BA.2*`,
             # label = `State Code`, 
             col = SVI, 
             size = pop))+
  geom_point(alpha = 0.5)+
  geom_smooth(aes(weight = pop),
              method = "lm",
              show.legend = F,
              alpha = 0.25)+
  geom_text(aes(label = `State Code`), 
            color = "white", 
            size = 2.5)+
  ggpubr::stat_cor(method = "pearson",  
                   cor.coef.name = "R", 
                   geom = "label",
                   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
                   r.digits = 2, 
                   r.accuracy = 0.01, size = 5, label.x = 0.50)+
  theme_minimal()+
  scale_x_continuous(labels = scales::percent)+
  scale_y_continuous(labels = scales::percent)+
  scale_color_met_c(palette = "Tam", 
                    # direction = -1, 
                    # show.limits = T,
                    guide = guide_colorbar(title = "SVI",
                                           title.position = "top",
                                           title.hjust = 0.5, 
                                           nbin = 10)
  )+
  scale_size_binned(labels = scales::comma,
                    range = c(1,14),
                    # nice.breaks = T,
                    # n.breaks = 7,
                    guide = guide_legend(title = "Population size",
                                         title.position = "top",
                                         title.hjust = 0.5))+
  labs(x = "Omicron BA.1* \n Attack Rate",
       y = "Omicron BA.2* \n Attack Rate")+
  theme(legend.position = "bottom", 
        legend.key.width = grid::unit(2.5, "cm"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))
ba1_ba2

ba1_ba5 <- states_ar_svi |> 
  ggplot(aes(x = `Omicron BA.1*`, 
             y = `Omicron BA.5*`,
             # label = `State Code`, 
             col = SVI, 
             size = pop))+
  geom_point(alpha = 0.5)+
  geom_smooth(aes(weight = pop),
              method = "lm",
              show.legend = F,
              alpha = 0.25)+
  geom_text(aes(label = `State Code`), 
            color = "white", 
            size = 2.5)+
  ggpubr::stat_cor(method = "pearson",  
                   cor.coef.name = "R", 
                   geom = "label",
                   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
                   r.digits = 2, 
                   r.accuracy = 0.01, 
                   size = 5)+
  theme_minimal()+
  scale_x_continuous(labels = scales::percent)+
  scale_y_continuous(labels = scales::percent)+
  scale_color_met_c(palette = "Tam", 
                    # direction = -1,
                    # show.limits = T,
                    guide = guide_colorbar(title = "SVI",
                                           title.position = "top",
                                           title.hjust = 0.5, 
                                           nbin = 10)
  )+
  scale_size_binned(labels = scales::comma,
                    range = c(1,14),
                    # nice.breaks = T,
                    # n.breaks = 7,
                    guide = guide_legend(title = "Population size",
                                         title.position = "top",
                                         title.hjust = 0.5))+
  labs(x = "Omicron BA.1* \n Attack Rate",
       y = "Omicron BA.5* \n Attack Rate")+
  theme(legend.position = "bottom", 
        legend.key.width = grid::unit(2.5, "cm"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))
ba1_ba5

library(patchwork)

patchwork_ar_svi <- (((ba1_ba2) / 
                        (ba1_ba5))|
                       (ar_svi))+
  plot_layout(guides = 'collect',
              tag_level = "new")+
  plot_annotation(tag_levels = 'A')&
  theme(legend.position = "bottom")
patchwork_ar_svi

ggsave(filename = "Output/Plots/patchwork_ar_svi.png",
       plot = patchwork_ar_svi, 
       width = 16, 
       height = 9,
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/patchwork_ar_svi.pdf",
       plot = patchwork_ar_svi, 
       width = 16, 
       height = 9,
       dpi = 100)

#