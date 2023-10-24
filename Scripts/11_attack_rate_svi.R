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
  # ## filtering to days before the peak of infections of each variant, all states
  # filter(days <= days[which.max(infections)],
  #        .by = c("name_states", "variant")) |>
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

# states_omicron_ba1_ar <- states_attack_rates|> 
#   filter(variant == "Omicron BA.1*") |> 
#   dplyr::select(name_states, attack_rate) |> 
#   rename("Omicron BA.1* Attack Rate" = attack_rate)
# 
# states_omicron_ba2_ar <- states_attack_rates|> 
#   filter(variant == "Omicron BA.2*") |> 
#   dplyr::select(name_states, attack_rate) |> 
#   rename("Omicron BA.2* Attack Rate" = attack_rate)
# 
# states_omicron_ba4_ar <- states_attack_rates|> 
#   filter(variant == "Omicron BA.4*") |> 
#   dplyr::select(name_states, attack_rate) |> 
#   rename("Omicron BA.4* Attack Rate" = attack_rate)
# 
# states_omicron_ba5_ar <- states_attack_rates|> 
#   filter(variant == "Omicron BA.5*") |> 
#   dplyr::select(name_states, attack_rate) |> 
#   rename("Omicron BA.5* Attack Rate" = attack_rate)
# 
# states_omicron_xbb_ar <- states_attack_rates|> 
#   filter(variant == "Omicron XBB*") |> 
#   dplyr::select(name_states, attack_rate) |> 
#   rename("Omicron XBB* Attack Rate" = attack_rate)
# 
# states_joined <- full_join(states_omicron_ba1_ar, 
#                            states_omicron_ba2_ar) |>
#   left_join(states_omicron_ba4_ar) |> 
#   left_join(states_omicron_ba5_ar) |> 
#   left_join(states_omicron_xbb_ar) |> 
#   left_join(states_abb)

## Needs to join with the SVI per states, binned into 3 or 5 categories
## low, medium, high SVI, or low, mid-low, medium, mid-high, high

# Loading functions
source("Scripts/Functions/functions.R")
source("Scripts/get_svi.R")

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
# st_write(states_svi, "Data/states_svi.shp")

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
  # ggrepel::geom_label_repel(aes(label = `State Code`, 
  #                               geometry = geometry),
  #                           size= 3,
  #                           stat = "sf_coordinates", 
  #                           max.overlaps = Inf,
  #                           color = "black")+
  theme_void()+
  scale_fill_met_c(palette = "Paquin", 
                   direction = -1,
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
  scale_color_met_d(palette = "Klimt")+
  scale_fill_met_d(palette = "Klimt")+
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
  scale_color_met_c(palette = "Paquin", 
                    direction = -1, 
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
  scale_color_met_c(palette = "Paquin", 
                    direction = -1,
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

# ba2_ba4 <-  states_ar_svi |>
#   ggplot(aes(x = `Omicron BA.4*`,
#              y = `Omicron BA.2*`,
#              label = `State Code`,
#              col = SVI_rank))+
#   geom_point()+
#   geom_text_repel(show.legend = F)+
#   theme_bw()+
#   scale_x_continuous(labels = scales::percent)+
#   scale_y_continuous(labels = scales::percent)+
#   scale_color_met_c(palette = "Johnson",
#                     direction = -1,
#                     # show.limits = T,
#                     guide = guide_colorbar(title = "Social Vulnerability Index",
#                                            title.position = "top",
#                                            title.hjust = 0.5,
#                                            nbin = 5)
#   )+
#   labs(x = "Omicron BA.4* \n Attack Rate",
#        y = "Omicron BA.2* \n Attack Rate")+
#   theme(legend.position = "bottom",
#         legend.key.width = grid::unit(3, "cm"),
#         axis.title = element_text(size = 14),
#         axis.text = element_text(size = 14))
# ba2_ba4
# 
# ggsave(filename = "Output/Plots/ba2_ba4_ar_svi.png",
#        plot = ba2_ba4,
#        width = 16,
#        height = 9,
#        dpi = 100)
# 
# ba2_ba5 <- states_ar_svi |>
#   ggplot(aes(x = `Omicron BA.5*`,
#              y = `Omicron BA.2*`,
#              label = `State Code`,
#              col = SVI_rank))+
#   geom_point()+
#   geom_text_repel(show.legend = F)+
#   theme_bw()+
#   scale_x_continuous(labels = scales::percent)+
#   scale_y_continuous(labels = scales::percent)+
#   scale_color_met_c(palette = "Johnson",
#                     direction = -1,
#                     # show.limits = T,
#                     guide = guide_colorbar(title = "Social Vulnerability Index",
#                                            title.position = "top",
#                                            title.hjust = 0.5,
#                                            nbin = 5)
#   )+
#   labs(x = "Omicron BA.5* \n Attack Rate",
#        y = "Omicron BA.2* \n Attack Rate")+
#   theme(legend.position = "bottom",
#         legend.key.width = grid::unit(3, "cm"),
#         axis.title = element_text(size = 14),
#         axis.text = element_text(size = 14))
# ba2_ba5
# 
# ggsave(filename = "Output/Plots/ba2_ba5_ar_svi.png",
#        plot = ba2_ba5,
#        width = 16,
#        height = 9,
#        dpi = 100)
# 
# ba4_ba5 <- states_ar_svi |>
#   ggplot(aes(x = `Omicron BA.4*`,
#              y = `Omicron BA.5*`,
#              label = `State Code`,
#              col = SVI_rank))+
#   geom_point()+
#   geom_text_repel(show.legend = F)+
#   theme_bw()+
#   scale_x_continuous(labels = scales::percent)+
#   scale_y_continuous(labels = scales::percent)+
#   scale_color_met_c(palette = "Johnson",
#                     direction = -1,
#                     # show.limits = T,
#                     guide = guide_colorbar(title = "Social Vulnerability Index",
#                                            title.position = "top",
#                                            title.hjust = 0.5,
#                                            nbin = 5)
#   )+
#   labs(x = "Omicron BA.4* \n Attack Rate",
#        y = "Omicron BA.5* \n Attack Rate")+
#   theme(legend.position = "bottom",
#         legend.key.width = grid::unit(3, "cm"),
#         axis.title = element_text(size = 14),
#         axis.text = element_text(size = 14))
# ba4_ba5
# 
# ggsave(filename = "Output/Plots/ba4_ba5_ar_svi.png",
#        plot = ba4_ba5,
#        width = 16,
#        height = 9,
#        dpi = 100)
# 
# ba5_xbb <- states_ar_svi |>
#   ggplot(aes(x = `Omicron XBB*`,
#              y = `Omicron BA.5*`,
#              label = `State Code`,
#              col = SVI_rank))+
#   geom_point()+
#   geom_text_repel(show.legend = F)+
#   theme_bw()+
#   scale_x_continuous(labels = scales::percent)+
#   scale_y_continuous(labels = scales::percent)+
#   scale_color_met_c(palette = "Johnson",
#                     direction = -1,
#                     # show.limits = T,
#                     guide = guide_colorbar(title = "Social Vulnerability Index",
#                                            title.position = "top",
#                                            title.hjust = 0.5,
#                                            nbin = 5)
#   )+
#   labs(x = "Omicron XBB* \n Attack Rate",
#        y = "Omicron BA.5* \n Attack Rate")+
#   theme(legend.position = "bottom",
#         legend.key.width = grid::unit(3, "cm"),
#         axis.title = element_text(size = 14),
#         axis.text = element_text(size = 14))
# ba5_xbb
# 
# ggsave(filename = "Output/Plots/ba5_xbb_ar_svi.png",
#        plot = ba5_xbb,
#        width = 16,
#        height = 9,
#        dpi = 100)
# 
# library(patchwork)
# 
# patchwork_ar <- (ba1_ba2 | ba2_ba4 | ba2_ba5 | ba4_ba5 | ba5_xbb)+
#   plot_layout(guides = 'collect')&
#   theme(legend.position = "bottom")
# patchwork_ar
# 
# ggsave(filename = "Output/Plots/patchwork_ar_svi.png",
#        plot = patchwork_ar,
#        width = 16,
#        height = 9,
#        dpi = 100)

# ## Bivariate maps to Attack Rates vs SVI
# library(multiscales)
# library(colorspace)
# library(MetBrewer)
#
# states_ar_svi_map <- states_ar_svi |>
#   pivot_longer(cols = c(`Omicron BA.1*`:`Omicron XBB*`),
#                names_to = "variants",
#                values_to = "attack_rate") |>
#   # bi_class(x = SVI, y = attack_rate, style = "quantile", dim = 4) |>
#   st_as_sf() |>
#   tigris::shift_geometry()
#
# ar_svi_map <-
#   ggplot(states_ar_svi_map |>
#            filter(variants == "Omicron BA.1*"))+
#   geom_sf(aes(fill = zip(attack_rate/100, SVI)),
#           color = "gray30",
#           size = 0.2) +
#   bivariate_scale("fill",
#                   pal_vsup(colorspace::divergingx_hcl(palette = "RdYlBu", 8)),
#                   name = c("Attack Rate", "SVI"),
#                   # limits = list(c(-40, 40), c(0, 1)),
#                   # breaks = list(c(-40, -20, 0, 20, 40), c(0, 0.25, 0.50, 0.75, 1.)),
#                   labels = list(scales::percent, waiver()),
#                   guide = "colourfan"
#   ) +
#   theme_void() +
#   theme(
#     legend.key.size = grid::unit(0.8, "cm"),
#     legend.title.align = 0.5,
#     plot.margin = margin(5.5, 20, 5.5, 5.5)
#   )
# ar_svi_map
#
# ## Labels
# labels1 <- bi_class_breaks(states_ar_svi_map,
#                            x = SVI,
#                            y = attack_rate,
#                            style = "quantile",
#                            dim = 4,
#                            dig_lab = 2,
#                            split = FALSE)
# ## Legend
# legend <- bi_legend(pal = custom_pal4_1,
#                     breaks = labels1,
#                     dim = 4,
#                     xlab = "Higher SVI",
#                     ylab = "Higher Attack Rate",
#                     size = 12)+
#   theme_minimal()+
#   # bi_theme()+
#   theme(axis.text.x = element_text(angle = 90))
#
