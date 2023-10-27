## Cleaning the ambient
rm(list = ls())
gc()

## Loading Libraries
packs = c("tidyverse", "vroom", "patchwork", "geofacet", "tigris", "ggforce", "ggthemes", "MetBrewer", "ggrepel", "gghighlight")
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

## Loading the frequencies of variants
variant_count<-vroom("Data/variant_counts_us.csv.xz") |>
  rename(days = epiweek,
         variant = voc_cdc) |>
  filter(!variant %in% c("Alpha*", "Beta*", "Gamma*", "Delta*", "Omicron BA.3*", "Recombinant", "Other")) |>
  mutate(variant = droplevels(factor(variant))) |>
  mutate(variant = factor(variant,
                          levels = c("Omicron BA.1*", "Omicron BA.2*", "Omicron BA.4*", "Omicron BA.5*",
                                     "Omicron XBB*")))

## Loading estimates_rt_incidence dataset
estimates_rt_incidence <- vroom("Data/state_full_data.csv.xz")


## Figure1
figure1a_data <- rt_estimates |> 
  filter(days <= "2023-03-01") |>
  ## Renaming infections
  rename(infections = I) |>
  group_by(name_states, days) |> 
  summarise(infections = sum(infections)) |> 
  left_join(pop_states, by = c("name_states" = "state")) |>
  ## Creating incidence per 100k
  mutate(incidence = (infections/pop)*1e5, 
         percentual_incidence = round(incidence/pop, 2)) 

figure1a_data_national <- estimates_rt_incidence |>
  filter(days <= "2023-03-01") |>
  reframe(mean_infections = sum(infections, na.rm = T),
          .by = c(days))

scaleFactor <- max(figure1a_data$infections, na.rm = T)/max(figure1a_data_national$mean_infections, na.rm = T)

figure1a<-figure1a_data |> 
  ggplot(aes(x = days, y = infections, 
             group = name_states))+
  geom_line(alpha = .1)+
  geom_line(data = figure1a_data_national,
            aes(x = days, y = mean_infections*scaleFactor, 
                col = "Mean infections per day"), 
            inherit.aes = F, 
            show.legend = F)+
  labs(x = NULL, 
       y = "Infections per days \n State-level", 
       title = "All variants")+
  theme_minimal()+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y.right = element_text(color = "firebrick1", size = 14),
        axis.title.y.left = element_text(size = 14),
        axis.text.y.right = element_text(colour = "firebrick1", size = 14))+
  scale_x_date(date_breaks = "2 months", date_labels = "%b '%y")+
  scale_y_continuous(labels = scales::label_comma(), 
                     breaks = pretty(figure1a_data$infections),
                     sec.axis = sec_axis(~./scaleFactor, 
                                         labels = scales::label_comma(),
                                         breaks = pretty(figure1a_data_national$mean_infections),
                                         name = "Infections per day \n Whole country"))+
  scale_color_manual(name = "", values = "firebrick1")
figure1a

ggsave(filename = "Output/Plots/figure1a_manuscript.png", 
       plot = figure1a, 
       width = 16, 
       height = 9, 
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figure1a_manuscript.pdf", 
       plot = figure1a, 
       width = 16,
       height = 9, 
       dpi = 100)

figure1b_data <- estimates_rt_incidence |> 
  filter(days <= "2023-03-01", incidence >= 10) |>
  reframe(infections = sum(infections, na.rm = T), 
          .by = c(days, name_states, variant))

figure1b_data_national <- estimates_rt_incidence |>
  filter(days <= "2023-03-01", incidence >= 10) |>
  reframe(mean_infections = sum(infections, na.rm = T),
          .by = c(days, variant))

scaleFactor <- max(figure1b_data$infections, 
                   na.rm = T)/max(figure1b_data_national$mean_infections, 
                                  na.rm = T)

figure1b <- figure1b_data |> 
  ggplot(aes(x = days, y = infections, 
             col = variant,
             group = name_states))+
  geom_line(alpha = .1, 
            show.legend = F)+
  geom_line(data = figure1b_data_national,
            aes(x = days, y = mean_infections*scaleFactor, 
                col = "National mean"),
            # show.legend = F,
            inherit.aes = F)+
  theme_minimal()+
  facet_wrap(variant~., ncol = 1, scales = "free_y")+
  labs(x = "Date", y = "Infections per day \n State-level")+
  scale_x_date(date_breaks = "2 months", date_labels = "%b '%y")+
  scale_y_continuous(labels = scales::label_comma(),
                     sec.axis = sec_axis(~./scaleFactor, 
                                         labels = scales::label_comma(),
                                         name = "Infections per day \n Whole country", ))+
  scale_color_manual(values = c("firebrick1", 
                                MetBrewer::met.brewer(palette_name = "Archambault", 
                                                      type = "discrete", 
                                                      n = 5)),
                     name = "")+
  theme(legend.position = "none", 
        strip.text.x.top = element_text(hjust = 0),
        axis.text.x = element_text(size = 14),
        axis.title.y.right = element_text(color = "firebrick1", size = 14),
        axis.title.y.left = element_text(size = 14),
        axis.text.y.right = element_text(colour = "firebrick1", size = 14), 
        axis.text = element_text(size = 14))
figure1b

ggsave(filename = "Output/Plots/figure1b_manuscript.png", 
       plot = figure1b, 
       width = 16, 
       height = 9, 
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figure1b_manuscript.pdf", 
       plot = figure1b, 
       width = 16,
       height = 9, 
       dpi = 100)

figure1_patchwork<-((figure1a+theme(axis.text.x = element_blank())) / figure1b)+
  plot_annotation(tag_levels = 'A')+
  plot_layout(guides = "collect", 
              heights = c(1.2,1.6))
figure1_patchwork

ggsave(filename = "Output/Plots/figure1_manuscript.png", 
       plot = figure1_patchwork, 
       width = 16,
       height = 9, 
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figure1_manuscript.pdf", 
       plot = figure1_patchwork, 
       width = 16,
       height = 9, 
       dpi = 100)

states_abb<-vroom("Data/states_abbrev_region_division.csv")

figure2a_data<-estimates_rt_incidence |> 
  # filter(days <= "2023-04-01") |> 
  reframe(total_incidence = sum(incidence, na.rm = T), 
          total_infections = sum(infections, na.rm = T),
          .by = c("name_states", "variant"))|> 
  left_join(states_abb, by = c("name_states" = "State")) |> 
  mutate(Region = factor(Region, 
                         levels = c("West", "Midwest", "South", "Northeast"))) |> 
  left_join(pop_states, by = c("name_states" = "state")) |> 
  mutate(attack_rate = (total_infections/pop))

vroom_write(x = figure2a_data, 
            file = "Data/state_attack_rate_variants.csv.xz")

# scaleFactor<-max(figure2a_data$total_incidence, na.rm = T)/max(figure2a_data$attack_rate)

figure2a<- figure2a_data |> 
  ggplot(aes(x = variant, 
             col = variant,
             y = attack_rate))+    
  geom_boxplot()+
  geom_point(position = position_jitterdodge(jitter.width = .5))+
  # geom_point(aes(y = attack_rate*scaleFactor))+
  geom_text(data = figure2a_data |> 
              filter(attack_rate == max(attack_rate, na.rm = T),
                     .by = c(variant)),
            aes(label = `State Code`), 
            vjust = -1.5, nudge_x = 0.2)+
  geom_text(data = figure2a_data |> 
              filter(attack_rate == min(attack_rate, na.rm = T),
                     .by = c(variant)),
            aes(label = `State Code`), 
            vjust = 2)+
  theme_minimal()+
  MetBrewer::scale_color_met_d(palette_name = "Archambault", name = "")+
  theme(legend.position = "none", 
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14))+
  labs(x = "", y = "Attack Rate \n (%) of pop. ever infected")+
  scale_y_continuous(labels = scales::label_percent(),
                     breaks = pretty(figure2a_data$attack_rate))
figure2a

ggsave(filename = "Output/Plots/figure2a_manuscript.png", 
       plot = figure2a, 
       width = 16, 
       height = 9, 
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figure2a_manuscript.pdf", 
       plot = figure2a, 
       width = 16,
       height = 9, 
       dpi = 100)

figure2_data<-estimates_rt_incidence |> 
  # filter(days <= "2023-04-01") |> 
  reframe(infections = sum(infections), 
          .by = c("name_states", "variant"))|>
  left_join(states_abb, by = c("name_states" = "State")) |> 
  left_join(pop_states, by = c("name_states" = "state")) |> 
  mutate(attack_rate = round(infections/pop*100, 0))

figure2b<-figure2_data |> 
  ggplot(aes(x = variant, 
             y = attack_rate, 
             fill = variant))+
  geom_col()+
  geom_text(aes(label = attack_rate, fill = NULL),
            size = 2.5,
            vjust = -0.5,
            color = "grey50", show.legend = F)+
  lims(y = c(NA, 65))+
  MetBrewer::scale_fill_met_d(palette_name = "Archambault", name = "")+
  facet_geo(`State Code`~., scales = "free", strip.position = "right")+
  labs(y = "Attack rate \n (%) of pop. ever infected")+
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "bottom", 
        panel.background = element_blank(), 
        panel.grid = element_blank())
figure2b  

ggsave(filename = "Output/Plots/figure2b_manuscript.png", 
       plot = figure2b, 
       width = 16, 
       height = 9, 
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figure2b_manuscript.pdf", 
       plot = figure2b, 
       width = 16,
       height = 9, 
       dpi = 100)

# figure2_patchwork<-(figure2a / figure2b)+
#   plot_layout(heights = c(1.2, 1.8))+
#   plot_annotation(tag_levels = 'a')
# figure2_patchwork
# 
# ggsave(filename = "Output/Plots/figure2_manuscript.png", 
#        plot = figure2_patchwork, 
#        width = 11, 
#        height = 9, 
#        dpi = 100)
# 
# ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figure2_manuscript.pdf", 
#        plot = figure2_patchwork, 
#        width = 11,
#        height = 9, 
#        dpi = 100)

usmap<-tigris::states(cb = TRUE) |> 
  tigris::shift_geometry()

map_attackrate<-usmap |> 
  left_join(figure2_data, 
            by = c("NAME" = "name_states")) |> 
  mutate(NAME = droplevels(factor(NAME))) |> 
  filter(!is.na(variant)) |> 
  mutate(attack_rate_binned = cut(attack_rate, 
                                  breaks = pretty(attack_rate, n = 10),
                                  right = T))

unique(map_attackrate$attack_rate_binned)

attack_rate_map<-function(x){
  # Remove plot axis
  no_axis <- theme(axis.title=element_blank(),
                   axis.text=element_blank(),
                   axis.ticks=element_blank())
  
  plt <- map_attackrate |> 
    filter(variant == x) |> 
    ggplot()+
    geom_sf(aes(fill = attack_rate_binned), 
            color = "grey50", 
            size = 0.2)+
    theme_map()+
    scale_fill_manual(values = met.brewer(palette_name = "Johnson", 
                                          n = 13, 
                                          type = "continuous", 
                                          direction = -1), 
                      drop = F, 
                      labels = waiver(),
                      name = "Attack rate \n (%) of pop. ever infected",
                      guide = guide_bins(keywidth = grid::unit(1.2, "cm"),
                                         title.position = "top", title.hjust = 0.5, 
                                         show.limits = T), 
                      aesthetics = c("color", "fill"))+
    theme(legend.position = "bottom", 
          legend.title = element_text(size = 14), 
          legend.text = element_text(size = 14))+
    labs(subtitle = x)
  # plt
  return(plt)
}

variants<- as.character(unique(map_attackrate$variant))

map_list<-lapply(variants, attack_rate_map)

patchwork_map<-(((map_list[[1]] / map_list[[4]]) | (map_list[[2]]/ map_list[[3]] / map_list[[5]]))/guide_area())+
  plot_annotation(title = 'B')+
  plot_layout(guides = "collect", 
              widths = c(5, .5),
              heights = c(5, .5))&
  theme(legend.position = "bottom")
patchwork_map

ggsave(filename = "Output/Plots/figure2_attack_rate_map.png", 
       plot = patchwork_map, 
       width = 16,
       height = 9, 
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figure2b_map_manuscript.pdf", 
       plot = patchwork_map, 
       width = 16,
       height = 9, 
       dpi = 100)

patchwork_figure2<-(figure2a + 
                      theme(legend.position = "none", 
                            axis.text.x = element_text(angle = 90, 
                                                       size = 14)) | 
                      (patchwork_map+theme(axis.title = element_text(size = 14))))+
  plot_layout(widths = c(1, 3), 
              heights = c(1,1), 
              tag_level = 'new')+
  plot_annotation(tag_levels = list(c('A', 'B'), '', '', ''))
patchwork_figure2

ggsave(filename = "Output/Plots/figure2_manuscript_map.png", 
       plot = patchwork_figure2, 
       width = 16,
       height = 9, 
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figure2_manuscript_map.pdf", 
       plot = patchwork_figure2, 
       width = 16,
       height = 9, 
       dpi = 100)

figure3_data<- estimates_rt_incidence |> 
  dplyr::select(days, variant, name_states, Rt, upper, lower) |> 
  drop_na()

figure3_data_national<- estimates_rt_incidence |> 
  reframe(median = mean(Rt, na.rm = T),
          lower = mean(lower, na.rm = T),
          upper = mean(upper, na.rm = T),
          .by = c(days, variant))

figure3<- figure3_data |> 
  ggplot(aes(x = days, 
             y = Rt, 
             # ymin = lower, 
             # ymax = upper,
             # fill = variant,
             group = name_states))+
  geom_hline(yintercept = 1, 
             show.legend = F, 
             col = "grey9")+
  # geom_line(aes(col = variant),
  #           alpha = .1,
  #           show.legend = F)+
  # geom_ribbon(alpha = .01,
  #             show.legend = F)+
  # geom_line(data = figure3_data_national,
  #           aes(x = days,
  #               y = median,
  #               col = "National average"),
  #           inherit.aes = F,
  #           alpha = .1,
#           lwd = .8)+
geom_smooth(aes(col = variant,
                fill = variant),
            alpha = .1,
            lwd = .1,
            method = "loess",
            show.legend = FALSE)+
  geom_smooth(data = figure3_data_national,
              aes(x = days,
                  y = median,
                  col = "National average"),
              inherit.aes = F,
              alpha = .1,
              lwd = .8,
              method = "loess")+
  theme_minimal()+
  labs(x = "Date", y = "Effective Reproduction \n (Rt)")+
  facet_wrap(variant~., 
             ncol = 1, 
             scales = "free_y", 
             strip.position = "top")+
  scale_x_date(date_breaks = "2 months",
               date_labels = "%b '%y"
  )+
  scale_color_manual(values = c("firebrick3",
                                met.brewer(palette_name = "Archambault",
                                           type = "discrete", 
                                           n = 5)))+
  scale_fill_manual(values = c("firebrick3",
                               met.brewer(palette_name = "Archambault", 
                                          type = "discrete", 
                                          n = 5)))+
  theme(legend.title = element_blank(), 
        legend.position = "bottom", 
        axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        strip.text = element_text(hjust = 0))+
  guides(color = "none", 
         fill = "none")
figure3

ggsave(filename = "Output/Plots/figure3_manuscript.png", 
       plot = figure3, 
       width = 16, 
       height = 9, 
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figure3_manuscript.pdf", 
       plot = figure3, 
       width = 16,
       height = 9, 
       dpi = 100)

## log-scale Rt figure
figure3_log <- figure3+scale_y_log10()+labs(x = "Date", y = "log Effective Reproduction \n (Rt)")
figure3_log

ggsave(filename = "Output/Plots/figure3_logscale_manuscript.png", 
       plot = figure3_log, 
       width = 16, 
       height = 9, 
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figure3_loscale_manuscript.pdf", 
       plot = figure3_log, 
       width = 16,
       height = 9, 
       dpi = 100)

rt_ratio<-vroom("Output/Tables/rt_ratios.csv.xz")|>
  filter(ratios != "XBB*/BA.2*", ratios != "BA.4*/BA.5*") |>
  mutate(ratios = factor(droplevels(factor(ratios)),
                         levels = c("BA.2*/BA.1*", "BA.4*/BA.2*", "BA.5*/BA.2*",
                                    "BA.5*/BA.4*", "XBB*/BA.5*")))

figure4a_data<-rt_ratio |> 
  reframe(median = mean(median, na.rm = T),
          lower = mean(lower, na.rm = T),
          upper = mean(upper, na.rm = T), 
          .by = c(days, ratios))

figure4a <- figure4a_data |> 
  mutate(ratios = factor(ratios, 
                         levels = c("BA.2*/BA.1*", "BA.4*/BA.2*", "BA.5*/BA.2*", "BA.5*/BA.4*", 
                                    "XBB*/BA.5*"))) |> 
  ggplot(aes(x = days, y = median, 
             ymin = lower, ymax = upper,
             color = ratios, fill = ratios))+
  geom_hline(yintercept = 1, show.legend = F, col = "grey9")+
  geom_line()+
  geom_ribbon(alpha = .1)+
  theme_minimal()+
  facet_wrap(ratios~., ncol = 1, scales = "fixed", strip.position = "top")+
  theme(legend.position = "none", 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14), 
        strip.text.x.top = element_text(hjust = 0))+
  scale_x_date(date_breaks = "2 months", date_labels = "%b '%y")+
  labs(y = "Rt ratio", x = "")+
  MetBrewer::scale_color_met_d(palette_name = "Johnson", name = "", direction = -1)+
  MetBrewer::scale_fill_met_d(palette_name = "Johnson", name = "", direction = -1)
figure4a

ggsave(filename = "Output/Plots/figure4a_manuscript.png", 
       plot = figure4a, 
       width = 16, 
       height = 9, 
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figure4a_manuscript.pdf", 
       plot = figure4a, 
       width = 16,
       height = 9, 
       dpi = 100)

figure4b_data <- rt_ratio |> 
  reframe(median = mean(median, na.rm = T),
          lower = mean(lower, na.rm = T),
          upper = mean(upper, na.rm = T), 
          .by = c(name_states, ratios)) |> 
  mutate(name_states2 = tidytext::reorder_within(name_states, median, within = ratios))

figure4b <- figure4b_data |> 
  # mutate(ratios = factor(ratios, 
  #                        levels = rev(c("BA.2*/BA.1*", "BA.4*/BA.2*", "BA.5*/BA.2*", "BA.5*/BA.4*", 
  #                                       "XBB*/BA.5*")))) |> 
  ggplot(aes(y = median, ymin = lower, ymax = upper, 
             x = ratios, 
             col = ratios))+
  geom_hline(yintercept = 1, show.legend = F, col = "grey9")+
  geom_boxplot(show.legend = F)+
  geom_point(position = position_dodge2(width = .1))+
  theme_minimal()+
  theme(legend.title = element_blank(), 
        legend.position = "none", 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))+
  MetBrewer::scale_color_met_d(palette_name = "Johnson", name = "", direction = -1)+
  MetBrewer::scale_fill_met_d(palette_name = "Johnson", name = "", direction = -1)+
  labs(y = 'Average (Rt) ratio', x = 'Ratios')
figure4b

ggsave(filename = "Output/Plots/figure4b_manuscript.png", 
       plot = figure4b, 
       width = 16, 
       height = 9, 
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figure4b_manuscript.pdf", 
       plot = figure4b, 
       width = 16,
       height = 9, 
       dpi = 100)

figure4c_data <- rt_ratio |> 
  reframe(median = mean(median, na.rm = T),
          lower = mean(lower, na.rm = T),
          upper = mean(upper, na.rm = T), 
          .by = c("name_states", "ratios")) |> 
  left_join(states_abb, by = c("name_states" = "State"))

figure4c <- figure4c_data |> 
  ggplot(aes(x = ratios, y = median, group = `State Code`)) +
  geom_hline(yintercept = 1, show.legend = F, col = "grey9")+
  geom_line(aes(color = ratios), alpha = .5, size = 1, show.legend = F) +
  geom_point(aes(color = ratios), alpha = .5, size = 5, show.legend = T) +
  geom_text_repel(data = figure4c_data %>%
                    filter(ratios == "BA.2*/BA.1*"
                           ,
                           `State Code` %in% c("NY", "CA", "CT", "KY", "SD", "DE", "MA", "ID")
                    ),
                  aes(label = `State Code`) ,
                  hjust = "left",
                  fontface = "bold",
                  size = 3,
                  nudge_x = -.45,
                  direction = "y",
                  show.legend = F) +
  geom_text_repel(data = figure4c_data %>%
                    filter(ratios == "XBB*/BA.5*"
                           ,
                           `State Code` %in% c("NY", "CA", "CT", "KY", "SD", "DE", "MA", "ID")
                    ),
                  aes(label = `State Code`) ,
                  hjust = "right",
                  fontface = "bold",
                  size = 3,
                  nudge_x = .5,
                  direction = "y",
                  show.legend = F) +
  # geom_label(aes(label = round(median, 2)),
  #            size = 2.5,
  #            # vjust = -0.5,
  #            label.padding = unit(0.05, "lines"),
  #            label.size = 0.0,
  #            show.legend = F)+
  MetBrewer::scale_color_met_d(palette_name = "Johnson", name = "", direction = -1)+
  MetBrewer::scale_fill_met_d(palette_name = "Johnson", name = "", direction = -1)+
  theme_minimal()+
  labs(x = "Ratios",
       y = "Average (Rt) ratio")+
  theme(legend.title = element_blank(),
        legend.position = "bottom", 
        axis.title = element_text(size = 14),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 14))

figure4c <- figure4c + 
  gghighlight(`State Code` %in% c("NY", "CA", "CT", "KY", "SD", "DE", "MA", "ID"),
              use_direct_label = F)
figure4c

ggsave(filename = "Output/Plots/figure4c_manuscript.png", 
       plot = figure4c, 
       width = 16, 
       height = 9, 
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figure4c_manuscript.pdf", 
       plot = figure4c, 
       width = 16,
       height = 9, 
       dpi = 100)

## Patchwork
patchwork_figure4 <- (figure4a + figure4b)+
  plot_annotation(tag_levels = 'A')
patchwork_figure4

ggsave(filename = "Output/Plots/figure4_manuscript.png", 
       plot = patchwork_figure4, 
       width = 16, 
       height = 9, 
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figure4_manuscript.pdf", 
       plot = patchwork_figure4, 
       width = 16,
       height = 9, 
       dpi = 100)

figureS2 <- figure4c_data |> 
  ggplot(aes(x = ratios, y = median, group = `State Code`)) +
  geom_hline(yintercept = 1, show.legend = F, col = "grey9")+
  geom_line(aes(color = ratios), alpha = .5, size = 1, show.legend = F) +
  geom_point(aes(color = ratios), alpha = .5, size = 5, show.legend = T) +
  # geom_label(aes(label = round(median, 2)),
  #            size = 2.5,
  #            # vjust = -0.5,
  #            label.padding = unit(0.05, "lines"),
  #            label.size = 0.0,
  #            show.legend = F)+
  MetBrewer::scale_color_met_d(palette_name = "Johnson", name = "", direction = -1)+
  MetBrewer::scale_fill_met_d(palette_name = "Johnson", name = "", direction = -1)+
  theme_minimal()+
  lims(y = c(NA, 1.7))+
  labs(x = "Ratios",
       y = "Average (Rt) ratio")+
  theme(legend.title = element_blank(),
        legend.position = "bottom", 
        axis.text.x = element_text(angle = 90))+
  facet_geo(`State Code`~.)
figureS2

ggsave(filename = "Output/Plots/figureS2_manuscript.png", 
       plot = figureS2, 
       width = 16, 
       height = 9, 
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figureS2_manuscript.pdf", 
       plot = figureS2, 
       width = 16,
       height = 9, 
       dpi = 100)

figureS1a <- variant_count |> 
  filter(days < "2023-05-01", days > "2021-12-01") |> 
  reframe(n = sum(n, na.rm = T), 
          .by = c("days", "variant")) |> 
  ggplot(aes(x = days, y = n, fill = variant))+
  geom_col(width = 6)+
  theme_minimal()+
  scale_fill_met_d(palette_name = "Archambault", 
                   name = "")+
  theme(legend.position = "bottom", 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))+
  scale_x_date(name = "Date", date_breaks = "2 months", date_labels = "%b '%y")+
  labs(y = "Number of genomic sequences")+
  guides(title.position = "top", title.vjust = .5)
figureS1a

ggsave(filename = "Output/Plots/figureS1a_manuscript.png", 
       plot = figureS1a, 
       width = 16, 
       height = 9, 
       dpi = 200)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figureS1a_manuscript.pdf", 
       plot = figureS1a, 
       width = 16,
       height = 9, 
       dpi = 200)

figureS1b <- variant_count |> 
  filter(days < "2023-05-01", days > "2021-12-01") |> 
  ggplot(aes(x = days, y = n, fill = variant, group = variant))+
  geom_col(width = 7)+
  theme_minimal()+
  scale_fill_met_d(palette_name = "Archambault", 
                   name = "")+
  theme(legend.position = "bottom")+
  scale_x_date(name = "Date", date_breaks = "2 months", date_labels = "%b %y")+
  labs(y = "Number of genomic sequences")+
  facet_geo(name_states~., scales = "free_y")+
  theme(axis.text.x = element_text(angle = 90))
figureS1b

ggsave(filename = "Output/Plots/figureS1b_manuscript.png", 
       plot = figureS1b, 
       width = 16, 
       height = 9, 
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figureS1b_manuscript.pdf", 
       plot = figureS1b, 
       width = 16,
       height = 9, 
       dpi = 100)

figureS1c <- variant_count |> 
  filter(days < "2023-04-01", days > "2021-12-01") |> 
  ggplot(aes(x = days, y = freq, fill = variant, group = variant))+
  geom_col(position = position_fill(), width = 7)+
  theme_minimal()+
  scale_fill_met_d(palette_name = "Archambault", 
                   name = "")+
  theme(legend.position = "bottom")+
  scale_x_date(name = "Date", date_breaks = "2 months", date_labels = "%b %y")+
  labs(y = "Frequency (%) \n from the raw number of genomic sequences")+
  facet_geo(name_states~.)+
  theme(axis.text.x = element_text(angle = 90))
figureS1c

ggsave(filename = "Output/Plots/figureS1c_manuscript.png", 
       plot = figureS1c, 
       width = 16, 
       height = 9, 
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figureS1c_manuscript.pdf", 
       plot = figureS1c, 
       width = 16,
       height = 9, 
       dpi = 100)

figureS1d <- estimates_rt_incidence |>
  ggplot(aes(x = days, y = infections, col = variant))+
  geom_line()+
  theme_minimal()+
  scale_color_met_d(palette_name = "Archambault", 
                    name = "")+
  theme(legend.position = "bottom")+
  scale_x_date(name = "Date", date_breaks = "2 months", date_labels = "%b %y")+
  labs(y = "Infections per day estimates")+
  facet_geo(name_states~., scales = "free_y")+
  theme(axis.text.x = element_text(angle = 90))
figureS1d

#