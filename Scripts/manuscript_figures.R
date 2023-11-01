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

## Figure 1A
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

## Figure 1B
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

## Figure patchwork
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

## Figure 2
states_abb<-vroom("Data/state_abbreviation.tsv")

## Figure 2A
figure2a_data<-estimates_rt_incidence |> 
  reframe(total_incidence = sum(incidence, na.rm = T), 
          total_infections = sum(infections, na.rm = T),
          .by = c("name_states", "variant"))|> 
  left_join(states_abb) |> 
  mutate(Region = factor(Region, 
                         levels = c("West", "Midwest", "South", "Northeast"))) |> 
  left_join(pop_states, by = c("name_states" = "state")) |> 
  mutate(attack_rate = (total_infections/pop))

## Saving the data to the figure
vroom_write(x = figure2a_data, 
            file = "Data/state_attack_rate_variants.csv.xz")

figure2a<- figure2a_data |> 
  ggplot(aes(x = variant, 
             col = variant,
             y = attack_rate))+    
  geom_boxplot()+
  geom_point(position = position_jitterdodge(jitter.width = .5))+
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

## Figure 2B
usmap<-tigris::states(cb = TRUE) |> 
  tigris::shift_geometry()

figure2_data<-estimates_rt_incidence |> 
  reframe(infections = sum(infections), 
          .by = c("name_states", "variant"))|>
  left_join(states_abb) |> 
  left_join(pop_states, by = c("name_states" = "state")) |> 
  mutate(attack_rate = round(infections/pop*100, 0))

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

## Figure 3
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
             group = name_states))+
  geom_hline(yintercept = 1, 
             show.legend = F, 
             col = "grey9")+
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
                                           type = "discrete", direction = 1,
                                           n = 5)))+
  scale_fill_manual(values = c("firebrick3",
                               met.brewer(palette_name = "Archambault", 
                                          type = "discrete", direction = 1,
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

## Figure 4
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

## Read directly from the saved .shp file
states_svi <- sf::st_read("Data/states_svi.shp") |> 
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

## Reading States_ar_svi_longer.csv
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

cor_mat <- data.frame(correlation::cor_test(data = states_ar_svi,
                                            x = "Omicron BA.1*",
                                            y = "SVI", 
                                            method = "Pearson"))

figure4c <-
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
                   # label.x.npc = 0,
                   # label.y.npc = 1,
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
figure4c

figure4c <- figure4c + 
  inset_element(p = SVI_states_map+
                  theme(legend.position = "none"), 
                left = 0.65, 
                bottom = 0.75, 
                right = 1, top = 1)
figure4c

figure4a <- states_ar_svi |> 
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
figure4a

figure4b <- states_ar_svi |> 
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
figure4b

library(patchwork)

figure4 <- (((figure4a) / 
               (figure4b))|
              (figure4c))+
  plot_layout(guides = 'collect',
              tag_level = "new")+
  plot_annotation(tag_levels = 'A')&
  theme(legend.position = "bottom")
figure4

ggsave(filename = "Output/Plots/figure4_manuscript.png",
       plot = figure4, 
       width = 16, 
       height = 9,
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figure4_manuscript.pdf",
       plot = figure4, 
       width = 16, 
       height = 9,
       dpi = 100)

## Supplementary Figures

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

## Figure S2
figureS2 <- variant_count |> 
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
figureS2

ggsave(filename = "Output/Plots/figureS2_manuscript.png", 
       plot = figureS2, 
       width = 16, 
       height = 9, 
       dpi = 200)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figureS2_manuscript.pdf", 
       plot = figureS2, 
       width = 16,
       height = 9, 
       dpi = 200)

## Figure S3
figureS3 <- variant_count |> 
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
figureS3

ggsave(filename = "Output/Plots/figureS3_manuscript.png", 
       plot = figureS3, 
       width = 16, 
       height = 9, 
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figureS3_manuscript.pdf", 
       plot = figureS3, 
       width = 16,
       height = 9, 
       dpi = 100)

## Figure S4
figureS4 <- variant_count |> 
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
figureS4

ggsave(filename = "Output/Plots/figureS4_manuscript.png", 
       plot = figureS4, 
       width = 16, 
       height = 9, 
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figureS4_manuscript.pdf", 
       plot = figureS4, 
       width = 16,
       height = 9, 
       dpi = 100)

## Figure S2B
figureS5_data<-estimates_rt_incidence |> 
  reframe(infections = sum(infections), 
          .by = c("name_states", "variant"))|>
  left_join(states_abb) |> 
  left_join(pop_states, by = c("name_states" = "state")) |> 
  mutate(attack_rate = round(infections/pop*100, 0))

figureS5<-figureS5_data |> 
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
figureS5

ggsave(filename = "Output/Plots/figureS5_manuscript.png", 
       plot = figureS5, 
       width = 16, 
       height = 9, 
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figureS5_manuscript.pdf", 
       plot = figureS5, 
       width = 16,
       height = 9, 
       dpi = 100)

## Figure S6 and S7 data
rt_ratio<-vroom("Output/Tables/rt_ratios.csv.xz")|>
  filter(ratios != "XBB*/BA.2*", ratios != "BA.4*/BA.5*") |>
  mutate(ratios = factor(droplevels(factor(ratios)),
                         levels = c("BA.2*/BA.1*", "BA.4*/BA.2*", "BA.5*/BA.2*",
                                    "BA.5*/BA.4*", "XBB*/BA.5*")))

figureS6_data <- rt_ratio |> 
  reframe(median = mean(median, na.rm = T),
          lower = mean(lower, na.rm = T),
          upper = mean(upper, na.rm = T), 
          .by = c(name_states, ratios)) |> 
  left_join(states_abb)

## Figure S6
figureS6 <- figureS6_data |> 
  ggplot(aes(x = ratios, y = median, group = `State Code`)) +
  geom_hline(yintercept = 1, show.legend = F, col = "grey9")+
  geom_line(aes(color = ratios), alpha = .5, size = 1, show.legend = F) +
  geom_point(aes(color = ratios), alpha = .5, size = 5, show.legend = T) +
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
figureS6

ggsave(filename = "Output/Plots/figureS6_manuscript.png", 
       plot = figureS6, 
       width = 16, 
       height = 9, 
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figureS6_manuscript.pdf", 
       plot = figureS6, 
       width = 16,
       height = 9, 
       dpi = 100)

## Figure S7
figureS7_data<-rt_ratio |> 
  reframe(median = mean(median, na.rm = T),
          lower = mean(lower, na.rm = T),
          upper = mean(upper, na.rm = T), 
          .by = c(days, ratios))

figureS7a <- figureS7_data |> 
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
figureS7a

ggsave(filename = "Output/Plots/figureS7a_manuscript.png", 
       plot = figureS7a, 
       width = 16, 
       height = 9, 
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figureS7a_manuscript.pdf", 
       plot = figureS7a, 
       width = 16,
       height = 9, 
       dpi = 100)

figure_dataS7b <- rt_ratio |> 
  reframe(median = mean(median, na.rm = T),
          lower = mean(lower, na.rm = T),
          upper = mean(upper, na.rm = T), 
          .by = c(name_states, ratios)) |> 
  mutate(name_states2 = tidytext::reorder_within(name_states, median, within = ratios))

figureS7b <- figure_dataS7b |>
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
figureS7b

ggsave(filename = "Output/Plots/figureS7b_manuscript.png", 
       plot = figureS7b, 
       width = 16, 
       height = 9, 
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figureS7b_manuscript.pdf", 
       plot = figureS7b, 
       width = 16,
       height = 9, 
       dpi = 100)

## Patchwork figure S7
patchwork_figureS7 <- (figureS7a + figureS7b)+
  plot_annotation(tag_levels = 'A')
patchwork_figureS7

ggsave(filename = "Output/Plots/figureS7_manuscript.png", 
       plot = patchwork_figureS7, 
       width = 16, 
       height = 9, 
       dpi = 100)

ggsave(filename = "~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/figureS7_manuscript.pdf", 
       plot = patchwork_figureS7, 
       width = 16,
       height = 9, 
       dpi = 100)

#