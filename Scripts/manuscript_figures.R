## Cleaning the ambient
rm(list = ls())
gc()

## Loading Libraries
packs = c("tidyverse", "vroom", "patchwork", "geofacet")
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
  ## Renaming infections
  rename(infections = I) |> 
  ## Creating incidence per 100k
  mutate(incidence = (infections/pop)*1e5, 
         percentual_incidence = round(incidence/pop, 2))
# |>
#   ## Filtering to greater incidence than 100/100k
#   filter(incidence >= 100)


## Figure1
figure1a <- rt_estimates |> 
  ## Renaming infections
  rename(infections = I) |>
  group_by(name_states, days) |> 
  summarise(infections = sum(infections)) |> 
  left_join(pop_states, by = c("name_states" = "state")) |>
  ## Creating incidence per 100k
  mutate(incidence = (infections/pop)*1e5, 
         percentual_incidence = round(incidence/pop, 2)) 

figure1a<-estimates_rt_incidence |> 
  filter(days <= "2023-04-01", incidence >= 10) |>
  reframe(infections = sum(infections, na.rm = T), 
          .by = c(name_states, days)) |> 
  ggplot(aes(x = days, y = infections, 
             group = name_states))+
  geom_line(alpha = .1)+
  geom_line(data = estimates_rt_incidence |>
              filter(days <= "2023-04-01", incidence >= 10) |>
              reframe(mean_infections = mean(infections, na.rm = T),
                      .by = c(days)),
            aes(x = days, y = mean_infections, 
                col = "Mean infections per day"), 
            inherit.aes = F)+
  labs(x = NULL, y = "Infections per days")+
  theme_minimal()+
  theme(legend.position = "bottom", 
        axis.text.x = element_blank())+
  scale_x_date(date_breaks = "2 months", date_labels = "%b %y")+
  scale_y_continuous(labels = scales::label_scientific())+
  scale_color_manual(name = "", values = "firebrick2")
figure1a

figure1b <- estimates_rt_incidence |> 
  filter(days <= "2023-04-01", incidence >= 10) |>
  reframe(infections = sum(infections, na.rm = T), 
          .by = c(days, name_states, variant)) |> 
  ggplot(aes(x = days, y = infections, 
             color = variant,
             group = name_states))+
  geom_line(alpha = .1, show.legend = F)+
  geom_line(data = estimates_rt_incidence |>
              filter(days <= "2023-04-01", incidence >= 10) |>
              reframe(mean_infections = mean(infections, na.rm = T),
                      .by = c(days, variant)),
            aes(x = days, y = mean_infections, 
                col = "Mean infections per day"),
            show.legend = F,
            inherit.aes = F)+
  theme_minimal()+
  facet_wrap(variant~., ncol = 1, scales = "free_y", strip.position = "right")+
  theme(legend.position = "none")+
  labs(x = "Date", y = "Infections per day")+
  scale_x_date(date_breaks = "2 months", date_labels = "%b %y")+
  scale_y_continuous(labels = scales::label_scientific())+
  scale_color_manual(values = c("firebrick1", MetBrewer::met.brewer(palette_name = "Archambault")))
figure1b

figure1_patchwork<-(figure1a / figure1b)+
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides = "collect", heights = c(1.2,1.6))&
  theme(legend.position = "bottom")
figure1_patchwork

ggsave(filename = "Output/Plots/figure1_manuscript.png", 
       plot = figure1_patchwork, 
       width = 11,
       height = 9, 
       dpi = 100)

states_abb<-vroom("Data/states_abbrev_region_division.csv")

figure2a<-estimates_rt_incidence |> 
  filter(days <= "2023-04-01") |> 
  reframe(total_incidence = sum(incidence, na.rm = T), 
          .by = c("name_states", "variant"))|> 
  left_join(states_abb, by = c("name_states" = "State")) |> 
  mutate(Region = factor(Region, 
                         levels = c("West", "Midwest", "South", "Northeast"))) |> 
  ggplot(aes(x = variant, 
             col = variant,
             y = total_incidence))+    
  geom_boxplot()+
  geom_point(position = position_jitterdodge(jitter.width = .9))+
  # stat_summary(fun.y = mean, 
  #              geom="text", 
  #              vjust=-20, 
  #              aes(label = round(..y.., digits = 0)))+
  theme_minimal()+
  MetBrewer::scale_color_met_d(palette_name = "Archambault", name = "")+
  theme(legend.position = "none")+
  labs(x = "", y = "Incidence \n infections/100K")+
  scale_y_continuous(labels = scales::label_scientific())
figure2a

figure2_data<-estimates_rt_incidence |> 
  filter(days <= "2023-04-01") |> 
  reframe(infections = sum(infections), 
          .by = c("name_states", "variant"))|>
  left_join(states_abb, by = c("name_states" = "State")) |> 
  left_join(pop_states, by = c("name_states" = "state")) |> 
  mutate(attack_rate = round(infections/pop*100, 1))

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
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "bottom", 
        panel.background = element_blank(), 
        panel.grid = element_blank())
figure2b  

figure2_patchwork<-(figure2a / figure2b)+
  plot_layout(heights = c(1.2, 1.8))+
  plot_annotation(tag_levels = 'a')
figure2_patchwork

ggsave(filename = "Output/Plots/figure2_manuscript.png", 
       plot = figure2_patchwork, 
       width = 11, 
       height = 9, 
       dpi = 100)

usmap<-tigris::states() |> 
  shift_geometry()

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
            lwd = .1)+
    theme_void()+
    scale_fill_manual(values = met.brewer(palette_name = "Hiroshige", 
                                          n = 13, 
                                          type = "continuous", 
                                          direction = -1), 
                      drop = F, 
                      name = "Attack rate \n (%) of pop. ever infected",
                      guide = guide_bins(keywidth = grid::unit(3, "cm"),
                                         title.position = "top", 
                                         show.limits = T), 
                      aesthetics = c("color", "fill"))+
    theme(legend.position = "bottom")+
    no_axis+
    labs(subtitle = x)
  # plt
  return(plt)
}

variants<- as.character(unique(map_attackrate$variant))

map_list<-lapply(variants, attack_rate_map)

patchwork_map<-((map_list[[1]] / map_list[[4]]) | (map_list[[2]]/ map_list[[3]] / map_list[[5]]))+
  plot_layout(guides = "collect")+
  plot_annotation(tag_levels = 'A')&
  theme(legend.position = "bottom")
patchwork_map

ggsave(filename = "Output/Plots/attack_rate_map.png", 
       plot = patchwork_map, 
       width = 11,
       height = 9, 
       dpi = 100)

figure3_data<- estimates_rt_incidence |> 
  filter(incidence >= 15) |>
  reframe(median = mean(Rt, na.rm = T),
          lower = mean(lower, na.rm = T),
          upper = mean(upper, na.rm = T),
          .by = c(days, variant, name_states))

figure3<- figure3_data |> 
  ggplot(aes(group = name_states))+
  geom_hline(yintercept = 1, show.legend = F, col = "grey9")+
  geom_line(aes(x = days, y = median, col = variant),
            alpha = .1)+
  geom_ribbon(aes(x = days, ymin = lower, ymax = upper, 
                  fill = variant),
              alpha = .01, show.legend = F)+
  theme_minimal()+
  labs(x = "Date", y = "Effective Reproduction \n (Rt)")+
  theme(legend.title = element_blank(), 
        legend.position = "bottom")+
  facet_wrap(variant~., ncol = 1, strip.position = "right")+
  scale_x_date(date_breaks = "2 months",
               date_labels = "%b %y"
  )+
  MetBrewer::scale_color_met_d(palette_name = "Archambault", name = "")+
  MetBrewer::scale_fill_met_d(palette_name = "Archambault", name = "")
figure3

ggsave(filename = "Output/Plots/figure3_manuscript.pdf", 
       plot = figure3, 
       width = 11, 
       height = 9, 
       dpi = 100)


rt_ratio<-vroom("Output/Tables/rt_ratios.csv.xz")|> 
  mutate(ratios = factor(ratios, 
                         levels = c("BA.2*/BA.1*", "BA.4*/BA.2*", "BA.5*/BA.2*", 
                                    "BA.4*/BA.5*", "XBB*/BA.5*")))

figure4a_data<-rt_ratio |> 
  reframe(median = mean(median, na.rm = T),
          lower = mean(lower, na.rm = T),
          upper = mean(upper, na.rm = T), 
          .by = c(days, ratios))
# |> 
#   filter(ratios %in% c("BA.2*/BA.1*", "BA.5*/BA.2*", "XBB*/BA.5*")) |> 
#   mutate(ratios = droplevels(factor(ratios)))

figure4a <- figure4a_data |> 
  ggplot(aes(x = days, y = median, 
             ymin = lower, ymax = upper,
             color = ratios, fill = ratios))+
  geom_hline(yintercept = 1, show.legend = F, col = "grey9")+
  geom_line()+
  geom_ribbon(alpha = .05)+
  theme_minimal()+
  facet_wrap(ratios~., ncol = 1, strip.position = "right")+
  theme(legend.position = "bottom")+
  scale_x_date(date_breaks = "2 months", date_labels = "%b %y")+
  labs(y = "Average (Rt) ratios", x = "Date")+
  MetBrewer::scale_color_met_d(palette_name = "Austria", name = "")+
  MetBrewer::scale_fill_met_d(palette_name = "Austria", name = "")
figure4a


