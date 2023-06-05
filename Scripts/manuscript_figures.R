## Cleaning the ambient
rm(list = ls())
gc()

## Loading Libraries
packs = c("tidyverse", "vroom", "patchwork")
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
  group_by(name_states, days) |> 
  summarise(infections = sum(infections, na.rm = T), 
            incidence = sum(incidence, na.rm = T)) |> 
  filter(days <= "2023-04-01") |> 
  ## Filtering to greater incidence than 100/100k
  filter(incidence >= 100) |> 
  ggplot(aes(x = days, y = incidence, 
             group = name_states))+
  geom_line(alpha = .1)+
  scale_color_viridis_d(option = "rocket")+
  labs(x = NULL, y = "Incidence \n infections/100K")+
  theme_minimal()+
  theme(legend.position = "none", 
        axis.text.x = element_blank())+
  scale_x_date(date_breaks = "2 months", date_labels = "%b %y")+
  scale_y_continuous(labels = scales::label_scientific())
figure1a

figure1b <- estimates_rt_incidence |> 
  filter(days <= "2023-04-01") |>  
  ## Filtering to greater incidence than 100/100k
  filter(incidence >= 200) |> 
  ggplot(aes(x = days, y = infections, 
             group = name_states))+
  geom_line(alpha = .1)+
  theme_minimal()+
  facet_wrap(variant~., ncol = 1, strip.position = "right", scales = "free_y")+
  theme(legend.position = "none")+
  labs(x = "Date", y = "Incidence \n infections/100K")+
  scale_x_date(date_breaks = "2 months", date_labels = "%b %y")+
  scale_y_continuous(labels = scales::label_scientific())
figure1b

figure1_patchwork<-(figure1a / figure1b)+
  patchwork::plot_annotation(title = "Incidence break by variant")
figure1_patchwork

states_abb<-vroom("Data/states_abbrev_region_division.csv")


figure2a<-estimates_rt_incidence |> 
  filter(days <= "2023-04-01") |> 
  reframe(total_incidence = sum(incidence, na.rm = T), 
          .by = c("name_states", "variant"))|> 
  left_join(states_abb, by = c("name_states" = "State")) |> 
  ggplot(aes(x = variant, 
             y = total_incidence, 
             col = Region))+    
  ggrepel::geom_text_repel(aes(label = `State Code`),
                           position = position_jitterdodge(jitter.width  = .5, 
                                                           dodge.width = .9))+
  geom_point(position = position_jitter(width = .9))+
  theme_minimal()+
  scale_color_viridis_d(option = "turbo")+
  theme(legend.position = "none")+
  labs(x = "VOCs", y = "Incidence \n infections/100K")+
  scale_y_continuous(labels = scales::label_scientific())
figure2a

figure2b<-estimates_rt_incidence |> 
  filter(days <= "2023-04-01") |> 
  reframe(total_incidence = sum(incidence, na.rm = T), 
          .by = c("name_states", "variant"))|>
  mutate(variant2 = tidytext::reorder_within(variant, -total_incidence, within = name_states)) |> 
  ggplot(aes(x = variant2, 
             y = total_incidence, 
             fill = variant))+
  geom_col()+
  geom_text(aes(label = round(total_incidence/1e5, 1)),
            size = 2.5,
            vjust = -0.5)+
  lims(y = c(NA,6.5e6))+
  scale_fill_discrete(name = "VOCs")+
  facet_geo(name_states~., scales = "free")+
  labs(y = "Attack rate \n (%) of pop. ever infected")+
  theme_minimal()+
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "bottom")+
  tidytext::scale_x_reordered()
figure2b  

figure2_patchwork<-(figure2a | figure2b)
figure2_patchwork

figure3<-estimates_rt_incidence |> 
  ggplot(aes(x = days, y = Rt, 
             ymin = lower, ymax = upper,
             col = variant, fill = variant))+
  geom_hline(yintercept = 1, show.legend = F, col = "grey9")+
  geom_line()+
  geom_ribbon(alpha = .05)+
  theme_minimal()+
  labs(x = "Date", y = "Median Rt")+
  theme(legend.title = element_blank(), 
        legend.position = "bottom")+
  facet_wrap(variant~., ncol = 1, strip.position = "right")+
  scale_x_date(date_breaks = "2 months",
               date_labels = "%b %y"
  )
figure3
