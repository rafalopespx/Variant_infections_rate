## CovidEstim State-level Estimates
url<-GET(paste('https://api2.covidestim.org/latest_runs?geo_type=eq.state&select=*%2Ctimeseries(*)'))

covidestim_state<-fromJSON(rawToChar(url$content))
name_states<-covidestim_state$geo_name
covidestim_state<-covidestim_state[[8]]
names(covidestim_state)<-name_states
covidestim_state<-covidestim_state |> 
  bind_rows(.id = "name_states") 


covidestim_state<-covidestim_state|> 
  select(name_states, date, infections, infections_p2_5, infections_p97_5, r_t, r_t_p2_5, r_t_p97_5) |> 
  mutate(epiweek = end.of.epiweek(as.Date(date))) |> 
  group_by(epiweek, name_states) |> 
  summarise(mean_rt = mean(r_t, na.rm = T), 
            mean_lower = mean(r_t_p2_5, na.rm = T), 
            mean_upper = mean(r_t_p97_5, na.rm = T)) |> 
  mutate(id = "CovidEstim") |> 
  filter(!is.nan(mean_rt) | !is.nan(mean_lower) | !is.nan(mean_upper))

variantestim_state<-rt_estimates |> 
  mutate(epiweek = end.of.epiweek(days)) |> 
  group_by(epiweek, name_states) |> 
  summarise(mean_rt = mean(Rt, na.rm = T), 
            mean_lower = mean(lower, na.rm = T), 
            mean_upper = mean(upper, na.rm = T)) |> 
  mutate(id = "VariantEstim") |> 
  filter(!is.nan(mean_rt) | !is.nan(mean_lower) | !is.nan(mean_upper))

comparison_df<-full_join(covidestim_state, variantestim_state)

comparison_df |> 
  filter(name_states == "Connecticut") |>
  ggplot(aes(x = epiweek, y = mean_rt, 
             ymin = mean_lower, ymax = mean_upper, 
             col = id, fill = id))+
  geom_hline(yintercept = 1,
             show.legend = F, 
             aes(col = "gray9"), alpha = .5)+
  geom_line(show.legend = T)+
  geom_ribbon(alpha = .5, 
              show.legend = F)+
  theme_minimal()+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90))+
  # facet_geo(name_states~.)+
  # facet_wrap(variant_reduced~.)+
  scale_x_date(date_breaks = "2 months", 
               date_labels = "%b %y")+
  labs(x = "Date", 
       y = "Instantenous Reproduction Number \n Rt(t)", 
       title = "Connecticut")
  
mean_rt_df|> 
  filter(name_states == "California") |>
  ggplot()+
  geom_hline(yintercept = 1,show.legend = F, aes(col = "gray9"), alpha = .5)+
  geom_line(aes(x = epiweek, y = mean_rt, 
                color = "firebrick1"),
            show.legend = F)+
  geom_ribbon(aes(x = epiweek, y = mean_rt, 
                  ymin = mean_lower, 
                  ymax = mean_upper, 
                  color = "firebrick1", 
                  fill = "firebrick1"),
              alpha = .5, 
              show.legend = F)+
  theme_minimal()+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90))+
  # facet_geo(name_states~.)+
  # facet_wrap(variant_reduced~.)+
  scale_x_date(date_breaks = "4 months", 
               date_labels = "%b %y")+
  labs(x = "Date", 
       y = "Instantenous Reproduction Number \n Rt(t)")

covidestim_state |> 
  # ungroup() |> 
  # filter(!is.nan(mean_rt) | !is.nan(mean_lower) | !is.nan(mean_upper)) |> 
  filter(name_states == "Connecticut") |>
  ggplot()+
  geom_hline(yintercept = 1,show.legend = F, aes(col = "gray9"), alpha = .5)+
  geom_line(aes(x = epiweek, y = mean_rt, 
                color = "firebrick1"),
            show.legend = F)+
  geom_ribbon(aes(x = epiweek, y = mean_rt, 
                  ymin = mean_lower, 
                  ymax = mean_upper, 
                  color = "firebrick1", 
                  fill = "firebrick1"),
              alpha = .5, 
              show.legend = F)+
  theme_minimal()+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90))+
  # facet_geo(name_states~.)+
  # facet_wrap(variant_reduced~.)+
  scale_x_date(date_breaks = "4 months", 
               date_labels = "%b %y")+
  labs(x = "Date", 
       y = "Instantenous Reproduction Number \n Rt(t)")

