test<-bind_rows(rt_list)
test_reduce<-bind_rows(rt_list) 

plot_rt_states<-test_reduce |> 
  filter(variant_reduced != "Other") |>
  ggplot(aes(x = days, y = Rt,
             ymin = lower, ymax = upper, 
             col = variant_reduced, fill = variant_reduced))+
  geom_line()+
  geom_ribbon(alpha = .5)+
  theme_minimal()+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90))+
  facet_geo(name_states~.)+
  # facet_wrap(variant_reduced~.)+
  scale_x_date(date_breaks = "4 months", 
               date_labels = "%b %y")+
  labs(x = "Date", 
       y = "Instantenous Reproduction Number \n Rt(t)")+
  scale_fill_brewer(aesthetics = c("color", "fill"), 
                    palette = "Spectral", 
                    name = "VOCs")
plot_rt_states

plot_rt<-test_reduce |> 
  filter(name_states == "Connecticut") |> 
  ggplot(aes(x = days, y = Rt,
             ymin = lower, ymax = upper,
             col = variant_reduced, fill = variant_reduced))+
  geom_hline(yintercept = 1,show.legend = F, aes(col = "gray9"), alpha = .5)+
  geom_line()+
  geom_ribbon(alpha = .5)+
  theme_minimal()+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90))+
  # facet_wrap(variant_reduced~.)+
  scale_x_date(date_breaks = "4 months", 
               date_labels = "%b %y")+
  labs(x = "Date", 
       y = "Instantenous Reproduction Number \n Rt(t)")+
  scale_fill_brewer(aesthetics = c("color", "fill"), 
                    palette = "Spectral", 
                    name = "VOCs")
plot_rt

test_split<-test_reduce %>% 
  split(list(.$name_states, .$variant_reduced))

mean_rt<-function(x, first_something_days){
  x<-x |> 
    filter(!is.na(Rt))
  
  if(!missing(first_something_days)){
    first_something_days<-as.numeric(first_something_days)
    x<-x[1:first_something_days,]
  }
  
  x<-x |> 
    group_by(name_states, variant_reduced) |> 
    summarise(mean_rt = mean(Rt, na.rm = T), 
              mean_lower = mean(lower, na.rm = T),
              mean_upper = mean(upper, na.rm = T))
  
  return(x)
}

## Calculating averages over time periods

## States names
states_division<-vroom("Data/states_abbrev_region_division.csv") |> 
  rename(name_states = State)

mean_rt_30<-lapply(test_split, mean_rt, 30) |> 
  bind_rows() |> 
  left_join(states_division)

mean_rt_60<-lapply(test_split, mean_rt, 60) |> 
  bind_rows()|> 
  left_join(states_division)

mean_rt_90<-lapply(test_split, mean_rt, 90) |> 
  bind_rows()|> 
  left_join(states_division)

mean_rt_120<-lapply(test_split, mean_rt, 120) |> 
  bind_rows()|> 
  left_join(states_division)

## Plotting
plt_mean_rt_30<-mean_rt_30 |> 
  filter(!is.na(Division)) |> 
  ggplot(aes(x = Division, 
             y = mean_rt, 
             col = variant_reduced, fill = variant_reduced))+
  geom_hline(yintercept = 1,show.legend = F, aes(col = "gray9"), alpha = .5)+
  geom_violin(position = position_dodge(width = .9),
              scale = "width",
              trim = FALSE,
              alpha = .5) +
  geom_point(position = position_jitterdodge(jitter.width = .1))+
  theme_minimal()+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, size = 9))+
  labs(x = "Division", 
       y = "Mean Reproduction Number \n Rt(t)", 
       subtitle = "Mean over 30 first days")+
  scale_fill_brewer(aesthetics = c("color", "fill"), 
                    palette = "Spectral", 
                    name = "VOCs")
plt_mean_rt_30

plt_mean_rt_60<-mean_rt_60 |> 
  filter(!is.na(Division)) |> 
  ggplot(aes(x = Division, 
             y = mean_rt, 
             col = variant_reduced, fill = variant_reduced))+
  geom_hline(yintercept = 1,show.legend = F, aes(col = "gray9"), alpha = .5)+
  geom_violin(position = position_dodge(width = 0.9), 
              scale = "width", 
              trim = FALSE, 
              alpha = .5) +
  geom_point(position = position_jitterdodge(jitter.width  = .1))+
  theme_minimal()+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, size = 9))+
  labs(x = "Division", 
       y = "Mean Reproduction Number \n Rt(t)", 
       subtitle = "Mean over 60 first days")+
  scale_fill_brewer(aesthetics = c("color", "fill"), 
                    palette = "Spectral", 
                    name = "VOCs")
plt_mean_rt_60

plt_mean_rt_90<-mean_rt_90 |> 
  filter(!is.na(Division)) |> 
  ggplot(aes(x = Division, 
             y = mean_rt, 
             col = variant_reduced, fill = variant_reduced))+
  geom_hline(yintercept = 1,show.legend = F, aes(col = "gray9"), alpha = .5)+
  geom_violin(position = position_dodge(width = 0.9), 
              scale = "width", 
              trim = FALSE, 
              alpha = .5) +
  geom_point(position = position_jitterdodge(jitter.width  = .1))+
  theme_minimal()+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, size = 9))+
  labs(x = "Division", 
       y = "Mean Reproduction Number \n Rt(t)", 
       subtitle = "Mean over 90 first days")+
  scale_fill_brewer(aesthetics = c("color", "fill"), 
                    palette = "Spectral", 
                    name = "VOCs")
plt_mean_rt_90

plt_mean_rt_120<-mean_rt_120 |> 
  filter(!is.na(Division)) |> 
  ggplot(aes(x = Division, 
             y = mean_rt, 
             col = variant_reduced, fill = variant_reduced))+
  geom_hline(yintercept = 1,
             show.legend = F, aes(col = "gray9"), alpha = .5)+
  geom_violin(position = position_dodge(width = 0.9), 
              scale = "width", 
              trim = FALSE, 
              alpha = .5) +
  geom_point(position = position_jitterdodge(jitter.width  = .1))+
  theme_minimal()+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, size = 9))+
  labs(x = "Division", 
       y = "Mean Reproduction Number \n Rt(t)", 
       subtitle = "Mean over 120 first days")+
  scale_fill_brewer(aesthetics = c("color", "fill"), 
                    palette = "Spectral", 
                    name = "VOCs")
plt_mean_rt_120


plt_mean_rt_30<-mean_rt_30 |> 
  filter(!is.na(Division)) |> 
  ggplot(aes(x = Division, 
             y = mean_rt, 
             col = variant_reduced, fill = variant_reduced))+
  geom_hline(yintercept = 1,show.legend = F, aes(col = "gray9"), alpha = .5)+
  geom_violin(position = position_dodge(width = 0.9), 
              scale = "width", 
              trim = FALSE, 
              alpha = .5) +
  geom_point(position = position_jitterdodge(jitter.width  = .1))+
  theme_minimal()+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, size = 9))+
  labs(x = "Division", 
       y = "Mean Reproduction Number \n Rt(t)", 
       subtitle = "Mean over 30 first days")+
  scale_fill_brewer(aesthetics = c("color", "fill"), 
                    palette = "Spectral", 
                    name = "VOCs")
plt_mean_rt_30

states_division<-vroom("Data/states_abbrev_region_division.csv") |> 
  rename(name_states = State)

mean_rt_df<-rt_estimates |> 
  mutate(week = end.of.epiweek(days)) |> 
  group_by(week, name_states) |> 
  summarise(mean_rt = mean(Rt, na.rm = T), 
            mean_lower = mean(lower, na.rm = T), 
            mean_upper = mean(upper, na.rm = T)) |> 
  filter(!is.nan(mean_rt) | !is.nan(mean_lower) | !is.nan(mean_upper)) |> 
  filter(name_states == "Connecticut") |>
  ggplot()+
  geom_hline(yintercept = 1,show.legend = F, aes(col = "gray9"), alpha = .5)+
  geom_line(aes(x = week, y = mean_rt, 
                color = "firebrick1"),
            show.legend = F)+
  geom_ribbon(aes(x = week, y = mean_rt, 
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
mean_rt_df




