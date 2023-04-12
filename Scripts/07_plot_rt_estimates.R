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
  geom_line()+
  geom_ribbon(alpha = .5)+
  theme_minimal()+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90))+
  facet_wrap(variant_reduced~.)+
  scale_x_date(date_breaks = "4 months", 
               date_labels = "%b %y")+
  labs(x = "Date", 
       y = "Instantenous Reproduction Number \n Rt(t)")+
  scale_fill_brewer(aesthetics = c("color", "fill"), 
                    palette = "Spectral", 
                    name = "VOCs")
plot_rt

mean_rt<-test_reduce |> 
  group_by(name_states, variant_reduced) |> 
  summarise(mean_rt = mean(Rt, na.rm = T), 
            mean_lower = mean(lower, na.rm = T),
            mean_upper = mean(upper, na.rm = T))

states_division<-vroom("Data/states_abbrev_region_division.csv") |> 
  rename(name_states = State)

mean_rt<-mean_rt |> 
  left_join(states_division)

mean_rt |> 
  ggplot(aes(x = Division, 
             y = mean_rt, 
             col = variant_reduced, fill = variant_reduced))+
  geom_violin(position = position_dodge(width = 0.9), 
              scale = "width", 
              trim = FALSE, 
              alpha = .5) +
  geom_point(position = position_jitterdodge(jitter.width  = .1))+
  theme_minimal()+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90))+
  # facet_wrap(Region~., ncol = 1, strip.position = "right")+
  labs(x = "Division", 
       y = "Instantenous Reproduction Number \n Rt(t)")+
  scale_fill_brewer(aesthetics = c("color", "fill"), 
                    palette = "Accent", 
                    name = "VOCs")





