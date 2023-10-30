## Cleaning the ambient
rm(list = ls())
gc()

## Loading Libraries
packs = c("tidyverse", "vroom", "jsonlite", "httr", "geofacet", "EpiEstim")
lapply(packs,require, character.only = TRUE)

## Loading functions
source("Scripts/Functions/functions.R")

## Loading the estimates
rt_estimates<-vroom("Output/Tables/rt_estimates_cori_method_daily.tsv.xz")

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
  mutate(incidence = (infections/pop)*1e5) 
# |> 
#   ## Filtering to greater incidence than 100/100k
#   filter(incidence >= 100, 
#          days <= "2023-03-01")

## Per states with facet per variant
states<-unique(rt_estimates$name_states)

## Variant vector
variant<-unique(rt_estimates$variant)

## Plotting Rt and Frequencies
## Making function to plot
plt_rt_freq<-function(data, x.title, y.title, title = NULL, sec.axis.name = NULL){
  
  ## Removing objects
  rm(scaleFactor)
  
  ## Scalling factor
  scaleFactor<-max(data$upper, na.rm = T)/max(data$freq, na.rm = T)
  
  rt_plot<-data |>
    ggplot(aes(x = days, y = Rt, 
               # ymin = lower, ymax = upper,
               color = variant, fill = variant))+
    geom_hline(yintercept = 1, show.legend = F, aes(col = "gray9"), alpha = .5)+
    geom_line()+
    geom_ribbon(aes(ymin = lower, 
                    ymax = upper),
                alpha = .15)+
    geom_col(aes(x = days, y = freq*scaleFactor), 
             col = NA,
             show.legend = F, 
             width = 5,
             # position = position_stack(),
             alpha = .15)+
    theme_minimal()+
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90))+
    scale_x_date(date_breaks = "1 month", 
                 date_labels = "%b %y")+
    labs(x = x.title, 
         y = y.title,
         title = title)+
    colorspace::scale_fill_discrete_divergingx(name = "VOCs",
                                               palette = "Zissou1",
                                               aesthetics = c("color", "fill"))+
    scale_y_continuous(sec.axis=sec_axis(~./scaleFactor, 
                                         name=sec.axis.name))+
    facet_wrap(~variant, ncol = 1)+
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    )
  # rt_plot
  
  return(rt_plot)
  
}


plot_rt_freq_list<-lapply(states, function(x){
  plt<-plt_rt_freq(data = estimates_rt_incidence |> 
                     filter(name_states == x) |>
                     filter(variant %in% c("Omicron BA.1*", "Omicron BA.2*", "Omicron BA.4*", "Omicron BA.5*", 
                                           "Omicron XBB*"), 
                            days >= "2021-12-01"),
                   x.title = "Date", 
                   y.title = "Instantenous Reproduction Number \n Rt(t)", 
                   title = x, 
                   sec.axis.name = "Frequency \n (%)")
  
  # ggsave(filename = paste0("Output/Plots/States_rt/rt_freq/plt_rt_estimates_frequence_", x, "_daily.png"),
  #        plot = plt,
  #        width = 15,
  #        height = 9,
  #        dpi = 100)
  
  return(plt)
})

names(plot_rt_freq_list)<-states

plot_rt_freq_list$California
plot_rt_freq_list$`New York`
plot_rt_freq_list$Connecticut

## Plotting Rt and Infections
## Joining Rt with variant counts
estimates_rt_incidence<-rt_estimates |>
  left_join(variant_count) |>
  mutate(variant = droplevels(factor(variant))) |> 
  left_join(pop_states, by = c("name_states" = "state")) |> 
  rename(infections = I) |> 
  mutate(incidence = (infections/pop)*10^5)

## Making function to plot
plt_rt_infec<-function(data, x.title, y.title, title = NULL, sec.axis.name = NULL){
  
  ## Removing objects
  rm(scaleFactor)
  
  ## Scalling factor
  scaleFactor<-max(data$upper, na.rm = T)/max(data$incidence, na.rm = T)
  
  rt_plot<-data |>
    ggplot(aes(x = days, y = Rt, 
               color = variant, fill = variant))+
    geom_hline(yintercept = 1, show.legend = F, aes(col = "gray9"), alpha = .5)+
    geom_line()+
    geom_ribbon(aes(ymin = lower, 
                    ymax = upper),
                alpha = .15)+
    geom_col(aes(x = days, y = incidence*scaleFactor), 
             col = NA,
             show.legend = F,
             width = 1,
             alpha = .15)+
    theme_minimal()+
    scale_x_date(date_breaks = "1 month", 
                 date_labels = "%b %y")+
    labs(x = x.title, 
         y = y.title,
         title = title)+
    colorspace::scale_fill_discrete_divergingx(name = "VOCs",
                                               palette = "Zissou1",
                                               aesthetics = c("color", "fill"))+
    scale_y_continuous(sec.axis=sec_axis(~./scaleFactor, 
                                         name=sec.axis.name))+
    facet_wrap(~variant, ncol = 1, scales = "free_y")+
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90),
          strip.background = element_blank(),
          strip.text.x = element_blank()
    )
  # rt_plot
  
  return(rt_plot)
  
}


plot_rt_infec_list<-lapply(states, function(x){
  plt<-plt_rt_infec(data = estimates_rt_incidence |> 
                      filter(name_states == x) |>
                      filter(variant %in% c("Omicron BA.1*", "Omicron BA.2*", "Omicron BA.4*", "Omicron BA.5*", 
                                            "Omicron XBB*"), 
                             days >= "2021-12-01"),
                    x.title = "Date", 
                    y.title = "Instantenous Reproduction Number \n Rt(t)", 
                    title = x, 
                    sec.axis.name = "Infections per 100k \n people per day")
  
  # ggsave(filename = paste0("Output/Plots/States_rt/rt_infec/plt_rt_estimates_infections_", x, "_daily.png"),
  #        plot = plt,
  #        width = 15,
  #        height = 9,
  #        dpi = 100)
  
  return(plt)
})

names(plot_rt_infec_list)<-states

plot_rt_infec_list$California
plot_rt_infec_list$`New York`
plot_rt_infec_list$Connecticut
plot_rt_infec_list$Texas

#