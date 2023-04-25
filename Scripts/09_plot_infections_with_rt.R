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

## Loading infections estimates
infections_estimates_daily<-vroom("Data/infections_estimates_variants_daily.csv.xz")

## Loading the frequencies of variants
variant_count<-vroom("Data/variant_counts_us.csv.xz") |> 
  rename(days = epiweek, 
         variant = voc_cdc) |> 
  filter(!variant %in% c("Alpha*", "Beta*", "Gamma*", "Delta*")) |> 
  mutate(variant = droplevels(factor(variant))) |> 
  mutate(variant = factor(variant,
                          levels = c("Omicron BA.1*", "Omicron BA.2*", 
                                     "Omicron BA.3*", "Omicron BA.4*", "Omicron BA.5*",
                                     "Omicron XBB*","Recombinant", "Other")))

## Joining Rt with variant counts
rt_estimates<-rt_estimates |> 
  full_join(variant_count) |> 
  mutate(variant = droplevels(factor(variant)))

## Making function to plot
plt_rt<-function(data, x.title, y.title, title = NULL, sec.axis.name = NULL){
  
  ## median Rt and ribbon for 95% IC
  scaleFactor<-max(data$Rt, na.rm = T)/max(data$freq, na.rm = T)
  
  rt_plot<-data |> 
    ggplot(aes(x = days, y = Rt, 
               ymin = lower, ymax = upper,
               color = variant, fill = variant))+
    geom_hline(yintercept = 1, show.legend = F, aes(col = "gray9"), alpha = .5)+
    geom_line()+
    geom_ribbon(alpha = .15)+
    geom_col(aes(x = days, y = freq*scaleFactor), 
             # col = "gray90", 
             show.legend = F, 
               position = position_stack(),
               alpha = .5)+
    theme_minimal()+
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90))+
    scale_x_date(date_breaks = "2 months", 
                 date_labels = "%b %y")+
    labs(x = x.title, 
         y = y.title,
         title = title)+
    scale_fill_brewer(aesthetics = c("color", "fill"), 
                      palette = "Paired", 
                      name = "VOCs")+
    scale_y_continuous(sec.axis=sec_axis(~./scaleFactor, name=sec.axis.name))
  
  return(rt_plot)
  
}

## Per states with facet per variant
states<-unique(rt_estimates$name_states)

variant<-unique(rt_estimates$variant)

plot_rt_list<-lapply(states, function(x){
  plt<-rt_estimates |> 
    filter(name_states == x) |>
    plt_rt(x.title = "Date", 
           y.title = "Instantenous Reproduction Number \n Rt(t)", 
           title = x, 
           sec.axis.name = "Frequency (%)")+
    facet_wrap(~variant, ncol = 1)+
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    )
  
  ggsave(filename = paste0("Output/Plots/States_rt/rt_freq/plt_rt_estimates_frequence_", x, "_daily.png"),
         plot = plt,
         width = 15,
         height = 9,
         dpi = 100)
  
  return(plt)
})

#