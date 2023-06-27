## Cleaning the workspace
rm(list = ls())
gc()

## Loading packages
packs<-c("tidyverse", "vroom", "geofacet", "colorspace")
lapply(packs, require, character.only = TRUE)

# Remove plot axis
no_axis <- theme(axis.title=element_blank(),
                 axis.text=element_blank(),
                 axis.ticks=element_blank())

## Loading predictions
pred_states_long<-vroom("Output/pred_multinomial_us_states.csv.xz") 
# |> 
#   mutate(voc_cdc = factor(voc_cdc,
#                           levels = c("Omicron BA.1*", "Omicron BA.2*", "Omicron BA.2.75*",
#                                      "Omicron BA.3*", "Omicron BA.4*",
#                                      "Omicron BA.5*", "Omicron BQ.1*", "Omicron BJ.1*",
#                                      "XBB.1*", "XBB.1.5*", "Recombinant",
#                                      "Other"))) 

pred_us_long<-vroom("Output/pred_multinomial_us.csv.xz") 
# |> 
#   mutate(voc_cdc = factor(voc_cdc,
#                           levels = c("Omicron BA.1*", "Omicron BA.2*", "Omicron BA.2.75*",
#                                      "Omicron BA.3*", "Omicron BA.4*",
#                                      "Omicron BA.5*", "Omicron BQ.1*", "Omicron BJ.1*",
#                                      "XBB.1*", "XBB.1.5*", "Recombinant",
#                                      "Other"))) 

## Plotting
### US Whole country
pred_states_long |> 
  # filter(epiweek >= "2022-01-01") |> 
  #data
  ggplot(
    #mapa
    aes(x = epiweek, y = median/51, 
        # ymax = upper, ymin = lower,
        col = voc_cdc, 
        fill = voc_cdc))+
  #geometria
  geom_col()+
  #formatações
  theme_minimal()+
  # facet_geo(name_states~., grid = "us_state_without_DC_grid3")+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90))+
  scale_fill_brewer(name = "VOCs", 
                    palette = "Paired", 
                    aesthetics = c("color", "fill"))+
  scale_x_date(date_breaks = "3 months", 
               date_labels = "%b %y")+
  labs(x = "Week of Collect", 
       y = "Frequency", 
       title = "Nowcast of Dominance of Variants",
       subtitle = "y ~ Multinomial(n, prob), \n g(Prob) = a+b*epiweek+c*state \n g is log",
       caption = "Source: GISAID")->
  plot_us_predictions
plot_us_predictions

### US States
pred_states_long |> 
  # filter(epiweek >= "2022-01-01") |> 
  #data
  ggplot(
    #mapa
    aes(x = epiweek, y = median, 
        ymax = upper, ymin = lower,
        col = voc_cdc, 
        fill = voc_cdc))+
  #geometria
  geom_col()+
  #formatações
  theme_minimal()+
  facet_geo(name_states~., grid = "us_state_grid1")+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90))+
  scale_fill_brewer(name = "VOCs", 
                    palette = "Paired", 
                    aesthetics = c("color", "fill"))+
  scale_x_date(date_breaks = "3 months", 
               date_labels = "%b %y")+
  labs(x = "Week of Collect", 
       y = "Frequency", 
       title = "Nowcast of Dominance of Variants",
       # subtitle = "Multinomial Logistic Regression", 
       caption = "Source: GISAID")->
  plot_states_predictions
plot_states_predictions


## CT
plot_ct_predictions<-
  pred_states_long |> 
  filter(name_states == "Connecticut") |> 
  #data
  ggplot(
    #mapa
    aes(x = epiweek, y = median, 
        ymax = upper, ymin = lower,
        col = voc_cdc, 
        fill = voc_cdc))+
  #geometria
  geom_col()+
  #formatações
  theme_minimal()+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90))+
  scale_fill_brewer(name = "VOCs", 
                    palette = "Paired", 
                    aesthetics = c("color", "fill"))+
  scale_x_date(date_breaks = "2 months", 
               date_labels = "%b %y")+
  labs(x = "Week of Collect", 
       y = "Frequency", 
       title = "Nowcast of Dominance of Variants, state: Connecticut",
       # subtitle = "Multinomial Logistic Regression", 
       caption = "Source: GISAID")
plot_ct_predictions
