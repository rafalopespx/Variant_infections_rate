## Cleaning the workspace
rm(list = ls())
gc()

## Loading packages
packs<-c("tidyverse", "vroom", "geofacet", "colorspace")
lapply(packs, require, character.only = TRUE)

## Loading the filtered data for the US
variant_count<-vroom("Data/variant_counts_us.csv.xz")

# Remove plot axis
no_axis <- theme(axis.title=element_blank(),
                 axis.text=element_blank(),
                 axis.ticks=element_blank())

## Plotting US counts
plt_us_variant_counts<-variant_count |>
  filter(epiweek >= "2021-09-01", epiweek <= "2023-03-01") |> 
  ggplot(aes(x = epiweek, y = n, fill = voc_cdc))+
  geom_col(width = 7)+
  theme_minimal()+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90))+
  scale_x_date(name = "Date of the end of week", 
               date_breaks = "2 months", 
               date_labels = "%b %Y")+
  labs(y = "Genomes counts \n per week", 
       title = "Raw data for US",
       caption = "Source data from: GISAID")+
  scale_fill_brewer(name = "VOCs", 
                    palette = "Paired")
plt_us_variant_counts

## Saving the plot
ggsave(plot = plt_us_variant_counts, 
       filename = "Output/plot_us_raw_count_genomes.png", 
       width = 11, 
       height = 9, 
       dpi = 100)

## Plotting States counts
plt_variant_counts<- variant_count |>
  filter(epiweek >= "2021-09-01", epiweek <= "2023-03-01") |>  
  ggplot(aes(x = epiweek, y = n, fill = voc_cdc))+
  geom_col(width = 7)+
  facet_geo(name_states~., grid = "us_state_grid1", scales = "free_y")+
  theme_minimal()+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90))+
  scale_x_date(name = "Date of the end of week", 
               date_breaks = "2 months", 
               date_labels = "%b %Y")+
  labs(y = "Genomes counts \n per week", 
       title = "Raw data from GISAID genomes", 
       subtitle = "per state", 
       caption = "Source data from: GISAID")+
  colorspace::scale_fill_discrete_divergingx(name = "VOCs", 
                                             palette = "Spectral")
plt_variant_counts

## ggsave the plots
ggsave(plot = plt_variant_counts, 
       filename = "Output/plot_raw_count_genomes_us.png", 
       width = 11, 
       height = 9, 
       dpi = 100)

## Plotting States percentage share
plt_variant_freq<- variant_count |>
  filter(epiweek >= "2021-09-01", epiweek <= "2023-03-01") |>  
  ggplot(aes(x = epiweek, y = n, fill = voc_cdc))+
  geom_col(position = position_fill(),
           width = 7)+
  facet_geo(name_states~., grid = "us_state_grid1", scales = "free_y")+
  theme_minimal()+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90))+
  scale_x_date(name = "Date of the end of week", 
               date_breaks = "2 months", 
               date_labels = "%b %Y")+
  labs(y = "Genomes counts \n per week", 
       title = "Raw data from GISAID genomes", 
       subtitle = "per state", 
       caption = "Source data from: GISAID")+
  colorspace::scale_fill_discrete_divergingx(name = "VOCs", 
                                             palette = "Spectral")
plt_variant_freq

## ggsave the plots
ggsave(plot = plt_variant_freq, 
       filename = "Output/plot_raw_freq_genomes_us.png", 
       width = 11, 
       height = 9, 
       dpi = 100)


## Exploratory analysis for CT
# plot_ct<-variant_count |> 
#   filter(name_states == "Connecticut", 
#          epiweek >= "2021-09-01", 
#          epiweek <= "2023-03-01")|>  
#   ggplot(aes(x = epiweek, y = n, fill = voc_cdc))+
#   geom_col(position = position_fill(),
#            width = 7) +
#   theme_minimal()+
#   theme(legend.position = "bottom", 
#         axis.text.x = element_text(angle = 90))+
#   scale_x_date(name = "Date of the end of week", 
#                date_breaks = "2 months", 
#                date_labels = "%b %Y")+
#   labs(y = "Genomes counts \n per week", 
#        title = "Raw data from GISAID genomes", 
#        subtitle = "Connecticut", 
#        caption = "Source data from: GISAID")+
#   colorspace::scale_fill_discrete_divergingx(name = "VOCs", 
#                                              palette = "Spectral")
# plot_ct

#