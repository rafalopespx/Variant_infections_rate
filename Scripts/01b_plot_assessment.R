## Cleaning workspace
rm(list = ls())
gc()

## Loading packages
packs<-c("tidyverse", "vroom")
lapply(packs, require, character.only = TRUE)

## Loading functions
source("Scripts/Functions/functions.R")

## Filtered metadata
metadata<-vroom("Data/metadata_us.csv.xz")

metadata |> 
  filter(division == "Connecticut", 
         voc_cdc == "Other") |> 
  ## Grouping by epiweek, state and voc_cdc
  summarise(n=n(), .by = c(epiweek, pango_lineage)) |> 
  mutate(freq = round(100*n/sum(n),2)) |> 
  ggplot(aes(x = epiweek, y = n, col = pango_lineage, fill = pango_lineage))+
  geom_col(position = position_stack())
