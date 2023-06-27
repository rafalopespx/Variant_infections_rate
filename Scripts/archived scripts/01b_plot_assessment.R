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

freq_fun<-function(x){
  table<- metadata |> 
    filter(voc_cdc == x, 
           !is.na(pango_lineage))|> 
    select(voc_cdc, pango_lineage) |> 
    summarise(n = n(), 
              .by = c(voc_cdc, pango_lineage)) |> 
    mutate(freq = round(100*n/sum(n),2))
}

variants<-unique(metadata$voc_cdc)

table_list<-sapply(variants, freq_fun)

table_list |> 
  filter(voc_cdc == "Omicron BA.5*") |> 
  pull(var = pango_lineage)
