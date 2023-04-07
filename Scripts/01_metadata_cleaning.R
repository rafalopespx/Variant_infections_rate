## Cleaning workspace
rm(list = ls())
gc()

## Loading packages
packs<-c("tidyverse", "vroom")
lapply(packs, require, character.only = TRUE)

## Loading functions
source("Scripts/functions.R")

## Loading the metadata and states
metadata<-vroom("Data/metadata_2023-03-07/metadata.tsv")

## States abbreviation
states<-vroom("Data/state_abbreviation.tsv")

## Filtering for the US, human only infections, and formatting date
metadata<-metadata |> 
  filter(country == "USA", 
         host == "Human") |> 
  mutate(date = as.Date(date, "%Y-%m-%d"), 
         ## Fixing DC
         division = case_when(division == "Washington DC" ~ "District of Columbia", 
                              TRUE ~ division)) |> 
  ## Filtering out any state name that is not matching 49 contiguous states plus Alaska and Hawaii
  filter(division %in% states$name_state)

## Formatting metadata in a more pratical form
metadata <- metadata |> 
  ## Dummy date variable
  mutate(copy_date = as.numeric(date)) |> 
  ## Choosen lineages to be analyzed
  mutate(voc_cdc = case_when(grepl(x = pango_lineage, 
                                   pattern = 'BA.1', 
                                   perl = TRUE) ~ 'Omicron BA.1*',
                             grepl(x = pango_lineage, 
                                   pattern = 'BA.2.75', 
                                   perl = TRUE) ~ 'Omicron BA.2.75*',
                             grepl(x = pango_lineage, 
                                   pattern = '^(?!.*BA.2.75)*BA.2.*$', 
                                   perl = TRUE) ~ 'Omicron BA.2*',
                             grepl(x = pango_lineage, 
                                   pattern = 'BA.3', 
                                   perl = TRUE) ~ 'Omicron BA.3*',
                             grepl(x = pango_lineage, 
                                   pattern = 'BA.4', 
                                   perl = TRUE) ~ 'Omicron BA.4*',
                             grepl(x = pango_lineage, 
                                   pattern = 'BA.5', 
                                   perl = TRUE) ~ 'Omicron BA.5*', 
                             grepl(x = pango_lineage, 
                                   pattern = 'BQ.1', 
                                   perl = TRUE) ~ 'Omicron BQ.1*',
                             grepl(x = pango_lineage, 
                                   pattern = 'BJ.1', 
                                   perl = TRUE) ~ 'Omicron BJ.1*',
                             grepl(x = pango_lineage, 
                                   pattern = 'XBB.1.5', 
                                   perl = TRUE) ~ 'XBB.1.5*',
                             grepl(x = pango_lineage, 
                                   pattern = 'XBB.1', 
                                   perl = TRUE) ~ 'XBB.1*',
                             grepl(x = pango_lineage, 
                                   pattern = '^(?!.*XBB.1 | XBB.1.5).*$', 
                                   perl = TRUE) ~ 'Recombinant',
                             TRUE ~ 'Other')) |> 
         ## Completing dates to avoid missing dates
         complete(date, nesting(voc_cdc), fill = list(copy_date = 0)) |> 
         ## Putting date at the ending date of the epiweek
         mutate(epiweek = end.of.epiweek(date))

## Saving as compressed file
vroom_write(x = metadata, 
            file = "Data/metadata_us.csv.xz")

## Variant sequences counts
variant_count<-metadata |> 
  ## First probable sequence of Omicron in the US
  filter(date >= "2021-09-01", !is.na(division)) |> 
  ## dummy date column
  mutate(copy_date = as.numeric(epiweek)) |> 
  ## Grouping by epiweek, state and voc_cdc
  group_by(epiweek, division, voc_cdc) |> 
  summarise_at(vars(copy_date), list(n = sum)) |>
  mutate(freq = round(100*n/sum(n),2)) |> 
  rename(name_states = division)

## Saving counts data
vroom_write(x = variant_count, file = "Data/variant_counts_us.csv.xz")

#