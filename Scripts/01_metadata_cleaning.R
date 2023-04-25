## Cleaning workspace
rm(list = ls())
gc()

## Loading packages
packs<-c("tidyverse", "vroom")
lapply(packs, require, character.only = TRUE)

## Loading functions
source("Scripts/Functions/functions.R")

## Loading the metadata and states
metadata<-vroom("Data/metadata/metadata_2023-04-09_23-54.tsv.gz")

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
  filter(division %in% states$name_states)
  
vroom_write(x = metadata, file = "Data/metadata_us_raw.csv.xz")

metadata<-vroom("Data/metadata_us_raw.csv.xz")

## Formatting metadata in a more pratical form
metadata <- metadata |> 
  ## Dummy date variable
  mutate(copy_date = as.numeric(date)) |> 
  ## Choosen lineages to be analyzed
  mutate(voc_cdc = case_when(grepl(x = pango_lineage, 
                                   pattern = '^(?=.*B.1.1.7)|(?=.*^Q)' , 
                                   perl = TRUE) ~ 'Alpha*',
                             grepl(x = pango_lineage, 
                                   pattern = '^(?=.*B.1.351)', 
                                   perl = TRUE) ~ 'Beta*',
                             grepl(x = pango_lineage, 
                                   pattern = '^(?=.*^P.1)', 
                                   perl = TRUE) ~ 'Gamma*',
                             grepl(x = pango_lineage, 
                                   pattern = '^(?=.*B.1.617.2)', 
                                   perl = TRUE) ~ 'Delta*',
                             grepl(x = pango_lineage, 
                                   pattern = '^(?=.*BA.1)|(?=.*B.1.1.529.1)|(?=.*BD.1)', 
                                   perl = TRUE) ~ 'Omicron BA.1*',
                             grepl(x = pango_lineage, 
                                   pattern = '^(?=.*BA.2)|(?=.*B.1.1.529.2)|(?=.*^B[JGHLMNSRY])|(?=.*^C[ABMJVH])|(?=.*^D[DVS])|(?=.*^E[PJ])|(?=.*^F[J])', 
                                   perl = TRUE) ~ 'Omicron BA.2*',
                             grepl(x = pango_lineage, 
                                   pattern = '^(?=.*BA.3)|(?=.*B.1.1.529.3)', 
                                   perl = TRUE) ~ 'Omicron BA.3*',
                             grepl(x = pango_lineage, 
                                   pattern = '^(?=.*BA.4)|(?=.*B.1.1.529.4)|(?=.*^C[S])|(?=.*^D[C])', 
                                   perl = TRUE) ~ 'Omicron BA.4*',
                             grepl(x = pango_lineage, 
                                   pattern = '^(?=.*BA.5)|(?=.*B.1.1.529.5)|(?=.*^B[EQKTFUVZW])|(?=.*^C[CDEKNULPYRFGTWZQ])|(?=.*^D[FAGBLHEJQYZURMNKPTW])|(?=.*^E[ACNEQFYDRHTWZSV])|(?=.*^F[A-CF])', 
                                   perl = TRUE) ~ 'Omicron BA.5*',
                             grepl(x = pango_lineage,
                                   pattern = '^(?=.*XBB)',
                                   perl = TRUE) ~ 'Omicron XBB*',
                             grepl(x = pango_lineage, 
                                   pattern = '^(?=.*X)(?!.*XBB)', 
                                   perl = TRUE) ~ 'Recombinant',
                             TRUE ~ 'Other')) |> 
  ## Completing dates to avoid missing dates
  complete(date, nesting(voc_cdc), fill = list(copy_date = 0)) |> 
  ## Putting date at the ending date of the epiweek
  mutate(epiweek = end.of.epiweek(date))

## Talvez tentar inverte a forma de como realizar a classifição, 
# fazer uma primeira estimativa de frequências nas pango lineages, e 
# filtrar somente aquelas linhagens que tem mais de 0.1% de frequencia, 
# com essas linhas realizar a estimativa de infecções

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