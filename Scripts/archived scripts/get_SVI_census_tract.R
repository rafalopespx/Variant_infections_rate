library(tidyverse)
library(tidycensus)
library(sf)
library(here)
source(here("utils", "get_svi.R"))

# primarily a wrapper to get census data using tidycensus into the structure of svi
# https://walker-data.com/tidycensus/articles/basic-usage.html

# set values
state_id <- "CT"
acs_year <- 2020
geo_unit <- "state"

# function to get svi
svi_df <- get_svi(geo_unit, acs_year, state_id)

# covert tibble to sf
svi_sf <- svi_df %>%
  st_as_sf()

# save data
# svi_df %>% write_rds(here("out", "svi_ct_2020_tract_df.rds"))
# 
# svi_sf %>% write_rds(here("out", "svi_ct_2020_tract_sf.rds"))

# vector of all states
us <- unique(fips_codes$state)[1:51]
svi_df_raw <- map_df(us, function(x) {
  get_svi(geo_unit, acs_year, x)
}) 
svi_df_raw<-svi_df_raw|> 
  select(GEOID, name_states, geometry, starts_with("EP_"))



