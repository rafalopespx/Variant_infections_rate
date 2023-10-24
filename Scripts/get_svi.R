#vars <- tidycensus::load_variables(2019, "acs5/subject")

# tidycensus call
get_svi <- function(geo_unit, acs_year, state_id) {
  
  # ACS variables used to calculate SVI
  acs_svi_vars <- c(
    
    ## socioeconomic variables
    "S0601_C01_001",  # population estimate
    "DP04_0001", # housing units estimate
    "DP02_0001", # households estimate
    "S1701_C01_040", # persons below 150% poverty (replaced B17001_002)
    "DP03_0005",  # civilian (age 16+) unemployed
    "S2503_C01_028", # occupied housing units with annual income of less than $20,000 with monthly housing costs of 30 percent or more of annual income (replaced B19301_001)
    "S2503_C01_032", # occupied housing units with annual income of $20,000 to $34,999 with monthly housing costs of 30 percent or more of annual income (replaced B19301_001)
    "S2503_C01_036", # occupied housing units with annual income of $35,000 to $49,999 with monthly housing costs of 30 percent or more of annual income (replaced B19301_001)
    "S2503_C01_040", # occupied housing units with annual income of $50,000 to $74,999 with monthly housing costs of 30 percent or more of annual income (replaced B19301_001)
    "B06009_002", # persons (age 25+) with no high school diploma
    "S2701_C04_001",  # uninsured in the total civilian non institutionalized population
    
    ##  household characteristics
    "S0101_C01_030", # persons aged 65 and older
    "B09001_001", # persons aged 17 and younger,
    "DP02_0072", # civilian non institutionalized population with a disability (replaced DP02_0071)
    "B11012_010", # male householder, no spouse or partner present, with own children under 18 years + (replaced DP02_0007)
    "B11012_015", # female householder, no spouse or partner present, with own children under 18 years + (replaced DP02_0007)
    # persons (age 5+) who speak English less than well
    "B16005_007",
    "B16005_008",
    "B16005_012",
    "B16005_013",
    "B16005_017",
    "B16005_018",
    "B16005_022",
    "B16005_023",
    "B16005_029",
    "B16005_030",
    "B16005_034",
    "B16005_035",
    "B16005_039",
    "B16005_040",
    "B16005_044",
    "B16005_045",
    
    ## racial & ethnic minority status
    # minority
    "DP05_0071", # Hispanic or Latino, Total Population
    "DP05_0078", # Black and African American Not Hispanic or Latino
    "DP05_0079", # American Indian and Alaska Native Not Hispanic or Latino
    "DP05_0080", # Asian Not Hispanic or Latino
    "DP05_0081", # Native Hawaiian and Other Pacific Islander Not Hispanic or Latino
    "DP05_0082", # Other Races Not Hispanic or Latino
    "DP05_0083", # Two or More Races Not Hispanic or Latino
    
    
    ## housing type/transportation
    "DP04_0012", # units in structure - total housing units - 10 to 19 units
    "DP04_0013", # units in structure - total housing units - 20 or more units
    "DP04_0014", # mobile homes
    "DP04_0078", # occupants per room, occupied housing units, 1.01 to 1.50 +
    "DP04_0079", # occupants per room, occupied housing units, 1.51 or more
    "DP04_0058", # households with no vehicle available
    "B26001_001", # persons in group quarters
    
    ## socioeconomic variables (percentages)
    "S1701_C01_001", # population for whom poverty status is determined
    "DP03_0009P", # unemployment rate
    "S2503_C01_001", # occupied housing units
    "S0601_C01_033", # percentage of persons with no high school diploma (age 25+)
    "S2701_C05_001", # percentage uninsured in the total civilian noninstitutionalized population
    
    ##  household characteristics (percentages)
    "S0101_C02_030", # percentage of persons aged 65 and older
    "DP02_0072P", # percentage of a civilian noninstitutionalized population with a disability
    "B16005_001", # population age 5 and older
    
    ## housing type/transportation (percentages)
    "DP04_0014P", # percentage of mobile homes
    "DP04_0002", # occupied housing units
    "DP04_0058P" # percentage of households with no vehicle available
  )
  
  svi_raw_vars <- get_acs(geography = geo_unit,
                          variables =  acs_svi_vars,
                          year = acs_year,
                          state =  state_id,
                          geometry = TRUE, 
                          keep_geo_vars = TRUE) %>% 
    tibble() %>% 
    dplyr::select(GEOID, 
           name_states = NAME.x,
           state = LSAD,
           ALAND, AWATER, geometry,
           variable, estimate, moe) %>% 
    rename(E = estimate,
           M = moe) %>%
    pivot_wider(names_from = variable,
                names_glue = "{variable}{.value}",
                values_from = c(E, M)) %>% 
    filter(!is.na(ALAND))
  
  svi_clean <- svi_raw_vars %>%
    mutate(
      ## socioeconomic variables
      E_TOTPOP = S0601_C01_001E,
      M_TOTPOP = S0601_C01_001M,
      E_HU = DP04_0001E,
      M_HU = DP04_0001M,
      E_HH = DP02_0001E,
      M_HH = DP02_0001M,
      E_POV150 = S1701_C01_040E,
      M_POV150 = S1701_C01_040M,
      E_UNEMP = DP03_0005E,
      M_UNEMP = DP03_0005M,
      E_HBURD = S2503_C01_028E + S2503_C01_032E + S2503_C01_036E + S2503_C01_040E,
      M_HBURD = sqrt(S2503_C01_028M ^ 2 + S2503_C01_032M ^ 2 + S2503_C01_036M ^ 2 + S2503_C01_040M ^ 2),
      E_NOHSDP = B06009_002E,
      M_NOHSDP = B06009_002M,
      E_UNINSUR = S2701_C04_001E,
      M_UNINSUR = S2701_C04_001M,
      
      ##  household characteristics
      E_AGE65 = S0101_C01_030E,
      M_AGE65 = S0101_C01_030M,
      E_AGE17 = B09001_001E,
      M_AGE17 = B09001_001M,
      E_DISABL = DP02_0072E,
      M_DISABL = DP02_0072M,
      E_SNGPNT = B11012_010E +
        B11012_015E,
      M_SNGPNT = sqrt(B11012_010M ^ 2 +
                        B11012_015M ^ 2),
      E_LIMENG = B16005_007E +
        B16005_008E +
        B16005_012E +
        B16005_013E +
        B16005_017E +
        B16005_018E +
        B16005_022E +
        B16005_023E +
        B16005_029E +
        B16005_030E +
        B16005_034E +
        B16005_035E +
        B16005_039E +
        B16005_040E +
        B16005_044E +
        B16005_045E,
      M_LIMENG = sqrt(B16005_007M ^ 2 +
                        B16005_008M ^ 2 +
                        B16005_012M ^ 2 +
                        B16005_013M ^ 2 +
                        B16005_017M ^ 2 +
                        B16005_018M ^ 2 +
                        B16005_022M ^ 2 +
                        B16005_023M ^ 2 +
                        B16005_029M ^ 2 +
                        B16005_030M ^ 2 +
                        B16005_034M ^ 2 +
                        B16005_035M ^ 2 +
                        B16005_039M ^ 2 +
                        B16005_040M ^ 2 +
                        B16005_044M ^ 2 +
                        B16005_045M ^ 2),
      
      ## racial & ethnic minority status
      E_MINRTY = DP05_0071E +
        DP05_0078E +
        DP05_0079E +
        DP05_0080E +
        DP05_0081E +
        DP05_0082E +
        DP05_0083E,
      M_MINRTY = sqrt(DP05_0071M ^ 2 +
                        DP05_0078M ^ 2 +
                        DP05_0079M ^ 2 +
                        DP05_0080M ^ 2 +
                        DP05_0081M ^ 2 +
                        DP05_0082M ^ 2 +
                        DP05_0083M ^ 2),
      
      ## housing type/transportation
      E_MUNIT = DP04_0012E + DP04_0013E,
      M_MUNIT = sqrt(DP04_0012M ^ 2 + DP04_0013M ^ 2),
      E_MOBILE = DP04_0014E,
      M_MOBILE = DP04_0014M,
      E_CROWD = DP04_0078E + DP04_0079E,
      M_CROWD = sqrt(DP04_0078M ^ 2 + DP04_0079M ^ 2),
      E_NOVEH = DP04_0058E,
      M_NOVEH = DP04_0058M,
      E_GROUPQ = B26001_001E,
      M_GROUPQ = B26001_001M) %>%
    # create percentages
    mutate(
      ## socioeconomic status
      EP_POV150 = (E_POV150/S1701_C01_001E) * 100,
      MP_POV150 = ((sqrt(M_POV150 ^ 2 - ((EP_POV150 / 100) ^ 2 * S1701_C01_001M ^ 2))) / S1701_C01_001E) * 100,
      EP_UNEMP = DP03_0009PE,
      MP_UNEMP = DP03_0009PM,
      EP_HBURD = (E_HBURD / S2503_C01_001E) * 100,
      MP_HBURD = ((sqrt(M_HBURD ^ 2 - ((EP_HBURD / 100) ^ 2 * S2503_C01_001M ^ 2))) / S2503_C01_001E) * 100,
      EP_NOHSDP = S0601_C01_033E,
      MP_NOHSDP = S0601_C01_033M,
      EP_UNINSUR = S2701_C05_001E,
      MP_UNINSUR = S2701_C05_001M,
      
      ## housing type/transportation percentage
      EP_AGE65 = S0101_C02_030E,
      MP_AGE65 = S0101_C02_030M,
      EP_AGE17 = (E_AGE17 / E_TOTPOP) * 100,
      MP_AGE17 = ((sqrt(M_AGE17 ^ 2 - ((EP_AGE17 / 100) ^ 2 * M_TOTPOP ^ 2))) / E_TOTPOP) * 100,
      EP_DISABL = DP02_0072PE,
      MP_DISABL = DP02_0072PM,
      EP_SNGPNT = (E_SNGPNT / E_HH) * 100,
      MP_SNGPNT = ((sqrt(M_SNGPNT ^ 2 - ((EP_SNGPNT / 100) ^ 2 * M_HH ^ 2))) / E_HH) * 100,
      EP_LIMENG = (E_LIMENG / B16005_001E) * 100,
      MP_LIMENG = ((sqrt(M_LIMENG ^ 2 - ((EP_LIMENG / 100) ^ 2 * B16005_001M ^ 2))) / B16005_001E) * 100,
      
      ## racial & ethnic minority status percentages
      EP_MINRTY = (E_MINRTY / E_TOTPOP) * 100,
      MP_MINRTY = ((sqrt(M_MINRTY ^ 2 - ((EP_MINRTY / 100) ^ 2 * M_TOTPOP ^ 2))) / E_TOTPOP) * 100,
      
      # ## breakdown of racial & ethnic minority status percentages
      # EP_LATINO = (DP05_0071E / E_TOTPOP) * 100,
      # MP_LATINO = ((sqrt(DP05_0071M ^ 2 - ((DP05_0071M / 100) ^ 2 * M_TOTPOP ^ 2))) / E_TOTPOP) * 100,
      # EP_BLACK = (DP05_0078E / E_TOTPOP) * 100,
      # MP_BLACK = ((sqrt(DP05_0078M ^ 2 - ((DP05_0078M / 100) ^ 2 * M_TOTPOP ^ 2))) / E_TOTPOP) * 100,
      # EP_ASIAN = (DP05_0080E / E_TOTPOP) * 100,
      # MP_ASIAN = ((sqrt(DP05_0080M ^ 2 - ((DP05_0080M / 100) ^ 2 * M_TOTPOP ^ 2))) / E_TOTPOP) * 100,
      # EP_NATIVE = ((DP05_0081E + DP05_0079E )/ E_TOTPOP) * 100,
      # MP_NATIVE = ((sqrt((DP05_0081E + DP05_0079E) ^ 2 - (((DP05_0081E + DP05_0079E) / 100) ^ 2 * M_TOTPOP ^ 2))) / E_TOTPOP) * 100,
      # EP_TWOPLUS = (DP05_0083E / E_TOTPOP) * 100,
      # MP_TWOPLUS = ((sqrt(DP05_0083M ^ 2 - ((DP05_0083M / 100) ^ 2 * M_TOTPOP ^ 2))) / E_TOTPOP) * 100,
      # EP_OTHERNL = (DP05_0082E / E_TOTPOP) * 100,
      # MP_OTHERNL = ((sqrt(DP05_0082M ^ 2 - ((DP05_0082M / 100) ^ 2 * M_TOTPOP ^ 2))) / E_TOTPOP) * 100,
      
      ## housing type/transportation percentages
      EP_MUNIT = (E_MUNIT / E_HU) * 100,
      MP_MUNIT = ((sqrt(M_MUNIT ^ 2 - ((EP_MUNIT / 100) ^ 2 * M_HU ^ 2))) / E_HU) * 100,
      EP_MOBILE = DP04_0014PE,
      MP_MOBILE = DP04_0014PM,
      EP_CROWD = (E_CROWD / DP04_0002E) * 100,
      MP_CROWD = ((sqrt(M_CROWD ^ 2 - ((EP_CROWD / 100) ^ 2 * DP04_0002M ^ 2))) / DP04_0002E) * 100,
      EP_NOVEH = DP04_0058PE,
      MP_NOVEH = DP04_0058PM,
      EP_GROUPQ = (E_GROUPQ / E_TOTPOP) * 100,
      MP_GROUPQ = ((sqrt(M_GROUPQ ^ 2 - ((EP_GROUPQ / 100) ^ 2 * M_TOTPOP ^ 2))) / E_TOTPOP) * 100) %>%
    # drop unneeded variables
    dplyr::select(!(B16005_001E:DP04_0058PM))
  
  svi_rank <- svi_clean  
  # |> 
  #   # group_by(ST) %>%
  #   mutate(
  #     EPL_POV150 = percent_rank(EP_POV150),
  #     EPL_UNEMP = percent_rank(EP_UNEMP),
  #     EPL_HBURD = percent_rank(EP_HBURD),
  #     EPL_NOHSDP = percent_rank(EP_NOHSDP),
  #     EPL_UNINSUR = percent_rank(EP_UNINSUR),
  #     SPL_THEME1 = EPL_POV150 + EPL_UNEMP +  EPL_HBURD + EPL_NOHSDP + EPL_UNINSUR,
  #     RPL_THEME1 = percent_rank(SPL_THEME1),
  #     
  #     EPL_AGE65 = percent_rank(EP_AGE65),
  #     EPL_AGE17 = percent_rank(EP_AGE17),
  #     EPL_DISABL = percent_rank(EP_DISABL),
  #     EPL_SNGPNT = percent_rank(EP_SNGPNT),
  #     EPL_LIMENG = percent_rank(EP_LIMENG),
  #     SPL_THEME2 = EPL_AGE65 + EPL_AGE17 + EPL_DISABL + EPL_SNGPNT + EPL_LIMENG,
  #     RPL_THEME2 = percent_rank(SPL_THEME2),
  #     
  #     EPL_MINRTY = percent_rank(EP_MINRTY),
  #     SPL_THEME3 = EPL_MINRTY,
  #     RPL_THEME3 = percent_rank(SPL_THEME3),
  #     
  #     EPL_MUNIT = percent_rank(EP_MUNIT),
  #     EPL_MOBILE = percent_rank(EP_MOBILE),
  #     EPL_CROWD = percent_rank(EP_CROWD),
  #     EPL_NOVEH = percent_rank(EP_NOVEH),
  #     EPL_GROUPQ = percent_rank(EP_GROUPQ),
  #     SPL_THEME4 = EPL_MUNIT + EPL_MOBILE + EPL_CROWD + EPL_NOVEH + EPL_GROUPQ,
  #     RPL_THEME4 = percent_rank(SPL_THEME4),
  #     
  #     SPL_THEMES = SPL_THEME1 + SPL_THEME2 + SPL_THEME3 + SPL_THEME4,
  #     RPL_THEMES = percent_rank(SPL_THEMES))
  
  return(svi_rank)
}