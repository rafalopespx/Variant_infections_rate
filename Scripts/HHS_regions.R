# HHS 10 regions

region_1 = c("Connecticut", "Maine", "Massachusetts", "New Hampshire", "Rhode Island", "Vermont")
region_2 = c("New Jersey", "New York")
region_3 = c("Delaware", "District of Columbia", "Maryland", "Pennsylvania", "Virginia", "West Virginia")
region_4 = c("Alabama", "Florida", "Georgia", "Kentucky", "Mississippi", "North Carolina", "South Carolina",  "Tennessee")
region_5 = c("Illinois", "Indiana", "Michigan", "Minnesota", "Ohio", "Wisconsin")
region_6 = c("Arkansas", "Louisiana", "New Mexico", "Oklahoma", "Texas")
region_7 = c("Iowa", "Kansas", "Missouri", "Nebraska")
region_8 = c("Colorado", "Montana", "North Dakota", "South Dakota", "Utah", "Wyoming")
region_9 = c("Arizona", "California", "Hawaii", "Nevada")
region_10 = c("Alaska", "Idaho", "Oregon", "Washington")

HHS_region = list()
for(i in 1:nrow(metadata_usa)) {
  if (metadata_usa$division[i] %in% region_1) {
    HHS_region[i] = "Region 1"
  } else if (metadata_usa$division[i] %in% region_2){
    HHS_region[i] = "Region 2"
  } else if (metadata_usa$division[i] %in% region_3){
    HHS_region[i] = "Region 3"
  }
  else if (metadata_usa$division[i] %in% region_4){
    HHS_region[i] = "Region 4"
  }
  else if (metadata_usa$division[i] %in% region_5){
    HHS_region[i] = "Region 5"
  }
  else if (metadata_usa$division[i] %in% region_6){
    HHS_region[i] = "Region 6"
  }
  else if (metadata_usa$division[i] %in% region_7){
    HHS_region[i] = "Region 7"
  }
  else if (metadata_usa$division[i] %in% region_8){
    HHS_region[i] = "Region 8"
  }
  else if (metadata_usa$division[i] %in% region_9){
    HHS_region[i] = "Region 9"
  }
  else if (metadata_usa$division[i] %in% region_10){
    HHS_region[i] = "Region 10"
  }
}
HHS_region = tibble(HHS_region)