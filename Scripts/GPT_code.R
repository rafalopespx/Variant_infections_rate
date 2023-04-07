library(tidyverse)
library(zoo)

# Load the dataset
data <- vroom("Data/variant_counts_us.csv.xz")

# Filter the time series of infections when it has over than 10 cases
filtered_data <- data %>%
  group_by(state, variant) %>%
  filter(sum(infections) >= 10) %>%
  ungroup()

# Function to estimate Rt for a given state-variant combination
estimate_Rt <- function(state, variant) {
  # Subset the data for the current state-variant combination
  subset_data <- filtered_data %>% filter(state == state, variant == variant)
  
  # Convert the data to a time series
  ts_data <- zoo(subset_data$infections, subset_data$dates)
  
  # Estimate the Rt using the EpiEstim package
  epi_result <- EpiEstim::estimate_R(zoo::na.omit(ts_data), method = "parametric_si",
                                     config = EpiEstim::make_config(list(mean_si = 4.7, std_si = 2.9)))
  
  # Return the estimated Rt as a data frame
  data.frame(state = state, variant = variant, dates = index(epi_result$R),
             Rt = coredata(epi_result$R))
}

# Estimate Rt for all state-variant combinations
Rt_results <- expand.grid(state = unique(filtered_data$state), variant = unique(filtered_data$variant)) %>%
  rowwise() %>%
  do.call(rbind, map2(.$state, .$variant, estimate_Rt))

# Plot the Rt for each state-variant combination
ggplot(Rt_results, aes(x = dates, y = Rt)) + 
  geom_line() + 
  facet_grid(rows = vars(variant), cols = vars(state)) + 
  labs(title = "Effective Reproduction Number (Rt) for Infections per Variant",
       x = "Date", y = "Rt")
