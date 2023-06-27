## Cleaning the ambient
rm(list = ls())
gc()

## Loading Libraries
packs = c("nnet", "tidyverse", "vroom", "marginaleffects", "doParallel")
lapply(packs,require, character.only = TRUE)

## Parallel
### Create cluster
cl<-makeCluster(detectCores() - 1)

### Register cluster
registerDoParallel(cl)

## Loading functions
source("Scripts/functions.R")
source("Scripts/fastmultinomHESS.R")

## Reading the database
variants_count<-vroom("Data/variant_counts_us.csv.xz")

## Multinomial Logistic Regression
# Logistic regression model for the growth rate of variants, 
# Y ~ epiweek + name_state, where Y is the probability of being a certian variant,
# or the frequency of variant over the totality of variants sequenced, 
# predicted by week of collection of the genome and per state of the US
# length(unique(variants_counts$epiweek))
# 79
# length(unique(variants_counts$name_states))
# 59
# 79 weeks*59 regions = 4661 predictors

multinomial_states<-nnet::multinom(data = variants_count, 
                                   formula = voc_cdc ~ epiweek + name_states, 
                                   weights = log(n),
                                   MaxNWts = 100000,
                                   Hess = FALSE, 
                                   maxit = 10000, 
                                   allowParallel = TRUE)

## Adding the Hessian
multinomial_states$Hessian <- fastmultinomHess(multinomial_states, 
                                               model.matrix(multinomial_states))

## Adding the VCOV
multinomial_states$vcov <- vcov(multinomial_states)

## Summary of model
summary(multinomial_states)

## Vectors of unique classes of each predictor, to construct predicts
name_states<-unique(variants_count$name_states)

epiweek<-variants_count |> 
  filter(epiweek >= "2021-09-01") |> 
  pull(var = epiweek) |> 
  unique()

voc_cdc<-unique(variants_count$voc_cdc)

## Newdata
newdata_states<-data.frame(name_states = rep(name_states, 
                                             length(epiweek)), 
                           epiweek = rep(epiweek, 
                                         length(name_states)))

## rebinding all together
pred_states<-cbind(newdata_states, 
                   predict(multinomial_states, 
                           newdata = newdata_states, 
                           type = "probs"))

## Predicts by States, 
## non-numbered are nowcast, until last date of submitted on the dataset
pred_states_long<-pred_states |>
  pivot_longer(cols = c(-epiweek, -name_states), 
               names_to = "voc_cdc", 
               values_to = "probability") |> 
  rename(median = probability) |> 
  mutate(upper = median + 1.96*sd(median), 
         lower = median - 1.96*sd(median))

pred_us_long<-pred_states |>
  pivot_longer(cols = c(-epiweek, -name_states), 
               names_to = "voc_cdc", 
               values_to = "probability") |> 
  group_by(epiweek, voc_cdc) |> 
  reframe(prob_quantile = quantile(probability, 
                                   probs = c(0.275, 0.5, 0.975), na.rm = T), 
          probs = c("lower", "median", "upper")) |> 
  pivot_wider(id_cols = c("epiweek", "voc_cdc"), 
              names_from = "probs", 
              values_from = "prob_quantile")

## Saving results
vroom_write(x = pred_us_long, 
            file = "Output/pred_multinomial_us.csv.xz")

vroom_write(x = pred_states_long, 
            file = "Output/pred_multinomial_us_states.csv.xz")


#