## Cleaning the ambient
rm(list = ls())
gc()

## Loading Libraries
packs = c("tidyverse", "vroom", "ggeffects", "marginaleffects", "broom.helpers", "metafor")
lapply(packs,require, character.only = TRUE)


states_sum <- vroom("Data/states_immune_attack_rate.csv.xz")

## meta-regression

metareg_data<-states_sum |> 
  filter(!is.na(type_protected)) |>
  pivot_wider(names_from = type_protected,
              values_from = values) |> 
  dplyr::select(name_states, variant, attack_rate, 
                starts_with("protectedInf"))
metareg_data$vi<-0

logit<-function(x){
  y<-log(x/(1-x))
  return(y)
}

variants <- unique(metareg_data$variant)
model_list<-list()

## Looping over the variants
for (i in variants) {
  dat<-metareg_data |>
    filter(variant == i) |> 
    mutate(attack_rate_logit = logit(attack_rate/100))
  
  var_names <- colnames(dat |> 
                          dplyr::select(-attack_rate, -attack_rate_logit, -variant, -name_states))
  
  model_list[[i]] <- rma(data = dat,
                         yi = attack_rate, 
                         vi = vi,
                         mods = reformulate(var_names, intercept = F), 
                         method = "REML")
}

ggcoef_compare(model_list, 
               # exponentiate = T,
               variable_labels = list("protectedInf_infection" ~ "By Infection",
                                      "protectedInf_vaccination" ~ "By Vaccination",
                                      "protectedInf_hybrid" ~ "Hybrid",
                                      "protectedInf" ~ "Total protection"))+
  scale_color_manual(values = MetBrewer::met.brewer(palette_name = "Archambault",
                                                    n = 5,
                                                    type = "discrete"))+
  labs(x = "(OR) \n Odds Ratio", 
       y = "Effectively Protected Estimates")

