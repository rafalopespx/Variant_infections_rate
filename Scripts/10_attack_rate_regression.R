## Cleaning the ambient
rm(list = ls())
gc()

## Loading Libraries
packs = c("tidyverse", "vroom", "MASS", "ggeffects", "marginaleffects", "broom.helpers", "sf", "tidycensus", "tidyverse")
lapply(packs,require, character.only = TRUE)

# Loading functions
source("Scripts/Functions/functions.R")
source("Scripts/get_svi.R")

## Loading data sources
states_fulldata<-vroom("Data/state_full_data.csv.xz")

states_attack_rates<-vroom("Data/state_attack_rate_variants.csv.xz")

## SVI variable vector of all states
us <- unique(fips_codes$state)[1:51]
acs_year <- 2020
geo_unit <- "state"
svi_df_raw <- map_df(us, function(x) {
  get_svi(geo_unit, acs_year, x)
}) 
svi_df_raw<-svi_df_raw|> 
  select(GEOID, name_states, geometry, starts_with("EP_"))

## Joining data streams
states_full_model<-left_join(states_attack_rates, svi_df_raw)

var_names<-colnames(states_full_model[,c(startsWith(colnames(states_full_model), "EP_"))])

formula <- reformulate(response = "attack_rate", termlabels = var_names)

variants<-unique(states_full_model$variant)

effect_size_fun<-function(x){
  plt<-x %>%
    # perform initial tidying of the model
    tidy_and_attach(exponentiate = TRUE, conf.int = TRUE) %>%
    # adding in the reference row for categorical variables
    tidy_add_reference_rows() %>%
    # adding a reference value to appear in plot
    tidy_add_estimate_to_reference_rows() %>%
    # adding the variable labels
    tidy_add_term_labels() %>%
    # removing intercept estimate from model
    tidy_remove_intercept()%>%
    mutate(
      plot_label = paste(label) %>%
        forcats::fct_inorder() %>%
        forcats::fct_rev()
    )
}

poisson_list<-list()
effect_size_list<-list()

for (i in variants) {
  data<-states_full_model |> 
    filter(variant == i)
  
  poisson_list[[i]]<-glm(data = data, 
                         formula = formula,
                         family = quasipoisson(link = "log"))
  
  effect_size_list[[i]]<-poisson_list[[i]] |> 
    effect_size_fun() |> 
    mutate(variant = i)
}

effect_size_list <- effect_size_list |> 
  bind_rows() |> 
  mutate(plt_label = case_when(label == 'EP_POV150' ~ "Below 150% poverty",
                               label == 'EP_UNEMP' ~ "Unemployed",
                               label == 'EP_HBURD' ~ "Housing Cost Burden",
                               label == 'EP_NOHSDP' ~ "No High School Diploma",
                               label == 'EP_UNINSUR' ~ "No Health Insurance",
                               label == 'EP_AGE65' ~ "Aged 65 and older",
                               label == 'EP_AGE17' ~ "Aged 17 and younger",
                               label == 'EP_DISABL' ~ "Civilian with a Disability",
                               label == 'EP_SNGPNT' ~ "Single-Parent Households",
                               label == 'EP_LIMENG' ~ "English Language Proficieny",
                               label == 'EP_MINRTY' ~ "Racial/Ethnic Minority",
                               label == 'EP_MUNIT' ~ "Multi-unit Structures",
                               label == 'EP_MOBILE' ~ "Mobile Homes",
                               label == 'EP_CROWD' ~ "Crowding",
                               label == 'EP_NOVEH' ~ "No Vehicle",
                               label == 'EP_GROUPQ' ~ "Group Quartes"))


effect_size_list|> 
  ggplot(aes(x = plt_label, y = estimate, 
             ymin = conf.low, ymax = conf.high, 
             color = variable)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_pointrange() +
  coord_flip() +
  theme(legend.position = "none") +
  labs(
    y = "IRR",
    x = "Social Vulnerability Index (SVI) components",
    title = title
  )+
  theme_minimal()+
  theme(legend.position = "bottom")

## Exploratory models
poi_model<-glm(data = states_full_model, 
               formula = formula,
               family = poisson(link = "log"))
summary(poi_model)

effect_size_poi<-poi_model|>
  effect_size_plt(title = "Poisson Model")
effect_size_poi

# quasipoi_model<-glm(data = states_model, 
#                     formula = formula, 
#                     family = quasipoisson(link = "log"))
# summary(quasipoi_model)
# 
# effect_size_quasipoi<-quasipoi_model|>
#   effect_size_plt(title = "Quasi-Poisson Model")
# effect_size_quasipoi
# 
# nb_model<-glm.nb(data = states_model, 
#                  formula = formula, 
#                  link = "log")
# summary(nb_model)
# 
# effect_size_nb <- nb_model |> 
#   effect_size_plt(title = "Negative Binomial Model")
# effect_size_nb

# 
# 
# 
# 
plot_predictions(poi_model,
                 condition = c("household_income", "variant"),
                 type = "response")+
  theme_minimal()
# 
plot_predictions(poi_model,
                 condition = c("household_size", "variant"),
                 type = "response")+
  theme_minimal()
# 
# plot_predictions(poi_model,
#                  condition = c("household_income", "variant", "household_size"),
#                  type = "response")+
#   theme_minimal()
# 
# newdata<-expand.grid(household_income = seq(min(states_model$household_income, na.rm = T), 
#                                             max(states_model$household_income, na.rm = T),
#                                             by = 1000), 
#                      household_size = seq(min(states_model$household_size, na.rm = T),
#                                           max(states_model$household_size, na.rm = T), 
#                                           by = 0.1),
#                      # name_states = unique(states_model$name_states),
#                      variant = unique(states_model$variant))
# 
# pred_poi_model<-data.frame(newdata,
#                            predict.glm(object = poi_model, 
#                                        newdata = newdata, 
#                                        type = "response", 
#                                        se.fit = T)) |> 
#   # mutate(attack_rate = round(attack_rate, 2)) |> 
#   mutate(model = "Poisson")
# 
# pred_quasipoi_model<-data.frame(newdata,
#                                 predict.glm(object = poi_model, 
#                                             newdata = newdata, 
#                                             type = "link", 
#                                             se.fit = T)) |> 
#   # mutate(attack_rate = round(attack_rate, 2)) |> 
#   mutate(model = "Quasi-Poisson")
# 
# pred_nb_model<-data.frame(newdata,
#                           predict.glm(object = poi_model, 
#                                       newdata = newdata, 
#                                       type = "link", 
#                                       se.fit = T)) |> 
#   # mutate(attack_rate = round(attack_rate, 2)) |> 
#   mutate(model = "Negative Binomial")
# 
# pred_nb_model<-data.frame(newdata,
#                           predict.glm(object = poi_model, 
#                                       newdata = newdata, 
#                                       type = "link", 
#                                       se.fit = T)) |> 
#   # mutate(attack_rate = round(attack_rate, 2)) |> 
#   mutate(model = "Negative Binomial")
# 
# 
# ## Exploratory plots
# # pred_all_models<-rbind(pred_poi_model, pred_quasipoi_model, pred_nb_model)
# 
# pred_poi_model |> 
#   ggplot(aes(x = household_size, y = exp(fit), 
#              # ymin = fit - 1.96*se.fit,
#              # ymax = fit + 1.96*se.fit,
#              # col = household_income
#              ))+
#   geom_line()+
#   # geom_ribbon()+
#   facet_wrap(variant~., scales = "free_y", nrow = 1)+
#   theme_minimal()+
#   labs(title = "Poisson model", 
#        y = "Attack Rate \n (%) of pop. ever infected", 
#        x = "Household Income \n ($)")+
#   MetBrewer::scale_color_met_c(name = "Average Household size",
#                                palette_name = "Demuth", 
#                                direction = -1,
#                                guide = guide_colorsteps(title.position = "top",ticks = T, even.steps = F))+
#   theme(legend.position = "bottom", 
#         legend.key.width = grid::unit(2, "cm"))
# 
# effect_size_poi<-ggpredict(poi_model, terms = c("household_size", "household_income", "variant"))
# effect_size_poi
# 
# pred_quasipoi_model |> 
#   ggplot(aes(x = household_size, y = attack_rate, 
#              col = household_income))+
#   geom_point()+
#   facet_wrap(variant~., scales = "free_y", nrow = 1)+
#   theme_minimal()+
#   labs(title = "Quasi-Poisson model", 
#        y = "Attack Rate \n (%) of pop. ever infected", 
#        x = "Household Income \n ($)")+
#   MetBrewer::scale_color_met_c(name = "Average Household size",
#                                palette_name = "Demuth", 
#                                direction = -1, 
#                                guide = guide_colorbar(keywidth = grid::unit(3, "cm"),
#                                                       title.position = "top", 
#                                                       show.limits = T))+
#   theme(legend.position = "bottom")
# 
# pred_nb_model|> 
#   ggplot(aes(x = household_size, y = attack_rate, 
#              col = household_income))+
#   geom_point()+
#   facet_wrap(variant~., scales = "free_y", nrow = 1)+
#   theme_minimal()+
#   labs(title = "Negative Binomial model", 
#        y = "Attack Rate \n (%) of pop. ever infected", 
#        x = "Household Income \n ($)")+
#   MetBrewer::scale_color_met_c(name = "Average Household size",
#                                palette_name = "Demuth", 
#                                direction = -1, 
#                                guide = guide_colorbar(keywidth = grid::unit(3, "cm"),
#                                                       title.position = "top", 
#                                                       show.limits = T))+
#   theme(legend.position = "bottom")
# 
# 
# # ## A try into GLMM
# # 
# # sample_data<-states_full_model |> 
# #   slice_sample(n = 5000)
# # 
# # poi_mm_model<-lme4::glmer(data = sample_data, 
# #                           formula = "incidence ~ 
# #                           household_size + household_income + variant + (days|name_states)", 
# #                           family = poisson(link = "log"))
# # summary(poi_mm_model)
# 
# 
