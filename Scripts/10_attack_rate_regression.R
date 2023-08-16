## Cleaning the ambient
rm(list = ls())
gc()

## Loading Libraries
packs = c("tidyverse", "vroom", "MASS", "ggeffects", "marginaleffects", "broom.helpers", "sf", "tidycensus", "tidyverse", "ggstats")
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
svi_df<-svi_df_raw|> 
  select(GEOID, name_states, geometry, starts_with("EP_"))

## Joining data streams
states_full_model<-left_join(states_attack_rates, svi_df)

var_names<-states_full_model |> 
  select(starts_with("EP_"), -EP_MINRTY, -EP_OTHERNL) |> 
  colnames()

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

poisson_list<-gaussian_list<-list()
effect_size_list<-list()

for (i in variants) {
  data<-states_full_model |> 
    filter(variant == i) |> 
    mutate(attack_rate_logit = Logit(attack_rate/100))
  
  poisson_list[[i]]<-glm(data = data, 
                         formula = reformulate(response = "attack_rate", termlabels = var_names), 
                         family = poisson(link = "log"))
  
  gaussian_list[[i]]<-glm(data = data, 
                          formula = reformulate(response = "attack_rate_logit", termlabels = var_names), 
                          family = gaussian(link = "identity"))
}

var_plt<-list('EP_POV150' ~ "Below 150% poverty",
              'EP_UNEMP' ~ "Unemployed",
              'EP_HBURD' ~ "Housing Cost Burden",
              'EP_NOHSDP' ~ "No High School Diploma",
              'EP_UNINSUR' ~ "No Health Insurance",
              'EP_AGE65' ~ "Aged 65 and older",
              'EP_AGE17' ~ "Aged 17 and younger",
              'EP_DISABL' ~ "Civilian with a Disability",
              'EP_SNGPNT' ~ "Single-Parent Households",
              'EP_LIMENG' ~ "English Language Proficieny",
              # 'EP_MINRTY' ~ "Racial/Ethnic Minority",
              'EP_LATINO' ~ "Hispanic or Latino (of any race)",
              'EP_BLACK' ~ "Black or African American, not Hispanic or Latino",
              'EP_ASIAN' ~ "Asian, not Hispanic or Latino",
              'EP_NATIVE' ~ "American, Alaskan, Hawaiian Native, and Pacific Islander, not Hispanic or Latino",
              'EP_TWOPLUS' ~ "Two or More Races, not Hispanic or Latino",
              # 'EP_OTHERNL' ~ "Other Races, not Hispanic or Latino",
              'EP_MUNIT' ~ "Multi-unit Structures",
              'EP_MOBILE' ~ "Mobile Homes",
              'EP_CROWD' ~ "Crowding",
              'EP_NOVEH' ~ "No Vehicle",
              'EP_GROUPQ' ~ "Group Quartes")

poisson_list |> 
  ggstats::ggcoef_compare(exponentiate = T, 
                          variable_labels = var_plt,
                          type = "faceted")+
  scale_color_manual(values = MetBrewer::met.brewer(palette_name = "Nizami",
                                                    n = 22,
                                                    type = "continuous"))+
  labs(x = "(IRR) \n Incidence Rate Ratio", y = "(SVI) \n Social Vulnerabilit Index components")+
  # xlim(c(NA,1.70))+
  guides(color = "none")

gaussian_list |> 
  ggstats::ggcoef_compare(exponentiate = T, 
                          variable_labels = var_plt,
                          type = "faceted")+
  scale_color_manual(values = MetBrewer::met.brewer(palette_name = "Nizami",
                                                    n = 22,
                                                    type = "continuous"))+
  labs(x = "Odds Ratio", y = "(SVI) \n Social Vulnerabilit Index components")+
  # xlim(c(NA,1.70))+
  guides(color = "none")

# ## Exploratory models
# poi_model<-glm(data = states_full_model, 
#                formula = formula,
#                family = poisson(link = "log"))
# summary(poi_model)
# 
# effect_size_poi<-poi_model|>
#   effect_size_plt(title = "Poisson Model")
# effect_size_poi

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
# plot_predictions(poi_model,
#                  condition = c("household_income", "variant"),
#                  type = "response")+
#   theme_minimal()
# # 
# plot_predictions(poi_model,
#                  condition = c("household_size", "variant"),
#                  type = "response")+
#   theme_minimal()
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
