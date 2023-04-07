regiao<- function(x){
  x<-as.data.frame(x)
  x$regiao <- NA
  x[x$name_state %in% c("Rio Grande do Sul", "Santa Catarina", "Paraná"), "regiao"] <- "Sul"
  x[x$name_state %in% c("Amazonas", "Acre", "Amapá", "Pará", "Rondônia", "Roraima", "Tocantins"), "regiao"] <- "Norte"
  x[x$name_state %in% c("Mato Grosso do Sul", "Mato Grosso", "Goiás", "Distrito Federal"), "regiao"] <- "Centro-Oeste"
  x[x$name_state %in% c("São Paulo", "Rio de Janeiro", "Minas Gerais", "Espírito Santo"), "regiao"] <- "Sudeste"
  x[x$name_state %in% c("Bahia", "Pernambuco", "Paraíba", "Maranhão", "Ceará", "Rio Grande do Norte", "Alagoas", "Sergipe", "Piauí"), "regiao"] <- "Nordeste"
  return(x)
}

end.of.epiweek <- function(x, end = 6) {
  offset <- (end - 4) %% 7
  num.x <- as.numeric(x)
  return(x - (num.x %% 7) + offset + ifelse(num.x %% 7 > offset, 7, 0))
}
plot.multinom<- function(data, 
                         is.tt = FALSE, 
                         log.OR.x = FALSE, 
                         names_sg_uf){
  if(!is.tt){
    tt <- broom::tidy(data,conf.int=TRUE)
    tt <- dplyr::filter(tt, term!="(Intercept)")
    tt_sp$term<-as.Date(substr(tt_sp$term, 4, 13))
  }
  if(missing(names_sg_uf)){
    names_sg_uf<-""
  }
  if (log.OR.x){
    plot<-data %>% ggplot(aes(x=estimate,y=term, col = y.level))+
      geom_pointrangeh(aes(xmin=conf.low,
                           xmax=conf.high),
                       position=position_dodgev(height = 0.75))+
      theme_minimal()+
      scale_y_date(name = "Mês", guide = guide_axis(angle = 0), date_breaks = "1 month", date_labels = "%b %Y")+
      labs(title = paste0("Estado ", names_sg_uf), y = "Log (OR)")+
      theme(legend.position = "bottom")
    return(plot)
  }
  else{
    plot<-data %>% ggplot(aes(y=estimate,x=term, col = y.level))+
      geom_pointrange(aes(ymin=conf.low,
                          ymax=conf.high),
                      position=position_dodge(width = 0.75))+
      theme_minimal()+
      scale_x_date(name = "Mês", guide = guide_axis(angle = 90), date_breaks = "1 month", date_labels = "%b %Y")+
      labs(title = paste0("Estado ", names_sg_uf), y = "Log (OR)")+
      theme(legend.position = "bottom")
    return(plot)
  }
}

plot_multinomCI<-function(data, estado = "SP", type.cases = "SRAG", type.info = "Casos"){
  data_ic <- data %>% 
    # filter(dt_sin_pri >= "2020-03-01" & dt_sin_pri < "2020-12-01") %>%
    filter(week >= "2020-03-01" & week < "2020-12-01") %>%
    group_by(age) %>%
    summarize(lower_95ic = quantile(est, probs = 0.025),
              upper_95ic = quantile(est, probs = 0.975),
              lower_50ic = quantile(est, probs = 0.25),
              upper_50ic = quantile(est, probs = 0.75)) %>%
    left_join(data, by = c("age"))
  # data_ic <- data_ic %>%
  #   group_by(age, week) %>%
  #   summarize(N = n()) %>%
  #   left_join(data_ic, by = c("age", "week"))
  
  p<-data_ic %>% 
    filter(week < "2021-05-01") %>%
    ggplot(aes(x = week, 
               y = est, ymin = lwr.ci, ymax = upr.ci, 
               group = age,
               col = age))+
    geom_pointrange(size = .5, 
                    fatten = .75
    )+
    # geom_col(aes(x = week, y = N, group = age, fill = 'darkgrey'), alpha = 0.15)+
    geom_ribbon(aes(ymin=lower_95ic,
                    ymax=upper_95ic),
                color = "grey", fill='darkgrey', linetype=0, alpha = 0.15) +
    geom_ribbon(aes(ymin=lower_50ic, ymax=upper_50ic),
                color = "grey", fill='darkgrey', linetype=0, alpha = 0.3) +
    theme_Publication()+
    labs(title = paste0("Estado de ", estado),
         subtitle = paste0(type.info," de ", type.cases),
         y = "Proporção",
         x = "Semana de Primeiro Sintomas", 
         color = "Faixa Etária", 
         fill = "Faixa Etária", 
         caption = "Parceria Observatório Covid19 BR e MAVE/FioCruz")+
    scale_color_colorblind()+
    scale_fill_colorblind()+
    # scale_fill_viridis_d(aesthetics = c("colour", "fill"),
    #                      option = "viridis",
    #                      direction = -1)+
    theme(legend.position = "none", 
          axis.text.x = element_text(angle = 90), 
          # strip.text = element_blank()
    )+
    # scale_x_date(date_labels = "%W %y", 
    #              date_breaks = "5 weeks")+
    facet_grid(age ~., scales = "free_y")
  return(p)
}

multinom_ci <- function(data){
  mCI<-DescTools::MultinomCI(data$N, conf.level = 0.95)
  mCI_merge<-cbind(data, estimate = mCI[,1], lower = mCI[,2], upper = mCI[,3])
  return(mCI_merge)
}

f1<- function(x, stratum, ...){
  m1 <- DescTools::MultinomCI(x, ...)
  cbind(stratum, data.frame(m1))
}
