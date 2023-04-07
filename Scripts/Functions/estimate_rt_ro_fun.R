################################################################################
## Source File for Cori_method.Rmd
################################################################################

################################################################################
## Carregar bibliotecas
################################################################################

library(zoo)
library(EpiEstim)
library(ggplot2)
library(coarseDataTools)
library(tidyverse)
library("ggpubr")
library(dplyr)
library(incidence)
library(coda)
library(plyr)
library(viridis)
library(lubridate)

################################################################################
## Parametros de formatacao comum aos plots
################################################################################
plot.formatos <- theme_bw()+
  theme(axis.text= element_text(size=10, face="bold"),
        axis.title = element_text(size=10, face="bold"),
        legend.text = element_text(size=12),
        title = element_text(size = 12),
        plot.margin = margin(5, 0, 0, 0, "pt"))

################################################################################
## Carregar Serial Interval (Nishiura et al 2020)
################################################################################
nd <- read.table("Data/serial_interval/nishi_si_table.txt", header = TRUE)
nishi_si <- read.table("Data/serial_interval/nish_si_posterior.txt", header = TRUE)

################################################################################
## Fun??es utilizadas
################################################################################

#### Fun??o para estimativa de Re baseada no m?todo de Cori et al. 2013

estimate.R0.cori <- function(novos.casos, day0 = NA, delay=7, method, parameter.table, 
                             p.distribution, bayes.control, si.data, si.sample, modified = FALSE, n2 = NA, ...) {
  
  if(is.na(day0)){
    
    mdif = mean(si.data$SL - si.data$EL) # empirical SI
    day0= min(which(cumsum(novos.casos) > 12 & 
                      1:length(novos.casos) > ceiling(mdif) 
                    & 1:length(novos.casos) > delay)) 
    
  }
  
  if(method == "uncertain_si")  {
    
    parameter.list  <- lapply(parameter.table, function(x) x=x)
    parameter.list$si_parametric_distr = p.distribution 
    parameter.list$t_start = seq(day0,length(novos.casos)-delay)
    parameter.list$t_end = seq(delay+day0,length(novos.casos))
    
    config <- make_config(parameter.list)
    return(estimate_R(novos.casos, method = method, config = config))
  }
  
  if(method == "parametric_si"){
    
    parameter.list  <- lapply(parameter.table, function(x) x=x)
    parameter.list$si_parametric_distr =  p.distribution
    parameter.list$t_start = seq(day0,length(novos.casos)-delay)
    parameter.list$t_end = seq(delay+day0,length(novos.casos))
    
    config <- make_config(parameter.list)
    return(estimate_R(novos.casos, method = "parametric_si", config = config))
  }
  
  if (method == "si_from_data") {
    config <- make_config(incid = novos.casos, 
                          method = method, 
                          list(si_parametric_distr = p.distribution,
                               mcmc_control = make_mcmc_control(burnin = bayes_control$burnin, 
                                                                thin = bayes_control$thin),
                               mean_prior = 5,
                               std_prior = 5,
                               n1 = bayes_control$n1, 
                               n2 = bayes_control$n2,
                               t_start = seq(day0,length(novos.casos)-delay),
                               t_end = seq(delay+day0,length(novos.casos))))
    
    return(estimate_R(novos.casos, method = method, 
                      si_data = si.data, config = config))
  }
  
  if(method == "si_from_sample"){
    
    parameter.list  <- list()
    parameter.list$t_start = seq(day0,length(novos.casos)-delay)
    parameter.list$t_end = seq(delay+day0,length(novos.casos))
    
    if(is.na(n2)){
      parameter.list$n2 = 50 } else {
        parameter.list$n2 = n2  
      }
    
    config <- make_config(parameter.list)
    
    if(modified == FALSE){
      return(estimate_R(novos.casos, method = method,
                        si_sample = si.sample, 
                        config = config)) } else {
                          return(estimate_R_modified(novos.casos, method = method,
                                                     si_sample = si.sample, 
                                                     config = config))
                        }
  }
  
}

estimate_R_modified <- function(incid,
                                method = c(
                                  "non_parametric_si", "parametric_si",
                                  "uncertain_si", "si_from_data",
                                  "si_from_sample"
                                ),
                                si_data = NULL,
                                si_sample = NULL,
                                config = make_config(incid = incid, method = method)) {
  
  method <- match.arg(method)
  config <- make_config(incid = incid, method = method, config = config)
  config <- EpiEstim:::process_config(config)
  EpiEstim:::check_config(config, method)
  
  if (method == "si_from_data") {
    ## Warning if the expected set of parameters is not adequate
    si_data <- EpiEstim:::process_si_data(si_data)
    config <- EpiEstim:::process_config_si_from_data(config, si_data)
    
    ## estimate serial interval from serial interval data first
    if (!is.null(config$mcmc_control$seed)) {
      cdt <- dic.fit.mcmc(
        dat = si_data,
        dist = config$si_parametric_distr,
        burnin = config$mcmc_control$burnin,
        n.samples = config$n1 * config$mcmc_control$thin,
        init.pars = config$mcmc_control$init_pars,
        seed = config$mcmc_control$seed
      )
    } else {
      cdt <- dic.fit.mcmc(
        dat = si_data,
        dist = config$si_parametric_distr,
        burnin = config$mcmc_control$burnin,
        n.samples = config$n1 * config$mcmc_control$thin,
        init.pars = config$mcmc_control$init_pars
      )
    }
    
    ## check convergence of the MCMC and print warning if not converged
    MCMC_conv <- check_cdt_samples_convergence(cdt@samples)
    
    ## thin the chain, and turn the two parameters of the SI distribution into a
    ## whole discrete distribution
    c2e <- coarse2estim(cdt, thin = config$mcmc_control$thin)
    
    cat(paste(
      "\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@",
      "\nEstimating the reproduction number for these serial interval",
      "estimates...\n",
      "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
    ))
    
    ## then estimate R for these serial intervals
    
    if (!is.null(config$seed)) {
      set.seed(config$seed)
    }
    
    out <- estimate_R_func(
      incid = incid,
      method = "si_from_data",
      si_sample = c2e$si_sample,
      config = config
    )
    out[["MCMC_converged"]] <- MCMC_conv
  } else {
    if (!is.null(config$seed)) {
      set.seed(config$seed)
    }
    
    out <- estimate_R_func(
      incid = incid, method = method, si_sample = si_sample,
      config = config
    )
  }
  return(out)
}

##########################################################
## estimate_R_func: Doing the heavy work in estimate_R  ##
##########################################################

#'
#' @importFrom stats median qgamma quantile rnorm sd
#'
#' @importFrom incidence as.incidence
#'
estimate_R_func <- function(incid,
                            si_sample,
                            method = c(
                              "non_parametric_si", "parametric_si",
                              "uncertain_si", "si_from_data", "si_from_sample"
                            ),
                            config) {
  
  #########################################################
  # Calculates the cumulative incidence over time steps   #
  #########################################################
  
  calc_incidence_per_time_step <- function(incid, t_start, t_end) {
    nb_time_periods <- length(t_start)
    incidence_per_time_step <- EpiEstim:::vnapply(seq_len(nb_time_periods), function(i) 
      sum(incid[seq(t_start[i], t_end[i]), c("local", "imported")]))
    return(incidence_per_time_step)
  }
  
  #########################################################
  # Calculates the parameters of the Gamma posterior      #
  # distribution from the discrete SI distribution        #
  #########################################################
  
  posterior_from_si_distr <- function(incid, si_distr, a_prior, b_prior,
                                      t_start, t_end) {
    nb_time_periods <- length(t_start)
    lambda <- overall_infectivity(incid, si_distr)
    final_mean_si <- sum(si_distr * (seq(0, length(si_distr) -
                                           1)))
    a_posterior <- vector()
    b_posterior <- vector()
    a_posterior <- EpiEstim:::vnapply(seq_len(nb_time_periods), function(t) if (t_end[t] >
                                                                                final_mean_si) {
      a_prior + sum(incid[seq(t_start[t], t_end[t]), "local"]) 
      ## only counting local cases on the "numerator"
    }
    else {
      NA
    })
    b_posterior <- EpiEstim:::vnapply(seq_len(nb_time_periods), function(t) if (t_end[t] >
                                                                                final_mean_si) {
      1 / (1 / b_prior + sum(lambda[seq(t_start[t], t_end[t])]))
    }
    else {
      NA
    })
    return(list(a_posterior, b_posterior))
  }
  
  #########################################################
  # Samples from the Gamma posterior distribution for a   #
  # given mean SI and std SI                              #
  #########################################################
  
  sample_from_posterior <- function(sample_size, incid, mean_si, std_si,
                                    si_distr = NULL,
                                    a_prior, b_prior, t_start, t_end) {
    nb_time_periods <- length(t_start)
    
    if (is.null(si_distr)) {
      si_distr <- discr_si(seq(0, T - 1), mean_si, std_si)
    }
    
    final_mean_si <- sum(si_distr * (seq(0, length(si_distr) -
                                           1)))
    lambda <- overall_infectivity(incid, si_distr)
    a_posterior <- vector()
    b_posterior <- vector()
    a_posterior <- EpiEstim:::vnapply(seq_len(nb_time_periods), function(t) if (t_end[t] >
                                                                                final_mean_si) {
      a_prior + sum(incid[seq(t_start[t], t_end[t]), "local"]) 
      ## only counting local cases on the "numerator"
    }
    else {
      NA
    })
    b_posterior <- EpiEstim:::vnapply(seq_len(nb_time_periods), function(t) if (t_end[t] >
                                                                                final_mean_si) {
      1 / (1 / b_prior + sum(lambda[seq(t_start[t], t_end[t])], na.rm = TRUE))
    }
    else {
      NA
    })
    sample_r_posterior <- vapply(seq_len(nb_time_periods), function(t) 
      if (!is.na(a_posterior[t])) {
        rgamma(sample_size,
               shape = unlist(a_posterior[t]),
               scale = unlist(b_posterior[t])
        )
      }
      else {
        rep(NA, sample_size)
      }, numeric(sample_size))
    if (sample_size == 1L) {
      sample_r_posterior <- matrix(sample_r_posterior, nrow = 1)
    }
    return(list(sample_r_posterior, si_distr))
  }
  
  method <- match.arg(method)
  
  incid <- EpiEstim:::process_I(incid)
  T <- nrow(incid)
  
  a_prior <- (config$mean_prior / config$std_prior)^2
  b_prior <- config$std_prior^2 / config$mean_prior
  
  EpiEstim:::check_times(config$t_start, config$t_end, T)
  nb_time_periods <- length(config$t_start)
  
  if (method == "si_from_sample") {
    if (is.null(config$n2)) {
      stop("method si_from_sample requires to specify the config$n2 argument.")
    }
    si_sample <- EpiEstim:::process_si_sample(si_sample)
  }
  
  min_nb_cases_per_time_period <- ceiling(1 / config$cv_posterior^2 - a_prior)
  incidence_per_time_step <- calc_incidence_per_time_step(
    incid, config$t_start,
    config$t_end
  )
  if (incidence_per_time_step[1] < min_nb_cases_per_time_period) {
    warning("You're estimating R too early in the epidemic to get the desired
            posterior CV.")
  }
  
  if (method == "non_parametric_si") {
    si_uncertainty <- "N"
    parametric_si <- "N"
  }
  if (method == "parametric_si") {
    si_uncertainty <- "N"
    parametric_si <- "Y"
  }
  if (method == "uncertain_si") {
    si_uncertainty <- "Y"
    parametric_si <- "Y"
  }
  if (method == "si_from_data" | method == "si_from_sample") {
    si_uncertainty <- "Y"
    parametric_si <- "N"
  }
  if (si_uncertainty == "Y") {
    if (parametric_si == "Y") {
      mean_si_sample <- rep(-1, config$n1)
      std_si_sample <- rep(-1, config$n1)
      for (k in seq_len(config$n1)) {
        while (mean_si_sample[k] < config$min_mean_si || mean_si_sample[k] >
               config$max_mean_si) {
          mean_si_sample[k] <- rnorm(1,
                                     mean = config$mean_si,
                                     sd = config$std_mean_si
          )
        }
        while (std_si_sample[k] < config$min_std_si || std_si_sample[k] >
               config$max_std_si) {
          std_si_sample[k] <- rnorm(1, mean = config$std_si,
                                    sd = config$std_std_si)
        }
      }
      temp <- lapply(seq_len(config$n1), function(k) sample_from_posterior(config$n2,
                                                                           incid, mean_si_sample[k], std_si_sample[k],
                                                                           si_distr = NULL, a_prior,
                                                                           b_prior, config$t_start, config$t_end
      ))
      config$si_distr <- cbind(
        t(vapply(seq_len(config$n1), function(k) (temp[[k]])[[2]], numeric(T))),
        rep(0, config$n1)
      )
      r_sample <- matrix(NA, config$n2 * config$n1, nb_time_periods)
      for (k in seq_len(config$n1)) {
        r_sample[seq((k - 1) * config$n2 + 1, k * config$n2), which(config$t_end >
                                                                      mean_si_sample[k])] <- (temp[[k]])[[1]][, which(config$t_end >
                                                                                                                        mean_si_sample[k])]
      }
      mean_posterior <- apply(r_sample, 2, mean, na.rm = TRUE)
      std_posterior <- apply(r_sample, 2, sd, na.rm = TRUE)
      quantile_0.025_posterior <- apply(r_sample, 2, quantile,
                                        0.025,
                                        na.rm = TRUE
      )
      quantile_0.05_posterior <- apply(r_sample, 2, quantile,
                                       0.05,
                                       na.rm = TRUE
      )
      quantile_0.25_posterior <- apply(r_sample, 2, quantile,
                                       0.25,
                                       na.rm = TRUE
      )
      median_posterior <- apply(r_sample, 2, median, na.rm = TRUE)
      quantile_0.25_posterior <- apply(r_sample, 2, quantile,
                                       0.75,
                                       na.rm = TRUE
      )
      quantile_0.25_posterior <- apply(r_sample, 2, quantile,
                                       0.95,
                                       na.rm = TRUE
      )
      quantile_0.975_posterior <- apply(r_sample, 2, quantile,
                                        0.975,
                                        na.rm = TRUE
      )
    }
    else {
      config$n1 <- dim(si_sample)[2]
      mean_si_sample <- rep(-1, config$n1)
      std_si_sample <- rep(-1, config$n1)
      for (k in seq_len(config$n1)) {
        mean_si_sample[k] <- sum((seq_len(dim(si_sample)[1]) - 1) * 
                                   si_sample[, k])
        std_si_sample[k] <- sqrt(sum(si_sample[, k] * 
                                       ((seq_len(dim(si_sample)[1]) - 1) - 
                                          mean_si_sample[k])^2))
      }
      temp <- lapply(seq_len(config$n1), function(k) sample_from_posterior(config$n2,
                                                                           incid,
                                                                           mean_si = NULL, std_si = NULL, si_sample[, k], a_prior,
                                                                           b_prior, config$t_start, config$t_end
      ))
      config$si_distr <- cbind(
        t(vapply(seq_len(config$n1), function(k) (temp[[k]])[[2]], 
                 numeric(nrow(si_sample)))),
        rep(0, config$n1)
      )
      r_sample <- matrix(NA, config$n2 * config$n1, nb_time_periods)
      for (k in seq_len(config$n1)) {
        r_sample[seq((k - 1) * config$n2 + 1,k * config$n2), which(config$t_end >
                                                                     mean_si_sample[k])] <- (temp[[k]])[[1]][, which(config$t_end >
                                                                                                                       mean_si_sample[k])]
      }
      mean_posterior <- apply(r_sample, 2, mean, na.rm = TRUE)
      std_posterior <- apply(r_sample, 2, sd, na.rm = TRUE)
      quantile_0.025_posterior <- apply(r_sample, 2, quantile,
                                        0.025,
                                        na.rm = TRUE
      )
      quantile_0.05_posterior <- apply(r_sample, 2, quantile,
                                       0.05,
                                       na.rm = TRUE
      )
      quantile_0.25_posterior <- apply(r_sample, 2, quantile,
                                       0.25,
                                       na.rm = TRUE
      )
      median_posterior <- apply(r_sample, 2, median, na.rm = TRUE)
      quantile_0.25_posterior <- apply(r_sample, 2, quantile,
                                       0.75,
                                       na.rm = TRUE
      )
      quantile_0.25_posterior <- apply(r_sample, 2, quantile,
                                       0.95,
                                       na.rm = TRUE
      )
      quantile_0.975_posterior <- apply(r_sample, 2, quantile,
                                        0.975,
                                        na.rm = TRUE
      )
    }
  } else {
    # CertainSI
    if (parametric_si == "Y") {
      config$si_distr <- discr_si(seq(0,T - 1), config$mean_si, config$std_si)
    }
    if (length(config$si_distr) < T + 1) {
      config$si_distr[seq(length(config$si_distr) + 1,T + 1)] <- 0
    }
    final_mean_si <- sum(config$si_distr * (seq(0,length(config$si_distr) -
                                                  1)))
    Finalstd_si <- sqrt(sum(config$si_distr * (seq(0,length(config$si_distr) -
                                                     1))^2) - final_mean_si^2)
    post <- posterior_from_si_distr(
      incid, config$si_distr, a_prior, b_prior,
      config$t_start, config$t_end
    )
    
    a_posterior <- unlist(post[[1]])
    b_posterior <- unlist(post[[2]])
    mean_posterior <- a_posterior * b_posterior
    std_posterior <- sqrt(a_posterior) * b_posterior
    quantile_0.025_posterior <- qgamma(0.025,
                                       shape = a_posterior,
                                       scale = b_posterior, lower.tail = TRUE, log.p = FALSE
    )
    quantile_0.05_posterior <- qgamma(0.05,
                                      shape = a_posterior,
                                      scale = b_posterior, lower.tail = TRUE, log.p = FALSE
    )
    quantile_0.25_posterior <- qgamma(0.25,
                                      shape = a_posterior,
                                      scale = b_posterior, lower.tail = TRUE, log.p = FALSE
    )
    median_posterior <- qgamma(0.5,
                               shape = a_posterior,
                               scale = b_posterior, lower.tail = TRUE, log.p = FALSE
    )
    quantile_0.25_posterior <- qgamma(0.75,
                                      shape = a_posterior,
                                      scale = b_posterior, lower.tail = TRUE, log.p = FALSE
    )
    quantile_0.25_posterior <- qgamma(0.95,
                                      shape = a_posterior,
                                      scale = b_posterior, lower.tail = TRUE, log.p = FALSE
    )
    quantile_0.975_posterior <- qgamma(0.975,
                                       shape = a_posterior,
                                       scale = b_posterior, lower.tail = TRUE, log.p = FALSE
    )
  }
  
  results <- list(R = as.data.frame(cbind(
    config$t_start, config$t_end, mean_posterior,
    std_posterior, quantile_0.025_posterior, quantile_0.05_posterior,
    quantile_0.25_posterior, median_posterior, quantile_0.25_posterior,
    quantile_0.25_posterior, quantile_0.975_posterior
  )))
  
  names(results$R) <- c(
    "t_start", "t_end", "Mean(R)", "Std(R)",
    "Quantile.0.025(R)", "Quantile.0.05(R)", "Quantile.0.25(R)",
    "Median(R)", "Quantile.0.75(R)", "Quantile.0.95(R)",
    "Quantile.0.975(R)"
  )
  results$method <- method
  results$si_distr <- config$si_distr
  if (is.matrix(results$si_distr)) {
    colnames(results$si_distr) <- paste0("t", seq(0,ncol(results$si_distr) - 1))
  } else {
    names(results$si_distr) <- paste0("t", seq(0,length(results$si_distr) - 1))
  }
  if (si_uncertainty == "Y") {
    results$SI.Moments <- as.data.frame(cbind(
      mean_si_sample,
      std_si_sample
    ))
  } else {
    results$SI.Moments <- as.data.frame(cbind(
      final_mean_si,
      Finalstd_si
    ))
  }
  names(results$SI.Moments) <- c("Mean", "Std")
  
  
  if (!is.null(incid$dates)) {
    results$dates <- EpiEstim:::check_dates(incid)
  } else {
    results$dates <- seq_len(T)
  }
  results$I <- rowSums(incid[, c("local", "imported")])
  results$I_local <- incid$local
  results$I_imported <- incid$imported
  
  results$r_sample <- r_sample
  
  class(results) <- "estimate_R"
  return(results)
}

default.R.cori <- partial(estimate.R0.cori,
                          delay = 7,
                          day0 = NA,
                          method = "si_from_sample",
                          si.data = nd, 
                          si.sample = nishi_si[,sample(1:ncol(nishi_si), 1)], # Samples 1 SI interval
                          p.distribution = "G",
                          modified = TRUE,
                          n2 = 50)

default.R.cori.k <- partial(estimate.R0.cori,
                            delay = k.partial,
                            day0 = NA,
                            method = "si_from_sample",
                            si.data = nd, 
                            si.sample = nishi_si[,sample(1:ncol(nishi_si), 1)], # Samples 1 SI interval
                            p.distribution = "G",
                            modified = TRUE,
                            n2 = 50)

default.R.cori.day0eq2 <- partial(estimate.R0.cori,
                                  delay = 7,
                                  day0 = 2,
                                  method = "si_from_sample",
                                  si.data = nd, 
                                  si.sample = nishi_si[,sample(1:ncol(nishi_si), 1)], # Samples 1 SI interval
                                  p.distribution = "G",
                                  modified = TRUE,
                                  n2 = 50)

default.R.cori.day0eq10 <- partial(estimate.R0.cori,
                                   delay = 7,
                                   day0 = 10,
                                   method = "si_from_sample",
                                   si.data = nd, 
                                   si.sample = nishi_si[,sample(1:ncol(nishi_si), 1)], # Samples 1 SI interval
                                   p.distribution = "G",
                                   modified = TRUE,
                                   n2 = 50)

default.R.cori.k2.day0 <- partial(estimate.R0.cori,
                                  delay = k.partial,
                                  day0 = day0,
                                  method = "si_from_sample",
                                  si.data = nd, 
                                  si.sample = nishi_si[,sample(1:ncol(nishi_si), 1)], # Samples 1 SI interval
                                  p.distribution = "G",
                                  modified = TRUE,
                                  n2 = 50)


#### Plotagem dos resultados da estimativa de Re

plot.estimate.R0 <- function(res.si, incidence.data, plot = TRUE){
  if(class(incidence.data) != "Date"){
    incidence.data <- as.Date(incidence.data, "%Y-%m-%d")
  }
  
  ## Converter para o formato zoo
  res.zoo <- zoo(res.si$R,  incidence.data[res.si$R$t_end])
  names(res.zoo) <- gsub("\\(R\\)", ".R", names(res.zoo))
  
  ## Plotar grafico
  ggplot.R0 <- 
    ggplot(data = res.zoo, aes(Index, Mean.R)) +
    geom_ribbon(aes(ymin = Quantile.0.025.R, ymax = Quantile.0.975.R), fill="lightgrey") +
    geom_line(size = 1.25, color="darkblue") +
    scale_x_date( date_labels = "%d/%b", name="") +
    ylim(0, max(res.zoo$Quantile.0.975.R)) +
    geom_hline(yintercept=1, linetype="dashed", col="red") +
    ylab("Numero de Reproducao") +
    plot.formatos
  
  if(plot == TRUE){
    print(ggplot.R0)
  } else
  {return(ggplot.R0)}
}

### Plotar distribui??es de SI estimados a partir dos dados emp?ricos:

plot.SI <- function(SItable, xlim = NA, ylim = NA, main = "") {
  if(is.na(xlim[1])){
    xlim = c(0,ncol(SItable))  
  }
  
  if(is.na(ylim[1])){
    ylim =  c(0,max(SItable)) 
  }  
  
  SItable = t(as.matrix(SItable) )
  
  plot(1,1, type = "n", xlim = xlim, ylim = ylim, ylab = "", xlab = "Tempo (dias)", main = main)
  
  for(i in 1:nrow(SItable)){
    lines(t(SItable)[i,], col = 2) #alpha("red", 1))
  }
}

### Retorna uma tabela de incid?ncia do n?mero de casos com as datas de forma cont?nua

return.incidence <- function(count, dates){
  time_sequence = seq(first(dates), last(dates), by="days")
  df = data.frame(row.names = time_sequence, onset = time_sequence, n.casos = rep(0,length(time_sequence)))
  df[as.character(dates),2] = count
  return(df)
}

# Create sample of SI estimates using bayesian estimation

sample.si <- function(si.table, dist, n1 = 50) {
  
  mcmc_control <- make_mcmc_control(burnin = 1000)
  config <- make_config(list(si_parametric_distr = dist,
                             mcmc_control = mcmc_control,
                             n1 = n1, 
                             n2 = 50))
  n_mcmc_samples <- config$n1 * mcmc_control$thin
  
  si_data = cbind(si.table, type = 0)
  
  SI_fit <- coarseDataTools::dic.fit.mcmc(dat = si_data,
                                          dist = dist,
                                          init.pars = init_mcmc_params(si_data, dist),
                                          burnin = mcmc_control$burnin,
                                          n.samples = n_mcmc_samples)
  
  si_sample <- coarse2estim(SI_fit, thin = mcmc_control$thin)$si_sample
  return(si_sample)
}

data[-which(1 %in% c("A", "B")),]


# Estimate posterior probabilities for nowcasting AND SI variation
# Estimate posterior probabilities for nowcasting AND SI variation
posteriors <- function(data){
  r_sample = data[-which(row.names(data) %in% c("t_start","t_end"))]
  mean_posterior <- apply(r_sample, 2, mean, na.rm = TRUE)
  std_posterior <- apply(r_sample, 2, sd, na.rm = TRUE)
  quantile_0.025_posterior <- apply(r_sample, 2, quantile,
                                    0.025,
                                    na.rm = TRUE
  )
  quantile_0.05_posterior <- apply(r_sample, 2, quantile,
                                   0.05,
                                   na.rm = TRUE
  )
  quantile_0.25_posterior <- apply(r_sample, 2, quantile,
                                   0.25,
                                   na.rm = TRUE
  )
  median_posterior <- apply(r_sample, 2, median, na.rm = TRUE)
  quantile_0.75_posterior <- apply(r_sample, 2, quantile,
                                   0.75,
                                   na.rm = TRUE
  )
  quantile_0.95_posterior <- apply(r_sample, 2, quantile,
                                   0.95,
                                   na.rm = TRUE
  )
  quantile_0.975_posterior <- apply(r_sample, 2, quantile,
                                    0.975,
                                    na.rm = TRUE
  )
  
  results <- list(R = as.data.frame(
    cbind(as.numeric(data["t_start", ]),
          as.numeric(data["t_end", ]),
          mean_posterior,
          std_posterior,
          quantile_0.025_posterior,
          quantile_0.05_posterior,
          quantile_0.25_posterior,
          median_posterior,
          quantile_0.75_posterior,
          quantile_0.95_posterior,
          quantile_0.975_posterior)))
  
  names(results$R) <- c(
    "t_start",
    "t_end",
    "Mean(R)",
    "Std(R)",
    "Quantile.0.025(R)",
    "Quantile.0.05(R)",
    "Quantile.0.25(R)",
    "Median(R)",
    "Quantile.0.75(R)",
    "Quantile.0.95(R)",
    "Quantile.0.975(R)"
  )
  
  return(results)  
}

####

fill.dates <- function(data, column){
  
  if(class(column) == "character"){cases = data[,column]}
  if(class(column) == "numeric"){cases = data[,column]}
  
  
  sequence = seq(min(data$onset), max(data$onset), by = 1)
  incidence = rep(0, length(sequence))
  
  for(i in 1:length(cases)){
    incidence[which(data$onset[i] == sequence)] = cases[i]
  }
  
  return(data.frame(dates = sequence, incidence = incidence))
}

#### Creates a data.frame where the notified number of daily cases is replaced by the nowcasting series
#### in their corresponding dates.
#### nowcasting: a data.frame containing multiple nowcasting series
#### notified: an incidence table with dates and number of cases per day

join_nowcasting <- function(nowcasting, notified){
  
  prev.dates <- notified$onset[!(notified$onset %in% nowcasting$date)]
  
  prev.cases = notified$n.casos[!(notified$onset %in% nowcasting$date)]
  
  prev.cases.table <- data.frame(matrix(rep(prev.cases,each = ncol(nowcasting)-1), ncol = ncol(nowcasting)-1, byrow=TRUE))
  
  previous <- cbind(date = prev.dates, prev.cases.table)
  
  colnames(previous)[-1] <- colnames(nowcasting)[-1]
  
  all <- rbind(nowcasting, previous)
  
  all <- all %>% arrange(date)
  
  all <- all %>% filter(all$date <= max(nowcasting$date))
  
  return(all)
}

#### Plotagem dos resultados da estimativa de Re

plot.estimate.Rt <- function(res.si, incidence.data, color.db="darkblue", xlim.days=80, plot = TRUE){
  
  if(class(incidence.data$onset) != "Date"){
    incidence.data$onset <- as.Date(incidence.data$onset, "%Y-%m-%d")
    incidence.data$dt_base <- as.Date(incidence.data$dt_base, "%Y-%m-%d")
  }
  
  res.si.f=data.frame()
  for (i in 1:length(names(res.si))){
    aux=cbind(res.si[[i]]$R,database=rep(names(res.si)[[i]],dim(res.si[[i]]$R)[1]))
    aux.ord=incidence.data[incidence.data$dt_base==names(res.si)[[i]],]
    aux=cbind(aux,onset=aux.ord$onset[aux$t_end])
    res.si.f=rbind(res.si.f, aux)
  }
  
  names(res.si.f) <- gsub("\\(R\\)", ".R", names(res.si.f))
  ggplot.Rt <- 
    ggplot(data = res.si.f, aes(x=onset,y=as.numeric(Mean.R), fill=factor(database),color=factor(database))) +
    geom_ribbon(aes(ymin = as.numeric(Quantile.0.025.R), ymax = as.numeric(Quantile.0.975.R)), alpha=0.3, color=NA) +
    geom_line() +
    scale_fill_viridis(discrete=TRUE) +
    scale_color_viridis(discrete=TRUE) +
    theme(legend.position="none")+
    geom_hline(yintercept=1, linetype="dashed", col="black") +
    scale_x_date( date_labels = "%d/%b", name="",limits = c(Sys.Date()-xlim.days, Sys.Date())) +
    ylab("Numero de Reproducao") 
  
  
  if(plot == TRUE){
    print(ggplot.Rt)
  } else
  {return(ggplot.Rt)}
}


#' Calcula R efetivo sobre as trajetórias de nowcasting retornadas pela função NobBS.posterior
#' @param R.method função parcialmente aplicada que calcula Re de uma trajetória
#' @param trajectories data.frame com trajetórias de nowcasting
#' @param Nsamples int número de amostrar das trajetórias de nowcasting
#' @param quantiles logical retorna quantis da distribuição posterior de R ou as osteriores completas?
#' @param ... argumentos extras são repassados pro ldply: .parallel=T é especialmente útil
Re.nowcasting <- function(R.method,
                          trajectories,
                          Nsamples,
                          quantiles = TRUE,
                          ...) {
  N <- ncol(trajectories)
  if (missing(Nsamples) | Nsamples > N - 1)
    Nsamples <- N - 1
  ## trajetórias são auto-correlacionas, melhor sampling é com maior distância
  ## possível entre os índices - usamos intervalos regulares
  samples <- round(seq(2, N, length.out = Nsamples))
  fun <- function(traj){
    casos <- fill.dates(data.frame(onset = trajectories$date, n.casos = traj), 2)
    res <- R.method(casos$incidence)$r_sample
  }
  re.df <- ldply(trajectories[, samples], fun, .id = NULL, ...)
  # run a single time to get full output
  casos <- fill.dates(data.frame(onset = trajectories$date, n.casos = trajectories[, 2]), 2)
  res <- R.method(casos$incidence)
  re.df["t_start", ] <- res$R$t_start
  re.df["t_end", ] <- res$R$t_end
  
  if (!quantiles)
    return(re.df)
  # Calculate quantiles and moments
  results <- posteriors(re.df)
}