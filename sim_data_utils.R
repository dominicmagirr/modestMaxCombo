library(mvtnorm)
library(nphRCT)
library(survival)
library(clustermq)
library(dplyr)

##---------------------------------------------------------------------------------------------

find_cor <- function(formula,
                     data,
                     s_stars = c(1, 0.5),
                     rhos = NULL,
                     gammas = NULL){
  
  
  ## at risk table
  df_events <- find_at_risk(formula = formula,
                            data = data,
                            include_cens = FALSE)
  
  ## extract info
  n_event_g <- as.matrix(df_events[, grep("n_event_",names(df_events))])
  n_event <- df_events$n_event
  n_risk_g <- as.matrix(df_events[, grep("n_risk_",names(df_events))])
  n_risk <- df_events$n_risk
  
  ## weights
  if (!is.null(s_stars)){
    
    w <- purrr::map(s_stars,
                    function(x) find_weights(formula = formula,
                                             data = data,
                                             method="mw",
                                             s_star = x,    
                                             include_cens = FALSE))
  }
  else {
    
    w <- purrr::map2(rhos, gammas,
                    function(x,y) find_weights(formula = formula,
                                             data = data,
                                             method="fh",
                                             rho = x, 
                                             gamma = y,
                                             include_cens = FALSE))
    
  }
  # test statistics
  observed <- purrr::map(w, function(x) x %*% n_event_g)
  expected <- purrr::map(w, function(x) x %*%  (n_risk_g * (n_event / n_risk)))
  u <- purrr::map2(observed, expected, function(x, y) (x - y)[,2])
  
  ## variance
  var_d <- apply(n_risk_g, 1, prod) * n_event * (n_risk - n_event) / (n_risk ^ 2 * (n_risk - 1))
  v_u <- purrr::map(w, function(x) sum(ifelse(n_risk == 1, 0, x^2 * var_d)))
  
  
  
  ## covariance
  n_w <- length(w)
  cov_mat <- matrix(NA, nrow = n_w, ncol = n_w)
  for (i in 1:n_w) {
    for (j in 1:n_w){
      if (i == j){ cov_mat[i, j] <- v_u[[i]]}
      else if (i < j){ cov_mat[i, j] <- sum(ifelse(n_risk == 1, 0, w[[i]] * w[[j]] * var_d))}
      else {cov_mat[i,j] <- cov_mat[j,i]}
    }
  }
  
  return(list(cov_mat = cov_mat,
              v_u = unlist(v_u),
              u = unlist(u)))
}

##---------------------------------------------------------------------------------------------

f_crit <- function(crit, k, alpha, rho_mat){
  
  (1 - alpha) - pmvnorm(lower = c(-Inf, -Inf),
                        upper = crit * qnorm(1 - k * alpha),
                        corr = rho_mat)[1]
  
}

##---------------------------------------------------------------------------------------------

find_crit <- function(rho_mat, k = c(0.5, 0.5), alpha = 0.025, crit_interval = c(0.8, 1)){
  
  if (dim(rho_mat)[1] != length(k)) stop("Dimensions of rho_mat and k must match")
  crit <- uniroot(f_crit, crit_interval, k = k, alpha = alpha, rho_mat = rho_mat)$root
  crit * qnorm(1 - k * alpha)
  
}

##---------------------------------------------------------------------------------------------

one_sim_result <- function(scenario,
                           s_stars = c(1, 0.5),
                           k = c(0.5, 0.5),
                           alpha = 0.025,
                           rhos = NULL,
                           gammas = NULL){
  
  
  
  sim_data <- sim_events_delay(event_model = scenario$event_model,
                               recruitment_model = scenario$recruitment_model,
                               n_c = scenario$n_c,
                               n_e = scenario$n_e,
                               max_cal_t = scenario$max_cal_t)
  
  ##### Calculate covariance
  cov_res <- find_cor(formula = Surv(event_time, event_status) ~ group, data = sim_data, s_stars = s_stars, rhos = rhos, gammas = gammas)
  z <- cov_res$u / sqrt(cov_res$v)
  df <- data.frame(t(z))
  colnames(df) <- paste0("z_", seq_along(z))
  df_2 <- data.frame(t(-find_crit(cov2cor(cov_res$cov_mat), k = k, alpha = alpha)))
  colnames(df_2) <- paste0("c_", seq_along(z))
  cbind(df, df_2, data.frame(scenario = scenario$scenario))
  
  
}

##---------------------------------------------------------------------------------------------



## Simulation wrapper for clustermq
run_sim_wrapper <- function(scenario, 
                            s_stars = c(1, 0.5),
                            k = c(0.5, 0.5),
                            alpha = 0.025,
                            rhos = NULL,
                            gammas = NULL){
  
  if(!exists('worker_setup_complete', mode='logical')) {
    cat("Calling setup function:\n")
    
    source("sim_data_utils.R")
    
    worker_setup_complete <<- TRUE
  } else {
    cat("Calling GC:\n")
    print(gc())
  }
  
  one_sim_result(scenario, s_stars, k, alpha, rhos, gammas)
  
}











