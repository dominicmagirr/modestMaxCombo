source("sim_data_utils.R")
##------------------------------------------------
common_recruitment_model <- list(rec_model="power",
                                 rec_period=12,
                                 rec_power=1)
##------------------------------------------------
## original scenarios based on an oncology study
##------------------------------------------------
event_model_1 <- list(duration_c = 24, duration_e = c(6 , 18), lambda_c =  log(2)/15, lambda_e =  c(log(2)/15, log(2)/24))
event_model_2 <- list(duration_c = 24, duration_e = c(6 , 18), lambda_c =  log(2)/15, lambda_e =  c(log(2)/19, log(2)/19))
event_model_3 <- list(duration_c = 24, duration_e = c(9, 9, 6), lambda_c =  log(2)/15, lambda_e =  c(log(2)/22, log(2)/17, log(2)/10))
event_model_4 <- list(duration_c = 24, duration_e = c(6 ,18), lambda_c =  log(2)/15, lambda_e =  c(log(2)/15, log(2)/15))
event_model_5 <- list(duration_c = c(2,4,18), duration_e = c(2 ,22), lambda_c =  c(log(2)/14, log(2)/10, log(2)/15) , lambda_e =  c(log(2)/7, log(2)/15))
##------------------------------------------------
scenarios <- list(
  list(scenario = 1, recruitment_model = common_recruitment_model, n_c = 500, n_e = 500, max_cal_t = 24, event_model = event_model_1),
  list(scenario = 2, recruitment_model = common_recruitment_model, n_c = 500, n_e = 500, max_cal_t = 24, event_model = event_model_2),
  list(scenario = 3, recruitment_model = common_recruitment_model, n_c = 500, n_e = 500, max_cal_t = 24, event_model = event_model_3),
  list(scenario = 4, recruitment_model = common_recruitment_model, n_c = 500, n_e = 500, max_cal_t = 24, event_model = event_model_4),
  list(scenario = 5, recruitment_model = common_recruitment_model, n_c = 500, n_e = 500, max_cal_t = 24, event_model = event_model_5)
)

##-----------------------------------------------
## scenarios more like in cardiovascular studies
##-----------------------------------------------
event_model_cv_1 <- list(duration_c = 36, duration_e = c(6 ,30), lambda_c = .1 * log(2)/15, lambda_e = .1 * c(log(2)/15, log(2)/19.7))
event_model_cv_2 <- list(duration_c = 36, duration_e = c(6 ,30), lambda_c = .1 * log(2)/15, lambda_e = .1 * c(log(2)/18.5, log(2)/18.5))
event_model_cv_3 <- list(duration_c = 36, duration_e = c(9, 9, 18), lambda_c = .1 * log(2)/15, lambda_e = .1 * c(log(2)/33, log(2)/24, log(2)/12))
event_model_cv_4 <- list(duration_c = 36, duration_e = c(6 ,30), lambda_c = .1 * log(2)/15, lambda_e = .1 * c(log(2)/15, log(2)/15))
event_model_cv_5 <- list(duration_c = c(4,9,23), duration_e = c(4 ,32), lambda_c = .1 * c(log(2)/18, log(2)/9, log(2)/15) , lambda_e = .1 * c(log(2)/6, log(2)/15))
##-----------------------------------------------
scenarios_cv <- list(
  list(scenario = 1, recruitment_model = common_recruitment_model, n_c = 3000, n_e = 3000, max_cal_t = 36, event_model = event_model_cv_1),
  list(scenario = 2, recruitment_model = common_recruitment_model, n_c = 3000, n_e = 3000, max_cal_t = 36, event_model = event_model_cv_2),
  list(scenario = 3, recruitment_model = common_recruitment_model, n_c = 3000, n_e = 3000, max_cal_t = 36, event_model = event_model_cv_3),
  list(scenario = 4, recruitment_model = common_recruitment_model, n_c = 3000, n_e = 3000, max_cal_t = 36, event_model = event_model_cv_4),
  list(scenario = 5, recruitment_model = common_recruitment_model, n_c = 3000, n_e = 3000, max_cal_t = 36, event_model = event_model_cv_5)
)
##-----------------------------------------------
## run simulations
##-----------------------------------------------
sim_func <- function(scenarios, N_runs, n_jobs, k, s_stars, rhos = NULL, gammas = NULL){
  
  sims <- Q(run_sim_wrapper, 
            scenario=rep(scenarios, each=N_runs), 
            n_jobs = n_jobs, 
            const = list(k = k,
                         s_stars = s_stars,
                         rhos = rhos,
                         gammas = gammas))
  
  sims <- do.call(rbind, sims)
  
  sims |> 
    group_by(scenario) |> 
    summarise(power_test_1 = mean(z_1 < qnorm(0.025)), 
              power_test_2 = mean(z_2 < qnorm(0.025)), 
              power_combo = mean(z_1 < c_1 | z_2 < c_2))
  
  
}
set.seed(3553)
sim_func(scenarios = scenarios, N_runs = 1e4, n_jobs = 200, k = c(0.5, 0.5), s_stars = c(1, 0.5))
sim_func(scenarios = scenarios, N_runs = 1e4, n_jobs = 200, k = c(0.5, 0.5), s_stars = NULL, rhos = c(0,0), gammas = c(0,0.5))
sim_func(scenarios = scenarios_cv, N_runs = 1e4, n_jobs = 200, k = c(0.5, 0.5), s_stars = c(1, 0.5))
sim_func(scenarios = scenarios_cv, N_runs = 1e4, n_jobs = 200, k = c(0.5, 0.5), s_stars = NULL, rhos = c(0,0), gammas = c(0,0.5))

set.seed(23524)
sim_func(scenarios = scenarios, N_runs = 1e4, n_jobs = 200, k = c(0.6, 0.4), s_stars = c(1, 0.5))
sim_func(scenarios = scenarios, N_runs = 1e4, n_jobs = 200, k = c(0.6, 0.4), s_stars = NULL, rhos = c(0,0), gammas = c(0,0.5))
sim_func(scenarios = scenarios_cv, N_runs = 1e4, n_jobs = 200, k = c(0.6, 0.4), s_stars = c(1, 0.5))
sim_func(scenarios = scenarios_cv, N_runs = 1e4, n_jobs = 200, k = c(0.6, 0.4), s_stars = NULL, rhos = c(0,0), gammas = c(0,0.5))
## -------------------------------------
## function to plot survival curves...
## -------------------------------------
plot_scenario <- function(event_model, max_cal_t, ylim = c(0,1), main = "", leg_pos = "topright"){
  
  get_t_pieces = function(t, ts){
    
    if (length(t) > 1) stop("t must be length 1")
    
    pmax(0, pmin(t - ts[-length(ts)], diff(ts)))
  }
  surv_pieces_simple = function(t, change_points, lambdas){
    
    if (length(t) > 1) stop("t must be length 1")
    if (length(change_points) != length(lambdas) - 1) stop("require one event rate per time period")
    ts = c(0, change_points, Inf)
    
    exp(-sum(get_t_pieces(t, ts) * lambdas))
    
  }
  t_seq <- seq(0, max_cal_t, length.out = 100)
  
  s_c <- purrr::map_dbl(t_seq, surv_pieces_simple, 
                        change_points = event_model$duration_c |> cumsum(), 
                        lambdas = c(event_model$lambda_c, event_model$lambda_c[length(event_model$lambda_c)]))
  
  s_e <- purrr::map_dbl(t_seq, surv_pieces_simple, 
                        change_points = event_model$duration_e |> cumsum(), 
                        lambdas = c(event_model$lambda_e, event_model$lambda_e[length(event_model$lambda_e)]))
  
  plot(t_seq, s_c, ylim = ylim, type = "l", 
       xlab = "Time (months)", ylab = "Survival",
       main = main)
  
  points(t_seq, s_e, lty = 2, type = "l")
  legend(leg_pos, c("Experimental", "Control"), lty = 2:1)

}


par(mfrow = c(2,3))
plot_scenario(event_model_1, 24, main = "Delayed effect, high event rate")
plot_scenario(event_model_2, 24, main = "Proportional hazards, high event rate")
plot_scenario(event_model_3, 24, main = "Diminishing effect, high event rate")
plot_scenario(event_model_4, 24, main = "Equal survival, high event rate")
plot_scenario(event_model_5, 24, main = "Early harm, high event rate")

par(mfrow = c(2,3))
plot_scenario(event_model_cv_1, 36, ylim = c(0.6, 1), main = "Delayed effect, low event rate", leg_pos = "bottomleft")
plot_scenario(event_model_cv_2, 36, ylim = c(0.6, 1), main = "Proportional hazards, low event rate", leg_pos = "bottomleft")
plot_scenario(event_model_cv_3, 36, ylim = c(0.6, 1), main = "Diminishing effect, low event rate", leg_pos = "bottomleft")
plot_scenario(event_model_cv_4, 36, ylim = c(0.6, 1), main = "Equal survival, low event rate", leg_pos = "bottomleft")
plot_scenario(event_model_cv_5, 36, ylim = c(0.6, 1), main = "Early harm, low event rate", leg_pos = "bottomleft")








