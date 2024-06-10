######################################################################
### get and plot data
template <- tempfile(fileext = ".xlsx")

httr::GET(url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-018-0134-3/MediaObjects/41591_2018_134_MOESM3_ESM.xlsx", 
          httr::write_disk(template))

dat <- readxl::read_excel(template,sheet = 2) %>% 
  select(PtID, ECOGGR, OS, OS.CNSR, TRT01P) %>%
  mutate(event = -1 * (OS.CNSR - 1),
         time = OS,
         arm = ifelse(TRT01P == "Docetaxel", "control", "experimental")) %>% 
  select(time, event, arm) %>%
  as.data.frame()

km <- survfit(Surv(time, event) ~ arm, data = dat)

survminer::ggsurvplot(km, 
                      data = dat, 
                      risk.table = TRUE, 
                      break.x.by = 6,
                      legend.title = "",
                      xlab = "Time (months)",
                      ylab = "Overall survival",
                      risk.table.fontsize = 4,
                      legend = c(0.8,0.8))


######################################################################
### perform analyses
nphRCT::wlrt(Surv(time, event) ~ arm,  data = dat, method = "mw", s_star = 1)$z |> pnorm()
nphRCT::wlrt(Surv(time, event) ~ arm,  data = dat, method = "mw", s_star = 0.5)$z |> pnorm()
nphRCT::wlrt(Surv(time, event) ~ arm,  data = dat, method = "fh", rho = 0, gamma = 0.5)$z |> pnorm()

### get correlations
find_cor(Surv(time, event) ~ arm,  data = dat, s_stars = c(1, 0.5))$cov_mat |> cov2cor() 
find_cor(Surv(time, event) ~ arm,  data = dat, s_stars = NULL, rhos = c(0, 0), gammas = c(0, 0.5))$cov_mat |> cov2cor() 


get_zcs <- function(formula, data, s_stars, rhos, gammas, k, alpha = 0.025, crit_interval = c(0.8, 1)){

  cov_res <- find_cor(formula = formula, data = data, s_stars = s_stars, rhos = rhos, gammas = gammas)
  cor_mat <- cov2cor(cov_res$cov_mat)
  z <- cov_res$u / sqrt(cov_res$v)
  df <- data.frame(t(z))
  colnames(df) <- paste0("z_", seq_along(z))
  df_2 <- data.frame(t(-find_crit(cor_mat, k = k, alpha = alpha, crit_interval = crit_interval)))
  colnames(df_2) <- paste0("c_", seq_along(z))
  cbind(df, df_2)
  
}
get_zcs(formula = Surv(time, event) ~ arm,  data = dat, alpha = 0.025, s_stars = c(1, 0.5), k = c(0.5, 0.5))
get_zcs(formula = Surv(time, event) ~ arm,  data = dat, alpha = 0.025, s_stars = NULL, rhos = c(0, 0), gammas = c(0, 0.5), k = c(0.5, 0.5))


## calculate p-value??
find_diff <- function(alpha, formula, data, s_stars, rhos, gammas, k){
  
  zcs <- get_zcs(formula = formula,  data = data, s_stars = s_stars,  rhos = rhos, gammas = gammas, k = k, alpha = alpha, crit_interval = c(-100, 100))
  
  min(zcs[stringr::str_detect(names(zcs), "z_")]) - min(zcs[stringr::str_detect(names(zcs), "c_")])
  
}
uniroot(find_diff, c(0.00001, 0.7), formula = Surv(time, event) ~ arm,  data = dat, s_stars = c(1, 0.5), k = c(0.5, 0.5))
uniroot(find_diff, c(0.00001, 0.7), formula = Surv(time, event) ~ arm,  data = dat, s_stars = NULL,  rhos = c(0, 0), gammas = c(0, 0.5), k = c(0.5, 0.5))
uniroot(find_diff, c(0.00001, 0.7), formula = Surv(time, event) ~ arm,  data = dat, s_stars = NULL,  rhos = c(0, 0), gammas = c(0, 0.5), k = c(0.6, 0.4))









