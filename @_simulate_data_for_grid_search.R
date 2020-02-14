# bound to search
cr_bound <- c(3, 30)
cf_bound <- c(0.4, 2.8)
bg_target_bound <- c(100, 150)

# exploring grid
cr_cf_bg_target_explo <- expand.grid(
  'cr' = seq(from = cr_bound[1], to = cr_bound[2], length.out = 10),
  'cf' = seq(from = cf_bound[1], to = cf_bound[2], length.out = 10),
  'bg_target' = seq(from = bg_target_bound[1], to = bg_target_bound[2], length.out = 3)
)

# check out environment
Sys.getenv()
## Sys.setenv(RETICULATE_PYTHON = '/usr/local/bin/python3')

# load *reticulate* library
library(reticulate)

# check out python configuration found (should be 3.6)
py_config()

# use_python(python = '/usr/local/bin/python3')
# use_condaenv(condaenv = 'r-diabetes', conda = "/anaconda3/bin/conda")
# use_virtualenv('~/.virtualenvs/r-reticulate')
# py_install(packages = 'simglucose')

# check availability of *simglucose* and *gym* modules
py_module_available(module = "simglucose")
py_module_available(module = 'gym')

# load gym module and register individual Adult #001
gym <- import(module = 'gym')
gym$envs$registration$register(
  id='simglucose-adult3-v0',
  entry_point='simglucose.envs:T1DSimEnv',
  kwargs = list('patient_name' = 'adult#003')
)

# simulate for every combination considered
for(index_simu in 1:nrow(cr_cf_bg_target_explo)){
  
  print(index_simu)
  cat("\n")
  
  # my new policy (updated parameters)
  my_baseline_policy <- function(cgm, tod){
    mg_per_dl_to_mmol_per_l <- 0.02586
    my_custom_scenario$meal$amount[tod] / cr_cf_bg_target_explo$cr[index_simu] +
      mg_per_dl_to_mmol_per_l * max(cgm - cr_cf_bg_target_explo$bg_target[index_simu], 0) / cr_cf_bg_target_explo$cf[index_simu]
  }
  
  # create environment
  env = gym$make(id = 'simglucose-adult3-v0')
  
  # choose custom scenario (meals + amount of carbohydrates in grams)
  my_custom_scenario <- list('meal' = list(
    'time' = list(
      list(6 * 60),
      list(12 * 60),
      list(15 * 60),
      list(20 * 60)
    ),
    'amount' = c(50, 60, 15, 80)
  ))
  
  # specify basal prior, fixed througout the experiment
  ## test[test$category == 'Adult' & test$id == 1,]
  basal_prior_daily_amount <- 0 # 1.37
  basal_prior_3mins_amount <- basal_prior_daily_amount / (24 * 60 / 3)
  
  # choose length of experiment
  nb_days_simulated <- 8
  nb_observations <- 24 * 60 / 3 * nb_days_simulated
  
  # get first observation (CGM) by resetting environment
  ### env$seed(seed = as.integer(886579619))
  observation = env$reset()
  observation = observation$CGM
  
  # set progress bar
  pb <- txtProgressBar(min = 0, max = nb_observations)
  
  # indicator of meal (grams of CHO)
  meal_level <- 0
  
  # simulate !!!
  for(time_index in 1:nb_observations){
    
    # update progress bar
    setTxtProgressBar(pb = pb, value = time_index)
    
    # making sure no randomness in scenario chosen ...
    env$env$scenario$reset()
    env$env$scenario$scenario <- my_custom_scenario
    
    # basal action,
    action = basal_prior_3mins_amount 
    
    # to which we add bolus action, slightly randomized
    if(meal_level > 0){
      
      # pick random bolus and formulate action
      corre_tod <- which(abs(my_custom_scenario$meal$amount - meal_level*3) < 1e-03)
      bolus_policy <- my_baseline_policy(cgm = observation, tod = corre_tod)
      ## bolus_policy <- my_policy(cgm = observation, tod = corre_tod)
      action <- basal_prior_3mins_amount + bolus_policy
      
    }
    
    # evaluate new state
    new_state <- env$step(action)
    
    # check if meal is taken and new observation of blood glucose
    meal_level <- new_state$info$meal
    observation <- new_state$observation$CGM
    
  }
  
  # frame results simulated
  res2 <- data.frame(
    'time' = 1,
    'BG' = unlist(env$env$BG_hist),
    'CGM' = unlist(env$env$CGM_hist),
    'CHO' = c(unlist(env$env$CHO_hist), NA),
    'HBGI' = unlist(env$env$HBGI_hist),
    'LBGI' = unlist(env$env$LBGI_hist),
    'insulin' = c(unlist(env$env$insulin_hist), NA)
  )
  res2$time = 1:nrow(res2)
  
  # save them
  filename <- paste0('./simu_grid/simulation_', index_simu, '_adult003_365days_basal0.csv')
  ## file.exists(filename)
  if(file.exists(filename) == FALSE){
    write.csv(x = res2, file = filename)
  }
  
}

##
####
####

# handle the case BG = 0 !

list.files('simu_grid/')
#
cr_cf_bg_target_explo$res <- NA
for(index_simu in 1:nrow(cr_cf_bg_target_explo)){
  simu_ <- read.csv(file = paste0('./simu_grid/simulation_', index_simu, '_adult003_365days_basal0.csv'))
  simu__ <- prepare_dataset(simu_adult001 = simu_)
  n1 = nrow(simu__)
  n0 = floor(n1 / 2)
  cr_cf_bg_target_explo$res[index_simu] <- mean(simu__$bg_reward[1:n1])
}

#
library(ggplot2)
ggplot(data = cr_cf_bg_target_explo) +
  geom_tile(mapping = aes(x = cr, y = cf, fill = res)) +
  facet_wrap(facets = ~ paste0('BG target = ', bg_target)) + 
  scale_fill_viridis_c() +
  labs(fill = '||V||') +
  xlab(label = 'Carbohydrate Ratio (CR)') +
  ylab(label = 'Correction Factor (CF)')

#
head(cr_cf_bg_target_explo[order(cr_cf_bg_target_explo$res, decreasing = TRUE),])
