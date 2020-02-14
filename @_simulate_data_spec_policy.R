#
my_baseline_policy <- function(cgm, tod, par = baseline_par_adult003){
  mg_per_dl_to_mmol_per_l <- 0.02586
  my_custom_scenario$meal$amount[tod] / par['cr'] +
    mg_per_dl_to_mmol_per_l * max(cgm - par['bg_target'], 0) / par['cf']
}

baseline_par_adult001 <- c('cr' = 15, 'cf' = 2.533333, 'bg_target' = 100)
baseline_par_adult002 <- c('cr' = 12, 'cf' = 2.533333, 'bg_target' = 125)
baseline_par_adult003 <- c('cr' = 15, 'cf' = 2.800000, 'bg_target' = 150)

# 
# my_policy <- function(cgm, tod){
#   i <- which.min( abs(res22$cgm - cgm) +
#                1000 * (as.numeric(res22$tod) != tod) )
#   opt <- optim(par = res22[i,'init_par'], fn = function(ins){
#     basis_index <- spline_basis_2(x = c(cgm, ins), i = 1:nb_splines)
#     - sum( set_alpha_used[,tod] * basis_index)
#   }, method = 'L-BFGS-B', lower = 0, upper = 10)
#   return(opt$par)
# }

pace <- (x_upper_bound[2] - x_lower_bound[2]) / 50

nonaggreg <- list(
  'alpha' = read.csv(file = './policy_found/adult003_gamma075_log_reward_alpha_nonaggreg.csv'),
  'res' = read.csv(file = './policy_found/adult003_gamma075_log_reward_res_nonaggreg.csv')
)
nonaggreg$res22 <- nonaggreg$res[nonaggreg$res$is_optimal,]
nonaggreg$res22$tod <- as.numeric(gsub(pattern = 'time_of_day_', replacement = '', nonaggreg$res22$time_of_day_label))
my_policy_notaggregated <- function(cgm, tod){
  i <- which.min( abs(nonaggreg$res22$cgm - cgm) +
                    1000 * (as.numeric(nonaggreg$res22$tod) != tod) )
  opt <- optim(par = nonaggreg$res22[i,'insulin'], fn = function(ins){
    basis_index <- spline_basis_2(x = c(cgm, ins), i = 1:nb_splines)
    - sum( nonaggreg$alpha[,tod] * basis_index)
  }, method = 'L-BFGS-B', lower = nonaggreg$res22[i,'insulin'] - 2 * pace, upper = nonaggreg$res22[i,'insulin'] + 2 * pace)
  return(opt$par)
}
###
aggreg <- list(
  'alpha' = read.csv(file = './policy_found/adult003_gamma075_log_reward_alpha.csv'),
  'res' = read.csv(file = './policy_found/adult003_gamma075_log_reward_res.csv')
)
aggreg$res22 <- aggreg$res[aggreg$res$is_optimal,]
aggreg$res22$tod <- as.numeric(gsub(pattern = 'time_of_day_', replacement = '', aggreg$res22$time_of_day_label))
my_policy_aggregated <- function(cgm, tod){
  i <- which.min( abs(aggreg$res22$cgm - cgm) +
                    1000 * (as.numeric(aggreg$res22$tod) != tod) )
  opt <- optim(par = aggreg$res22[i,'insulin'], fn = function(ins){
    basis_index <- spline_basis_2(x = c(cgm, ins), i = 1:nb_splines)
    - sum( aggreg$alpha[,tod] * basis_index)
  }, method = 'L-BFGS-B', lower = aggreg$res22[i,'insulin'] - 2 * pace, upper = aggreg$res22[i,'insulin'] + 2 * pace)
  return(opt$par)
}
my_policy(cgm = 250, tod = 1)

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

# reading patient data
test <- read.table('./data_patients.txt')
colnames(test) <- c('category', 'pid', 'body_weight', 'age', 'optimal_basal', 'optimal_CF', 'optimal_CR')
test$id <- test$pid
while(TRUE){
  test_id_sup_10 <- (test$id > 10)
  if(sum(test_id_sup_10)==0) break
  test$id[test_id_sup_10] <- test$id[test_id_sup_10] - 10
}
test[test$category == 'Adult' & test$id == 2,]

# load gym module and register individual Adult #001
gym <- import(module = 'gym')
gym$envs$registration$register(
  id='simglucose-adult3-v0',
  entry_point='simglucose.envs:T1DSimEnv',
  kwargs = list('patient_name' = 'adult#003')
)

for(index_simu in 1){
  
  print(index_simu)
  cat("\n")
  
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
test[test$category == 'Adult' & test$id == 2,]
basal_prior_daily_amount <- 0
basal_prior_3mins_amount <- basal_prior_daily_amount / (24 * 60 / 3)

# choose length of experiment
nb_days_simulated <- ceiling(365/8)
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
    ## bolus_policy <- my_policy_aggregated(cgm = observation, tod = corre_tod)
    ## bolus_policy <- my_policy_notaggregated(cgm = observation, tod = corre_tod)
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
# filename <- paste0('./simu_noexplo/simulation_policy_suggested_gamma_075_adult003_log_reward_basal000_', index_simu, '.csv')
# filename <- paste0('./simu_noexplo/simulation_policy_suggested_nonaggregated_gamma_075_adult003_log_reward_basal000_', index_simu, '.csv')
filename <- paste0('./simu_noexplo/simulation_policy_baseline_cr_15_cf_2_8_target_150_gamma_075_adult003_log_reward_basal000_', index_simu, '.csv')
file.exists(filename)
if(file.exists(filename) == FALSE){
  write.csv(x = res2, file = filename)
}

}
