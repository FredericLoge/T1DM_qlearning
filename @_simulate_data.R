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
test[test$category == 'Adult' & test$id == 3,]

# load gym module and register individual Adult #001
gym <- import(module = 'gym')
gym$envs$registration$register(
  id='simglucose-adult3-v0',
  entry_point='simglucose.envs:T1DSimEnv',
  kwargs = list('patient_name' = 'adult#003')
)

for(index_simu in 17:30){
  
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
test[test$category == 'Adult' & test$id == 3,]
basal_prior_daily_amount <- 0 ## 0.89
basal_prior_3mins_amount <- basal_prior_daily_amount / (24 * 60 / 3)

# choose length of experiment
nb_days_simulated <- 365
nb_observations <- 24 * 60 / 3 * nb_days_simulated

# get first observation (CGM) by resetting environment
### env$seed(seed = as.integer(886579619))
observation = env$reset()
observation = observation$CGM

# set progress bar
pb <- txtProgressBar(min = 0, max = nb_observations)

# indicator of meal (grams of CHO)
meal_level <- 0

# compute bolus 
compute_bolus <- function(obs, bg_target = 150, cf, cr, meal_level){
  ## mgperdL_to_mmolperL <- 1/100 * 1/0.18
  1/100 * pmax(observation - bg_target, 0) / cf + meal_level / cr 
}

#
pick_bolus <- function(n = 1, upper_level, lower_level){
  p_range <- 0.30
  p_belong <- 0.75
  m <- (upper_level + lower_level) / 2
  s <- ((upper_level - lower_level) * p_range) / qnorm(p = 1 - (1 - p_belong)/2)
  sam <- rnorm(n = n, mean = m, sd = s)
  sam <- pmin(pmax(sam, lower_level), upper_level)
  return(sam)
}
# test <- pick_bolus(n = 1000, upper_level = 100, lower_level = -100)
# hist(test)
# mean(abs(test) < 200*0.25)

count_too_low = 0
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
    
    # bounds of parameters
    cr_bound <- c(3, 30)
    cf_bound <- c(0.4, 2.8)
    bg_target_bound <- c(100, 150)

    # upper and lower bound of insulin level to prescribe
    upper_bound <- compute_bolus(obs = observation, bg_target = bg_target_bound[1], cf = cf_bound[1], cr = cr_bound[1], meal_level = meal_level)
    upper_bound <- upper_bound / 2
    lower_bound <- compute_bolus(obs = observation, bg_target = bg_target_bound[2], cf = cf_bound[2], cr = cr_bound[2], meal_level = meal_level)
    
    # pick random bolus and formulate action
    bolus_rand <- pick_bolus(n = 1, upper_level = upper_bound, lower_level = lower_bound)
    action <- basal_prior_3mins_amount + bolus_rand
    
  }
  
  # evaluate new state
  new_state <- env$step(action)

  # check if meal is taken and new observation of blood glucose
  meal_level <- new_state$info$meal
  observation <- new_state$observation$CGM
  if(abs(observation - 39) < 1e-10){
    count_too_low = count_too_low + 1
  }else{
    count_too_low = 0
  }
  if(count_too_low > 24 * 60 / 3) break 
}

# plot.ts(unlist(env$env$BG_hist))
# lines(unlist(env$env$insulin_hist) * 50, col = 'red')

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
filename <- paste0('./simu/simulation_adult003_365days_basal000_00', index_simu, '.csv')
file.exists(filename)
if(file.exists(filename) == FALSE){
  write.csv(x = res2, file = filename)
}

}

# EXTRA STUFF TO CLEAN UP
#
#
#

res2$prevCHO <- factor(c(NA, res2$CHO[-nrow(res2)]))
ggplot(data = res2) +
  geom_point(mapping = aes(x = CGM, y = insulin)) +
  facet_grid(~prevCHO)
plot.ts(res2$CGM, res2$insulin)

res2$reward = 1*(res2$BG >= 70 & res2$BG < 150)
plot.ts(res2)
library(ggplot2)
ggplot(data = res2, aes(x = time, y = reward)) + 
  geom_point() + geom_smooth(method = 'loess')
plot.ts(res2$BG, ylim = c(0, 300))
lines(res2$insulin*10, col = 'red')

###
meal_time_indexes <- which(res2$CHO > 0)
total_indexes <- c(1, meal_time_indexes, nrow(res2))
res3 <- array(data = NA, dim = c(length(meal_time_indexes), 6))
for(index in 2:length(meal_time_indexes)){
  preceding_meal_time_index = meal_time_indexes[index-1]
  meal_time_index = meal_time_indexes[index]
  new_data <- c(
      res2$time[meal_time_index],
      res2$CGM[meal_time_index],
      res2$CHO[meal_time_index],
      mean(res2$HBGI[preceding_meal_time_index:meal_time_index] + res2$LBGI[preceding_meal_time_index:meal_time_index]),
      mean(res2$reward[preceding_meal_time_index:meal_time_index]),
      res2$insulin[meal_time_index+1]
    )
    res3[index,] <- new_data
}
colnames(res3) <- c('time', 'noisy_bg', 'CHO', 'cost', 'reward', 'insulin')
plot.ts(res3)

res3 <- data.frame(res3)
res3 <- res3[-1,]
res3$time_of_day_index = rep(1:4, 20)[1:39]
ggplot(data = res3, mapping = aes(x = noisy_bg, y = reward, color = insulin)) +
  geom_point() +
  facet_grid(~factor(time_of_day_index))
noisy_bg <- function(x){
  x + 0
}

plot(res3$insulin + c(0, res3$insulin[-27]), res3$reward)
# colnames(res) <- c('bg', 'insulin', 'reward', 'meal')
# plot.ts(res)
# 
# #
# plot.ts(res[,'bg'])
# abline(v = which(res[,"meal"]>0), col = 'red')

# percentage of time spent in the different categories of BG levels
breaks_bg <- c(0, 50, 70, 150, 300, Inf)
hist(res2$BG)
abline(v = breaks_bg, col = 'red')
library(arules)
table(discretize(x = res2$BG, method = 'fixed', breaks = breaks_bg))/nrow(res2)

# the meals seem to be taken very much randomly ... although they are specified in the scenario ... is it drawn at random at every time step ?
env$env$scenario$scenario
env$env$BG_hist

plot(x = res2$CGM, res2$insulin)

diff(unlist(env$env$scenario$scenario$meal$time))/3

write.csv(x = res2, file = "simulation_adult_001_random_policy.csv")
