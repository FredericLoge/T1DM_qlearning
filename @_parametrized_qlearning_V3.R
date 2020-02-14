#### LOAD LIBRARIES, SOURCE FUNCTIONS
####
####

# load libraries
library(tidyverse)
library(ggplot2)
library(Rcpp)

# source functions
source('@_pql_V3_foo.R')

#### RECALL THE CUSTOM SCENARIO
####
####

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

#### READ SIMULATIONS SAVED IN .CSV FILES
####
####

# read bunch of cvs files and compile results
list_csv_files <- list.files(path = './simu/', pattern = 'simulation_adult003_365days_basal000_*', full.names = TRUE)
##list_csv_files <- list.files(path = './simu/', pattern = 'simulation_adult001_365days_basal127_*', full.names = TRUE)
## list_csv_files <- list.files(path = './simu/', pattern = 'simulation_adult002_365days_basal137_*', full.names = TRUE)
simu <- lapply(X = list_csv_files, FUN = function(fn){
  simu_fn <- read.csv(file = fn)
  simm <- prepare_dataset(simu_adult001 = simu_fn, reward_aggregated = FALSE)
  simm <- simm[1:(nrow(simm)-4),]
  return(simm)
})
simu <- do.call(rbind.data.frame, simu)

#
summary(simu[,'cgm'])
summary(simu[,'bg'])

# change CGM by BG + Noise
## simu$cgm <- simu$bg + noise

# update to matrix to transfer easily afterwards
simu_mat <- as.matrix(simu[,c('time_of_day', 'cgm', 'insulin', 'bg_reward')])

#### SOME DESCRIPTIVE STATISTICS
####
####

# observe correlation between BG and CGM measures
cor(simu$bg, simu$cgm)

# enum of ToD
time_of_day_enum <- sort(unique(simu$time_of_day))
nb_tod <- length(time_of_day_enum)

# represent reward plot
ggplot(data = simu) +
  geom_point(mapping = aes(x = bg, y = bg_reward, col = 'Between-meals reward'), alpha = 0.5) +
  geom_point(mapping = aes(x = bg, y = bg_reward(bg), col = 'Single-measure reward')) +
  xlab('Blood glucose level in mg/L') +
  ylab('Reward associated') +
  facet_wrap(~paste0('Time of Day : ', time_of_day)) +
  labs(colour = 'Reward')
ggplot(data = simu) +
  geom_point(mapping = aes(y = bg_reward(bg), x = bg_reward, col = 'Between-meals reward'), alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~time_of_day) 

#### SET BASIS FUNCTIONS
####
####

# sqrt(nb knots) 
sqrt_nb_knots <- 8

# support(x) bound
x_lower_bound <- c(min(simu$cgm), min(simu$insulin))
x_upper_bound <- c(max(simu$cgm), max(simu$insulin))

# grid for basis
ix <- expand.grid(
  'cgm' = seq(from = x_lower_bound[1], to = x_upper_bound[1], length.out = sqrt_nb_knots),
  'insulin' = seq(from = x_lower_bound[2], to = x_upper_bound[2], length.out = sqrt_nb_knots)
)
nb_splines <- nrow(ix) # sqrt_nb_knots^2

# "appropriate" bandwidth
prob <- 0.2
h_cgm <- (x_upper_bound[1] - x_lower_bound[1]) / sqrt_nb_knots * (- log(2 * prob))^(-1/2)
h_insulin <- (x_upper_bound[2] - x_lower_bound[2]) / sqrt_nb_knots * (- log(2 * prob))^(-1/2)
# h <- 0.4
# h_cgm <- h * sd(simu$cgm)
# h_insulin <- h * sd(simu$insulin)

# highest level attained by a basis function
niveau_max <- 1 / (2 * pi * h_cgm * h_insulin) 

# basis function, takes x a two-column matrix as input, i is scalar
# this is just for representation purposes 
spline_basis_1 <- function(x, i){
  dnorm(x = x[1,], mean = ix$cgm[i], sd = h_cgm) *
    dnorm(x = x[2,], mean = ix$insulin[i], sd = h_insulin)
}
spline_basis_2 <- function(x, i){
  dnorm(x = x[1], mean = ix$cgm[i], sd = h_cgm) *
    dnorm(x = x[2], mean = ix$insulin[i], sd = h_insulin)
}

### VIZ IS EXPECTED ...

#### OPTIM PARAMETERS
####
####

# learning rate
learning_rate <- 0.1

# number of iterations of algorithm
nb_turns <- 5e05 
nb_global_turns <- 10

# gamma parameter of state(-action) value functions
gamma <- 0.75

# representing cumulative effect of discount
time_index <- 1:10
value <- cumsum(gamma^(time_index-1))
plot(x = c(0, time_index), y = c(0, value), ylim = c(0, 1 / (1-gamma)), type = 'o', pch = 20)
abline(h = 1 / (1 - gamma), col = 'red')
abline(h = 0.8 * 1 / (1 - gamma), col = 'red', lty = 2)

# initialize set_alpha and memory 
set_alpha <- array(data = 0, dim = c(nb_splines, nb_tod))
set_alpha_mem <- list()

# record start time
start_time <- Sys.time()

# consider discrete grid for insulin level tested in gradient update (max_{a'} Q(s',a'))
insulin_test <- seq(from = x_lower_bound[2], to = x_upper_bound[2], length.out = 100)
insulin_test[2] - insulin_test[1]

# global turns of gradient descent (Rcpp implementation does not support 1e06 observations, so we are splitting the work)
for(index in 1:nb_global_turns){
  
  # sample row indexes to analyze
  sample_indexes <- sample.int(n = nrow(simu)-1, size = nb_turns, replace = TRUE)
  
  # perform gradient descent algorithm (max supported: 500 000)
  set_alpha_mem_temp <- gradient_descent(nb_saves = 100, 
                                         set_alpha = set_alpha, 
                                         sample_index = sample_indexes, 
                                         length_sample_index = length(sample_indexes),
                                         simu_tod = as.numeric(simu_mat[,'time_of_day']),
                                         simu_cgm = simu_mat[,'cgm'], 
                                         simu_insulin = simu_mat[, 'insulin'], 
                                         simu_reward = simu_mat[, 'bg_reward'],
                                         cgmRef = ix$cgm, 
                                         insulinRef = ix$insulin, 
                                         h_cgm = h_cgm, h_insulin = h_insulin, gamma = gamma,
                                         insulin_test = insulin_test, 
                                         final_learning_rate = learning_rate / niveau_max,
                                         nb_splines = nb_splines, 
                                         length_insulin_test = length(insulin_test), 
                                         previous_nb_iterations = (index-1)*nb_turns)
  
  # save !
  set_alpha_mem[[index]] <- set_alpha_mem_temp
  set_alpha <- set_alpha_mem_temp[[length(set_alpha_mem_temp)]]
  if(max(abs(set_alpha)) > 1e5) break
}

# record end time and compare with start time
end_time <- Sys.time()
difftime(time2 = start_time, time1 = end_time)

# 2.23 minutes pour 1 millions d'iterations ...

#### CHECK CONVERGENCE OF ESTIMATE OF (\alpha) SET
####
####

# prep set_alpha_mem list
set_alpha_mem_prep <- NULL
for(i in 1:length(set_alpha_mem)){
  test <- do.call(rbind.data.frame, set_alpha_mem[[i]])
  set_alpha_mem_prep <- rbind(set_alpha_mem_prep, test)
}
colnames(set_alpha_mem_prep) <- paste0('time_of_day_', 1:4)
set_alpha_mem_prep$spline_index <- factor(paste0('alpha_', 1:nb_splines))
set_alpha_mem_prep$time_index <- rep(1:(nrow(set_alpha_mem_prep)/nb_splines), each = nb_splines)
set_alpha_mem_prep <- set_alpha_mem_prep %>% 
    gather(key = param_key, value = param_value, -time_index, -spline_index, factor_key = TRUE)
str(set_alpha_mem_prep)
summary(set_alpha_mem_prep[[length(set_alpha_mem_prep)]])

# ggplot of spline basis
nb_modulo <- (nrow(set_alpha_mem_prep)/nb_splines) / 100
##ggplot(data = set_alpha_mem_prep) + 
ggplot(data = set_alpha_mem_prep %>% filter(time_index %% nb_modulo == 1)) +
  geom_line(mapping = aes(x = time_index, y = param_value, color = spline_index)) +
  facet_grid(~param_key) + 
  theme(legend.position = 'none')

#### REPRESENT "FINAL" Q(s,a) ESTIMATE
####
####

# reach for selection
set_alpha_used <- set_alpha_mem[[1]][[length(set_alpha_mem[[1]])]]
summary(set_alpha_used)
set_alpha_used <- set_alpha_mem[[nb_global_turns]][[length(set_alpha_mem[[nb_global_turns]])]]

# grid for basis
nb <- 50
ix2 <- expand.grid(
  'cgm' = seq(from = x_lower_bound[1], to = x_upper_bound[1], length.out = nb),
  'insulin' = seq(from = x_lower_bound[2], to = x_upper_bound[2], length.out = nb)
)

# compute basis
res <- array(NA, dim = c(nrow(ix2), nb_tod))
for(ix2_row_index in 1:nrow(ix2)){
  basis_index <- spline_basis_2(x = as.numeric(ix2[ix2_row_index, c('cgm', 'insulin')]), i = 1:nb_splines)
  for(tod_index in 1:nb_tod){
    res[ix2_row_index, tod_index] <- sum( set_alpha_used[,tod_index] * basis_index)
  }
}
colnames(res) <- paste0('time_of_day_', 1:nb_tod)
res <- data.frame(res)

# find optimal policy
res2 <- cbind.data.frame(ix2, res)
res2 <- res2 %>% gather(key = time_of_day_label, value = qsa_estimate, -cgm, -insulin, factor_key = TRUE) 
res2$is_optimal <- FALSE
for(tod in 1:nb_tod){
  for(cgm_level in unique(res2$cgm)){
    cond <- (res2$time_of_day_label == paste0('time_of_day_', tod)) &
      (res2$cgm == cgm_level)
    index_max <- which.max(res2$qsa_estimate[cond])
    res2$is_optimal[cond][index_max] <- TRUE
  }
}

#
res22 <- res2[res2$is_optimal,c('cgm', 'time_of_day_label', 'insulin')]
colnames(res22)[3] <- 'init_par'
res22$tod <- as.numeric(gsub(pattern = 'time_of_day_', replacement = '', res22$time_of_day_label))
res22$opt_par <- NA
pace <- (x_upper_bound[2] - x_lower_bound[2]) / 50
for(i in 1:nrow(res22)){
  opt <- optim(par = res22[i,'init_par'], fn = function(ins){
    basis_index <- spline_basis_2(x = c(res22[i,'cgm'], ins), i = 1:nb_splines)
    - sum( set_alpha_used[,res22[i,'tod']] * basis_index)
  }, method = 'L-BFGS-B', lower = res22[i,'init_par'] - 2 * pace, upper = res22[i,'init_par'] + 2 * pace)
  res22$opt_par[i] <- opt$par
}

# represent policy derived
ggplot(data = res2, mapping = aes(x = cgm, y = insulin, fill = qsa_estimate)) +
  geom_tile() + 
  geom_point(inherit.aes = FALSE, data = ix, mapping = aes(x = cgm, y = insulin), col = 'red', alpha = 0.3) +
  geom_point(inherit.aes = FALSE, data = simu, aes(x = cgm, y = insulin), alpha = 0.1) +
  geom_line(data = res2[res2$is_optimal,], col = 'white') +
  geom_line(inherit.aes = FALSE, data = res22, mapping = aes(x = cgm, y = opt_par), col = 'black') +
  facet_wrap(~time_of_day_label) +
  scale_fill_viridis_c()

summary(res2$qsa_estimate[res2$is_optimal])
1 / (1 - gamma)

res2_my_policy <- res2
fn <- './policy_found/adult003_gamma075_log_reward_res_nonaggreg.csv'
if(file.exists(fn) == FALSE){
  write.csv(x = res2_my_policy, file = fn)
}
fn <- './policy_found/adult003_gamma075_log_reward_alpha_nonaggreg.csv'
if(file.exists(fn) == FALSE){
  write.csv(x = set_alpha_used, file = fn, row.names = FALSE)
}

res2 <- read.csv(file = './policy_found/adult002_gamma075.csv')

my_policy <- function(cgm, tod){
  res2_tod_cond <- (res2_my_policy$time_of_day_label == paste0('time_of_day_', tod))
  res2_subset <- res2_my_policy[res2_tod_cond & res2_my_policy$is_optimal,]
  cond <- which.min(abs(res2_subset$cgm - cgm))
  res2_subset$insulin[cond]
}
my_policy(cgm = 0, tod = 1)

#### DERIVE SMOOTHED POLICY FUNCTION
####
####

# policy1
res3 <- res2 %>% filter(time_of_day_label == 'time_of_day_1')
ggplot(data = res3, mapping = aes(x = cgm, y = insulin, fill = qsa_estimate)) +
  geom_tile() + 
  geom_point(inherit.aes = FALSE, data = ix, mapping = aes(x = cgm, y = insulin), col = 'red', alpha = 0.3) +
  geom_point(inherit.aes = FALSE, data = simu, aes(x = cgm, y = insulin), alpha = 0.1) +
  geom_line(data = res3[res3$is_optimal,], col = 'white') +
  geom_smooth(data = res3[res3$is_optimal & (res3$cgm < 80 | (res3$cgm > 120 & res3$cgm < 220)),], method = 'gam') +
  scale_fill_viridis_c()
lm(formula = insulin ~ cgm, data = res3[res3$is_optimal & (res3$cgm < 80 | (res3$cgm > 120 & res3$cgm < 220)),])
policy_tod_1 <- function(cgm){
  0.022 * max(100, cgm) + 2.371
}

# policy2
res3 <- res2 %>% filter(time_of_day_label == 'time_of_day_2')
ggplot(data = res3, mapping = aes(x = cgm, y = insulin, fill = qsa_estimate)) +
  geom_tile() + 
  geom_point(inherit.aes = FALSE, data = ix, mapping = aes(x = cgm, y = insulin), col = 'red', alpha = 0.3) +
  geom_point(inherit.aes = FALSE, data = simu, aes(x = cgm, y = insulin), alpha = 0.1) +
  geom_line(data = res3[res3$is_optimal,], col = 'white') +
  geom_smooth(data = res3[res3$is_optimal & (res3$cgm < 220),], method = 'gam') +
  scale_fill_viridis_c()
lm(formula = insulin ~ cgm, data = res3[res3$is_optimal & res3$cgm < 220,])
policy_tod_2 <- function(cgm){
  0.02487 * min(240, cgm) + 3.49282
}

# policy3
res3 <- res2 %>% filter(time_of_day_label == 'time_of_day_3')
ggplot(data = res3, mapping = aes(x = cgm, y = insulin, fill = qsa_estimate)) +
  geom_tile() + 
  geom_point(inherit.aes = FALSE, data = ix, mapping = aes(x = cgm, y = insulin), col = 'red', alpha = 0.3) +
  geom_point(inherit.aes = FALSE, data = simu, aes(x = cgm, y = insulin), alpha = 0.1) +
  geom_line(data = res3[res3$is_optimal,], col = 'white') +
  geom_smooth(data = res3[res3$is_optimal & (res3$cgm < 150),], method = 'gam') +
  scale_fill_viridis_c()
lm(formula = insulin ~ cgm, data = res3[res3$is_optimal & res3$cgm < 150,])
lm(formula = insulin ~ cgm, data = res3[res3$is_optimal & res3$cgm > 270,])
policy_tod_3 <- function(cgm){
  thresh <- 265
  (cgm < thresh) * 1.06 + (cgm >= thresh) * 7.57
}

# policy4
res3 <- res2 %>% filter(time_of_day_label == 'time_of_day_4')
ggplot(data = res3, mapping = aes(x = cgm, y = insulin, fill = qsa_estimate)) +
  geom_tile() + 
  geom_point(inherit.aes = FALSE, data = ix, mapping = aes(x = cgm, y = insulin), col = 'red', alpha = 0.3) +
  geom_point(inherit.aes = FALSE, data = simu, aes(x = cgm, y = insulin), alpha = 0.1) +
  geom_line(data = res3[res3$is_optimal,], col = 'white') +
  geom_smooth(data = res3[res3$is_optimal & (res3$cgm < 150),], method = 'gam') +
  scale_fill_viridis_c()
lm(formula = insulin ~ cgm, data = res3[res3$is_optimal & res3$cgm < 180,])
lm(formula = insulin ~ cgm, data = res3[res3$is_optimal & res3$cgm > 220,])
policy_tod_4 <- function(cgm){
  thresh <- 200
  (cgm < thresh) * 2.3 + (cgm >= thresh) * 6
}

# establish policy
my_policy <- function(cgm, tod){
  if(tod == 1){
    policy_tod_1(cgm)
  }else if(tod == 2){
    policy_tod_2(cgm)
  }else if(tod == 3){
    policy_tod_3(cgm)
  }else{
    policy_tod_4(cgm) 
  }
}

###
### <!> nous avons considéré avoir accès au blood glucose sans
### erreur de mesure, cf la construction de l'objet "simu" plus haut
### dans le script.
###

## read patient data
test <- read.table('./data_patients.txt')
colnames(test) <- c('category', 'pid', 'body_weight', 'age', 'optimal_basal', 'optimal_CF', 'optimal_CR')
test$id <- test$pid
while(TRUE){
  test_id_sup_10 <- (test$id > 10)
  if(sum(test_id_sup_10)==0) break
  test$id[test_id_sup_10] <- test$id[test_id_sup_10] - 10
}
test[test$category == 'Adult' & test$id == 1,]

# linear policy, as baseline (to work on !)
my_baseline_policy <- function(cgm, tod){
  mg_per_dl_to_mmol_per_l <- 0.02586
  my_custom_scenario$meal$amount[tod] / 15 + mg_per_dl_to_mmol_per_l * max(cgm - 125, 0) / 2
}
