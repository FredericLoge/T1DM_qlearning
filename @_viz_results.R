#
my_baseline_policy <- function(cgm, tod, par = baseline_par_adult003){
  mg_per_dl_to_mmol_per_l <- 0.02586
  my_custom_scenario$meal$amount[tod] / par['cr'] +
    mg_per_dl_to_mmol_per_l * max(cgm - par['bg_target'], 0) / par['cf']
}

# load libraries
library(tidyverse)
library(ggplot2)
library(Rcpp)

# source functions
source('@_pql_V3_foo.R')

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

# read bunch of cvs files and compile results
list_csv_files <- list.files(path = './simu/', pattern = 'simulation_adult001_365days_basal127_*', full.names = TRUE)
## list_csv_files <- list.files(path = './simu/', pattern = 'simulation_adult002_365days_basal137_*', full.names = TRUE)
## list_csv_files <- list.files(path = './simu/', pattern = 'simulation_adult003_365days_basal000_*', full.names = TRUE)
simu <- lapply(X = list_csv_files, FUN = function(fn){
  simu_fn <- read.csv(file = fn)
  simm <- prepare_dataset(simu_adult001 = simu_fn, reward_aggregated = TRUE)
  simm <- simm[1:(nrow(simm)-4),]
  return(simm)
})
simu <- do.call(rbind.data.frame, simu)

# update to matrix to transfer easily afterwards
simu_mat <- as.matrix(simu[,c('time_of_day', 'cgm', 'insulin', 'bg_reward')])

# enum of ToD
time_of_day_enum <- sort(unique(simu$time_of_day))
nb_tod <- length(time_of_day_enum)

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

###
aggreg <- list(
  'alpha' = read.csv(file = './policy_found/adult001_gamma075_log_reward_alpha.csv'),
  'res' = read.csv(file = './policy_found/adult001_gamma075_log_reward_res.csv')
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
###
nonaggreg <- list(
  'alpha' = read.csv(file = './policy_found/adult001_gamma075_log_reward_alpha_nonaggreg.csv'),
  'res' = read.csv(file = './policy_found/adult001_gamma075_log_reward_res_nonaggreg.csv')
)
nonaggreg$res22 <- nonaggreg$res[nonaggreg$res$is_optimal,]
nonaggreg$res22$tod <- as.numeric(gsub(pattern = 'time_of_day_', replacement = '', nonaggreg$res22$time_of_day_label))
my_policy_notaggregated <- function(cgm, tod){
  i <- which.min( abs(nonaggreg$res22$cgm - cgm) +
                    1000 * (as.numeric(nonaggreg$res22$tod) != tod) )
  opt <- optim(par = nonaggreg$res22[i,'insulin'], fn = function(ins){
    basis_index <- spline_basis_2(x = c(cgm, ins), i = 1:nb_splines)
    - sum( nonaggreg$alpha[,tod] * basis_index)
  }, method = 'L-BFGS-B', lower = nonaggreg$res22[i,'insulin'] - 20 * pace, upper = nonaggreg$res22[i,'insulin'] + 20 * pace)
  return(opt$par)
}

#
library(ggplot2)
library(gridExtra)

##
##
##
df_plot_policy <- expand.grid(
  'bg_observed' = seq(from = x_lower_bound[1], to = x_upper_bound[1], length.out = 1000),
  'tod' = 1:4
)
df_plot_policy$baseline = NA
df_plot_policy$policy_aggreg = NA
df_plot_policy$policy_nonaggreg = NA
for(i in 1:nrow(df_plot_policy)){
  df_plot_policy$baseline[i] <- my_baseline_policy(cgm = df_plot_policy$bg_observed[i], tod = df_plot_policy$tod[i], par = baseline_par_adult001)
  df_plot_policy$policy_aggreg[i] <- my_policy_aggregated(cgm = df_plot_policy$bg_observed[i], tod = df_plot_policy$tod[i])
  df_plot_policy$policy_nonaggreg[i] <- my_policy_notaggregated(cgm = df_plot_policy$bg_observed[i], tod = df_plot_policy$tod[i])
}
temp <- gather(data = df_plot_policy, key = policy_index, value = bolus, -bg_observed, -tod)  
temp$tod_label <- ''
temp$tod_label[temp$tod == 1] <- '[1] - Breakfast (6 am), CHO = 50g'
temp$tod_label[temp$tod == 2] <- '[2] - Lunch (12 am), CHO = 60g'
temp$tod_label[temp$tod == 3] <- '[3] - Snack (3pm), CHO = 15g'
temp$tod_label[temp$tod == 4] <- '[4] - Dinner (8 pm), CHO = 80g'
temp$policy_label <- ''
temp$policy_label[temp$policy_index == 'baseline'] <- 'Baseline'
temp$policy_label[temp$policy_index == 'policy_aggreg'] <- 'Q-Learning, rewards between meals'
temp$policy_label[temp$policy_index == 'policy_nonaggreg'] <- 'Q-Learning, rewards at meals'
ggplot(data = temp[temp$policy_index %in% c('baseline', 'policy_aggreg'),]) +
  geom_point(inherit.aes = FALSE, data = simu, mapping = aes(x = cgm, y = insulin), col  ='black', alpha = 0.1) +
  geom_line(mapping = aes(x = bg_observed, y = bolus, col = tod_label), lwd = 1) +
  facet_grid(~ policy_label) +
  ggtitle('Optimized bolus advisors') + 
  xlab(label = 'Blood Glucose observed (mg/dL)') +
  ylab(label = 'Bolus Insulin prescribed (U)') +
  labs(colour = 'Time of Day')
ggplot(data = temp) +
  geom_point(inherit.aes = FALSE, data = simu, mapping = aes(x = cgm, y = insulin), col  ='black', alpha = 0.1) +
  geom_line(mapping = aes(x = bg_observed, y = bolus, col = policy_label), lwd = 1) +
  facet_grid(~ tod_label) +
  ggtitle('Optimized bolus advisors') + 
  xlab(label = 'Blood Glucose observed (mg/dL)') +
  ylab(label = 'Bolus Insulin prescribed (U)') +
  labs(colour = 'Time of Day')

g0_policy_baseline <- ggplot(data = df_plot_policy) +
  geom_line(mapping = aes(x = bg_observed, y = baseline, col = factor(tod)), lwd = 2, alpha = 0.5) +
  ggtitle('Baseline') + 
  xlab(label = 'Blood Glucose observed in mg/dL') +
  ylab(label = 'Insulin') +
  labs(colour = 'Time of Day')
g0_policy_found <- ggplot(data = df_plot_policy) +
  geom_line(mapping = aes(x = bg_observed, y = policy_found, col = factor(tod)), lwd = 2, alpha = 0.5) +
  ggtitle('Policy found') + 
  xlab(label = 'Blood Glucose observed in mg/dL') +
  ylab(label = 'Insulin') +
  labs(colour = 'Time of Day')
grid.arrange(g0_policy_baseline, g0_policy_found)

G0 <- ggplot(data = temp[temp$policy_index %in% c('baseline', 'policy_aggreg'),]) +
  geom_line(mapping = aes(x = bg_observed, y = bolus, col = tod_label), lwd = 1) +
  facet_grid(~ policy_label) +
  ggtitle('Optimized bolus advisors') + 
  xlab(label = 'Blood Glucose observed (mg/dL)') +
  ylab(label = 'Bolus Insulin prescribed (U)') +
  labs(colour = 'Time of Day') +
  theme(legend.position = 'bottom', legend.text = element_text(size = 10))
G0
temp$heuristic <- NA
for(i in 1:nrow(temp)){
  index_simu <- which(simu$time_of_day == temp$tod[i])
  distances <- (temp$bolus[i] - simu$insulin[index_simu])^2 / var(simu$insulin[index_simu]) +
    (temp$bg_observed[i] - simu$cgm[index_simu])^2 / var(simu$cgm[index_simu]) 
  temp$heuristic[i] <- mean(sort(distances, decreasing = TRUE)[1:50])
}
summary(temp$heuristic)
temp$heuristic_normalized <- 1 - (temp$heuristic - min(temp$heuristic)) / (max(temp$heuristic) - min(temp$heuristic))
G1 <- ggplot(data = temp[temp$policy_index == 'policy_aggreg',]) +
  geom_line(mapping = aes(x = bg_observed, y = bolus, col = tod_label, alpha = heuristic_normalized^2), lwd = 2) +
  xlab(label = 'Blood Glucose observed (mg/dL)') +
  ylab(label = 'Bolus Insulin prescribed (U)') +
  labs(colour = 'Time of Day', alpha = NULL) + 
  scale_alpha_continuous(guide = FALSE)
mysimu <- simu
mysimu$tod_label <- sort(unique(temp$tod_label))[mysimu$time_of_day]
G2 <- ggplot(data = temp[temp$policy_index == 'policy_aggreg',]) +
  geom_point(inherit.aes = FALSE, data = mysimu, mapping = aes(x = cgm, y = insulin), col  ='black', alpha = 0.1) +
  geom_line(mapping = aes(x = bg_observed, y = bolus, col = tod_label, alpha = heuristic_normalized^2), lwd = 2) +
  facet_wrap(~tod_label) +
  xlab(label = 'Blood Glucose observed (mg/dL)') +
  ylab(label = 'Bolus Insulin prescribed (U)') +
  labs(colour = 'Time of Day', alpha = NULL) + 
  scale_alpha_continuous(guide = FALSE)
G0
G1
G2
show_guide = FALSE
grid.arrange(G0, G1, ncol = 2)

###
###
compute_time_day <- function(CHO){
  some_index <- which(CHO > 0)[1]
  first_meal <- abs(my_custom_scenario$meal$amount - (CHO[some_index] * 3)) < 1e-03
  true_time <- rep(NA, length(CHO))
  true_time[some_index] <- my_custom_scenario$meal$time[[which(first_meal)]][[1]]
  if(some_index > 1){
    true_time[1:(some_index-1)] <- true_time[some_index] - (((some_index-1):1) *3)
  }
  true_time[(some_index+1):length(CHO)] <- true_time[some_index] + 3 * (1 : (length(CHO) - (some_index+1) + 1))
  return(true_time %% (24 * 60))
}

# read simulations
simu_baseline_raw <- read.csv(file = 'simu_noexplo/simulation_policy_baseline_cr_15_cf_2_53333_target_100_gamma_075_adult001_log_reward_basal127_1.csv')
simu_baseline_raw$time_day <- compute_time_day(CHO = simu_baseline_raw$CHO)
#plot(x = simu_baseline_raw$time_day, y = simu_baseline_raw$BG)
#plot.ts(simu_baseline_raw$BG[1:4000])

simu_policy_agg_raw <- read.csv(file = 'simu_noexplo/simulation_policy_suggested_gamma_075_adult001_log_reward_basal127_1.csv')
simu_policy_agg_raw$time_day <- compute_time_day(CHO = simu_policy_agg_raw$CHO)
# plot(x = simu_policy_agg_raw$time_day, y = simu_policy_agg_raw$BG)

simu_policy_nonagg_raw <- read.csv(file = 'simu_noexplo/simulation_policy_suggested_nonaggregated_gamma_075_adult001_log_reward_basal127_1.csv')
simu_policy_nonagg_raw$time_day <- compute_time_day(CHO = simu_policy_nonagg_raw$CHO)
# plot(x = simu_policy_nonagg_raw$time_day, y = simu_policy_nonagg_raw$BG)

# simu_my_policy_raw <- read.csv(file = 'simu_noexplo/simulation_noexplo_adult001_365days_basal127_001.csv')
# simu_my_policy_raw <- read.csv(file = 'simu_noexplo/simulation_my_policy_2_noexplo_adult001_365days_basal127_001.csv1.csv')
# simu_my_policy_raw <- read.csv(file = 'simu_noexplo/simulation_my_policy_3_gamma09_noexplo_adult001_365days_basal127_001.csv1.csv')
# simu_my_policy_raw <- read.csv(file = 'simu_noexplo/simulation_my_policy_4_gamma075_noexplo_adult001_365days_basal127_001.csv1.csv')
# simu_baseline_raw <- read.csv(file = 'simu_noexplo/simulation_baseline_noexplo_adult001_365days_basal127_001.csv')
# simu_baseline_raw <- read.csv(file = 'simu_noexplo/simulation_baseline_gridtuned_noexplo_adult001_365days_basal127_001.csv')
# simu_baseline_raw <- read.csv(file = 'simu_noexplo/simulation_baseline_gridtuned_reward2_noexplo_adult001_365days_basal127_1.csv')

simu_raw <- rbind.data.frame(
  cbind('policy_index' = 'baseline', simu_baseline_raw),
  cbind('policy_index' = 'policy_aggreg', simu_policy_agg_raw),  
  cbind('policy_index' = 'policy_nonaggreg', simu_policy_nonagg_raw)
)
str(simu_raw)
simu_raw$policy_label <- ''
simu_raw$policy_label[simu_raw$policy_index == 'baseline'] <- 'Baseline'
simu_raw$policy_label[simu_raw$policy_index == 'policy_aggreg'] <- 'Q-Learning, rewards between meals'
simu_raw$policy_label[simu_raw$policy_index == 'policy_nonaggreg'] <- 'Q-Learning, rewards at meals'
rowsel <- sample(x = which(simu_raw$policy_index %in% c("baseline", "policy_aggreg")), size = 5000, replace = FALSE)
rowsel_ok <- (simu_raw$time[rowsel] > (24 * 60 / 3))
rowsel <- rowsel[rowsel_ok]
ggplot(data = simu_raw[rowsel,], mapping = aes(x = time_day, y = BG)) +
  geom_hline(yintercept = c(70, 180), col = 'black', lty = 2, lwd = 1) +
  geom_hline(yintercept = c(40, 350), col = 'black', lty = 1, lwd = 1) +
  geom_hline(yintercept = 112.5, col = 'black', lty = 3, lwd = 1) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = 'loess', span = 0.1, col = 'blue', lwd = 1.5) +
  facet_grid(~ policy_label) +
  ggtitle('Blood glucose daily profiles') + 
  ylab(label = 'Average Blood Glucose observed (mg/dL)') +
  xlab(label = 'Time of day (mins)') +
  ylim(40, 350) 

G3 <- ggplot(data = simu_raw[rowsel,], mapping = aes(x = time_day, y = BG, col = factor(policy_label))) +
  geom_hline(yintercept = c(70, 180), col = 'black', lty = 2, lwd = 1) +
  geom_hline(yintercept = c(40, 350), col = 'black', lty = 1, lwd = 1) +
  geom_hline(yintercept = 112.5, col = 'black', lty = 3, lwd = 1) +
  geom_smooth(method = 'loess', span = 0.12, lwd = 2, se = TRUE) +
  ## geom_quantile(method = "rqss", quantiles = c(0.95, 0.75, 0.5, 0.25, 0.05), lambda = 100) + 
  ggtitle('Blood glucose daily profiles') + 
  labs(col = 'Policy') +
  ylab(label = 'Average Blood Glucose observed (mg/dL)') +
  xlab(label = 'Time of day (mins)') +
  ylim(40, 350) +
  theme(legend.position = 'bottom')
G3

# compare distribution
bg_reference_levels <- c(0, 40, 70, 112.5, 180, 350, 600)
library(arules)
tapply(X = simu_raw$BG, INDEX = simu_raw$policy_index, FUN = function(x){
  round(c(mean(x[x < 70]),   mean(x[x > 180])),2)
})
tapply(X = simu_raw$BG, INDEX = simu_raw$policy_index, FUN = function(x){
  xx <- discretize(x = x, method = 'fixed', breaks = bg_reference_levels)
  round(prop.table(table(xx)),2)
})
tapply(X = simu_raw$BG, INDEX = simu_raw$policy_index, FUN = function(x){
  round(summary(bg_reward(x)),2)
})



simu_raw$BG_reward <- bg_reward(simu_raw$BG)
summary(simu_raw$BG_reward)
vec <- tapply(X = simu_raw$BG_reward, INDEX = simu_raw$policy_label, FUN = mean)
vec <- data.frame('key' = names(vec), 'value' = as.numeric(vec))
ggplot(data = simu_raw, mapping = aes(x = factor(policy_label), y = BG_reward)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(inherit.aes = FALSE, data = vec, mapping = aes(x = key, y = value), col = 'red', pch = 'X', cex = 3) +
  ggtitle('Reward distribution') + 
  xlab(label = 'Policy followed') +
  ylab(label = 'Reward perceived')
ggplot(data = simu_raw, mapping = aes(fill = factor(policy_label), x = BG_reward)) +
  geom_density(alpha = 0.2) +
  ggtitle('Reward distribution') + 
  xlab(label = 'Policy followed') +
  ylab(label = 'Reward perceived')


# prepare datasets
simu_policy_nonagg <- prepare_dataset(simu_adult001 =  simu_policy_nonagg_raw)
simu_policy_agg <- prepare_dataset(simu_adult001 =  simu_policy_agg_raw)
simu_baseline <- prepare_dataset(simu_adult001 =  simu_baseline_raw)

### BLOOD GLUCOSE GRAPHICS
###
###

summary(bg_reward(simu_policy_nonagg_raw$BG))
summary(bg_reward(simu_policy_agg_raw$BG))
summary(bg_reward(simu_baseline_raw$BG))

summary(simu_policy_nonagg_raw$BG)
summary(simu_policy_agg_raw$BG)
summary(simu_baseline_raw$BG)

boxplot(simu_policy_nonagg_raw$BG, simu_policy_agg_raw$BG, simu_baseline_raw$BG)

### distribution of glucose (global)
#
g0_my_policy <- ggplot(data = simu_my_policy_raw, mapping = aes(x = BG)) +
  geom_histogram(col = 'white') +
  xlim(x_lower_bound[1], x_upper_bound[1]) +
  ggtitle('Policy found - BG histogram')
#
g0_baseline <- ggplot(data = simu_baseline_raw, mapping = aes(x = BG)) +
  geom_histogram(col = 'white') +
  xlim(x_lower_bound[1], x_upper_bound[1]) +
  ggtitle('Baseline - BG histogram')
#
grid.arrange(g0_baseline, g0_my_policy)

### much less hypo-glycemia -45% !
hypo_threshold <- 80
mean(simu_my_policy_raw$BG < hypo_threshold)
mean(simu_baseline_raw$BG < hypo_threshold)
mean(simu_my_policy_raw$BG < hypo_threshold) - mean(simu_baseline_raw$BG < hypo_threshold)

### traded for hyper-glycemia +44%
hyper_threshold <- 120
mean(simu_my_policy_raw$BG > hyper_threshold)
mean(simu_baseline_raw$BG > hyper_threshold)
mean(simu_my_policy_raw$BG > hyper_threshold) - mean(simu_baseline_raw$BG > hyper_threshold)

### distribution of glucose (globally)
#
some_index <- which(simu_my_policy_raw$CHO > 0)[1]
first_meal <- abs(my_custom_scenario$meal$amount - (simu_my_policy_raw$CHO[some_index] * 3)) < 1e-03
simu_my_policy_raw$true_time <- NA
simu_my_policy_raw$true_time[some_index] <- my_custom_scenario$meal$time[[which(first_meal)]][[1]]
simu_my_policy_raw$true_time[1:(some_index-1)] <- simu_my_policy_raw$true_time[some_index] - (((some_index-1):1) *3)
simu_my_policy_raw$true_time[(some_index+1):nrow(simu_my_policy_raw)] <- simu_my_policy_raw$true_time[some_index] + 3 * (1 : (nrow(simu_my_policy_raw) - (some_index+1) + 1))

simu_my_policy_raw$time_day <- simu_my_policy_raw$time %% (24 * 60 / 3)
sample_times_my_policy_raw <- sample.int(n = nrow(simu_my_policy_raw), size = nrow(simu_my_policy_raw)/10)
sample_times_my_policy_raw <- sample_times_my_policy_raw[sample_times_my_policy_raw > 24 * 60 / 3]
g1_my_policy_raw <- ggplot(data = simu_my_policy_raw[sample_times_my_policy_raw,], mapping = aes(x = time_day, y = BG)) +
  geom_hline(yintercept = c(hypo_threshold, hyper_threshold), col = 'red', lwd = 2) +
  geom_point(alpha = 0.2, lwd = 2) +
  geom_smooth(method = 'loess', span = 0.1) +
  geom_rug(data = simu_baseline_raw[simu_baseline_raw$CHO > 0,], mapping = aes(x = time_day, y = NULL), lwd = 2) +
  ylim(x_lower_bound[1], 250) +
  ggtitle('Policy found') + 
  xlab(label = 'Time of Day')
#
simu_baseline_raw$time_day <- (simu_baseline_raw$time + 100) %% (24 * 60 / 3)
sample_times_baseline_raw <- sample.int(n = nrow(simu_baseline_raw), size = nrow(simu_baseline_raw)/10)
sample_times_baseline_raw <- sample_times_baseline_raw[sample_times_baseline_raw > 24 * 60 / 3]
g1_baseline_raw <- ggplot(data = simu_baseline_raw[sample_times_baseline_raw,], mapping = aes(x = time_day, y = BG)) +
  geom_hline(yintercept = c(hypo_threshold, hyper_threshold), col = 'red', alpha = 0.5, lwd = 2) +
  geom_point(alpha = 0.2, lwd = 2) +
  geom_smooth(method = 'loess', span = 0.1) +
  geom_rug(data = simu_baseline_raw[simu_baseline_raw$CHO > 0,], mapping = aes(x = time_day, y = NULL), lwd = 2) +
  ylim(x_lower_bound[1], 250) +
  ggtitle('Baseline') + 
  xlab(label = 'Time of Day')
#
grid.arrange(g1_my_policy_raw, g1_baseline_raw)

### distribution of glucose (at measured times)
#
g1_my_policy <- ggplot(data = simu_my_policy, mapping = aes(x = bg, fill = time_of_day_label)) +
  geom_histogram() +
  xlim(x_lower_bound[1], x_upper_bound[1]) +
  ggtitle('Policy found - BG by Time of Day')
#
g1_baseline <- ggplot(data = simu_baseline, mapping = aes(x = bg, fill = time_of_day_label)) +
  geom_histogram() +
  xlim(x_lower_bound[1], x_upper_bound[1]) +
  ggtitle('Baseline - BG by Time of Day')
#
grid.arrange(g1_my_policy, g1_baseline)

### REWARD RESULTS
###
###

# global distribution
summary(bg_reward(simu_my_policy_raw$BG))
summary(bg_reward(simu_baseline_raw$BG))
summary(simu_my_policy$bg_reward)
summary(simu_baseline$bg_reward)

# representation as function of time of day
g2_my_policy <- ggplot(data = simu_my_policy, mapping = aes(x = bg_reward, fill = time_of_day_label)) +
  geom_histogram() +
  xlim(0, 1) +
  ggtitle('Policy found - Reward by Time of Day')
g2_my_baseline <- ggplot(data = simu_baseline, mapping = aes(x = bg_reward, fill = time_of_day_label)) +
  geom_histogram() +
  xlim(0, 1) +
  ggtitle('Baseline - Reward by Time of Day')
grid.arrange(g2_my_policy, g2_my_baseline)

### CLARKE GRID
###
###

# ega
library(ega)

### From one meal to another
###

#
simu_my_policy$next_bg <- c(simu_my_policy$bg[-1], NA)
g4_my_policy <- plotClarkeGrid(simu_my_policy$bg, simu_my_policy$next_bg, factor(simu_my_policy$time_of_day_label)) + 
  xlab('Glucose Concentration (mg/dL) at ToD t') +
  ylab('Glucose Concentration (mg/dL) at ToD t+1') + 
  ggtitle('Policy found - Clarke Error Grid - From meal to another')
#
simu_baseline$next_bg <- c(simu_baseline$bg[-1], NA)
g4_baseline <- plotClarkeGrid(simu_baseline$bg, simu_baseline$next_bg, factor(simu_baseline$time_of_day_label)) + 
  xlab('Glucose Concentration (mg/dL) at ToD t') +
  ylab('Glucose Concentration (mg/dL) at ToD t+1') + 
  ggtitle('Baseline - Clarke Error Grid - From meal to another')
#
grid.arrange(g4_my_policy, g4_baseline)

g4_my_policy <- ggplot(data = simu_my_policy, mapping = aes(x = bg, y = next_bg, col = factor(time_of_day_label))) + 
  geom_point() +
  xlab('Glucose Concentration (mg/dL) at ToD t') +
  ylab('Glucose Concentration (mg/dL) at ToD t+1') + 
  ggtitle('Policy found - From one meal to another') +
  labs(colour = 'Time of Day') +
  xlim(20, 150) +
  ylim(20, 150)
#
simu_baseline$next_bg <- c(simu_baseline$bg[-1], NA)
g4_baseline <- ggplot(data = simu_baseline, mapping = aes(x = bg, y = next_bg, col = factor(time_of_day_label))) + 
  geom_point() +
  xlab('Glucose Concentration (mg/dL) at ToD t') +
  ylab('Glucose Concentration (mg/dL) at ToD t+1') + 
  ggtitle('Baseline - From one meal to another') +
  labs(colour = 'Time of Day') + 
  xlim(20, 150) +
  ylim(20, 150) 
  
#
grid.arrange(g4_my_policy, g4_baseline)


### Globally
###

#
simu_my_policy_raw$NEXT_BG <- c(simu_my_policy_raw$BG[-1], NA)
g5_my_policy <- plotClarkeGrid(simu_my_policy_raw$BG[sample_times_my_policy_raw], simu_my_policy_raw$NEXT_BG[sample_times_my_policy_raw], "NULL") + 
  xlab('Glucose Concentration (mg/dL) at ToD t') +
  ylab('Glucose Concentration (mg/dL) at ToD t+1') + 
  ggtitle('Policy found - Clarke Error Grid - 3 mins intervals')
#
simu_baseline_raw$NEXT_BG <- c(simu_baseline_raw$BG[-1], NA)
g5_baseline <- plotClarkeGrid(simu_baseline_raw$BG[sample_times_baseline_raw], simu_baseline_raw$NEXT_BG[sample_times_baseline_raw], "NULL") + 
  xlab('Glucose Concentration (mg/dL) at ToD t') +
  ylab('Glucose Concentration (mg/dL) at ToD t+1') +  
  ggtitle('Baseline - Clarke Error Grid - 3 mins intervals')
#
grid.arrange(g5_my_policy, g5_baseline)

