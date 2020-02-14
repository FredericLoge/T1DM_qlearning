### reward functions
# log_transform_bg <- function(bg){
#     gPerDL_to_mM <- 1/18
#     alpha <- 1.026 # 1.084
#     beta <- 1.861 # 5.381
#     gamma <- 1.794 # 1.509
#     gamma * (log(bg * gPerDL_to_mM)^alpha - beta)
# }
# bg_reward <- function(bg){
#   bg <- pmax(20, pmin(600, bg))
#   f_bg <- log_transform_bg(bg = bg / 100)
#   return(1 - f_bg^2 / 10)
# }

# calibrated in @_own_reward_function.R# BG reference levels, in mg/dl
# bg_reference_levels <- c('low' = 40, 'hypo' = 70, 'center' = 112.5, 'hyper' = 180, 'high' = 350)
log_transform_bg <- function(bg){
  gamma = 4.0955978
  beta = 2.1495359
  alpha = 0.4937448
  gamma * (log(bg)^alpha - beta)
}
bg_reward <- function(bg){
  bg <- pmax(40, pmin(350, bg))
  f_bg <- log_transform_bg(bg = bg)
  return(1 - f_bg^2)
}
## round(bg_reward(bg = bg_reference_levels), 2)

# par_bg_reward <- list("a" = 80, "b" = 120, "c" = 30, "d" = 60)
# bg_reward <- function(bg, par = par_bg_reward){
#   rew <- 0 * bg + 1
#   rew[bg < par$a] <- exp( - (bg[bg < par$a] - par$a)^2 / par$c^2 )
#   rew[bg > par$b] <- exp( - (bg[bg > par$b] - par$b)^2 / par$d^2 )
#   return(rew)
# }

if(FALSE){
  reward_equiv <- data.frame(
    'reward' = seq(from = 1, to = 0, length.out = 100),
    'hypo' = NA,
    'hyper' = NA
  )
  reward_equiv$hypo[1] <- par_bg_reward$a
  reward_equiv$hyper[1] <- par_bg_reward$b
  for(i in 2:nrow(reward_equiv)){
    reward_equiv$hypo[i] <- par_bg_reward$a - sqrt(log(1/reward_equiv$reward[i])*par_bg_reward$c^2)
    reward_equiv$hyper[i] <- par_bg_reward$b + sqrt(log(1/reward_equiv$reward[i])*par_bg_reward$d^2)
  }
  g <- ggplot(data = reward_equiv) +
    geom_line(mapping = aes(x = reward, y = hypo, col = 'hypo value')) +
    geom_line(mapping = aes(x = reward, y = hyper, col = 'hyper value')) 
  library(plotly)
  ggplotly(g)
  bg <- seq(from = 0, to = 4, length.out = 100)
  plot.ts(log_transform_bg(bg = bg))
  plot(x = bg, y = bg_reward(bg = bg * 100), type = 'l')
  abline(h = 0.5, lty = 2)
  80 - sqrt(log(1/0.5)*30^2)
  plot(x = bg, y = bg_reward_0(bg = bg * 100), type = 'l')
  lines(x = bg, y = bg_reward(bg = bg * 100, par = list("a" = 80, "b" = 120, "c" = 70, "d" = 300)), type = 'l', col = 'red')
}

### transform to an appropriate data format
prepare_dataset <- function(simu_adult001, reward_aggregated = TRUE){
  
  # identify meal times
  meal_time_indexes <- which(simu_adult001$CHO > 0)
  nb_meal_time <- length(meal_time_indexes)
  
  # prepare container
  simu <- array(data = NA, dim = c(nb_meal_time-1, 7))
  
  # iterate over meal indexes
  for(index in 1:(nb_meal_time-1)){
    ## preceding_meal_time_index = meal_time_indexes[index-1]
    meal_time_index = meal_time_indexes[index]
    next_meal_time_index = meal_time_indexes[index+1]
    between_meal_indexes <- (meal_time_index+2):(next_meal_time_index+1)
    if(reward_aggregated){
      reward_value <- with(simu_adult001, mean(bg_reward(BG[between_meal_indexes])))
    }else{
      reward_value <- with(simu_adult001, bg_reward(BG[next_meal_time_index+1]))
    }
    new_data <- with(simu_adult001, 
                     c( time[meal_time_index], 
                        BG[meal_time_index], 
                        CGM[meal_time_index], 
                        CHO[meal_time_index],
                        mean(HBGI[between_meal_indexes] + LBGI[between_meal_indexes]),
                        reward_value,
                        insulin[meal_time_index+1]
                     )
    )
    simu[index,] <- new_data
  }
  
  ## add column names and switch to dataframe
  colnames(simu) <- c('time', 'bg', 'cgm', 'cho', 'bgi', 'bg_reward', 'insulin')
  simu <- data.frame(simu)
  
  ## temporal difference between meal 1 and 2, 2 and 3, etc
  temp <- unlist(my_custom_scenario$meal$time)
  temp <- c(temp, temp[1] + 24 * 60)
  
  ## deduce initial time index
  initial_index <- which(diff( simu$time[1:10] )[1] == diff(temp/3))
  
  ## time of day
  temp <- (1:4 + initial_index - 1) 
  temp[temp > 4] <- temp[temp > 4] - 4
  simu$time_of_day <- rep(x = temp, times = ceiling(nrow(simu)/4))[1:nrow(simu)]
  simu$time_of_day_label <- paste0('time_of_day_', simu$time_of_day)
  
  ## return "simu" dataset
  return(simu)
  
}

### my gradient descent function
cppFunction('
std::vector<NumericMatrix> gradient_descent(int nb_saves, int previous_nb_iterations, NumericMatrix set_alpha, NumericVector sample_index, int length_sample_index, NumericVector simu_tod, NumericVector simu_cgm, NumericVector simu_reward, NumericVector simu_insulin, NumericVector insulin_test, NumericVector cgmRef, NumericVector insulinRef, double h_cgm, double h_insulin, double gamma, double final_learning_rate, int nb_splines, int length_insulin_test){

  std::vector<NumericMatrix> set_alpha_memory;
  NumericMatrix set_alpha_frozen = set_alpha ; 
  NumericVector basis_index = nb_splines ;
  double target_index, val0, val1, val2, new_val, next_qsa, prev_qsa ;
  int i, tod, next_tod, j, k ;
  int memory_index = 0 ;
  double percent_update = 0.05 ; 
  int threshold_for_update = length_sample_index * percent_update, counter_since_update = threshold_for_update ; 
  double pct = 0, val = 0 ;

  // ITERATE OVER sample_index ELEMENTS
  for(int ii = 0 ; ii < length_sample_index ; ii++){
    
    // UPDATE STATUS, EVERY percent_update %
    counter_since_update += 1 ; 
    if(counter_since_update > threshold_for_update){
      val = ii + 1 ;
      pct = val / length_sample_index * 100 ;
      Rcout << pct << "%" << std::endl ;
      counter_since_update = 0 ;
    }

    // IDENTIFY ROW INDEX OF DATASET
    i = sample_index[ii] ;
    
    // IDENTIFY TIME-OF-DAY FOR INDEX i and i+1
    next_tod = simu_tod[i+1]-1 ;
    tod = simu_tod[i]-1 ; 
    
    // COMPUTING TARGET VALUE
    next_qsa = 0.00 ;
    for(j = 0 ; j < length_insulin_test ; j++){
      new_val = 0.00 ;
      for(k = 0 ; k < nb_splines ; k++){
        val0 = 2 * PI * h_cgm * h_insulin ;
        val1 = exp(- pow(simu_cgm[i+1] - cgmRef[k], 2) / (2 * pow(h_cgm, 2))) ;
        val2 = exp(- pow(insulin_test[j] - insulinRef[k], 2) / (2 * pow(h_insulin, 2))) ;
        new_val += set_alpha_frozen(k, next_tod) * val1 * val2 / val0 ;
      }
      if(new_val > next_qsa){
        next_qsa = new_val ;
      }
    }
    target_index = simu_reward[i] + gamma * next_qsa ; 
    
    // COMPUTING CURRENT QSA ESTIMATE
    prev_qsa = 0.00 ;
    for(k = 0 ; k < nb_splines ; k++){
      val0 = 2 * PI * h_cgm * h_insulin ;
      val1 = exp(- pow(simu_cgm[i] - cgmRef[k], 2) / (2 * pow(h_cgm, 2))) ;
      val2 = exp(- pow(simu_insulin[i] - insulinRef[k], 2) / (2 * pow(h_insulin, 2))) ;
      basis_index[k] = val1 * val2 / val0 ;
      prev_qsa += set_alpha(k, tod) * basis_index[k] ;
    }

    // PERFORM GRADIENT DESCENT
    for(k = 0 ; k < nb_splines ; k++){
      set_alpha(k, tod) = set_alpha(k, tod) - 2 * final_learning_rate * basis_index[k] * (prev_qsa - target_index) ;
    }
    
    // PERIODICALLY UPDATE FROZEN APPROXIMATION
    if(ii%100 == 1){
      set_alpha_memory.push_back(clone(set_alpha)) ; 
      // for(k = 0 ; k < nb_splines ; k++){
      //  for(int l = 0; l < 4 ; l++){
      //    set_alpha_memory(k, l, memory_index) = set_alpha(k,l) ;
      //  }
      // }
      memory_index += 1 ; 
    }
    
  }
  return(set_alpha_memory) ;
}')
