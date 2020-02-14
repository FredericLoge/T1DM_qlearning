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
bg_reference_levels <- c('low' = 40, 'hypo' = 70, 'center' = 112.5, 'hyper' = 180, 'high' = 350)
bg <- seq(from = bg_reference_levels['low'], to = bg_reference_levels['high'], length.out = 1000)
to_plot <- data.frame(bg, bg_reward = bg_reward(bg))
library(ggplot2)
ggplot(data = to_plot) +
  geom_line(mapping = aes(x = bg, y = bg_reward), col = 'blue', lwd = 1.5) +
  geom_vline(xintercept = bg_reference_levels['center'], lty = 3) +  
  geom_vline(xintercept = bg_reference_levels[c('hypo', 'hyper')], lty = 2) + 
  geom_vline(xintercept = bg_reference_levels[c('low', 'high')], lty = 1) +
  xlab('Blood glucose observed (mg/dL)') +
  ylab('') +
  ggtitle('Reward function associated to blood glucose readings')
  
  
  
  
  
  