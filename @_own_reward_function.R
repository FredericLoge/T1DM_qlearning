# BG reference levels, in mg/dl
bg_reference_levels <- c('low' = 40, 'hypo' = 70, 'center' = 112.5, 'hyper' = 180, 'high' = 350)

#
foo___ <- function(val, vec){
  (log(val)^vec[1] + vec[2]) * vec[3]
}
foo___2 <- function(val, vec12){
  (log(val)^vec12[1] + vec12[2]) * 80
}
opt <- optim(par = c(2, 0), fn = function(vec){
  vec[3] <- 80
  (foo___2(bg_reference_levels['center'], vec))^2 +
    (foo___2(bg_reference_levels['hypo'], vec) + foo___2(bg_reference_levels['hyper'], vec))^2 +
    (foo___2(bg_reference_levels['low'], vec) + foo___2(bg_reference_levels['high'], vec))^2 +
    1e08 * (vec[1] < 0)
}, control = list('maxit' = 1000, 'abstol' = 1e-10), method = 'BFGS')  
opt

# normalization
vec_sel <- c(opt$par, 80)
vec_sel[3] <- vec_sel[3] / foo___(bg_reference_levels['high'], vec_sel)

#
bg_val <- seq(from = 10, to = 600, length.out = 1000)
plot(x = bg_val, y = foo___(val = bg_val, vec = vec_sel), type = 'l')
abline(h = 0, v = bg_reference_levels['center'])
abline(h = 1, v = bg_reference_levels['high'], lty = 2)
abline(h = -1, v = bg_reference_levels['low'], lty = 2)
val_hyp_thresholds <- foo___(val = bg_reference_levels['hyper'], vec = vec_sel)
abline(h = val_hyp_thresholds, v = bg_reference_levels['hyper'], lty = 2, col = 'red')
abline(h = -val_hyp_thresholds, v = bg_reference_levels['hypo'], lty = 2, col = 'red')

#
new_reward___ <- function(val, vec){
  v <- foo___(val = val, vec = vec)
  v[v > 1] <- 1
  v[v < -1] <- -1
  1 - v^2
}
plot(x = bg_val, y = new_reward___(val = bg_val, vec = vec_sel), type = 'l')
abline(v = bg_reference_levels['center'])
abline(v = bg_reference_levels['high'], lty = 2)
abline(v = bg_reference_levels['low'], lty = 2)
abline(v = bg_reference_levels['hyper'], lty = 2, col = 'red')
abline(v = bg_reference_levels['hypo'], lty = 2, col = 'red')

new_reward___(val = bg_reference_levels, vec = vec_sel)
## new_reward___(val = bg_reference_levels, vec = c(0.4937448, -2.1495359, 4.0955978))
bg_reward <- function(bg){
  bg <- pmax(40, pmin(350, bg))
  new_reward___(val = bg, vec = vec_sel)
}
