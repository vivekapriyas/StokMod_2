

theta_eval = c(0.30,0.35,0.39,0.41,0.45)
y_eval = c(0.5,0.32,0.40,0.35,0.60)

EY <- function(theta){0.5}
VarY <- function(theta){0.5**2}
CorrY <- function(theta1,theta2){(1+15*abs(theta1-theta2))*exp(-15*abs(theta1-theta2))}

mu_eval = lapply(theta_eval,EY)


theta_grid = seq(from = 0.25, to = 0.5, by = 0.005 )
