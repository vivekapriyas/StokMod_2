
#observated theta and y
theta_eval <- c(0.30,0.35,0.39,0.41,0.45)
y_eval <- c(0.5,0.32,0.40,0.35,0.60)

#mean, variance and correlation functions
EY <- function(theta){0.5}
VarY <- function(theta){0.5**2}
CorrY <- function(theta1,theta2){(1+15*abs(theta1-theta2))*exp(-15*abs(theta1-theta2))}

#function to create sigma matrix from theta vector
sigma <- function(x){
  n <- length(x)
  S <- matrix(0,nrow = n,ncol = n,byrow = TRUE)
  
  for (i in 1:n) {
    for (j in 1:n){
      if (i == j){
        S[i,j] <- VarY(x[i])
      }else{
        S[i,j] <- CorrY(x[i],x[j])
      }
    }    
  }
  return(S)
}

#mu and sigma for observed data points
mu_eval <- sapply(theta_eval,EY)
sigma_eval <- sigma(theta_eval)

#gridpoints
theta_grid <- seq(from = 0.25, to = 0.5, by = 0.005 )
mu_grid <- sapply(theta_grid,EY)
sigma_grid <- sigma(theta_grid)

