
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


#functions to create conditional mean vector and covariance matrix
mu_c <- function(mu_A, sigma_AB, sigma_BB, x_B, mu_B){
  return(mu_A + sigma_AB%*%solve(sigma_BB)%*%(x_B - mu_B))
}

sigma_c <- function(sigma_AA, sigma_AB, sigma_BB, sigma_BA){
  return(sigma_AA - sigma_AB%*%solve(sigma_BB)%*%sigma_BA)
}
  

#function to create cross-covariance matrix
corr <- function(x_A, x_B){
  n <- length(x_A)
  m <- length(x_B)
  S <- matrix(0,nrow = n,ncol = m,byrow = TRUE)
  
  for (i in 1:n){
    for (j in 1:m){
      S[i,j] = CorrY(x_A[i],x_B[j])
    }
  }
  return(S)
}

  
#sigma_AB, A = grid, B = eval
sigma_AB <- corr(theta_grid, theta_eval)
sigma_BA <- corr(theta_eval, theta_grid)

#mu_grid_c: conditional mean vector of the process at the 51 grid points conditional on the 5 evaluation points
mu_grid_c <- mu_c(mu_grid, sigma_AB, sigma_eval, theta_eval, mu_eval)

#sigma_grid_c: conditional covariance matrix of the process at the 51 grid points conditional on the 5 evaluation points
sigma_grid_c <- sigma_c(sigma_grid, sigma_AB, sigma_eval, sigma_BA)

var_grid_c <- diag(sigma_grid_c)
print(solve(sigma_eval))
print(sigma_eval)
print(var_grid_c)

print(mu_grid_c)

#Nytt forsøk
mu_grid_cNew <- mu_c(mu_grid, sigma_AB, sigma_eval, y_eval, mu_eval)
sigma_grid_cNew <- sigma_c(sigma_grid, sigma_AB, sigma_eval, sigma_BA)
var_grid_cNew <- diag(sigma_grid_cNew)
print(var_grid_cNew)
