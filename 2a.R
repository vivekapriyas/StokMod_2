
#observated theta and y
theta_eval <- c(0.30,0.35,0.39,0.41,0.45)
y_eval <- c(0.5,0.32,0.40,0.35,0.60)

#mean, variance and correlation functions
EY <- function(theta){0.5}
VarY <- function(theta){0.5^2}
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
        S[i,j] <- CorrY(x[i],x[j])*0.5^2
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
cross_cov <- function(x_A, x_B){
  n <- length(x_A)
  m <- length(x_B)
  S <- matrix(0,nrow = n,ncol = m,byrow = TRUE)
  
  for (i in 1:n){
    for (j in 1:m){
      S[i,j] = CorrY(x_A[i],x_B[j])*0.5^2
    }
  }
  return(S)
}

  
#sigma_AB, A = grid, B = eval
sigma_AB <- cross_cov(theta_grid, theta_eval)
sigma_BA <- cross_cov(theta_eval, theta_grid)

#mu_grid_c: conditional mean vector of the process at the 51 grid points conditional on the 5 evaluation points
mu_grid_c <- mu_c(mu_grid, sigma_AB, sigma_eval, y_eval, mu_eval)

#sigma_grid_c: conditional covariance matrix of the process at the 51 grid points conditional on the 5 evaluation points
sigma_grid_c <- sigma_c(sigma_grid, sigma_AB, sigma_eval, sigma_BA)
var_grid_c <- diag(sigma_grid_c)
sd_grid_c <- (var_grid_c)^(1/2)


plot(NULL, NULL, xlim = c(0.25, 0.5), ylim = c(0,1))
points(theta_grid, mu_grid_c)


#b) 
#Y(theta) < 0.30
q_grid = seq(from = 0, to = 1, by = 0.005)


cdf_grid_c <- c()

for (i in 1:51){
  Q_grid <- c(Q_grid, q_grid)
  cdf_grid_c <- c(cdf_grid_c,pnorm(q = q_grid , mean = mu_grid_c[i], sd = sd_grid_c[i]))
}


len_q = length(q_grid)

plot(NULL, NULL, xlim = c(0, 1), ylim = c(0,1))
for (i in seq(from = 1, to = len_q*51, by  = len_q)){
  
  lines(q_grid, cdf_grid_c[i:(i+len_q-1)], col = i)
}


Q_grid[max(which(cdf_grid_c < 0.3))]

