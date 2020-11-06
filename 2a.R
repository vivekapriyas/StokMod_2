
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



plot(NULL, NULL, xlim = c(0.25, 0.5), ylim = c(0,1.2), main = "Prediction of Y(??) conditional on 5 evauation points", cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5, xlab = "??", ylab = "Y(??)")
legend("top", lty = 1:1, bty = "n", cex = 1, legend = c("E[Y(??)|??_E]", "Prediction interval"), col = c("black","#66CCFF"))
legend(x = 0.3515, y = 1.115, pch = 4, bty = "n", cex = 1, legend = c("Evaluation points"), col = c("red"))
polygon(c(theta_grid, rev(theta_grid)), c(mu_grid_c - 1.64*sd_grid_c, rev(mu_grid_c + 1.64*sd_grid_c)), col = "#CCFFFF", border = NA)
lines(theta_grid, mu_grid_c)
lines(theta_grid, mu_grid_c + 1.64*sd_grid_c, col = "#66CCFF", lwd = 2)
lines(theta_grid, mu_grid_c - 1.64*sd_grid_c, col = "#66CCFF", lwd = 2, legend = "Prediction interval")
points(theta_eval, y_eval, pch = 4, col = "red", lwd = 3)





#b)
a <- 0.3
cdf_t <- c()

for (i in 1:51){
  cdf_t <- c(cdf_t, pnorm(a, mean = mu_grid_c[i], sd = sd_grid_c[i]))
  
}


plot(NULL, NULL, xlim = c(0.25, 0.5), ylim = c(0,0.3), xlab = "??", ylab = "Pr{Y(??)|??_E < 0.30}", main = "Probability of Y(??)<0.30 conditional on 5 evaluation points", cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
lines(theta_grid, cdf_t, col = "#3333CC", lwd = 2)
points(theta_grid, cdf_t, pch = 19, col = "#3333CC", lwd = 2)

#c) 
theta_eval2 <- c(theta_eval, 0.33)
y_eval2 <- c(y_eval, 0.40)

#mu and sigma for observed data points
mu_eval2 <- sapply(theta_eval2,EY)
sigma_eval2 <- sigma(theta_eval2)

#sigma_AB, A = grid, B = eval
sigma_AB2 <- cross_cov(theta_grid, theta_eval2)
sigma_BA2 <- cross_cov(theta_eval2, theta_grid)

#mu_grid_c: conditional mean vector of the process at the 51 grid points conditional on the 5 evaluation points
mu_grid_c2 <- mu_c(mu_grid, sigma_AB2, sigma_eval2, y_eval2, mu_eval2)

#sigma_grid_c: conditional covariance matrix of the process at the 51 grid points conditional on the 6 evaluation points
sigma_grid_c2 <- sigma_c(sigma_grid, sigma_AB2, sigma_eval2, sigma_BA2)
var_grid_c2 <- diag(sigma_grid_c2)+1e-10
sd_grid_c2 <- (var_grid_c2)^(1/2)


plot(NULL, NULL, xlim = c(0.25, 0.5), ylim = c(0,1.2), main = "Prediction of Y(??) conditional on 6 evauation points", cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5, xlab = "??", ylab = "Y(??)")
legend("top", lty = 1:1, bty = "n", cex = 1, legend = c("E[Y(??)|??_E]", "Prediction interval"), col = c("black","#66CCFF"))
legend(x = 0.3515, y = 1.115, pch = 4, bty = "n", cex = 1, legend = c("Evaluation points"), col = c("red"))
polygon(c(theta_grid, rev(theta_grid)), c(mu_grid_c2 - 1.64*sd_grid_c2, rev(mu_grid_c2 + 1.64*sd_grid_c2)), col = "#CCFFFF", border = NA)
lines(theta_grid, mu_grid_c2)
lines(theta_grid, mu_grid_c2 + 1.64*sd_grid_c2, col = "#66CCFF", lwd = 2)
lines(theta_grid, mu_grid_c2 - 1.64*sd_grid_c2, col = "#66CCFF", lwd = 2, legend = "Prediction interval")
points(theta_eval2, y_eval2, pch = 4, col = "red", lwd = 3)



a <- 0.3
cdf_t2 <- c()

for (i in 1:51){
  cdf_t2 <- c(cdf_t2, pnorm(a, mean = mu_grid_c2[i], sd = sd_grid_c2[i]))
  
}


plot(NULL, NULL, xlim = c(0.25, 0.5), ylim = c(0,0.3), xlab = "??", ylab = "Pr{Y(??)|??_E < 0.30}", main = "Probability of Y(??)<0.30 conditional on 6 evaluation points", cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
lines(theta_grid, cdf_t2, col = "#3333CC", lwd = 2)
points(theta_grid, cdf_t2, pch = 19, col = "#3333CC", lwd = 2)


print(theta_grid[which(cdf_t2 == max(cdf_t2))])




