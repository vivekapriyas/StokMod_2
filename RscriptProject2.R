#Problem 1 ---- 
library(latexexp)
#(S,I_L,I_H) = (1,2,3)
#parameters
lambda = 0.01
mu_L = 1/7
mu_H = 1/20
alpha = 0.1

ST = c(lambda,mu_L,mu_H)

#time
Yrs = 1000
P = 365*Yrs
f = 365*5 #time for plotting

#starting conditions
i = 1.0
s = 0.0
X = c(i)
S = c(s)

#long-run mean fractions of time
timeS = 0.0
timeI_L = 0.0
timeI_H = 0.0

I_H_times = c()
lastI_H = 0.0

while(s<P){
  #finds sojourn time
  this_sojourn_time =rexp(1,rate=ST[i])
  
  if (s + this_sojourn_time > P){
    this_sojourn_time = P - s
  }
  
  s = s + this_sojourn_time
  
  
  S = c(S,s)
  
  #finds next state
  p = runif(1)
  if (i==1){
    timeS = timeS + this_sojourn_time
    if (p < alpha){
      i = 3
    }else {
      i = 2
    }
  }else if (i == 2){
    timeI_L = timeI_L + this_sojourn_time
    i = 1
  }else if (i == 3){
    I_H_times = c(I_H_times, s-lastI_H)
    
    lastI_H = s
    
    timeI_H = timeI_H + this_sojourn_time
    i = 1
  }
  
  
  X = c(X,i)
}

#c)
Splot = S[S<=f]
N = length(Splot)
Xplot = X[1:N]



plot(NULL, NULL, xlim = c(0,f), ylim = c(0, 4), lwd = 2, xlab = "Time [days]", ylab = "State", cex.axis = 1.5,
     main = TeX("$X(t)$ for t $\\in \\[0, 1825\\]$"), axes = FALSE, frame.plot = TRUE)

axis(2, at = c(0,1,2,3,4), labels = c("", "S", TeX("I_L"), TeX("I_H"), ""), las=0, cex.axis = 1.5)
axis(1, xlim= c(0,f), cex.axis = 1.5)



for(i in 1:N){
  lines(Splot[i:(i+1)], Xplot[i]*c(1,1), lwd = 3)
}
lines(c(Splot[N],f), Xplot[N]*c(1,1), lwd = 3)



#d)
print((timeI_L/P+timeI_H/P)*365)

#e)
print(mean(I_H_times))

#Problem 2 ----

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


plot(NULL, NULL, xlim = c(0.25, 0.5), ylim = c(0,1), main = "Prediction of Y[theta] conditional on observation points", cex.main = 1.0, xlab = "theta", ylab = "Y[theta]")
legend("top", lty = 1:1, bty = "n", cex = 0.8, legend = c("E[Y(theta)|theta_E]", "Prediction interval"), col = c("black","orange"))
legend(x = 0.355, y = 0.95, pch = 4, bty = "n", cex = 0.8, legend = c("Evaluation points"), col = c("red"))
lines(theta_grid, mu_grid_c)
lines(theta_grid, mu_grid_c + 1.64*sd_grid_c, col = "orange", lwd = 2)
lines(theta_grid, mu_grid_c - 1.64*sd_grid_c, col = "orange", lwd = 2, legend = "Prediction interval")
points(theta_eval, y_eval, pch = 4, col = "red", lwd = 3)




#b)
a <- 0.3
cdf_t <- c()

for (i in 1:51){
  cdf_t <- c(cdf_t, pnorm(a, mean = mu_grid_c[i], sd = sd_grid_c[i]))
  
}


plot(NULL, NULL, xlim = c(0.25, 0.5), ylim = c(0,0.3))
#lines(theta_grid, cdf_t)
points(theta_grid, cdf_t)

print(theta_grid[which(cdf_t == max(cdf_t))])

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

#sigma_grid_c: conditional covariance matrix of the process at the 51 grid points conditional on the 5 evaluation points
sigma_grid_c2 <- sigma_c(sigma_grid, sigma_AB2, sigma_eval2, sigma_BA2)
var_grid_c2 <- diag(sigma_grid_c2)+1e-10
sd_grid_c2 <- (var_grid_c2)^(1/2)


plot(NULL, NULL, xlim = c(0.25, 0.5), ylim = c(0,1.2), main = "Prediction of Y[theta] conditional on observation points", cex.main = 1.0, xlab = "theta", ylab = "Y[theta]")
legend("top", lty = 1:1, bty = "n", cex = 0.8, legend = c("E[Y(theta)|theta_E]", "Prediction interval"), col = c("black","orange"))
legend(x = 0.355, y = 1.15, pch = 4, bty = "n", cex = 0.8, legend = c("Evaluation points"), col = c("red"))
lines(theta_grid, mu_grid_c2)
lines(theta_grid, mu_grid_c2 + 1.64*sd_grid_c2, col = "orange", lwd = 2)
lines(theta_grid, mu_grid_c2 - 1.64*sd_grid_c2, col = "orange", lwd = 2, legend = "Prediction interval")
points(theta_eval, y_eval, pch = 4, col = "red", lwd = 3)



a <- 0.3
cdf_t2 <- c()

for (i in 1:51){
  cdf_t2 <- c(cdf_t2, pnorm(a, mean = mu_grid_c2[i], sd = sd_grid_c2[i]))
  
}


plot(NULL, NULL, xlim = c(0.25, 0.5), ylim = c(0,0.3))
points(theta_grid, cdf_t2)


print(theta_grid[which(cdf_t2 == max(cdf_t2))])

