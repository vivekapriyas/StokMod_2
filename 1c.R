
#parameters
lambda = 0.01
mu_L = 1/7
mu_H = 1/20
alpha = 0.1

#time
Yrs = 1000
P = 365*Yrs
f = 365*5 #time for plotting

#starting conditions
i = 1.0
s = 0.0
X = c(i)
S = c(t)

while(s<P){
  #finds sojourn time
  s = s + rexp(1,rate=lambda)
  S = c(S,s)
    
  #finds next state
  p = runif(1)
  if (i==1){
    if (p < lambda*alpha){
      i = 3
    }else {
      i = 2
    }
  }else{
    i = 1
  }
  X = c(X,i)
}

"""
plot(NULL, NULL, xlim = c(0,f), ylim = c(0, 4), lwd = 2, cex.lab = 1.5, xlab = "Time", ylab = "State", cex.axis = 1.5,
     main = "X(t)")
for(i in 1:N){
  lines(cT[i:(i+1)], x[i]*c(1,1), lwd = 3)
}
"""

