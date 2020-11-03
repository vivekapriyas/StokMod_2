
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

while(s<P){
  #finds sojourn time
  s = s + rexp(1,rate=ST[i])
  S = c(S,s)
    
  #finds next state
  p = runif(1)
  if (i==1){
    if (p < alpha){
      i = 3
    }else {
      i = 2
    }
  }else{
    i = 1
  }
  X = c(X,i)
}

Splot = S[S<f]
N = length(Splot)
Xplot = X[1:N]


plot(NULL, NULL, xlim = c(0,f), ylim = c(0, 4), lwd = 2, cex.lab = 1.5, xlab = "Time", ylab = "State", cex.axis = 1.5,
     main = "X(t)")
for(i in 1:N){
  lines(Splot[i:(i+1)], Xplot[i]*c(1,1), lwd = 3)
}


