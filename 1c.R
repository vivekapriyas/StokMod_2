
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

Splot = S[S<=f]
N = length(Splot)
Xplot = X[1:N]


plot(NULL, NULL, xlim = c(0,f), ylim = c(0, 4), lwd = 2, cex.lab = 1.5, xlab = "Time", ylab = "State", cex.axis = 1.5,
     main = "X(t)")
for(i in 1:N){
  lines(Splot[i:(i+1)], Xplot[i]*c(1,1), lwd = 3)
}
lines(c(Splot[N],f), Xplot[N]*c(1,1), lwd = 3)



#d)

print(timeS/P)
print(timeI_L/P)
print(timeI_H/P)

print((timeI_L/P+timeI_H/P)*365)

#e)
print(I_H_times)
print(mean(I_H_times))
print(sum(I_H_times))
print(P)
