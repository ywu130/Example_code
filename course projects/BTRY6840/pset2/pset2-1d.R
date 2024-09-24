#pset2-1d

pi = c(1/2,1/4,1/8,1/16,1/32,1/32)


h_Y <- function(n, pi) {
  #***YOUR CODE FOR CALCULATING THE PROBABILITIES OF h_Yn HERE***
  sum_n_f = 1:6
  prob_f = pi
  for(i in 2:n){
    sum_n = i:(6*i)
    prob = rep(0,5*i+1)
    for(j in sum_n_f){ #iterate values from former value vector
      for(x in 1:6){ #iterate 1 to 6
        index = which((j+x)==sum_n)
        prob[index] = prob[index] + prob_f[which(j==sum_n_f)]*pi[x] 
      }
    }
    sum_n_f = sum_n
    prob_f = prob
  }
  return(prob)
}

compute_stats <- function(n, probs) {
  # initialize
  expect <- 0.0
  variance <- 0.0
  #***YOUR CODE FOR COMPUTING THE EXPECTATION AND VARIANCE FROM probs HERE***
  value = rep(0, 5*n+1) #record weighted value
  index_var = seq(n,6*n)
  for(i in 1:(5*n+1)){
    value[i] = (i+n-1)*probs[i]
  }
  expect = sum(value)
  variance = sum((index_var - expect)^2*probs)
  return(c(expect, variance))
}


h_Y_50 = h_Y(50,pi)

h_mu_var = compute_stats(50,h_Y_50)
print(h_mu_var)


#in geometrical model:
#E(Xi) = sum( xi*pi)
Ex = sum(pi*seq(1,6))
Var.X = sum((seq(1,6)-Ex)^2*pi)
n = 50
#1/phi = Ex
phi = 1/Ex

g_Y = function(n,phi){
  g_probs = rep(0,5*n+1)
  for(y in n:(6*n)){
    #g_probs[y-n+1] = dbinom(n-1, y-1, phi)*phi
    #g_probs[y-n+1] = choose(y-1,n-1)*(1-phi)^(y-n)*(phi^n)
    g_probs[y-n+1] = dnbinom(y-n, n, phi, log = FALSE)
    
  }
  return(g_probs)
}


g_mu = n /phi
g_var = n*(1-phi)/phi^2
print(c(g_mu,g_var))
g_Y_50 = g_Y(50,phi)
g_mu2 = sum(g_Y_50 * seq(50,300))
g_var2 = sum((seq(50,300)-g_mu2)^2*g_Y_50)





#in normal distribution model
f_mu = n*Ex
var_f = n * Var.X

print(c(f_mu,var_f))
f_Y = function(n,mu,sd){
  f_Y_probs = rep(0,5*n+1)
  for(y in n:(6*n)){
    f_Y_probs[y-n+1] = dnorm(y,mean = mu, sd =sd)
  }
  return(f_Y_probs)
}
f_Y_50 = f_Y(50,f_mu,sqrt(var_f))






#plot

plot(seq(50,300),h_Y_50, ylab='probability')
par(new=TRUE)
plot(seq(50,300),g_Y_50,type = 'l',col = 'red',xaxt ='n',yaxt = 'n',ylab='',xlab='')
par(new=TRUE)
plot(seq(50,300),f_Y_50,type = 'l',col = 'blue',xaxt ='n',yaxt = 'n',ylab='',xlab='')
legend('topright' ,legend = c('h_Y','g_Y','f_Y'),col = c('black','red','blue') )

print(pnbinom(249,50,phi,lower.tail = F))
print(pnorm(299,f_mu,sqrt(var_f),lower.tail = F))
#par(new=TRUE)
#plot(seq(50,300),f_Y_50_2,type = 'l',col = 'green',xaxt ='n',yaxt = 'n',ylab='',xlab='')






