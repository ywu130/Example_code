#!/usr/bin/env Rscript

# Script for computing h_{Y_n} given n and the probabilities of observing 1-6.
# Arguments:
#    n - integer representing the number of observations
#    p - sequence of 6 float values summing to 1 representing pi_i, i in [1,6]
#
# Outputs:
#    h_dp_<n>.csv - a file the probabilities h_{Y_n}(y), y in [n,6n] in CSV
#                   format, for use in later parts of question 1.
#
# Example usage with n = 50, pi = [.1, .2, .2, .2, .1, .2]:
#    Rscript 1c.R  50  .1 .2 .2 .2 .1 .2
#
#    This would output a file called h_dp_50.csv.


# Computes h_{Y_n}(y) for all y in {n, ..., 6n} for a given n and pi.
# Arguments:
#    n - the number of random variables being summed
#    pi - the probabilities of 1-6 for a single observation
#
# Returns:
#    a vector of probabilities for y in {n, ..., 6n}.
#        index 1 should correspond to y = n, index 2 to n+1, etc.
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


# Returns the minimum ten probabilities of a given array.
# Arguments:
#    n - the number of random variables being summed
#    probs - the probabilities of [n, 6n]
# Returns:
#    a vector of indices indicating the minimum 10 probabilities in probs
#        priority is given to higher indices in case of ties
#        should be a list of values in the range [n, 6n]
min10 <- function(n, probs) {
  #***YOUR CODE FOR FINDING AND RETURNING THE INDICES, FROM probs, OF THE
  #   10 LOWEST PROBABILITY VALUES HERE***
  sorted_low = sort(probs)
  #print(sorted_low)
  low_index = rep(0,10)
  for(i in 1:10){
    low_index[i] = which(sorted_low[i]==probs)
  }
  return(low_index+n-1)
}


# Returns the maximum ten probabilities of a given array.
# Arguments:
#    n - the number of random variables being summed
#    probs - the probabilities of [n, 6n]
# Returns:
#    a vector of indices indicating the maximum 10 probabilities in probs
#        priority is given to higher indices in case of ties
#        should be a list of values in the range [n, 6n]
max10 <- function(n, probs) {
  #***YOUR CODE FOR FINDING AND RETURNING THE INDICES, FROM probs, OF THE
  #   10 HIGHEST PROBABILITY VALUES HERE***
  sorted_high = sort(probs, decreasing = TRUE)
  #print(sorted_high)
  high_index = rep(0,10)
  for(i in 1:10){
    high_index[i] = which(sorted_high[i]==probs)
  }
  return(high_index+n-1)
}
#pi = c(1/2, 1/4, 1/8, 1/16, 1/32, 1/32)
#probs.test = h_Y(2,pi)
#print(probs.test)
#print(min10(50,probs.test))
#print(max10(50,probs.test))


# Computes the mean and variance of an array of probabilities for [n, 6n].
# Arguments:
#    n - the value of n defining the range [n, 6n]
#    probs - the array of probabilities
# Returns tuple with two elements:
#    expect - the expectation
#    variance - the variance
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

#pi = c(1/2, 1/4, 1/8, 1/16, 1/32, 1/32)
#probs.test = h_Y(2,pi)
#print(probs.test)
#print(compute_stats(50,probs.test))

# Prints a given value of y together with its probability.
# Arguments:
#    y - the observed y value to print the probability of
#    n - the number of random variables being summed
#    probs - the probabilities of [n, 6n]
# Returns:
#    nothing
print_y_prob <- function(y, n, probs) {
  cat(sprintf("%d %.4g\n", y, probs[y - n + 1]))
}


main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 7) {
    cat("usage: Rscript 1c.R N p1 p2 p3 p4 p5 p6\n")
    cat("  where N is the number of random variables being summed,\n")
    cat("        p1 ... p6 are the probabilities of observing 1 through 6\n")
    cat("\n")
    stop()
  }
  n <- as.integer(args[1])
  pi <- as.double(args[2:7])
  print(pi)
  print(sum(pi))
  if (sum(pi) != 1.0) {
    stop("pi values do not sum to 1")
  }
  
  h <- h_Y(n, pi)
  err <- 1.0 - sum(h)
  if (err * err > 10 ^ -10) {
    stop("h probabilities do not sum to 1")
  }
  min10_indices <- min10(n, h)
  max10_indices <- max10(n, h)
  
  cat("Min probabilities:\n")
  dontcare <- lapply(min10_indices, print_y_prob, n, h)
  cat("\nMax probabilities:\n")
  dontcare <- lapply(max10_indices, print_y_prob, n, h)
  write.table(format(h, digits=5), file = paste("h_", n, "_dp.csv", sep = ""),
              col.names = FALSE, row.names=seq(n, 6*n), quote = FALSE,
              sep=",")
  mean_var <- compute_stats(n, h)
  cat(sprintf("\nMean is %.4g\n", mean_var[1]))
  cat(sprintf("Variance is %.4g\n", mean_var[2]))
}

main()