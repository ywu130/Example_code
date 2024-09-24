#!/usr/bin/env Rscript

# Script for computing posterior probabilities of hidden states at each
#   position of a given sequence.
# Arguments:
#   f: file containing the sequence (fasta file)
#   mu: the probability of switching states
#
# Outputs:
#    posteriors.csv - a KxN matrix outputted as a CSV with the posterior
#                     probability of each state at each position
#
# Example Usage:
#    Rscript 1b.R hmm-sequence.fa 0.01

# Reads the fasta file and outputs the sequence to analyze.
# Arguments:
#    filename: name of the fasta file
# Returns:
#    s: string with relevant sequence
read_fasta <- function(filename) {
  con <- file(filename, "r")
  s_name <- readLines(con, n = 1)
  s <- readLines(con, n = 1)
  close(con)
  return(s)
}


sumLogProb= function(a,b){
  if(a>b){return (a + log1p(exp(b-a)))}else{
    return(b + log1p(exp(a-b)))
  }
}

# Outputs the forward and backward probabilities of a given observation.
# Arguments:
#        obs: observed sequence of emitted states (vector of emissions)
#        trans_probs: transition log-probabilities (KxK matrix)
#        emiss_probs: emission log-probabilities (KxE matrix, E: # of emissions)
#        init_probs: initial log-probabilities for each hidden state (vector)
# Returns:
#        Fwd: matrix of forward probabilities
#        likelihood_f: P(obs) calculated using the forward algorithm
#        Back: matrix of backward probabilities
#        likelihood_b: P(obs) calculated using the backward algorithm
#        Posteriors: matrix of posterior probabilities
forward_backward <- function(obs, trans_probs, emiss_probs, init_probs) {
  # COMPLETE THIS FUNCTION
  states <- c('h', 'l')
  ncols = nchar(obs)
  nrows = nrow(trans_probs)
  Fwd = matrix(rep(0,ncols*nrows),nrow = nrows,ncol=ncols)
  likelihood_f = 0
  Back = matrix(rep(0,ncols*nrows),nrow = nrows,ncol=ncols)
  likelihood_b = 0
  posterior = matrix(rep(0,ncols*nrows),nrow = nrows,ncol=ncols)
  #forward initiazation:
  for(i in 1:nrows){
    Fwd[i,1] = init_probs[states[i]] + emiss_probs[states[i],substr(obs,1,1)]
  }
  #iteration
  for(j in 2:ncols){
    for(i in 1:nrows){
      pre_sum = Fwd[1,j-1] + trans_probs[states[1],states[i]]
      for(i2 in 2:nrows){
        pre_sum = sumLogProb(pre_sum,Fwd[i2,j-1]+trans_probs[states[i2],states[i]]) 
      }
      Fwd[i,j] = emiss_probs[states[i],substr(obs,j,j)] + pre_sum
    }
    
  }
  #calculate forward P(obs)
  likelihood_f = Fwd[1,ncols]
  for(i in 2:nrows){
    likelihood_f = sumLogProb(likelihood_f,Fwd[i,ncols])
  }
  
  #backward initiazation:
  #backward start with bk_N = 1, with log prob 0
  #has been initiated
  #iteration for backward
  for(j in (ncols-1):1){
    for(i in 1:nrows){
      Back[i,j] = trans_probs[states[1],states[i]] + emiss_probs[states[1],substr(obs,j+1,j+1)] + Back[1,j+1] 
      for(i2 in 2:nrows){
        Back[i,j] = sumLogProb(Back[i,j],trans_probs[states[i2],states[i]] + emiss_probs[states[i2],substr(obs,j+1,j+1)] + Back[i2,j+1])
      }
    }
    
  }
  likelihood_b = init_probs[states[1]] + emiss_probs[states[1],substr(obs,1,1)] + Back[1,1]
  for(i in 2:nrows){
    likelihood_b = sumLogProb(likelihood_b,init_probs[states[i]] + emiss_probs[states[i],substr(obs,1,1)] + Back[i,1])
  }
  
  
  #posterior
  denominator = rep(0,ncols)
  for(j in 1:ncols){
    denominator[j] = Fwd[1,j] + Back[1,j]
    for(i in 2:nrows){
      denominator[j] = sumLogProb(denominator[j],Fwd[i,j] + Back[i,j])
    }
    for(i in 1:nrows){
      posterior[i,j] = exp(Fwd[i,j] + Back[i,j]-denominator[j])
    }
    
  }
  print(posterior)
  
  return(list(Fwd = Fwd, likelihood_f = likelihood_f,Back = Back, likelihood_b = likelihood_b, posterior = posterior))
  
}

main <- function() {

  obs_sequence <- 'TTTAGCACCGGATGCGGTATCAATCCTGGTATCGTTAAAGCCTAGTGTTTCAAAAGTTCGAAAAACGGGCCGCGCAGCCCCGGGGCGTTCCCGGTAGGCTCGCGGCGGCTGGAGCACAGCGCCGCCGCGCCAGGGGACCCGCCGACTCCGTGGGCGTTCGGACCGCCTGTGTCCCTTCCGGAAGGCCTGGCAGGGTCGCCTCTTCCCGGAGCGGCCGCCCACACCGCGGCGCATGCAGGGCCCCGGTACGTTGACGTCTCTGACCCGCCGCTGCGGTGCGCCGGCGCGCCCGGCCGGGCGCCCCGGCGCGGGGCTCCCACGAGGCCCGCATCACGCGCGCCACCCAGGTTACCGGGCTCGGCCGGTGGACCATCGGGAGGGGCGGGGGCGGACGCCCGACCCCGCGGGCCACTCCAGTTTCTTAACTTGATAACGAAGCACTGAATACAGAAACAAGTTAAATCTCCCGGGTGTCGACGCGGTCCTGGTAACTTTTCAAGTAGGGTCGGTCAATGAAGGAGTGTTAGGCCTGGCCCGCCGGAGGGTTGAAGGACGGGTACGCGGTGTCAACGCCCGCCCCGGTCCGCCGAAGTCCCCCTTCACCGAGGGCCGAGCCGTCCTGCCCCCGGAGTTTCGGCGCGACCCCTACCGGGCTGGCCGCCGGGCGGGAGTCCCGGCTGTGGGCGGCTGGTCACCGGGTTTGGACCTCGGCAGGACAGCAGTTTGCATTCATATACAGGAAAAACACTCCCCACATGTGGTAAATACGTTGAACAAAGTATGTTACCATGAATAAGTAGCACAGCTAGTATTTTGTTGTTAGCAAGGAAAATAAATTCGTAAATTATTATTAATACGTATTTATTGAGTACTTTATTAAACTAAATATGTAGTTGATATCAGGCTGACACAACAACATGGTCTAGCAAAGTTCCTAACATTGGTATTTACCGCCACTTGCAGGCTCGAGGCGCCAGAGGGGCCCTGCGCTAGGTGCGGC'
  mu <- 0.01
  states <- c('h', 'l')
  transition_probabilities <- matrix(runif(4), nrow=2, ncol=2)
  rownames(transition_probabilities) <- states
  colnames(transition_probabilities) <- states
  transition_probabilities['h', 'h'] <- transition_probabilities['l', 'l'] <-
    log(1-mu)
  transition_probabilities['h', 'l'] <- transition_probabilities['l', 'h'] <-
    log(mu)
  emission_probabilities <- matrix(runif(8), nrow=2, ncol=4)
  rownames(emission_probabilities) <- states
  colnames(emission_probabilities) <- c('A', 'C', 'G', 'T')
  emission_probabilities['h', 'A'] <- emission_probabilities['h', 'T'] <-
    log(0.13)
  emission_probabilities['h', 'C'] <- emission_probabilities['h', 'G'] <-
    log(0.37)
  emission_probabilities['l', 'A'] <- emission_probabilities['l', 'T'] <-
    log(0.32)
  emission_probabilities['l', 'C'] <- emission_probabilities['l', 'G'] <-
    log(0.18)
  initial_probabilities <- log(c(0.5, 0.5))
  names(initial_probabilities) <- states
  
  results <- forward_backward(obs_sequence, transition_probabilities,
                              emission_probabilities, initial_probabilities)
  Fwd <- results$Fwd
  like_f <- results$likelihood_f
  print(c('like_f',like_f))
  Back <- results$Back
  like_b <- results$likelihood_b
  Posteriors <- results$posterior
  plot(seq(1,nchar(obs_sequence)),Posteriors[1,],type = 'l')

  
}

main()
