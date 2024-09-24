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


#sample state only

samplestate <- function(obs, trans_probs, emiss_probs, init_probs) {
  # COMPLETE THIS FUNCTION
  states <- c('h', 'l')
  ncols = nchar(obs)
  nrows = nrow(trans_probs)
  Fwd = matrix(rep(0,ncols*nrows),nrow = nrows,ncol=ncols)
  #likelihood_f = 0
  numtrans = rep(0,50)
  samplestate = matrix(,nrow = 50,ncol = 1000) #matrix to store state 
  samplikelihood = rep(0,50)
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
  #likelihood_f = Fwd[1,ncols]
  #for(i in 2:nrows){
    #likelihood_f = sumLogProb(likelihood_f,Fwd[i,ncols])
  #}
  
  
  #initialize the last colomun
  sumkN = sumLogProb(Fwd[1,1000],Fwd[2,1000])
  samplestate[,1000] = sample(c('h','l'),50,replace = TRUE,prob = c(exp(Fwd[1,1000]-sumkN),exp(Fwd[2,1000]-sumkN))) #exp the log probability

  
  for(j in 999:1){
    for(m in 1:50){
      fhj = Fwd[1,j] + trans_probs[states[1],samplestate[m,j+1]]
      flj = Fwd[2,j] + trans_probs[states[2],samplestate[m,j+1]]
      sumkj = sumLogProb(fhj,flj)
      samplestate[m,j] = sample(c('h','l'),1,prob = c(exp(fhj-sumkj),exp(flj-sumkj)) )
      if(samplestate[m,j]!=samplestate[m,j+1]){
        numtrans[m] = numtrans[m]+1
      }
      if(samplestate[m,j] =='h'){
        samplikelihood[m] = samplikelihood[m] +  trans_probs[states[1],samplestate[m,j+1]] + emiss_probs[samplestate[m,j+1],substr(obs,j+1,j+1)]
      }else{
        samplikelihood[m] = samplikelihood[m] +  trans_probs[states[2],samplestate[m,j+1]] + emiss_probs[samplestate[m,j+1],substr(obs,j+1,j+1)]
      }
      }
  }
  for(i in 1:50){
    samplikelihood[i] = samplikelihood[i] + init_probs[samplestate[i,1]] + emiss_probs[samplestate[i,1],substr(obs,1,1)]
  }
  
  return(list(likelihood = samplikelihood,states = samplestate))
  
}

find_intervals <- function(sequence) {
  # COMPLETE THIS FUNCTION
  h_state = c()
  
  for(i in 2:length(sequence)){
    if(sequence[i-1]=='l'&&sequence[i]=='h'){
      h_state = c(h_state,i)
    }else if(sequence[i-1]=='h'&&sequence[i] == 'l'){
      h_state = c(h_state,i-1)
    }
  }
  if(length(h_state)%%2!=0){
    h_state = c(h_state,length(sequence)) #add last number to the vector in case of no transition at last
  }
  #print(h_state)
  return(h_state)
}


main <- function() {
  #set.seed(1234)
  
  #obs_sequence <- 'TTTAGCACCGGATGCGGTATCAATCCTGGTATCGTTAAAGCCTAGTGTTTCAAAAGTTCGAAAAACGGGCCGCGCAGCCCCGGGGCGTTCCCGGTAGGCTCGCGGCGGCTGGAGCACAGCGCCGCCGCGCCAGGGGACCCGCCGACTCCGTGGGCGTTCGGACCGCCTGTGTCCCTTCCGGAAGGCCTGGCAGGGTCGCCTCTTCCCGGAGCGGCCGCCCACACCGCGGCGCATGCAGGGCCCCGGTACGTTGACGTCTCTGACCCGCCGCTGCGGTGCGCCGGCGCGCCCGGCCGGGCGCCCCGGCGCGGGGCTCCCACGAGGCCCGCATCACGCGCGCCACCCAGGTTACCGGGCTCGGCCGGTGGACCATCGGGAGGGGCGGGGGCGGACGCCCGACCCCGCGGGCCACTCCAGTTTCTTAACTTGATAACGAAGCACTGAATACAGAAACAAGTTAAATCTCCCGGGTGTCGACGCGGTCCTGGTAACTTTTCAAGTAGGGTCGGTCAATGAAGGAGTGTTAGGCCTGGCCCGCCGGAGGGTTGAAGGACGGGTACGCGGTGTCAACGCCCGCCCCGGTCCGCCGAAGTCCCCCTTCACCGAGGGCCGAGCCGTCCTGCCCCCGGAGTTTCGGCGCGACCCCTACCGGGCTGGCCGCCGGGCGGGAGTCCCGGCTGTGGGCGGCTGGTCACCGGGTTTGGACCTCGGCAGGACAGCAGTTTGCATTCATATACAGGAAAAACACTCCCCACATGTGGTAAATACGTTGAACAAAGTATGTTACCATGAATAAGTAGCACAGCTAGTATTTTGTTGTTAGCAAGGAAAATAAATTCGTAAATTATTATTAATACGTATTTATTGAGTACTTTATTAAACTAAATATGTAGTTGATATCAGGCTGACACAACAACATGGTCTAGCAAAGTTCCTAACATTGGTATTTACCGCCACTTGCAGGCTCGAGGCGCCAGAGGGGCCCTGCGCTAGGTGCGGC'
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 3) {
    cat("usage: Rscript 3b.R file.fasta mu output.txt\n")
    stop()
  }
  fasta_file <- args[1]
  mu <- as.numeric(args[2])
  intervals_file <- args[3]
  
  obs_sequence <- read_fasta(fasta_file)
  states <- c('h', 'l')
  transition_probabilities <- matrix(runif(4), nrow=2, ncol=2)
  rownames(transition_probabilities) <- states
  colnames(transition_probabilities) <- states
  
  
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
  
  mu <- 0.01

  transition_probabilities['h', 'h'] <- transition_probabilities['l', 'l'] <-
    log(1-mu)
  transition_probabilities['h', 'l'] <- transition_probabilities['l', 'h'] <-
    log(mu)
  
  results <- samplestate(obs_sequence, transition_probabilities,
                              emission_probabilities, initial_probabilities)
  #print(results$likelihood)

  teststate = results$states
  for(i in 1:50){
    find_intervals(teststate[i,])
  }
  count = rep(0,1000)
  for(i in 1:50){
    for(j in 1:1000){
      if(teststate[i,j]=='h'){
        count[j] = count[j] + 1
      }
    }
  }
  #plot(seq(1,1000),count,type = 'l')
  plot(seq(1,50),results$likelihood)
  #count those better than the viterbi
  count = 0
  for(i in 1:50){
    if(results$likelihood[i]>-1303.642){
      count = count + 1
    }
  }
  #print(count)
  y = matrix('',nrow = 50, byrow = TRUE)
  
  for(i in 1:50){
    for(j in 1:1000){
      y[i] = paste0(y[i],teststate[i,j])
    }
  }
  
  file_con <- file(intervals_file, "w")
  for (index in seq(1, length(y), 1)) {
    cat(sprintf("%s\n", y[index]),
        file = file_con)
  }
  close(file_con)

}

main()
