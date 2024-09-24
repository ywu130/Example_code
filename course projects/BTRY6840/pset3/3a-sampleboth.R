read_fasta <- function(filename) {
  con <- file(filename, "r")
  s_name <- readLines(con, n = 1)
  s <- readLines(con, n = 1)
  close(con)
  return(s)
}

main <- function(){

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2) {
    cat("usage: Rscript 3a.R mu output.txt\n")
    stop()
  }
  
  mu <- as.numeric(args[1])
  out_file <- args[2]

  #obs_sequence <- read_fasta(fasta_file)
  states <- c('h', 'l')
  transition_probabilities <- matrix(runif(4), nrow=2, ncol=2)
  rownames(transition_probabilities) <- states
  colnames(transition_probabilities) <- states
  transition_probabilities['h', 'h'] <- transition_probabilities['l', 'l'] <-log(1 - mu)
  transition_probabilities['h', 'l'] <- transition_probabilities['l', 'h'] <-log(mu)
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
  
  samplestate = matrix(,nrow = 50,ncol = 1000) #matrix to store state 
  sampleemiss = matrix(,nrow = 50,ncol = 1000) #matrix to store smission
  #initialization
  samplestate[,1] = sample(c('h','l'),size=50,replace = TRUE, prob = c(0.5,0.5))
  for(i in 1:50){
    sampleemiss[i,1] = sample(c('A', 'C', 'G', 'T'),1,prob = exp(c(emission_probabilities[samplestate[i,1],'A'],emission_probabilities[samplestate[i,1],'C'],emission_probabilities[samplestate[i,1],'G'],emission_probabilities[samplestate[i,1],'T'])))
  }
  
  f_likelihood = rep(0,50)
  for(i in 1:50){
    f_likelihood[i] = log(0.5) + emission_probabilities[samplestate[i,1],sampleemiss[i,1]]
  }
  
  
  numtrans = rep(0,50)
  for(i in 1:50){
    for(j in 2:1000){
      #print(sample(c('h','l'),size=1,prob =exp(c(transition_probabilities[samplestate[i,j-1],'h'],transition_probabilities[samplestate[i,j-1],'l']))))
      samplestate[i,j] = sample(c('h','l'),1,prob =exp(c(transition_probabilities[samplestate[i,j-1],'h'],transition_probabilities[samplestate[i,j-1],'l'])))
      sampleemiss[i,j] = sample(c('A', 'C', 'G', 'T'),1,prob = exp(c(emission_probabilities[samplestate[i,j],'A'],emission_probabilities[samplestate[i,j],'C'],emission_probabilities[samplestate[i,j],'G'],emission_probabilities[samplestate[i,j],'T'])) )
      f_likelihood[i] = f_likelihood[i] + emission_probabilities[samplestate[i,j],sampleemiss[i,j]]
      if(samplestate[i,j]!=samplestate[i,j-1]){
        numtrans[i] = numtrans[i]+1
      }
      
    }
  }
  outformat_e = matrix('',nrow = 50,ncol = 1,byrow = TRUE)
  outformat_s = matrix('',nrow = 50,ncol = 1,byrow = TRUE)
  
  for(i in 1:50){
    for(j in 1:1000){
      outformat_e[i] = paste0(outformat_e[i],sampleemiss[i,j])
      outformat_s[i] = paste0(outformat_s[i],samplestate[i,j])
    }
  }
  
  file_con <- file(out_file, "w")
  for (index in seq(1, 50, 1)) {
    cat(sprintf("%d %s\n%s\n%s\n%s\n", index, 'Sample States:',outformat_s[index],'Sample emissioin:',outformat_e[index]),
        file = file_con)
  }
  close(file_con)
  

  #print(numtrans)
  #print(mean(numtrans))
  #print(var(numtrans))
  #print(f_likelihood)
  #plot(seq(1,50,1),f_likelihood)
  #abline(h = -1293.1416634001)
  #plot(seq(1,50,1),f_likelihood+1293.1416634001)
  #abline(h = 0)
  
}
main()
