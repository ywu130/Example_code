#coal-HMM
#input sequence for 4 species
#estimate parameters with Bam-Walch
#"usage: Rscript 1a.R hmm-sequence.fa s u v a b a2\n"

#argument: 
#f: sequence to read in
#transition prob is determined by s,u,v
#emission prob is determined by a,b,a2
#7 free parameters to estimate from the sequence f

#to read in fasta file
library(seqinr)


states = c('HC1','HC2','HG','CG')
#species in the tree
species = c('H','C','G','O')
u = 0.0001
#longest length:
#timescale = a+b+c = a2+b2+c2
#unit: 18mya
timescale = 36

#sum from a and b in lof space
sumLogProb= function(a,b){
  if(a>b){return (a + log1p(exp(b-a)))}else{
    return(b + log1p(exp(a-b)))
  }
}


# Reads data from filename in fasta format.
#
# Arguments:
#    filename: name of fasta file to read
#    9158 nts including '-'
# Returns:
#    sequences: vector indexed by sequence id containing sequence data (string)
read_fasta <- function(filename) {
  con <- file(filename, "r")
  labels <- c()
  sequences <- c()
  while( TRUE ) {
    label <- substring(readLines(con, n = 1), 9)
    if ( length(label) == 0 ) break
    seq <- readLines(con, n = 1)
    labels <- c(labels, label)
    sequences <- c(sequences, seq)
  }
  close(con)
  names(sequences) <- labels
  return(sequences)
}

#filter out the gaps in the aligned multiple sequneces
gap_filter = function(sequences){
  seqlen = nchar(sequences[1])
  nseq = length(sequences)
  #record the index of non-gaps
  #start from first nucleotide,since the first nucleotide is not gap
  gap_mark = c(1)
  for(i in 1:nseq){
    for(j in 1:(seqlen-1)){
      if(substr(sequences[i],j,j)!='-'&&substr(sequences[i],j+1,j+1)=='-'){
        gap_mark = c(gap_mark,j)
      }else if(substr(sequences[i],j,j)=='-'&&substr(sequences[i],j+1,j+1)!='-'){
        gap_mark = c(gap_mark,j+1)
      }
    }
  }
  gap_mark = c(gap_mark,seqlen)#add the last number
  updated_seq = c()
  while(length(gap_mark)>0){
    non_gap = substr(sequences,gap_mark[1],gap_mark[2])
    updated_seq = paste0(updated_seq,non_gap)
    gap_mark = gap_mark[-c(1,2)] #remove first two elements
  }
  return(updated_seq)
}


#input aligned sequences for 4 speciees
#output into 0111 form
#1 means mismatch from the mutation
#sequnced inupt in the order of species(HCGO, etc)
seqto01 = function(sequences){
  seqlen = nchar(sequences[1])
  nseq = length(sequences)
  seqin01 = rep(0,seqlen)
  for(i in 1:seqlen){
    for(j in 2:nseq ){
      if(substr(sequences[j],i,i)==substr(sequences[1],i,i)){
        seqin01[i] = paste0(seqin01[i],'0')
      }else{
        #means mismatch and 1 mutation happened
        seqin01[i] = paste0(seqin01[i],'1')
      }
    }
  }
  return(seqin01)
}

#jukes cantor:
#return in log scale
jcm <- function(b, a, t, u) {
  if (missing(u)) u <- 0.01
  #t comes from in mya as unit
  t = t*100000
  prob = 0
  if(b == a){
    prob = 1/4*(1+3*exp(-4*u*t/3))
  }else{
    prob = 3/4*(1-exp(-4*u*t/3))
  }
  
  return(log(prob))
}


#build branch length matrix for jcm
#given branch length
#a,b,a2,b2
get_t = function(a,b,a2){
  t = matrix(runif(12), nrow=4, ncol=3)
  rownames(t) = states
  colnames(t) = species[-1]
  #constraints from coalescent recombination model
  b2 = 1.5*(a+b-a2)
  t[1,1] = 2*a
  t[1,2] = 2*(a+b)
  t[1,3] = timescale
  t[2,1] = t[3,2] = 2*a2
  t[2,2] = t[3,1] = t[4,1] = t[4,2] = 2*(a2+b2)
  t[2,3] = t[3,3] = t[4,3] = timescale
  return(t)
}


#esmission probability matrix:
#given state, s
#branch length matrix, t
#sequence alignment string
emiss_P = function(states,t){
  emiss = matrix(0, nrow=4, ncol=8)
  rownames(emiss) = states
  #0 means nt a=b, 1 means mutation detected
  #8 emissions
  colnames(emiss) = c('0000','0001','0011','0111',
                      '0110','0101','0100','0010')
  for(s in states){
    emiss[s,1] = jcm(0,0,t[s,1])+jcm(0,0,t[s,2])+jcm(0,0,t[s,3])
    emiss[s,2] = jcm(0,0,t[s,1])+jcm(0,0,t[s,2])+jcm(0,1,t[s,3])
    emiss[s,3] = jcm(0,0,t[s,1])+jcm(0,1,t[s,2])+jcm(0,1,t[s,3])
    emiss[s,4] = jcm(0,1,t[s,1])+jcm(0,1,t[s,2])+jcm(0,1,t[s,3]) 
    emiss[s,5] = jcm(0,1,t[s,1])+jcm(0,1,t[s,2])+jcm(0,0,t[s,3])
    emiss[s,6] = jcm(0,1,t[s,1])+jcm(0,0,t[s,2])+jcm(0,1,t[s,3])
    emiss[s,7] = jcm(0,1,t[s,1])+jcm(0,0,t[s,2])+jcm(0,0,t[s,3])
    emiss[s,8] = jcm(0,0,t[s,1])+jcm(0,1,t[s,2])+jcm(0,0,t[s,3])  
  }
  #return the emission prob matrix in log scale
  print(exp(emiss))
  return(emiss)
}



#argument: given ekb from EM
#output estimated a,b,a2,b2
decode_aba2b2 = function(ekb,u){
  if (missing(u)) u <- 0.01
  est_t = matrix(0,nrow = 4,ncol = 3)
  rownames(est_t) = states
  colnames(est_t) = species[-1]
  for(i in 1:4){
    #calculate 00s1/01s1
    est_t[i,1] = ekb[i,1]/ekb[i,7] + ekb[i,2]/ekb[i,6] + 
                 ekb[i,3]/ekb[i,4] + ekb[i,8]/ekb[i,5]
    #calculate 00s2/01s2
    est_t[i,2] = ekb[i,1]/ekb[i,8] + ekb[i,2]/ekb[i,3] + 
                 ekb[i,6]/ekb[i,4] + ekb[i,7]/ekb[i,5]
    #calculate 00s3-01s3
    est_t[i,3] = ekb[i,1]/ekb[i,2] + ekb[i,8]/ekb[i,3] + 
                 ekb[i,5]/ekb[i,4] + ekb[i,7]/ekb[i,6]
  }
  est_t = est_t/4 #average of 4 groups of ratio
  #00s1/01s1 = (1+3exp(-4ut/3))/(1-exp(-4ut/3))
  est_t[,3] = mean(est_t[,3])
  for(i in 1:4){
    est_t[i,1] = log((est_t[i,1]-1)/(est_t[i,1]+3))*3/(-4*u)
    est_t[i,2] = log((est_t[i,2]-1)/(est_t[i,2]+3))*3/(-4*u)
    est_t[i,3] = log((est_t[i,3]-1)/(est_t[i,3]+3))*3/(-4*u)
  }
  ratio = timescale/est_t[1,3]
  #corrected to timescale
  est_t = ratio*est_t
  return(est_t)
}



#transition probability matrix
#input: states, s,u,v
trans_P = function(states,s,u,v){
  trans = matrix(0,nrow=4,ncol=4)
  rownames(trans) = colnames(trans) =states
  trans[1,2:4] = s
  trans[2,3:4] = trans[3,2] = trans[3,4] = trans[4,2:3] = v
  trans[2:4,1] = u
  trans[1,1] = 1-3*s
  trans[2,2] = trans[3,3] = trans[4,4] = 1-u-2*v
  #return trans in log space
  #print(trans)
  return(log(trans))
}

# Outputs the forward and backward probabilities of a given observation.
# Arguments:
#        obs: observed sequence of emitted states (vector of emissions)
#        trans_probs: transition log probabilities (KxK matrix)
#        emiss_probs: emission log probabilities (KxE matrix, E: # of emissions)
#        init_probs: initial log probabilities for each hidden state (vector)
# Returns:
#        Fwd: matrix of forward probabilities
#        likelihood_f: P(obs) calculated using the forward algorithm
#        Back: matrix of backward probabilities
#        likelihood_b: P(obs) calculated using the backward algorithm
#        Posteriors: matrix of posterior probabilities
forward_backward <- function(obs, trans_probs, emiss_probs, init_probs) {
  #obs are the aligned seq in 01 form

  ncols = length(obs)
  nrows = nrow(trans_probs)
  Fwd = matrix(0,nrow = nrows,ncol=ncols)
  likelihood_f = 0
  Back = matrix(rep(0,ncols*nrows),nrow = nrows,ncol=ncols)
  likelihood_b = 0
  posterior = matrix(rep(0,ncols*nrows),nrow = nrows,ncol=ncols)
  #forward initiazation:
  for(i in 1:nrows){
    Fwd[i,1] = init_probs[states[i]] + emiss_probs[states[i],obs[1]]
  }
  
  for(j in 2:ncols){
    for(i in 1:nrows){
      pre_sum = Fwd[1,j-1] + trans_probs[states[1],states[i]]
      for(i2 in 2:nrows){
        pre_sum = sumLogProb(pre_sum,Fwd[i2,j-1]+trans_probs[states[i2],states[i]]) 
      }
      Fwd[i,j] = emiss_probs[states[i],obs[j]] + pre_sum
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
      Back[i,j] = trans_probs[states[1],states[i]] + emiss_probs[states[1],obs[j+1]] + Back[1,j+1] 
      for(i2 in 2:nrows){
        Back[i,j] = sumLogProb(Back[i,j],trans_probs[states[i2],states[i]] + emiss_probs[states[i2],obs[j+1]] + Back[i2,j+1])
      }
    }
    
  }
  likelihood_b = init_probs[states[1]] + emiss_probs[states[1],obs[1]] + Back[1,1]
  for(i in 2:nrows){
    likelihood_b = sumLogProb(likelihood_b,init_probs[states[i]] + emiss_probs[states[i],obs[1]] + Back[i,1])
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
  test = matrix(0,nrow = 8,ncol = ncols)
  
  
  print(c(likelihood_f,likelihood_b))

  return(list(Fwd = Fwd, likelihood_f = likelihood_f,Back = Back, likelihood_b = likelihood_b, posterior = posterior))
  
}



# Performs 1 EM step.
# Arguments:
#    fb_values: relevant output variables from forward-backward
#    obs: the sequence under analysis
#    tp: transition probabilities in log space
#    ep: emission probabilities in log space
#    ip: initalization probabilities in log space
# Returns:
#    tp: updated transition probabilities, in log space
#    ep: updated emission probabilities, in log space
em <- function(fb_output, obs, tp, ep) {

  P_x = fb_output$likelihood_f
  fwd = fb_output$Fwd
  back = fb_output$Back
  #initialize expected transtions as a K*K matrix
  Akt=matrix(0,nrow = 4, ncol = 4)
  rownames(Akt)=colnames(Akt) = states
  for(i in 1:(length(obs)-1)){
    for(j in 1:length(states)){
      for(j2 in 1:length(states)){
        #print(exp(fwd[j,i]+tp[j,j2]+ep[j2,obs[i+1]]+back[j2,i+1]-P_x))
        Akt[j,j2] = sumLogProb(Akt[j,j2],fwd[j,i]+tp[j,j2]+ep[j2,obs[i+1]]+back[j2,i+1]-P_x)
      }
    }
    }
  
  #expected transitions divided by P(x|theta)
  Akt = exp(Akt)-1#due to initialization at 0
  print(sum(Akt))
  a_kt = Akt/rowSums(Akt)
  #print(a_kt)
  
  #mu = mean(c(a_kt['h','l'],a_kt['l','h']))
  #a_kt[1,1] = 1-3*s, a_kt[1,2:4] = s
  s = (1/3-a_kt[1,1]/3 + a_kt[1,2] + a_kt[1,3] + a_kt[1,4])/4
  u = mean(a_kt[2:4,1])
  v = (a_kt[2,3]+a_kt[2,4]+a_kt[3,2]+a_kt[3,4]+a_kt[4,2]+a_kt[4,3])/6
  #v = a_kt[2,3]
  #print(c(s,u,v))
  #initialize Ekb
  Ekb = matrix(0,nrow = 4,ncol = 8)
  rownames(Ekb) = states
  colnames(Ekb) = c('0000','0001','0011','0111',
                    '0110','0101','0100','0010')
  rownames(fwd) = rownames(back) = states
  
  for(i in 1:length(obs)){
    #Ekb[states[1],obs[i]] = sumLogProb(Ekb[states[1],obs[i]],fwd[states[1],i]+back[states[1],i]-P_x)
    for(st in states){
      Ekb[st,obs[i]] = sumLogProb(Ekb[st,obs[i]],fwd[st,i]+back[st,i]-P_x)
    }
  }
  
  Ekb = exp(Ekb)-1 #due to initializatoin at 0
  print(Ekb)
  print(sum(Ekb))
  ekb = Ekb/rowSums(Ekb)
  #print(ekb)
  #print((ekb[1,1]/ekb[1,7] + ekb[1,2]/ekb[1,6] + 
          #ekb[1,3]/ekb[1,4] + ekb[1,8]/ekb[1,5])/4)
  #branch length t matrix 
  #t = decode_aba2b2(ekb)
  #print(t)
  #correct t with length of a+b+c
  #ratio = mean(timescale/t[,3])
  #t = t*ratio
  #t[,3] = mean(t[,3])
  #print(t)
  #t[2,1] = t[3,2] = mean(c(t[2,1],t[3,2]))
  ekb[2,1] = ekb[3,1] = mean(c(ekb[2,1],ekb[3,1]))
  ekb[2,2] = ekb[3,2] = mean(c(ekb[2,2], ekb[3,2]))
  ekb[2,3] = ekb[3,6] = mean(c(ekb[2,3],ekb[3,6]))
  ekb[2,4] = ekb[3,4] = mean(c(ekb[2,4], ekb[3,4]))
  ekb[2,5] = ekb[3,5] = mean(c(ekb[2,5],ekb[3,5]))
  ekb[2,6] = ekb[3,3] = mean(c(ekb[2,6],ekb[3,3]))
  ekb[2,7] = ekb[3,8] = mean(c(ekb[2,7],ekb[3,8]))
  ekb[2,8] = ekb[3,7] = mean(c(ekb[2,8],ekb[3,7]))
  ekb[4,3] = ekb[4,6] = mean(c(ekb[4,3],ekb[4,6]))
  ekb[4,7] = ekb[4,8] = mean(c(ekb[4,7],ekb[4,8]))
  
  ekb = ekb/rowSums(ekb)
  #print(ekb)
  #t[2,2] = t[3,1] = t[4,1] = t[4,2] = mean(c(t[2,2],t[3,1],t[4,1],t[4,2]))
  
  #ratio = mean(timescale/t[,3])
  t = decode_aba2b2(ekb)
  #print(t)
  a = t[1,1]/2
  b = t[1,2]/2-a
  a2 = (t[2,1]/2 + t[3,2]/2)/2
  #corrected with 3(a+b)-a2
  #r = (t[2,2]+t[3,1]+t[4,1]+t[4,2])/4/(3*a+3*b-a2)
  #a = r*a
  #b = r*b
  #a2 = r*a2
  #print(c(s,u,v,a,b,a2))
  #corrected t branch length with means of a,b,a2,b2 above
  #t = get_t(a,b,a2)

  #print(t)
  return(list(tp = trans_P(states,s,u,v),ep = emiss_P(states,get_t(a,b,a2)),s = s,u = u,v=v,a = a,b = b,a2 =a2))
}

# Helper function to save plot of log likelihoods over iterations to file for
#    visualization.
# Arguments:
#    log_likelihoods: vector of log likelihoods over iterations
#    init_mu, init_theta_h, init_theta_l: the initial values of parameters used
#        (for naming the file containing the plot)
# Outputs:
#    plot of log likelihoods to file
saveplot <- function(log_likelihoods, s,u,v,a,b,a2) {
  png(sprintf("em_%.2f_%.2f_%.2f%.2f_%.2f_%.2f.png", s, u, v,a,b,a2),
      height=480, width=640)
  plot(seq(1,length(log_likelihoods)), log_likelihoods, type='l', col='red',
       main=sprintf("EM log likelihoods with initialization %.2f, %.2f, %.2f", s, u, v),
       xlab="Iteration", ylab="Log likelihood")
  dev.off()
}

saveplot2 = function(posterior,s,u,v,a,b,a2){
  png(sprintf("em_%.2f_%.2f_%.2f%.2f_%.2f_%.2f.png", s, u, v,a,b,a2),
      height=480, width=640)
  matplot(t(posterior),type = 'l',col = rainbow(ncol(posterior)),lend = states,
          main = sprintf("Posterior likelihoods with %.2f, %.2f, %.2f, %.2f, %.2f, %.2f", s, u, v,a,b,a2))
  dev.off()
}

# Uses EM to infer the parameters ```s,u,v,a,b,a2,b2```, iterating until
#    a valid stopping condition is reached.
# Arguments:
#    sequence: sequence data to train on
#    s,u,v: the value used for initializing the transition probabilities
#    theta_h, theta_l: parameters of the emission probability distributions
# Returns:
#    mu: parameter of trained transition probability distribution
#    theta_h, theta_l: parameters of trained emission probability distribution
train <- function(sequence,s,u,v,a,b,a2,stop_diff = 0.0001) {
  init_s = s
  init_u = u
  init_v = v
  init_a = a
  init_b = b
  init_a2 = a2
  
  t = get_t(a,b,a2)
  trans_probs <- trans_P(states,s,u,v)
  emiss_probs <- emiss_P(states,t)
 
  init_probs <- log(c(0.7,0.1,0.1,0.1))
  names(init_probs) <- states
  
  log_likelihoods = c() # vector of log likelihoods from each iteration
  
  flag = 1
  log_likelihoods = c(log_likelihoods,-200000)#initialize with a impossible value
  while(flag){
    fb_output = forward_backward(sequence, trans_probs, emiss_probs, init_probs)
    log_likelihoods = c(log_likelihoods,fb_output$likelihood_f)
    para_t = em(fb_output, sequence, trans_probs, emiss_probs)

    trans_probs = para_t$tp
    emiss_probs = para_t$ep

    n = length(log_likelihoods)
    if(n == 500){flag = 0}

    
    if(log_likelihoods[n]-log_likelihoods[n-1] < stop_diff){
      flag = 0
    }
    
  }
  

  s = para_t$s
  u = para_t$u
  v = para_t$v

  
  #t[1,1] = 2*a
  a = para_t$a
  #t[1,2] = 2*a + 2*b
  b = para_t$b
  #t[2,1] = 2*a2
  a2 = para_t$a2
  #t[2,2] = 2*(a2+b2)
  b2 = para_t$b2
  posterior = forward_backward(sequence, trans_probs, emiss_probs, init_probs)$posterior
  #print(posterior)

  saveplot(log_likelihoods[-1],init_s,init_u,init_v,init_a,init_b,init_a2)
  #saveplot2(posterior,s,u,v,a,b,a2)
  return(c(s,u,v,a,b,a2,b2))
}


main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 7) {
    cat("usage: Rscript 1a.R hmm-sequence.fa s u v a b a2\n")
    stop()
  }
  
  fasta_file <- args[1]
  
  s <- as.numeric(args[2])
  u <- as.numeric(args[3])
  v <- as.numeric(args[4])
  a <- as.numeric(args[5])
  b <- as.numeric(args[6])
  a2 <- as.numeric(args[7])

  aligned_sequence <- read.fasta(file = fasta_file,as.string = TRUE,seqtype = 'DNA')
  #gap_del = gap_filter(aligned_sequence)
  obs_sequence = seqto01(aligned_sequence)
  
  count = rep(0,8)
  names(count) = c('0000','0001','0011','0111', '0110','0101','0100','0010')
  for(x in obs_sequence){
    count[x] = count[x]+1
  }
  print(count)
  print(sum(count))
  
  
  
  em_results <- train(obs_sequence,s,u,v,a,b,a2)
  
  s <- em_results[1]
  u <- em_results[2]
  v <- em_results[3]
  a <- em_results[4]
  b <- em_results[5]
  a2<- em_results[6]

  
  cat(sprintf("s: %.5f\nu: %.5f\nv: %.5f\na: %.5f\nb: %.5f\na2: %.5f\n",
              s,u,v,a,b,a2))
  
}

main()

