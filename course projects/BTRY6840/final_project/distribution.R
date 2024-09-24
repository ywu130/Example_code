library(seqinr)
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


saveplot <- function(log_likelihoods) {
  png(sprintf("disttest"),
      height=480, width=640)
  plot(seq(1,length(log_likelihoods)), log_likelihoods, type='l', col='green',
       main=sprintf("Distribution of '0000'among sequnces"),
       xlab="1-1000bp as example", ylab="Y/N '0000")
  dev.off()
}



main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 1) {
      cat("usage: Rscript 1a.R hmm-sequence.fa s u v a b a2\n")
      stop()
    }
    
    fasta_file <- args[1]
    aligned_sequence <- read.fasta(file = fasta_file,as.string = TRUE,seqtype = 'DNA')
    obs_sequence = seqto01(aligned_sequence)
    distr = matrix(0,nrow = 8,ncol = length(obs_sequence))
    rownames(distr) = c('0000','0001','0011','0111',
                        '0110','0101','0100','0010')
    for(i in 1:length(obs_sequence)){
      distr[obs_sequence[i],i] = 1
    }
    print(distr)
    saveplot(distr[1,1:1000])
}

main()
