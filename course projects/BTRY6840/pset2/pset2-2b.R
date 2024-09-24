#!/usr/bin/env Rscript

# Script for computing sequence alignments using Needleman-Wunsch with
#   linear gap penalties.
# Arguments:
#    f - FASTA file with sequences in FASTA format.
#    s - JSON with the score matrix for alignment.
#    d - The gap penalty for the alignment.
#
# Outputs:
#    Prints alignment to console.
#
# Example usage with d = 100:
#    Rscript 2b.R sequences.fasta 100


# Computes the actual string alignments given the traceback matrix.
# Arguments:
#    x: the first string we're aligning
#    y: the second string we're aligning
#    t: the traceback matrix
# Returns:
#    a: two element vector with first and second aligned sequences as strings
traceback <- function(x, y, t) {
  # COMPLETE THIS FUNCTION
  a = c('','')
  m = nchar(x)+1
  n = nchar(y)+1
  while(!is.na(t[m,n])){
    if(t[m,n]=='g'){
      a[1] = paste(substr(x,m-1,m-1),a[1],sep='')
      a[2] = paste(substr(y,n-1,n-1),a[2],sep='')
      m = m-1
      n = n-1
    } else if(t[m,n]=='u'){
      a[1] =  paste(substr(x,m-1,m-1),a[1],sep='')
      a[2] = paste('-',a[2],sep='')
      m = m-1

    } else if(t[m,n]=='l'){
      a[1] = paste('-',a[1],sep='')
      a[2] = paste(substr(y,n-1,n-1),a[2],sep='')
      n = n-1

    } 
    
  }
  return(a)
}


# Computes the score and alignment of two strings.
# Arguments:
#    x: the first string we're aligning
#    y: the second string we're aligning
#    s: the score matrix
#    d: the gap opening/extension penalty
# Returns a vector with two elements:
#    score: the score of the optimal sequence alignment
#    a: two element vector with first and second aligned sequences as strings
# The latter two are computed using the above traceback method.
sequence_alignment <- function(x, y, s, d) {
  m <- matrix(nrow = nchar(x)+1, ncol = nchar(y)+1) # Recurrence matrix
  t <- matrix(nrow = nchar(x)+1, ncol=nchar(y)+1) # Traceback matrix
  
  m[1,1] = 0 #base case
  for(i in 2:(nchar(x)+1)){
    m[i,1] = m[i-1,1] - d
    t[i,1] = 'u'
  }
  for(j in 2:(nchar(y)+1)){
    m[1,j] = m[1,j-1] - d
    t[1,j] = 'l'
  }
  
  for(i in 2:(nchar(x)+1)){
    for(j in 2:(nchar(y)+1)){
      vg = m[i-1,j-1]+s[substr(x,i-1,i-1),substr(y,j-1,j-1)] 
      vu = m[i-1,j]-d
      vl = m[i,j-1]-d

      m[i,j] = max(c(vg,vu,vl))
      if(m[i,j]==vg){
        t[i,j]= 'g'
      }else if(m[i,j]==vu) {
        t[i,j]='u'
        }else if(m[i,j]==vl){
          t[i,j]='l'
          }
    }
  }
  a <- traceback(x, y, t)
  # COMPLETE THIS FUNCTION
  return(c((m[nchar(x)+1,nchar(y)+1]),a[1],a[2])) #has to split 'a' in this way, do not know how to
  #contain 2 elements in one element to return
  
}


# Prints two aligned sequences formatted for convenient inspection.
# Arguments:
#    a: two element vector with first and second aligned sequences as strings
# Outputs:
#    Prints aligned sequences (80 characters per line) to console
print_alignment <- function(a) {
  a_x <- a[1]
  a_y <- a[2]

  if(nchar(a_x)!=nchar(a_y)) {
    stop("Sequence alignment lengths must be the same.")
  }
  nlines <- as.integer(1 + nchar(a_x) / 80)
  for(i in 1:nlines) {
    start <- (i - 1) * 80 + 1
    end <- i * 80
    cat(paste0( substr(a_x, start, end), "\n" ))
    cat(paste0( substr(a_y, start, end), "\n" ))
    cat("\n")
  }
}

# Reads in a fasta format file that is assumed to contain two sequences, with
# the complete sequence on single lines following the sequence labels
read_fasta <- function(fasta_file) {
  con <- file(fasta_file, "r")
  x_name <- readLines(con, n = 1)
  x <- readLines(con, n = 1)
  y_name <- readLines(con, n = 1)
  y <- readLines(con, n = 1)
  close(con)
  
  return(c(x,y))
}

# Generates the score matrix that can be indexed using the characters
# "A","C","G","T"
make_score_matrix <- function() {
  mat <- matrix(runif(16), nrow=4, ncol=4)
  rownames(mat) <- c("A", "C", "G", "T")
  colnames(mat) <- rownames(mat)
  
  mat["A", "A"] <- mat["T", "T"] <- 91
  mat["C", "C"] <- mat["G", "G"] <- 100
  mat["A", "C"] <- mat["C", "A"] <- mat["G", "T"] <- mat["T", "G"] <- -114
  mat["A", "G"] <- mat["G", "A"] <- mat["C", "T"] <- mat["T", "C"] <- -31
  mat["A", "T"] <- mat["T", "A"] <- -123
  mat["C", "G"] <- mat["G", "C"] <- -125
  
  return(mat)
}



main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2) {
    cat("usage: Rscript 2b.R file.fasta d\n")
    stop()
  }
  fasta_file <- args[1]
  d <- as.integer(args[2])
  
  seqs <- read_fasta(fasta_file)
  x <- seqs[1]
  y <- seqs[2]
  s <- make_score_matrix()
  
  results <- sequence_alignment(x, y, s, d)
  score <- results[1]
  a = c('','')
  a[1] <- results[2]
  a[2]<- results[3]
  cat("Alignment:\n")
  dontcare <- print_alignment(a)
  cat(sprintf("Score: %s\n", score))
}

main()