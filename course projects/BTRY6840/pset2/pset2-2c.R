#!/usr/bin/env Rscript

# Script for computing sequence alignments using Needleman-Wunsch with
#   affine gap penalties.
# Arguments:
#    f - FASTA file with sequences in FASTA format.
#    s - JSON with the score matrix for alignment.
#    d - The gap penalty for the alignment.
#    e - The gap extension penalty for the alignment.
#
# Outputs:
#    Prints alignment to console.
#
# Example usage with d = 430, e = 30:
#    Rscript 2c.R sequences.fasta 430 30


# Computes the actual string alignments given the traceback matrix.
# Arguments:
#    x: the first string we're aligning
#    y: the second string we're aligning
#    t: the traceback matrix, which stores values that point to which
#       prior matrix was used to reach a given location in each of the
#       3 matrices.
#    start: value indicating the starting matrix (that had the optimal value)
# Returns:
#    a: two element vector with first and second aligned sequences as strings
traceback <- function(x, y, t, start) {
  # COMPLETE THIS FUNCTION
  a = c('','')
  m = nchar(x)+1
  n = nchar(y)+1
  while(m!=1 | n!=1){
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


# Computes the score and alignment of two strings using an affine gap penalty.
# Arguments:
#    x: the first string we're aligning
#    y: the second string we're aligning
#    s: the score matrix
#    d: the gap opening penalty
#    e: the gap extension penalty
# Returns a vector with two elements:
#    score: the score of the optimal sequence alignment
#    a: two element vector with first and second aligned sequences as strings
# The latter two are computed using the above traceback method.
affine_sequence_alignment <- function(x, y, s, d, e) {
  mrows = nchar(x)+1
  mcols = nchar(y)+1
  m <- matrix(nrow = mrows, ncol=mcols) # Recurrence matrix
  i_x <- matrix(nrow = mrows, ncol=mcols) # Recurrence matrix
  i_y <- matrix(nrow = mrows, ncol=mcols) # Recurrence matrix
  t1 <- matrix(nrow = mrows, ncol=mcols) # Traceback matrix for m
  t2 <- matrix(nrow = mrows, ncol=mcols) # Traceback matrix for i_x
  t3 <- matrix(nrow = mrows, ncol=mcols) # Traceback matrix for i_y
  m[1,1] = i_x[1,1] = i_y[1,1] = 0 #initializations
  #boundaries
  for(i in 2:mrows){ #x as rows
    i_x[i,1] = i_x[i-1,1] -e
    i_y[i,1] = -Inf
    m[i,1] = -Inf
  } 
  for(j in 2:mcols){#y as colomuns
    i_y[1,j] = i_y[1,j-1]-e
    i_x[1,j] = -Inf
    m[1,j] = -Inf
  }
    for(i in 2:mrows){
      for(j in 2:mcols){
        m1 = m[i-1,j-1]+ s[substr(x,i-1,i-1),substr(y,j-1,j-1)]
        m2 = i_x[i-1,j-1] + s[substr(x,i-1,i-1),substr(y,j-1,j-1)]
        m3 = i_y[i-1,j-1] + s[substr(x,i-1,i-1),substr(y,j-1,j-1)]
        m[i,j] = max(m1,m2,m3)
        #traceback matrix 1
        if(m1==m[i,j]){
          t1[i,j] = 'g'
          t1[i-1,j-1] = 'g'
        }else if(m2 ==m[i,j]){
          t1[i,j] = 'g'
          t1[i-1,j-1] = 'u' #u/l???
        }else if(m3==m[i,j]){
          t1[i,j] = 'g'
          t1[i-1,j-1] = 'l'
        }
        
        #i_x
        i_x[i,j] = max(m[i-1,j]-d, i_x[i-1,j]-e)
        #traceback matrix 2
        if(m[i-1,j]-d==i_x[i,j]){
          t2[i,j] = 'u'
          t2[i-1,j] = 'g'
        }else if(i_x[i-1,j]-e ==i_x[i,j]){
          t2[i,j] = 'u'
          t2[i-1,j] = 'u' #u/l???
        }
        
        #i_y
        i_y[i,j] = max(m[i,j-1]-d, i_y[i,j-1]-e)
        #traceback matrix 3
        if(m[i,j-1]-d==i_y[i,j]){
          t2[i,j] = 'l'
          t2[i,j-1] = 'g'
        }else if(i_y[i,j-1]-e ==i_y[i,j]){
          t2[i,j] = 'l'
          t2[i,j-1] = 'l' #u/l???
        }

      }

    }
  
  score = max(m[mrows,mcols],i_x[mrows,mcols],i_y[mrows,mcols]) 

  start <- NULL # Indicator of the starting point (which matrix has optimal value).
  if(m[mrows,mcols]==score){
    #start = c(mrows,mcols)
    t = t1
  } else if(i_x[mrows,mcols]==score){
    t = t2#c(mrows-1,mcols)
  } else if(i_y[mrows,mcols]==score){
    t = t3#c(mrows,mcols-1)
  }
  #print(score)
  a <- traceback(x, y, t, start)
  return(c(score,a[1],a[2]))
  
}


# Prints two aligned sequences formatted for convenient inspection.
# Arguments:
#    a: two element vector with first and second aligned sequences as strings
# Outputs:
#    Prints aligned sequences (80 characters per line) to console
print_alignment <- function(a) {
  a_x <- a[1]
  a_y <- a[2]
  if(nchar(a_x) != nchar(a_y)) {
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
  if (length(args) != 3) {
    cat("usage: Rscript 2b.R file.fasta d e\n")
    stop()
  }
  fasta_file <- args[1]
  d <- as.integer(args[2])
  e <- as.integer(args[3])
  
  seqs <- read_fasta(fasta_file)
  x <- seqs[1]
  y <- seqs[2]
  s <- make_score_matrix()
  
  results <- affine_sequence_alignment(x, y, s, d, e)
  score <- results[1]
  a = c('','')
  a[1]<- results[2]
  a[2]<- results[3]
  cat("Alignment:\n")
  dontcare <- print_alignment(a)
  cat(sprintf("Score: %s\n", score))
}

main()