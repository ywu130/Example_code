#!/bin/bash
#SBATCH --job-name=jPair1.sh # Job name
#SBATCH --mail-type=BEGIN,END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL) 
#SBATCH --mail-user=yw867@cornell.edu # Where to send mail
#SBATCH --ntasks=1 # Run a single task
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH --mem=20gb # Job memory request
#SBATCH --time=3:30:00 # Time limit hrs:min:sec
#SBATCH --output=jPair1_piped.log # Standard output and error log
pwd; hostname; date
module load hisat2
module load samtools
hisat2 --summary-file Pair1_align.log --rna-strandness RF -x SL_index -1 SRR2239886_1.fastq -2 SRR2239886_2.fastq | tee >(samtools flagstat - > Pair1.flagstat) | samtools sort -O BAM | tee Pair1.bam | samtools index - Pair1.bam.bai
date
