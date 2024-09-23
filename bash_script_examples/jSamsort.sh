#!/bin/bash
#SBATCH --job-name=jSamsort.sh # Job name
#SBATCH --mail-type=BEGIN,END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL) 
#SBATCH --mail-user=yw867@cornell.edu # Where to send mail
#SBATCH --ntasks=1 # Run a single task
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH --mem=10gb # Job memory request
#SBATCH --time=3:30:00 # Time limit hrs:min:sec
#SBATCH --output=jSamsort.log # Standard output and error log
pwd; hostname; date
module load samtools
samtools view -Su 4140trimalign.sam | samtools sort -o SRR7774140_sorted.bam
samtools view -Su 4141trimalign.sam | samtools sort -o SRR7774141_sorted.bam
samtools view -Su 4142trimalign.sam | samtools sort -o SRR7774142_sorted.bam
samtools view -Su 4148trimalign.sam | samtools sort -o SRR7774148_sorted.bam
samtools view -Su 4149trimalign.sam | samtools sort -o SRR7774149_sorted.bam
samtools view -Su 4150trimalign.sam | samtools sort -o SRR7774150_sorted.bam
date
