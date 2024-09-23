#!/bin/bash
#SBATCH --job-name=jhisat2.sh # Job name
#SBATCH --mail-type=BEGIN,END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL) 
#SBATCH --mail-user=yw867@cornell.edu # Where to send mail
#SBATCH --ntasks=1 # Run a single task
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH --mem=5gb # Job memory request
#SBATCH --time=3:30:00 # Time limit hrs:min:sec
#SBATCH --output=jhisat2.log # Standard output and error log
pwd; hostname; date
module load hisat2
hisat2 --summary-file alignSRR7774140.log -x at_genome -U SRR7774140trim.fastq -S 4140trimalign.sam
date
