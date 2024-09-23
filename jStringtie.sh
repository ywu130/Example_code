#!/bin/bash
#SBATCH --job-name=jStringtie.sh # Job name
#SBATCH --mail-type=BEGIN,END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL) 
#SBATCH --mail-user=yw867@cornell.edu # Where to send mail
#SBATCH --ntasks=1 # Run a single task
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH --mem=10gb # Job memory request
#SBATCH --time=3:30:00 # Time limit hrs:min:sec
#SBATCH --output=jStringtie.log # Standard output and error log
pwd; hostname; date
stringtie SRR7774140_sorted.bam -G TAIR10_GFF3_genes.gff -o SRR7774140.gtf -p 10 -A SRR7774140gene.tab
stringtie SRR7774141_sorted.bam -G TAIR10_GFF3_genes.gff -o SRR7774141.gtf -p 10 -A SRR7774141gene.tab
stringtie SRR7774142_sorted.bam -G TAIR10_GFF3_genes.gff -o SRR7774142.gtf -p 10 -A SRR7774142gene.tab
stringtie SRR7774148_sorted.bam -G TAIR10_GFF3_genes.gff -o SRR7774148.gtf -p 10 -A SRR7774148gene.tab
stringtie SRR7774149_sorted.bam -G TAIR10_GFF3_genes.gff -o SRR7774149.gtf -p 10 -A SRR7774149gene.tab
stringtie SRR7774150_sorted.bam -G TAIR10_GFF3_genes.gff -o SRR7774150.gtf -p 10 -A SRR7774150gene.tab
date
