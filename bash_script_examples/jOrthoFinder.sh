#!/bin/bash
#Indicates this is a shell script
#SBATCH --job-name=jOrthoFinder.sh #JobName #SBATCH --mail-type=END,FAIL
#Mail events (NONE, BEGIN, END, FAIL, ALL) - When an email will be sent to the following address
#SBATCH --mail-user=<yw867@cornell.edu>
#SBATCH -p RM #Submission queue in Regular Memory. You don't need Large Memory for everything we will do in this class
#SBATCH --ntasks=1
# Unless you use MPI computing, leave this as 1. Almost always 1. #SBATCH --cpus-per-task=<XX> #How many threads/processors
#SBATCH --mem=<Xgb> #Memory you need
#SBATCH --time=<HH:MM:SS> #Time you want to reserve
#SBATCH --output=<jobname>.joblog pwd;hostname;date
<Commands>
date
