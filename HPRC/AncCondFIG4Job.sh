#!/bin/bash

##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE                #Do not propagate environment
#SBATCH --get-user-env=L             #Replicate login environment

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=AncCondFig4Job     #Set the job name to "JobExample1"
#SBATCH --time=50:00:00              #Set the wall clock limit to 1hr and 30min 
#SBATCH --nodes=1                    #Request 1 node
#SBATCH --ntasks-per-node=28          #Request 8 tasks/cores per node
#SBATCH --mem=25600M                  #Request 2560MB (2.5GB) per node
#SBATCH --output=AncCondFig4Job.%j      #Send stdout/err to "Example1Out.[jobID]"

##OPTIONAL JOB SPECIFICATIONS
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=nwanderson@tamu.edu    #Send all emails to email_address

#First Executable Line
module load GSL/2.4-GCCcore-6.4.0
module load FFTW/3.3.7-iompi-2017b
ml R_tamu/3.5.0-iomkl-2017b-recommended-mt

Rscript HPRCAncCondTestScriptFIG4.R
#submit command to call script



