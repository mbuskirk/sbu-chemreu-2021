#!/bin/bash
#
#SBATCH --job-name=Lab4starter-HF
#SBATCH --output=Lab4starter-HF.txt
#SBATCH --ntasks-per-node=28
#SBATCH --nodes=1
#SBATCH --time=60:00
#SBATCH -p short-28core

# This is a batch script, meaning that the below linux 
# commands will be run in order.  # symbols indicate comments.
# however, comments that begin with #SBATCH will be read and 
# processed by the SLURM job scheduler to process the job.  A 
# useful tutorial for using SLURM can be found at:
# https://it.stonybrook.edu/help/kb/using-the-slurm-workload-manager  
 
# As you normally would, load the psifour module 
module load psifour

# By default, this script will begin in your home directory,
# so you need to cd into the directory that contains your job.
# Change the below line so that it will cd to the director
# that contains your Lab4starter-HF.py script.
cd /gpfs/home/bglevine/Lab4/

# Run Lab4starter-HF.py
python Lab4starter-HF.py

# Once the script is done, the job will end. 
