#!/bin/bash

#SBATCH --job-name=Alisim
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --account=gely021262
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --partition=Innovation

#directory
job_dir=$( pwd )
cd $job_dir

source ~/initMamba.sh
conda activate Physalia

mpirun  -d FcC_supermatrix.fas -cat -gtr Mollusca_chain1