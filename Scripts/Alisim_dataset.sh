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

iqtree --alisim LG_C60_mollusca -t mollusca_tree1.tree --length 2000 --num-alignments 30 --seqtype AA -m LG+C60+G -af fasta
