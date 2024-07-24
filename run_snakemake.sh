#!/bin/bash

#BSUB -W 10:00                 # for 120 hours of wall clock time
#BSUB -J run_snakemake_e2g_benchmarking          # Job title
#BSUB -o run_snakemake_e2g_benchmarking.out   # Output file name
#BSUB -e run_snakemake_e2g_benchmarking.err   # Error file name

# activate snakemake mamba environment
source $HOME/.bashrc
mamba activate snakemake

# run snakemake for whole pipeline (ending with volcano plot)
snakemake --use-conda --jobs 1 --cluster 'bsub -W 10:00 -n 4 -R "rusage[mem=16G]" -o out.%J.txt -e err.%J.txt' build_e2g_modules

# Random comments pertaining to various components of the job submission.
# .%J adds job ID number to output files
# run using bsub < run_snakemake.sh. '<'
