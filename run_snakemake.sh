#!/bin/bash

#BSUB -W 10:00                 # for 120 hours of wall clock time
#BSUB -J e2g_benchmark          # Job title
#BSUB -o e2g_benchmark.out   # Output file name
#BSUB -e e2g_benchmark.err   # Error file name

# activate snakemake mamba environment
source $HOME/.bashrc
mamba activate snakemake

# run snakemake for whole pipeline (ending with volcano plot)
snakemake --use-conda --jobs 32 --cluster 'bsub -W 6:00 -n 1 -R "rusage[mem=32G]" -o out.%J.txt -e err.%J.txt' all

# Random comments pertaining to various components of the job submission.
# .%J adds job ID number to output files
# run using bsub < run_snakemake.sh. '<'
