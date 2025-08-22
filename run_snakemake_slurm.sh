#!/bin/bash
#SBATCH --job-name=e2g_gwas_pops_benchmark
#SBATCH --time=2:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/e2g_gwas_pops_benchmark.out
#SBATCH --error=logs/e2g_gwas_pops_benchmark.err
#SBATCH --partition=deyk,cpushort,cpu

# activate snakemake conda environment
source $HOME/.bashrc
conda activate snakemake

# # run snakemake
# snakemake --use-conda --jobs 32 --slurm --default-resources slurm_account=guruvak slurm_partition=deyk,cpushort,cpu

snakemake --use-conda all
