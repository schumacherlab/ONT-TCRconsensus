#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=275GB
#SBATCH --time=1-00:00:00
#SBATCH --output=./logs/run_nanpore_tcr_consensus_%A.out
#SBATCH --error=./logs/run_nanopore_tcr_consensus_%A.err

mkdir -p ./logs

source ~/miniconda3/etc/profile.d/conda.sh
conda activate ont_tcr_consensus_env

tcr_consensus $1
