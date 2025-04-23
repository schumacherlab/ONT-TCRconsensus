#!/usr/bin/env bash
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --time=02:00:00
#SBATCH --output=./logs/rsync_fastq_to_nano_tcr_cpu_processing_dir%A.out
#SBATCH --error=./logs/rsync_fastq_to_nano_tcr_cpu_processing_dir%A.err

######### HOW TO USE THIS SCRIPT #########
### 0. FIRST RUN THE run_basecall_pipeline_multi-gpu.sh SCRIPT
### 1. COPY THIS SCRIPT TO THE TARGET DIRECTORY 
### 2. ADAPT WITH CORRECT PROJECT_NAME
### 3. run with sbatch rsync_fastq_to_nano_tcr_cpu_processing_dir.sh
##########################################

PROJECT_NAME=$(basename "$PWD") 
INDIR="./in"
OUTDIR="./out"
CPU_NANO_DIR="/processing/${USER}/nanopore/${PROJECT_NAME}"

# check if a file called barcodes_run.txt exists, else exit with error "barcodes_run.txt does not exist"
if [ ! -f barcodes_run.txt ]; then
  echo "barcodes_run.txt does not exist"
  exit 1
fi
 
# if they don't exist, create the directories
mkdir -p /processing/${USER}/nanopore
mkdir -p $CPU_NANO_DIR
mkdir -p ${CPU_NANO_DIR}/fastq_pass

echo 'START concatenation'
barcodes=$(<"barcodes_run.txt") # Collect all barcodes into a variable

for barcode in $barcodes; do
  mkdir -p "${CPU_NANO_DIR}/fastq_pass/${barcode}"
  rsync -avzP "${OUTDIR}/concatenated_${barcode}.fastq.gz" "${CPU_NANO_DIR}/fastq_pass/${barcode}/concatenated_${barcode}.fastq.gz"
done
