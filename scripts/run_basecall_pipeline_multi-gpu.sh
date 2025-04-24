#!/usr/bin/env bash
#SBATCH --job-name nanopore_basecall
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH --gpus 4
#SBATCH --gpus-per-task 4
#SBATCH --cpus-per-task 16
#SBATCH --partition a100
#SBATCH --time=5-00:00:00
#SBATCH --output=./logs/nano_basecall_pipeline_multi-gpu%A.out
#SBATCH --error=./logs/nano_basecall_pipeline_multi-gpu%A.err

######### HOW TO USE THIS SCRIPT #########
### 0. SETUP:
###   0.1: MAKE THE TARGET DIR WITH A logs/ DIRECTORY and
###   0.2: get correct environment for dorado and nanopore_env
###   0.3: make sure to have the correct barcodes_run.txt file in the directory
### 1. COPY THIS SCRIPT TO THE TARGET DIRECTORY
### 2. ADAPT WITH CORRECT PROJECT_NAME
### 3. RUN WITH `sbatch run_basecall_pipeline_multi-gpu.sh`
##########################################

source ~/miniconda3/etc/profile.d/conda.sh
conda activate nanopore_env

PROJECT_NAME=$(basename "$PWD")
INDIR="./in"
OUTDIR="./out"  
mkdir -p $OUTDIR

if [ ! -f barcodes_run.txt ]; then
  echo "barcodes_run.txt does not exist"
  exit 1
fi

number_of_pods=$(ls -d $INDIR/pod*/ | wc -l)

date
echo $number_of_pods
echo "runnung on $number_of_pods pods, in multi-gpu mode (4 gpus)"

#####basecalling#####
echo 'START basecalling'
for SEQ in $(seq 1 $number_of_pods); do 
  date

  mkdir -p $OUTDIR/bam_$SEQ
  ~/dorado-0.9.0-linux-x64/bin/dorado basecaller sup@v5.0.0 $INDIR/pod$SEQ/ --min-qscore 10 | \
      ~/dorado-0.9.0-linux-x64/bin/dorado demux --kit-name SQK-NBD114-24 -o $OUTDIR/bam_$SEQ/

  echo 'START .fastq.gz conversion'
  for BAM in $OUTDIR/bam_$SEQ/*.bam
  do
    BAM_basename="${BAM%.*}"
    echo $BAM_basename
    echo  "${BAM_basename}.fastq.gz"
    echo  "$BAM"
    samtools fastq "$BAM" | pigz > "${BAM_basename}.fastq.gz"
  done
  rm $OUTDIR/bam_$SEQ/*.bam  # remove all .bam files
  rm -rf $INDIR/pod$SEQ # remove all pod files
done


#####concatenate#####
echo 'START concatenation'
barcodes=$(<"barcodes_run.txt") # Collect all barcodes into a variable

for barcode in $barcodes; do
  echo "$barcode"
  cat ${OUTDIR}/bam_*/*${barcode}.fastq.gz > $OUTDIR/concatenated_$barcode.fastq.gz
  rm ${OUTDIR}/bam_*/*${barcode}.fastq.gz
done

echo 'END'
