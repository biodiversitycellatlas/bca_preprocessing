#!/bin/bash

##################
# slurm settings #
##################
#SBATCH --job-name=fastq-dump
#SBATCH --output=../logs/%x_%A_%a.out       # %x=job-name, %A=jobID, %a=array task ID
#SBATCH --error=../logs/%x_%A_%a.err
#SBATCH --qos=short        
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G                        # per task
#SBATCH --time=05:00:00
#SBATCH --array=0-2                      # <-- adjust “3” to (#accessions−1)

# ───────────────── ENVIRONMENT SETUP ─────────────────
# Load conda environment
# conda activate ncbi_datasets

# ───────────────── ACCESSION LIST ─────────────────
accessions=(
  SRR30604769
  SRR30604770
  SRR30604771
  # SRR30604772
  # SRR30604773
  # SRR30604774
  # SRR30604775
  # SRR30604776
  # SRR30604777
  # SRR30604778
  # SRR30604779
  # SRR30604780
  # SRR30604781
  # SRR30604782
  # SRR30604783
  # SRR30604784
  # SRR30604785
  # SRR30604786
  # SRR30604787
)

# ───────────────── derive THIS TASK’S accession ─────────────────
idx=$SLURM_ARRAY_TASK_ID
acc=${accessions[$idx]}
echo "[$(date +%F\ %T)] Task $idx downloading $acc …"

# ───────────────── DOWNLOAD OPTIONS ─────────────────
outdir="/no_backup/asebe/bvanwaardenburg/external_data/Hs_1046517_oakseq_v2/fastq/"
mkdir -p "$outdir"

fasterq-dump "$acc" \
    --threads "$SLURM_CPUS_PER_TASK" \
    -O "$outdir"

gzip "$outdir/$acc"*.fastq

echo "[$(date +%F\ %T)] Finished $acc"

