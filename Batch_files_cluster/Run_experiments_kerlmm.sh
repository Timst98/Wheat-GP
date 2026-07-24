#!/bin/bash
#SBATCH --account=gratis
#SBATCH --mail-type=end,fail
#SBATCH --job-name="ker1prodLMM"
#SBATCH --time=00:21:00
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --array=1-30

export KERNEL=1

export NLOPTR=0
export ADDONLY=0
export PRODONLY=1

export GA_KERNEL=0
export LONG_ADAM=0
export GS_PARAMS=0

module load R/4.4.2-gfbf-2024a

echo "Starting batch job on $(date)"  


export MODEL='kerlmm'
echo "Running kerlmm"
Rscript ./GP_modeling/LMM_BGLR.R


echo "Batch job completed on $(date)"  
