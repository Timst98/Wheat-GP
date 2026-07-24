#!/bin/bash
#SBATCH --mail-type=end,fail
#SBATCH --job-name="1"
#SBATCH --time=15:00:00
#SBATCH --mem=30G
#SBATCH --ntasks=16  # Number of workers
#SBATCH --cpus-per-task=1
#SBATCH --array=1-30

export KERNEL=1

export NLOPTR=0
export ADDONLY=0
export PRODONLY=0

export GA_KERNEL=0
export LONG_ADAM=0
export GS_PARAMS=0

module load R/4.4.2-gfbf-2024a

echo "Starting batch job on $(date)"  
echo "Running short Adam"

Rscript ./GP_modeling/GP.R

export MODEL=1
echo "Running BGLR 1"
Rscript ./GP_modeling/LMM_BGLR.R

export MODEL=2

echo "Running BGLR 2"
Rscript ./GP_modeling/LMM_BGLR.R

export MODEL=3

echo "Running BGLR 3"
Rscript ./GP_modeling/LMM_BGLR.R

export LONG_ADAM=1

echo "Running long  Adam"  

Rscript ./GP_modeling/GP.R

export MODEL=1


echo "Running BGLR 1"
Rscript ./GP_modeling/LMM_BGLR.R

export MODEL=2

echo "Running BGLR 2"
Rscript ./GP_modeling/LMM_BGLR.R

export MODEL=3

echo "Running BGLR 3"
Rscript ./GP_modeling/LMM_BGLR.R

export NLOPTR=1
export LONG_ADAM=0

echo "Running LBFGS"  

Rscript ./GP_modeling/GP.R

export MODEL=1

echo "Running BGLR 1"
Rscript ./GP_modeling/LMM_BGLR.R

export MODEL=2

echo "Running BGLR 2"
Rscript ./GP_modeling/LMM_BGLR.R

export MODEL=3

echo "Running BGLR 3"
Rscript ./GP_modeling/LMM_BGLR.R

export GS_PARAMS=1

echo "Running BGlR  GRID SEARCH"
export MODEL=1

echo "Running BGLR 1"
Rscript ./GP_modeling/LMM_BGLR.R

export MODEL=2

echo "Running BGLR 2"
Rscript ./GP_modeling/LMM_BGLR.R

export MODEL=3

echo "Running BGLR 3"
Rscript ./GP_modeling/LMM_BGLR.R

export GS_PARAMS=0

export ADDONLY=1
export NLOPTR=0

echo "Running short Adam"

Rscript ./GP_modeling/GP.R

export MODEL=1
echo "Running BGLR 1"
Rscript ./GP_modeling/LMM_BGLR.R

export MODEL=2

echo "Running BGLR 2"
Rscript ./GP_modeling/LMM_BGLR.R

export MODEL=3

echo "Running BGLR 3"
Rscript ./GP_modeling/LMM_BGLR.R

export LONG_ADAM=1

echo "Running long  Adam"  

Rscript ./GP_modeling/GP.R

export MODEL=1


echo "Running BGLR 1"
Rscript ./GP_modeling/LMM_BGLR.R

export MODEL=2

echo "Running BGLR 2"
Rscript ./GP_modeling/LMM_BGLR.R

export MODEL=3

echo "Running BGLR 3"
Rscript ./GP_modeling/LMM_BGLR.R

export NLOPTR=1
export LONG_ADAM=0


echo "Running LBFGS"  

Rscript ./GP_modeling/GP.R

export MODEL=1

echo "Running BGLR 1"
Rscript ./GP_modeling/LMM_BGLR.R

export MODEL=2

echo "Running BGLR 2"
Rscript ./GP_modeling/LMM_BGLR.R

export MODEL=3

echo "Running BGLR 3"
Rscript ./GP_modeling/LMM_BGLR.R

export GS_PARAMS=1

echo "Running BGlR  GRID SEARCH"
export MODEL=1

echo "Running BGLR 1"
Rscript ./GP_modeling/LMM_BGLR.R

export MODEL=2

echo "Running BGLR 2"
Rscript ./GP_modeling/LMM_BGLR.R

export MODEL=3

echo "Running BGLR 3"
Rscript ./GP_modeling/LMM_BGLR.R

export GS_PARAMS=0
echo "Batch job completed on $(date)"  
