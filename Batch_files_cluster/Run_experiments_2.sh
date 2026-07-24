#!/bin/bash
#SBATCH --account=paygo
#SBATCH --wckey=imsv_biggp
#SBATCH --mail-type=end,fail
#SBATCH --job-name="GP bglr lmm"
#SBATCH --time=40:00:00
#SBATCH --nodes=1
#SBATCH --mem=30G
#SBATCH --ntasks=16  # Number of workers
#SBATCH --cpus-per-task=1
#SBATCH --array=1-30

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export VECLIB_MAXIMUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

export KERNEL=1

export NLOPTR=0
export ADDONLY=0
export PRODONLY=0

export GA_KERNEL=0
export LONG_ADAM=1
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

export MODEL='kerlmm'

echo "Running kerlmm"  

Rscript ./GP_modeling/LMM_BGLR.R

echo "Running Kernel 5"

export KERNEL=5

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

export MODEL='kerlmm'

echo "Running kerlmm"  

Rscript ./GP_modeling/LMM_BGLR.R


echo "Running additive only setting"

export ADDONLY=1
export KERNEL=1

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

echo "Running Kernel 5"

export KERNEL=5

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





echo "Running product only setting"

export ADDONLY=0
export PRODONLY=1
export KERNEL=1

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



echo "Running Kernel 5"

export KERNEL=5

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




echo "Batch job completed on $(date)"  
