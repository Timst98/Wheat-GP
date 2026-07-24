#!/bin/bash
#SBATCH --account=paygo
#SBATCH --wckey=imsv_biggp
#SBATCH --mail-type=end,fail
#SBATCH --job-name="GP with noise"
#SBATCH --time=60:00:00
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
echo "Running Adam kernel 1"

Rscript ./GP_modeling/GP.R

echo "Running Adam kernel 2"
export KERNEL=2

Rscript ./GP_modeling/GP.R

echo "Running Adam kernel 4"
export KERNEL=4

Rscript ./GP_modeling/GP.R

echo "Running Adam kernel 5"
export KERNEL=5

Rscript ./GP_modeling/GP.R


echo "Running product only setting"

export ADDONLY=0
export PRODONLY=1

Rscript ./GP_modeling/GP.R


echo "Running additive only setting"

export ADDONLY=1
export PRODONLY=0

Rscript ./GP_modeling/GP.R




echo "Batch job completed on $(date)"  
