#!/bin/bash
#SBATCH --account=paygo
#SBATCH --wckey=imsv_biggp
#SBATCH --mail-type=end,fail
#SBATCH --job-name="Env only new limits"
#SBATCH --time=10:00:00
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

export KERNEL=5

export NLOPTR=0
export ADDONLY=0
export PRODONLY=0

export GA_KERNEL=0
export LONG_ADAM=1
export GS_PARAMS=0
export SETUP=1
export LEAKAGE=0
module load R/4.4.2-gfbf-2024a

echo "Starting batch job on $(date)"  
echo "Running Adam no,1"

Rscript ./GP_modeling/Env_only/Wheat_envonly.R
export SETUP=2

echo "Running Adam no,2"


Rscript ./GP_modeling/Env_only/Wheat_envonly.R

export LEAKAGE=1
export SETUP=1

echo "Running Adam yes,1"

Rscript ./GP_modeling/Env_only/Wheat_envonly.R


export LEAKAGE=1
export SETUP=2

echo "Running Adam yes,2"

Rscript ./GP_modeling/Env_only/Wheat_envonly.R



echo "Batch job completed on $(date)"  
