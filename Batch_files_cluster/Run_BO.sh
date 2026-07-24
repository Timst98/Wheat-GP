#!/bin/bash
#SBATCH --account=paygo
#SBATCH --wckey=imsv_biggp
#SBATCH --mail-type=end,fail
#SBATCH --job-name="BO 20 tp "
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

export KERNEL=5

export NLOPTR=0
export ADDONLY=0
export PRODONLY=0

export GA_KERNEL=0
export LONG_ADAM=1
export GS_PARAMS=0
export train_prop=0.2
module load R/4.4.2-gfbf-2024a

echo "Starting batch job on $(date)"  
echo "Running short Adam"

Rscript ./GP_modeling/BO.R




echo "Batch job completed on $(date)"  
