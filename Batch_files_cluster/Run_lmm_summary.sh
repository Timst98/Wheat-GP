#!/bin/bash
#SBATCH --account=paygo
#SBATCH --wckey=imsv_biggp
#SBATCH --mail-type=end,fail
#SBATCH --job-name="BGLR3 summary"
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --ntasks=1 # Number of workers
#SBATCH --cpus-per-task=1
#SBATCH --array=1-30

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export VECLIB_MAXIMUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load R/4.4.2-gfbf-2024a
Rscript ./GP_modeling/Lmm_Summary.R


echo "Batch job completed on $(date)"  
