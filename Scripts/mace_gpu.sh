#!/bin/bash
#
#SBATCH -J testjob_gpu
#SBATCH -A liu-gpu-2025-13
#SBATCH --reservation devel
#SBATCH -t 00:05:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --exclusive
#SBATCH --gpus-per-task=v100:1
#
export NSC_MODULE_SILENT=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

source /proj/liu-compute-2025-38/software/init_gpu.sh
time  python3 <program to run>

echo "job completed"
