#!/bin/bash
#
#SBATCH -J testjob
#SBATCH -A liu-compute-2025-38
#SBATCH --reservation devel
#SBATCH -t 00:05:00
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --exclusive
#
export OPENBLAS_NUM_THREADS=32
export MKL_NUM_THREADS=32
export NUMEXPR_NUM_THREADS=32
export OMP_NUM_THREADS=32

source /proj/liu-compute-2025-38/software/init.sh
time mpprun python3 <program to run>

echo "job completed"
