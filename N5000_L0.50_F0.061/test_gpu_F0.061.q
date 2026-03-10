#!/bin/bash
#
#SBATCH -J testjob
#SBATCH -A liu-gpu-2025-14
#SBATCH --reservation gpu
#SBATCH -t 10:00:00
#SBATCH -n 1
#SBATCH --gpus-per-task=v100:1
#
export NSC_MODULE_SILENT=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1
source /proj/liu-gpu-2025-14/venv/bin/activate
source /proj/liu-gpu-2025-14/users/x_ishag/httk/init.shell

STEPS=5000
DT=0.50
FRIC=0.061

# --- ISOLATED EXECUTION SETUP ---

# 1. Store the main workspace path so we can securely point back to our inputs
ROOT_DIR=$PWD

# 2. Define the nested directory paths
PARENT_DIR="N${STEPS}_L${DT}"
CHILD_DIR="N${STEPS}_L${DT}_F${FRIC}"
DIR_NAME="${PARENT_DIR}/${CHILD_DIR}"

# 3. Create the unique folder and step inside it BEFORE running the physics
mkdir -p ${DIR_NAME}
cd ${DIR_NAME}

# 4. (Optional but safe) Copy the tracking files into the isolated folder so Python can read them
cp ${ROOT_DIR}/key . 2>/dev/null
cp ${ROOT_DIR}/defect_info . 2>/dev/null

# --- SIMULATION ---

# 5. Execute the simulation
# Because we are inside the nested folder, we use ${ROOT_DIR} to point exactly to where the files live
time python3 ${ROOT_DIR}/Isacks-Co/SourceCode/HTP_F0.061.py ${ROOT_DIR}/MoS2_C_int0.cif /proj/liu-gpu-2025-14/scripts/mace-mpa-0-medium.model

# --- CLEANUP ---

# 6. Pull the specific duplicate scripts and the Slurm log down into this folder
mv ${ROOT_DIR}/Isacks-Co/SourceCode/HTP_F${FRIC}.py . 2>/dev/null
mv ${ROOT_DIR}/test_gpu_F${FRIC}.q . 2>/dev/null
mv ${ROOT_DIR}/slurm-${SLURM_JOB_ID}.out . 2>/dev/null

echo "job completed"
