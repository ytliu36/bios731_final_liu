#!/bin/env bash
#SBATCH --account=yli2635
#SBATCH --job-name=mem_1k1k
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yutong.liu@emory.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=10gb  # Memory for each job
#SBATCH --time=6:00:00  # Time for each job
#SBATCH --output=mem_%A_%a.out  # Output file for each array task
#SBATCH --partition=su_lab
#SBATCH --array=1-1000  # This specifies the job array, here for 100 jobs

# activate conda environment
source /sulab/users/yli2635/miniconda3/etc/profile.d/conda.sh
conda activate memento_test

# Set directories
simulation_dir="/sulab/users/yli2635/memento/data/sim1k_null_400/"
py_script="/sulab/users/yli2635/memento/data/mem_1k1k.py"

# Get all unique base names (without _Exp.txt part)
base_names=($(ls "${simulation_dir}" | sed -n 's/_Exp.txt//p' | sort | uniq))

# Use the SLURM_ARRAY_TASK_ID to select the base name for this job
base_name="${base_names[$SLURM_ARRAY_TASK_ID-1]}"  # SLURM_ARRAY_TASK_ID starts at 1, so subtract 1 for zero-indexed arrays

# Run R script
python ${py_script} --base ${base_name} --input ${simulation_dir}
