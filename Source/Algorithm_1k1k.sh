#!/bin/env bash
#SBATCH --account=yli2635
#SBATCH --job-name=alg_sim_li
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yutong.liu@emory.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=5gb  # Memory for each job
#SBATCH --time=1:00:00  # Time for each job
#SBATCH --output=Alg_%A_%a.out  # Output file for each array task
#SBATCH --partition=su_lab
#SBATCH --array=1-1000  # This specifies the job array, here for 100 jobs

# activate conda environment
source /sulab/users/yli2635/miniconda3/etc/profile.d/conda.sh
conda activate coeQTL

# Set directories
simulation_dir="/sulab/users/yli2635/memento/data/sim1k_beta_100/"
output_dir="/sulab/users/yli2635/memento/data/sim1k_beta_100/"
r_script="/sulab/users/yli2635/memento/data/Algorithm_1k_Li.R"

# Check if simulation_dir is the same as output_dir
if [ "$simulation_dir" != "$output_dir" ]; then
  mkdir -p "${output_dir}"
else
  echo "Output directory is the same as simulation directory, skipping mkdir."
fi

# Get all unique base names (without _Exp.txt part)
base_names=($(ls "${simulation_dir}" | sed -n 's/_Exp.txt//p' | sort | uniq))

# Use the SLURM_ARRAY_TASK_ID to select the base name for this job
base_name="${base_names[$SLURM_ARRAY_TASK_ID-1]}"  # SLURM_ARRAY_TASK_ID starts at 1, so subtract 1 for zero-indexed arrays

# Run R script
Rscript ${r_script} --base ${base_name} --input ${simulation_dir} --out ${output_dir}
