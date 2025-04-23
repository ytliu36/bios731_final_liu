#!/bin/env bash
#SBATCH --account=yli2635
#SBATCH --job-name=1k1k_sim
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yutong.liu@emory.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=5gb  # Memory for each job
#SBATCH --time=1:00:00  # Time for each job
#SBATCH --output=sim_%A_%a.out
#SBATCH --partition=su_lab
#SBATCH --array=1-1000

# activate conda environment
source /sulab/users/yli2635/miniconda3/etc/profile.d/conda.sh
conda activate coeQTL

# Set working directory
workDir="/sulab/users/yli2635/memento/data/"
outDir="/sulab/users/yli2635/memento/data/sim1k_beta_100/"

# Use the array task ID to specify the sim run
n_sim=${SLURM_ARRAY_TASK_ID}

cd "${workDir}"
mkdir -p "${outDir}"

Rscript Simulation_1k1k.R --sim ${n_sim}  --prop_gene 0.6 --mu 100 --out "${outDir}"
