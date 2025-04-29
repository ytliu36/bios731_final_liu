import scanpy as sc
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import scipy.sparse as sp
from scipy.io import mmread
import itertools
import memento

import statsmodels.formula.api as smf
import statsmodels.api as sm

import os
import pickle as pkl
import anndata
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Process input and output files.")

# Add multiple arguments
parser.add_argument("--input", required=True, help="Input file dir")
parser.add_argument("--base", required=True, help="Input file name")

# Parse arguments
args = parser.parse_args()
sim_path = args.input
base_name = args.base

print(sim_path)
print(base_name)

# Read the raw counts matrix (cells as rows, genes as columns)
count_path = "/sulab/users/yli2635/memento/data/OneK1K_CD4/OneK1K_CD4_count.mtx"
counts = mmread(count_path).tocsr()

# Read barcodes metadata
bar_path = "/sulab/users/yli2635/memento/data/OneK1K_CD4/Barcodes.txt"
barcodes = pd.read_csv(bar_path,  sep=" ", index_col=1)

geno_path = sim_path + base_name+ "_Geno.txt"
genotypes = pd.read_csv(geno_path, sep=" ", header=0)

corr_dat = pd.read_csv(sim_path+base_name +"_Exp.txt", sep=" ", header=0)
corr_dat = corr_dat.loc[barcodes.index]

corr_dat_sparse = sp.csr_matrix(corr_dat) 
counts_appended = sp.hstack([corr_dat_sparse, counts], format='csr')

obs = barcodes.merge(genotypes, left_on="ind_ID", right_index=True, how="left")
var = pd.DataFrame(index=[f"Gene_{i}" for i in range(counts_appended.shape[1])])
adata = anndata.AnnData(X=counts_appended, obs=obs, var = var)

# Show basic information
print(adata)
type(adata.X) == sp.csr_matrix
print(adata.var)

adata.obs['capture_rate'] = 0.07
memento.setup_memento(adata, q_column='capture_rate')
memento.create_groups(adata, label_columns=['SNP1', 'ind_ID'])
memento.compute_1d_moments(adata,min_perc_group=.1)

sample_meta = memento.get_groups(adata)
sample_meta['ind_ID'] = sample_meta['ind_ID'].astype('category')
treatment_df = sample_meta[['SNP1']]
cov_df = pd.get_dummies(sample_meta['ind_ID'].astype('category'))

gene_pairs = list(itertools.product(['Gene_0'], ['Gene_1']))
memento.compute_2d_moments(adata, gene_pairs)
np.random.seed(36)
print(np.random.rand())
memento.ht_2d_moments(
    adata, 
    treatment=treatment_df,
    covariate=cov_df,
    num_boot=5000, 
    verbose=1,
    num_cpus=1)
print(np.random.rand())
result_2d = memento.get_2d_ht_result(adata)
result_2d.index = [base_name.split("_")[0]]
# Save to TXT file while keeping row and column names
result_2d.to_csv(sim_path+base_name +"_Mem.txt", sep=" ", index=True)
