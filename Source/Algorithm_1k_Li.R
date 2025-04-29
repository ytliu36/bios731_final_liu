# Algorithm Implementation
library(optparse)

# Define command-line arguments
option_list <- list(
  make_option(c("--input"), type = "character", default = "/projects/chang/yutong/coeQTL/sim_no_cov/",
              help = "Input direction"),
  make_option(c("--base"), type = "character", default = "sim1_ind1000_cell2000_rho0_var0_5_gene0_2_cov_0_5_0_5",
              help = "Simulation parameters"),
  make_option(c("--out"), type = "character", default = "/projects/chang/yutong/coeQTL/sim_no_cov/",
              help = "Output name")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print(opt)

# Validate input arguments
if (is.null(opt$base) || is.null(opt$out) || is.null(opt$input)) {
  stop("Arguments --input --base --out are required.")
}

# Simulated Data Loading
X<-read.table(paste0(opt$input, opt$base , "_Exp.txt"), header = T, row.names = 1)
# X as the sc count matrix is inputted as n_genes*n_cells
Bar<-read.table("/sulab/users/yli2635/memento/data/OneK1K_CD4/Barcodes.txt", header = TRUE)
# barcode as the ind vs cell matrix is inputted as ind_ID:cell_ID
Geno<-read.table(paste0(opt$input, opt$base , "_Geno.txt"), header = T, row.names = 1)
# geno as the genotype matrix is inputted as n_ind*n_SNPs HERE I only have 1
seq_depth <-as.vector(read.table(paste0(opt$input, opt$base , "_Seq.txt"), row.names = 1)[,1])
# seq_depth is a vector with all n_cells sequencing depth, could also be generated afterwards


# if (is.null(seq_depth)) {
#   seq_depth = apply(X, 1, sum, na.rm = T)
# }
# dimension check
if(nrow(X) != length(seq_depth)){
  stop('The length of the sequencing depth must match the number of cells.')
}
if(!identical(rownames(X), Bar$cell_ID)){
  X <- X[rownames(X) %in% Bar$cell_ID, , drop = FALSE]
  Bar <- Bar[Bar$cell_ID %in% rownames(X), , drop = FALSE]
} # match the cells in barcode and X if they are not matched
if(length(unique(Bar$ind_ID))!=nrow(Geno)){
  stop('# of subject of genotype need to numch # of subject in scRNA-seq.')
}
### QUESTION HERE: how to set the input format of data in real data?
start <-proc.time()
n_cell = nrow(Bar) # total number of cells
n_gene = ncol(X)
n_ind = length(unique(Bar[,1]))
# match Geno with Bar$ind_ID
ind_match <- match(unique(Bar$ind_ID), rownames(Geno))
Geno <- Geno[ind_match, ]
# group X from the same individual together and reorder seq_depth
cell_match <- match(Bar$cell_ID, rownames(X))
X <- X[cell_match, ]
seq_depth <- seq_depth[cell_match]
# get a matching list of indexes and order it by indexes
ind_cell<- split(seq_along(cell_match), Bar$ind_ID)
ind_cell <- ind_cell[order(sapply(ind_cell, min))]
n_ind_cell <- sapply(ind_cell, length)

# normalization
X<-log1p(X)
# co-eQTL filtering
z_stat<-c()
for (i in 1:n_ind){
  spearman_corr <- cor(X[ind_cell[[i]],1], X[ind_cell[[i]],2], method = "spearman")
  t<-spearman_corr * sqrt((n_ind_cell[i] - 2) / (1 - spearman_corr^2))
  z<-qnorm(pt(abs(t), df = (n_ind_cell[i] - 2), lower.tail = FALSE)) * sign(-t)
  z_stat<-c(z_stat, z)
}

alpha<-0.05
sig_prop<-sum(abs(z_stat) > qnorm(1-(alpha/2)), na.rm = T)/ length(z_stat)

# Remove NA values from z_stat and corresponding values in Geno
valid_idx <- !is.na(z_stat)  # Identify non-NA indices

# Compute Spearman correlation on filtered values
cor_QTL <- cor(z_stat[valid_idx], Geno[valid_idx], method = "spearman")
t_QTL<-cor_QTL * sqrt((sum(valid_idx) - 2) / (1 - cor_QTL^2))
z_QTL<-qnorm(pt(abs(t_QTL), df = (sum(valid_idx) - 2), lower.tail = FALSE)) * sign(-t_QTL)

end<-proc.time()
end-start


## result check
# estimated value
result <- data.frame(
  Corr = cor_QTL,
  z_stat = z_QTL,
  sig_prop = sig_prop,
  sprase = sum(valid_idx)!=n_ind
)
row.names(result) = sub("_.*", "", opt$base)
write.table(result, file = paste0(opt$out, opt$base, "_Li.txt"), row.names = T, col.names = T)
