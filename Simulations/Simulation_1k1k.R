start<-proc.time()
library(dplyr)
library(truncnorm)
library(MASS)
library(optparse)
option_list = list(
  make_option(c("--sim"), type="integer", default=10,
              help="times of simulation"),
  make_option(c("--mu"), type="double", default=0,
              help="assume a rank of mu level"),
  make_option(c("--rho"), type="double", default=0.3,
              help="assume a rho for total variance"),
  make_option(c("--prop_var"), type="double", default=0.5,
              help="the proportion of ind variance"),
  make_option(c("--prop_cov1"), type="double", default=0.5,
              help="proportion of cov beta0 in ind"),
  make_option(c("--prop_cov2"), type="double", default=0.5,
              help="proportion of cov beta1 in ind"),
  make_option(c("--prop_gene"), type="double", default=0.2,
              help="assume a proportion beta1 to beta0 (proportion of covariance related to gene)"),
  make_option(c("--out"), type="character", default="/Users/yutong/Downloads/",
              help="saving dictionary")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Validate input correlation value
if (opt$rho+2*opt$prop_gene*opt$rho>1) {
  stop("Total correlation larger than 1, please use smaller --rho or --prop_gene.")
}

# Data Simulation Setting
set.seed(22)
n_sim = 1
bar<-read.table("/sulab/users/yli2635/memento/data/OneK1K_CD4/Barcodes.txt", header = TRUE)
n_sub = length(unique(bar$ind_ID))
n_cell=nrow(bar)

## Simulate g_is to run more efficiently
g_i.mat<-data.frame(matrix(rbinom(n_sub*n_sim, 2, 0.5), nrow = n_sub, ncol = n_sim))
rownames(g_i.mat)<-unique(bar$ind_ID)
Geno<-apply(g_i.mat, 2, sort)

# reorder barcode and further X to match the reordered genotypes
bar$ind_ID <- factor(bar$ind_ID, levels = rownames(Geno))
bar<- bar[order(bar$ind_ID), ]
rownames(bar) <- NULL
n_cell_ind = table(bar$ind_ID)
print(n_cell==sum(n_cell_ind))

## Simulate s_ics with truncated normal(takes time 28.758)
s_ic.mat<-read.table("/sulab/users/yli2635/memento/data/OneK1K_CD4/Seq_depth.txt", row.names = 1)
# Ensure row names of s_ic.mat match cell_ID in reordered bar
s_ic.mat<- s_ic.mat[match(bar$cell_ID, rownames(s_ic.mat)), , drop = FALSE]
print(nrow(s_ic.mat)==n_cell)

## values can be changed
rho <- opt$rho # assume a rho for total variance
prop_var<-opt$prop_var # the proportion of ind variance
prop_cov<-c(opt$prop_cov1,opt$prop_cov2) # proportion of cov betas in ind
prop_gene <-opt$prop_gene # assume a proportion beta1 to beta0 (proportion of variance related to gene)

## values based on data
theta<-c(7.24, 9.81) ### could be tunned -- highly expressed genes are highly overdispersed
if (opt$mu ==100){
  mu_ij<-c(0.001745296, 0.0009766125)
} else if (opt$mu ==200){
  mu_ij<-c(0.0003340761,0.0003072830)
} else if (opt$mu ==400){
  mu_ij<-c(0.0001527515,0.0001480471)
} else if (opt$mu ==800){
  mu_ij<-c(8.410669e-05,8.288212e-05)
} else if (opt$mu ==1600){
  mu_ij<-c(4.350453e-05,4.321084e-05)
} else if (opt$mu ==3200){
  mu_ij<-c(2.080791e-05,2.067313e-05)
} else {
  mu_ij<-c(0.008078490, 0.007963918)
} # set mu_ij based on rank fo mu, 100,200,400,800,1600,3200

## values processed
sigma_ij<- mu_ij^2 /theta
cov_max<-sqrt(sigma_ij[1]*sigma_ij[2])# maximum covariance
beta<-c(rho*cov_max,rho*prop_gene*cov_max)
var_ind<-prop_var*sigma_ij
var_cell<-(1-prop_var)*sigma_ij
# formula for covariance
# cov_ij<-beta[0]+beta[1]*g_i
shape_1 <-mu_ij^2/var_ind
scale_1 <-var_ind/ mu_ij

####################### Set simulation results storage #####################
z_ij.mat<-data.frame(matrix(NA, nrow=2*n_sub,ncol = n_sim))
z_icj.mat<-data.frame(matrix(NA, nrow=2*n_cell,ncol = n_sim))

############################ Simulation ###################################
# set seed to ensure diffrent simulation result
seeds <- sample(1:1000000,opt$sim, replace = F)
set.seed(seeds[opt$sim])

for (i in 1:n_sim){
  counts<-table(factor(g_i.mat[,i], levels = c(0, 1, 2)))
  cov_ind<-prop_cov[1]*beta[1]+prop_cov[2]*beta[2] *c(0,1,2)
  cov_cell<-(1-prop_cov[1])*beta[1]+(1-prop_cov[2])*beta[2] *c(0,1,2)
  # Please make sure rho<1
  rho_ind<-cov_ind/prop_var/cov_max
  rho_cell<-cov_cell/(1-prop_var)/cov_max
  ## Simulate z_ij
  u_ij_0<-if (counts["0"] == 0) c() else mvrnorm(counts["0"], mu = c(0, 0), Sigma = matrix(c(1, rho_ind[1], rho_ind[1], 1), nrow = 2))
  u_ij_1<-if (counts["1"] == 0) c() else mvrnorm(counts["1"], mu = c(0, 0), Sigma = matrix(c(1, rho_ind[2], rho_ind[2], 1), nrow = 2))
  u_ij_2<-if (counts["2"] == 0) c() else mvrnorm(counts["2"], mu = c(0, 0), Sigma = matrix(c(1, rho_ind[3], rho_ind[3], 1), nrow = 2))
  u_ij<-rbind(u_ij_0,u_ij_1,u_ij_2)
  z_i1<-qgamma(pnorm(u_ij[,1]), shape = shape_1[1], scale = scale_1[1])
  z_i2<-qgamma(pnorm(u_ij[,2]), shape = shape_1[2], scale = scale_1[2])
  z_ij.mat[,i]<-c(z_i1,z_i2)
  ## Simulate z_icj
  u_icj_0<-if (counts["0"] == 0) c() else mvrnorm(sum(n_cell_ind[1:counts["0"]]), mu = c(0, 0), Sigma = matrix(c(1, rho_cell[1], rho_cell[1], 1), nrow = 2))
  u_icj_1<-if (counts["1"] == 0) c() else mvrnorm(sum(n_cell_ind[(counts["0"]+1):(counts["0"]+counts["1"])]), mu = c(0, 0), Sigma = matrix(c(1, rho_cell[2], rho_cell[2], 1), nrow = 2))
  u_icj_2<-if (counts["2"] == 0) c() else mvrnorm(sum(n_cell_ind[(counts["0"]+counts["1"]+1):n_sub]), mu = c(0, 0), Sigma = matrix(c(1, rho_cell[3], rho_cell[3], 1), nrow = 2))
  u_icj<-rbind(u_icj_0,u_icj_1,u_icj_2)
  shape_2_1 <-z_i1^2/var_cell[1]
  shape_2_2 <-z_i2^2/var_cell[2]
  scale_2_1 <-var_cell[1]/z_i1
  scale_2_2 <-var_cell[2]/z_i2
  z_ic1<-qgamma(pnorm(u_icj[,1]), shape = rep(shape_2_1,times = n_cell_ind), scale = rep(scale_2_1,times = n_cell_ind))
  z_ic2<-qgamma(pnorm(u_icj[,2]), shape = rep(shape_2_2,times = n_cell_ind), scale = rep(scale_2_2,times = n_cell_ind))
  z_icj.mat[,i]<-c(z_ic1,z_ic2)
}
## Simulate x_icj
s_ic.calc<-rbind(s_ic.mat,s_ic.mat)
x_icj.mat <- data.frame(matrix(rpois(n = 2*n_cell*n_sim, lambda=as.vector(as.matrix(z_icj.mat*s_ic.calc))),nrow=2*n_cell,ncol = n_sim))

# format optput name
form_rho<-gsub("\\.", "_", as.character(rho))
form_var<-gsub("\\.", "_", as.character(prop_var))
form_gene<-gsub("\\.", "_", as.character(prop_gene))
form_cov<-gsub("\\.", "_", as.character(prop_cov))

output_filename <- sprintf("sim%d_ind%d_cell%d_mu%d_rho%s_var%s_gene%s_cov_%s_%s", opt$sim, n_sub, n_cell, opt$mu, form_rho,form_var, 
                           form_gene, form_cov[1], form_cov[2])

# format output table
# Simulated Data cleaning
X<-cbind(x_icj.mat[1:n_cell,1],x_icj.mat[(n_cell+1):(2*n_cell),1])
colnames(X)<-c("gene1","gene2")
rownames(X)<-bar$cell_ID
write.table(X, file = paste0(opt$out, output_filename, "_Exp.txt"), row.names = T, col.names = T)
# barcode as the ind vs cell matrix is inputted as ind_ID:cell_ID
colnames(Geno)<-"SNP1"
write.table(Geno, file = paste0(opt$out, output_filename, "_Geno.txt"), row.names = T, col.names = T)
# geno as the genotype matrix is inputted as n_ind*n_SNPs HERE I only have 1
s_ic.mat[,1]<-s_ic.mat[,1]+X[,1]+X[,2]
write.table(s_ic.mat, file = paste0(opt$out, output_filename, "_Seq.txt"), row.names = T,col.names = F )

end<-proc.time()
end-start
