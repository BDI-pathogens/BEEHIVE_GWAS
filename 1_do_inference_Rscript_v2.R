args = commandArgs(trailingOnly=TRUE)
# argument has 2 components:
# (1) a filename parsed into
#             (i)   the index of the criterion for the G matrix definition (criteria_table.csv below)
#             (ii)  the phenotype
#             (iii) the alignment
# and (2) the name of a folder containing the data (".../Gmatrix_data")
# this code infers maximum likelihood parameters of the mixed effect model without any SNP effect included
library(MASS)
library(Matrix)
library(methods)

source("~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/function_mixed_model.R")
source("~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/functions_2_draft_GWAS_v2.R")
#crit_tab <- read.csv("criteria_table.csv")

# parse arguments:
tmp <- process_args(args)
ll <- tmp$best_filtering_idx
aln_name <- tmp$aln_name
pheno_name <- tmp$pheno_name
filename <- tmp$RDatafile_name
Gmatrixdata_folder <- tmp$Gmatrixdata_folder
stopifnot(file.exists(Gmatrixdata_folder))

result_filename <- paste0(Gmatrixdata_folder, "prepared_", pheno_name, "_", aln_name, "_", ll, ".csv")

load(paste0(Gmatrixdata_folder, filename))
cat("loading successful for conditions ", ll, pheno_name, aln_name, "\n")

# stopifnot(isDiagonal(Xrand[[2]]))
# if(length(Xrand) == 3) stopifnot(isDiagonal(Xrand[[3]]))
# X <- as(X, "dgCMatrix") # sparse matrix class from Matrix package

print("now testing likelihood functions...")
npar2 <- length(Xrand) + 2
npar1 <- length(Xrand)
Xrand <- lapply(Xrand, as.matrix)
X <- as.matrix(X)
#Xrand <- lapply(Xrand, function(x) as(x, "dgeMatrix")) # transform in dense matrix class of Matrix package

# define initialising function
min_pheno <- quantile(suba[, pheno_name], probs = c(0.05), na.rm = T)
max_pheno <- quantile(suba[, pheno_name], probs = c(0.5), na.rm = T)
init_fun2 <- function(){
  return(c(runif(npar2 - 1, 0, 1), runif(1, min_pheno, max_pheno)))
}
init_fun1 <- function(){
  return(c(runif(npar1, 0, 1)))
}
lower2 <- c(rep(0, npar2-1), min_pheno)
upper2 <- c(rep(5, npar2 - 1), max_pheno)

lower1 <- c(rep(0, npar1))
upper1 <- c(rep(1, npar1))

t0 <- Sys.time()
lik_test1 <- mynegativeloglik1(par = init_fun1(), data = list(phenotype = suba[, pheno_name]), X = X, n = nind, Xrand = Xrand, npar = npar1)
print(lik_test1)
Sys.time() - t0

t0 <- Sys.time()
lik_test2 <- mynegativeloglik2(par = init_fun2(), data = list(phenotype = suba[, pheno_name]), X = X, n = nind, Xrand = Xrand, npar = npar2)
print(lik_test2)
Sys.time() - t0


# optimise the likelihood properly
if(infer_two_error_variances <- F){ # this takes a lot of time and is not used
  print("inferring maximum likelihood random effect parameters for model with two error variances...")
  t0 <- Sys.time()
  #cl <- makeCluster(3) ; clusterExport(cl, c("nmkb", "mynegativeloglik2", "mynegativeloglik1", "solve.for.lik", "generate_sigma2", "generate_sigma1", "X", "Xrand", "nind", "suba", "pheno_name", "npar1", "npar2", "lower1", "upper1", "lower2", "upper2","min_pheno", "max_pheno", "init_fun1", "init_fun2", "get.inv.det.sigma", "get.beta.lik"))
  all_opt2 <- optim.fun.repeated(n.repeats = 3, lik.fun = mynegativeloglik2,
                                  init.fun = init_fun2, verbose = T, lower = lower2, upper = upper2, control = list(tol = 0.001),
                                  data = list(phenotype = suba[, pheno_name]), X = X, n = nind, Xrand = Xrand, npar = npar2)
  # stopCluster(cl)
  Sys.time() - t0
}
print("inferring maximum likelihood random effect parameters for main model...")
all_opt1 <- optim.fun.repeated(n.repeats = 3, lik.fun = mynegativeloglik1,
                               init.fun = init_fun1, verbose = T, lower = lower1, upper = upper1, control = list(tol = 0.001),
                               data = list(phenotype = suba[, pheno_name]), X = X, n = nind, Xrand = Xrand, npar = npar1)

opt1 <- clean_lik(all_opt1, npar1)
print("computing heritability and saving optimal parameters...")
h2_1 <- get_h2(par = opt1$pars, data = list(phenotype = suba[, pheno_name]), X = X, n = nind, Xrand = Xrand, npar = npar1, lik_fun_name = "mynegativeloglik1")
print(h2_1)

if(infer_two_error_variances){
  opt2 <- clean_lik(all_opt2, npar2)
  h2_2 <- get_h2(par = opt2$pars, data = list(phenotype = suba[, pheno_name]), X = X, n = nind, Xrand = Xrand, npar = npar2, lik_fun_name = "mynegativeloglik2")
  print(h2_2)
} else {
  opt2 <- list(pars = rep(NA, length(opt1$pars) + 2), lik = NA); h2_2 <- NA
}

results_tab <- matrix(c(opt1$pars, NA, NA, opt1$lik, h2_1, opt2$pars, opt2$lik, h2_2), ncol = 2)
write.csv(results_tab, file = result_filename, row.names = F)
quit()

# look at predicted distributions
# if(F){
#   h2_list1 <- get_h2(par = opt1$pars, data = list(phenotype = suba[, pheno_name]), X = X, n = nind, Xrand = Xrand, npar = npar1, lik_fun = "mynegativeloglik1", nsim = 2, return_full = T)
#   h2_list2 <- get_h2(par = opt2$pars, data = list(phenotype = suba[, pheno_name]), X = X, n = nind, Xrand = Xrand, npar = npar2, lik_fun = "mynegativeloglik2", nsim = 2, return_full = T)
#   
#   
#   with(h2_list1,
#        plot(density(e[1,] + g[1,]), type = "l", col = "red", xlim = c(-5, 5), ylim = c(0, 3))
#   )
#   points(density(suba[, pheno_name]), type = "l", col = "blue")
#   
#   with(h2_list2,
#        plot(density(e[1,] + g[1,]), type = "l", col = "red", xlim = c(-5, 5), ylim = c(0, 3))
#   )
#   points(density(suba[, pheno_name]), type = "l", col = "blue")
# }





