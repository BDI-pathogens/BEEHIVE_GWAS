 args = commandArgs(trailingOnly=TRUE)
# argument is a suffix describing to what GWAS the matrix corresponds
if(length(args)==0){
  args <- c("CD4_slope", "freq", "0") 
  wd <- "~/DID/BEEHIVE_Hackathon/"
  setwd(paste0(wd, "Code/DevelopMethods/LMM/SolvingLMMinR/bonferroni"))
}
library(Rcpp)
library(Matrix)

sourceCpp(code = '
          #include <RcppArmadillo.h>
          // [[Rcpp::depends(RcppArmadillo)]]
          // [[Rcpp::export]]
          arma::vec getEigenValues(arma::mat M) {
          return arma::eig_sym(M);
          }
          ')


#### SECOND STEP: COMPUTE EIGENVALUES ####
output_filename <- paste0("dist_table_", args[1], "_", args[2], "_", args[3], ".csv")
print(getwd())
print(output_filename)
mat <- read.csv(output_filename, header = F)
n_tests <- max(mat$V2)
#stopifnot(max(mat$V1)==n_tests-1)
#stopifnot(nrow(mat) == n_tests*(n_tests-1)/2)

# set infinite elements to 0
mat$V3[which(is.infinite(mat$V3))] <- 0

# set NA elements to 0
mat$V3[which(is.na(mat$V3))] <- 0

# transform into a matrix
mm <- sparseMatrix(i = mat$V1, j = mat$V2, x = mat$V3, dims = c(n_tests, n_tests)) # this matrix is all 0 except at positions specified by i and j
mm <- as.matrix(mm)
#stopifnot(all(diag(mm)==0))
#stopifnot(all(mm[lower.tri(mm, diag = T)]==0))

# the case where we did not compute distance for diagonal and lower triangular -> set diagonal to 1
if(all(diag(mm)==0) & all(mm[lower.tri(mm, diag = T)]==0)){
  diag(mm) <- 1 # set diagonal elements to 1
  mm[lower.tri(mm)] <- t(mm)[lower.tri(mm)] # copy uper triangular to lower triangular
  stopifnot(isSymmetric.matrix(mm)) # check symmetry of matrix
} else { # otherwise (this is for the distance measure corresponding to within-host frequency) we computed all elements
  stopifnot(all(diag(mm)==1)) # check diagonal is all 1
  stopifnot(any(mm[lower.tri(mm, diag = T)]!=0)) # check lower triangular is not all 0
}

t0 <- Sys.time()
tmp <- getEigenValues(mm)
Sys.time()-t0

sorted_ev <- sort(tmp, decreasing = T)
sorted_ev <- sorted_ev / sum(sorted_ev)
n_tests_effective <- which(cumsum(sorted_ev) > 0.995)[1]

n_tests_effective
plot(log10(sorted_ev), type = "o", pch = 20)

write.csv(file = paste0("n_effective_tests_", args[1], "_", args[2], "_", args[3], ".csv"),
          x = list(n_effective_tests = n_tests_effective, n_tests = n_tests), row.names = F)


