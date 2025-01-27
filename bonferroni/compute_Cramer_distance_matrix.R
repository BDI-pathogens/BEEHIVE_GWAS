args = commandArgs(trailingOnly=TRUE)
# args is a vector of length 3
# 1 phenotype name
# 2 type of GWAS = prot, freq, kmer, len or len3
# 3 do we take subsets of the positions or not?
#           0 NO
#           1 YES, only 1-mers or positions 1:5000
#           2 YES, 2-mers      or positions 5001:10000
#           3 YES, 3-mers      or positions 10001:15000
#           ...
#           7 YES, 7-mers      or positions 30001:35000
#           8 YES, 8-mers      or positions 35001:40000

if(length(args)==0){
  args <- c("BEEHIVE_LVL", "freq", 0) # argument = phenotype and type of GWAS
  wd <- "~/DID/BEEHIVE_Hackathon/"
  #wd <- "/Users/Christophe/Dropbox (Infectious Disease)/PROJECTS/HIV/BEEHIVE_Hackathon/"
  setwd(paste0(wd, "Code/DevelopMethods/LMM/SolvingLMMinR/"))
}

library(MASS)
library(Matrix)
library(methods)
#library(rcompanion)

# HERE WE COMPUTE A MATRIX OF DISTANCE AS CRAMER'S V BETWEEN PAIRS OF POSITIONS

# compute cramer's v between pairs of positions
mycramersV <- function(tab){
  chi2 <- chisq.test(tab, correct = F)$statistic # use correct = F because with the Yates correction done for 2x2 table, the cramer's v is not 1 even for perfect correlation
  k <- min(dim(tab))
  n <- sum(tab)
  V <- sqrt(chi2/(n * (k - 1)))
  V
}
mycorr <- function(partition1, partition2){
  return(cramerV(table(partition1, partition2), bias.correct = T))
  # rcompanion cramer'V with bias correction 
}
cramerV <- function (x, y = NULL, digits = 4, bias.correct = FALSE, ...) # function from rcompanion package
{
  CV = NULL
  if (is.factor(x)) {
    x = as.vector(x)
  }
  if (is.factor(y)) {
    x = as.vector(y)
  }
  if (is.vector(x) & is.vector(y)) {
    N = length(x)
    Chi.sq = suppressWarnings(chisq.test(x, y, correct = FALSE, 
                                         ...)$statistic)
    Phi = Chi.sq/N
    R = length(unique(x))
    C = length(unique(y))
    CV = sqrt(Phi/min(R - 1, C - 1))
  }
  if (is.matrix(x)) {
    x = as.table(x)
  }
  if (is.table(x)) {
    N = sum(x)
    Chi.sq = suppressWarnings(chisq.test(x, correct = FALSE, 
                                         ...)$statistic)
    Phi = Chi.sq/N
    R = nrow(x)
    C = ncol(x)
    CV = sqrt(Phi/min(R - 1, C - 1))
  }
  if (bias.correct == TRUE) {
    Phi = max(0, Phi - ((R - 1) * (C - 1)/(N - 1)))
    CC = C - ((C - 1)^2/(N - 1))
    RR = R - ((R - 1)^2/(N - 1))
    CV = sqrt(Phi/min(RR - 1, CC - 1))
  }
  CV = signif(as.numeric(CV), digits = digits)
  names(CV) = "Cramer V"
  return(CV)
}
matrix_correlation <- function(x1, y1){
  
  fit <- lm(y1 ~ x1 + 0) 
  predicted.y <- x1 %*% coef(fit)
  predicted.y <- predicted.y + mean(y1 - predicted.y)
  
  sum.norm.sqr.residuals <- sum((y1 - predicted.y)^2)
  sum.norm.sqr.diff.from.means <-sum((y1 - scale(y1, scale = FALSE))^2)
  r2 <- 1 - sum.norm.sqr.residuals / sum.norm.sqr.diff.from.means
  r2 <- ifelse(r2 < 0, 0, r2)
  return(r2)
}
get_idx_in_design <- function(idx_in_tab, mytab, mydf){
  if(idx_in_tab==1){
    accumulated_dfs <- 0
  } else {
    accumulated_dfs <- sum(mytab$df[1:(idx_in_tab-1)]) # dfs accumulated till focal position
  }
  return(
    (accumulated_dfs + 1):(accumulated_dfs + mydf)
  )
}

stopifnot(length(args) == 3)
RDatafile_name <- paste0("firstGWAS_", args[1], "_full_17.RData")
gwas_name <- args[2]
subset_only <- as.numeric(args[3])


source("function_mixed_model.R")

# LOAD FILE
print("loading data...")
load(paste0("GWAS_results/", RDatafile_name)) # load the prepared data

output_filename <- paste0("dist_table_", args[1], "_", args[2], "_", args[3], ".csv")

my_GWAS <- get(paste0("GWAS_", args[2]))               # result-of-GWAS matrix
my_dm   <- get(paste0("full_design_matrix_", args[2])) # design matrix for each position
my_ind  <- get(paste0("full_ind_tokeep_", args[2]))    # individuals to keep at each position

to_keep <- c("output_filename", "args", "my_GWAS", "my_dm", "my_ind", "mycorr", "cramerV", "matrix_correlation", "subset_only")

rm(list = ls()[!ls() %in% to_keep]) # ditch the rest to free memory space

stopifnot(sum(my_GWAS$df)==ncol(my_dm))

# number of each type of test (1, 2, ..., 6-mer)
# table(GWAS_kmer$end_position_alignment - GWAS_kmer$start_position_alignment + 1)

n_tests_total <- nrow(my_GWAS)
nind <- nrow(my_ind)

# first re-construct multiple alleles per position if categorical alleles (i.e. not a frequency)
if(args[2] != "freq"){
  
  multialleles <- matrix("NA", nrow = nind, ncol = n_tests_total)
  k <- 1
  for(i in 1:n_tests_total){
    if(i/100==floor(i/100)) print(i)
    my_df <- my_GWAS$df[i]
    my_design <- my_dm[ , k:(k+my_df-1), drop = F]
    multialleles[, i] <- apply(my_design, 1, function(vec) paste(vec, collapse = ""))
    k <- k + my_df
  }
  stopifnot(k - my_df==ncol(my_dm)) # check that k is as intended
  
  # replace the NAs
  multialleles[grepl(pattern = "NA", multialleles)] <- NA
  multialleles <- as.matrix(multialleles, nrow = nind)
  
}
index_in_dm <- cumsum(my_GWAS$df) # index of each test in design matrix

if(subset_only == 0){ # do not subset (only for freq)
  subset_tests <- seq(n_tests_total)
} else {
  # subset the positions if we want
  if(args[2] != "kmer"){  # for len, len3 variants, we want to subset by blocks of 5000 positions
    subset_tests <- (5000*(subset_only - 1) + 1):(5000*subset_only)
    if(min(subset_tests) > n_tests_total) stop('index of positions to consider greater than total number of positions')
    if(max(subset_tests) > n_tests_total) subset_tests <- (5000*(subset_only - 1) + 1):n_tests_total
  } else {                # for kmer variants, we want to subset by k-mer length
    subset_tests <- which(my_GWAS$end_position_alignment - my_GWAS$start_position_alignment + 1 == subset_only)
  }
}

n_tests <- length(subset_tests)

# NOW COMPUTING MATRIX OF DISTANCES BETWEEN POSITIONS
print("computing all distances... (takes ~ 10 hours)")
#all_dist <- matrix(NA, nrow = n_tests, ncol = n_tests) # would be wise to NOT define such a big matrix
#all_dist <- matrix(NA, ncol = n_tests, nrow = n_tests)

if(args[2] == "freq") end_idx <- n_tests else end_idx <- n_tests-1
for(i in 1:end_idx){
  print(i)
  if(args[2] == "freq") begin_idx <- 1 else begin_idx <- i + 1
  
  if(subset_tests[i] == 1){
    x <- my_dm[, 1:my_GWAS$df[1], drop = F]
  } else {
    x <- my_dm[ , (index_in_dm[subset_tests[i] - 1] + 1):index_in_dm[subset_tests[i]], drop = F]
  }
  
  for(j in begin_idx:n_tests){
    
    my_subset <- which(my_ind[, subset_tests[i]] & my_ind[ , subset_tests[j]]) # subset of non-NA individuals for this position
    
    if(args[2] %in% c("prot", "kmer", "len", "len3")){
      my_distance_measure <- mycorr(multialleles[my_subset, subset_tests[i]], multialleles[my_subset, subset_tests[j]])
    }
    if(subset_tests[j] == 1){
      y <- my_dm[, 1:my_GWAS$df[1], drop = F]
    } else {
      y <- my_dm[ , (index_in_dm[subset_tests[j] - 1] + 1):index_in_dm[subset_tests[j]], drop = F]
    }
    if(args[2] == "freq"){
      my_distance_measure <- matrix_correlation(x[my_subset, , drop = F], y[my_subset, , drop = F]) # distance for frequency
    }
    
    #if(ncol(x)==1)stopifnot(all(table(x[my_subset, ] == multialleles[my_subset, subset_tests[i]])))
    #if(ncol(y)==1)stopifnot(all(table(y[my_subset, ] == multialleles[my_subset, subset_tests[j]])))
    
    # write directly in a file (faster than write.table)
    if(!is.na(my_distance_measure)){
      if(my_distance_measure != 0){ # write only non-zero
        data.table::fwrite(
          x = list(i, j, my_distance_measure),
          file = output_filename, append = TRUE, sep = "," , quote = TRUE
        )
      }
    }
  }
}
#diag(all_dist) <- 1

quit()