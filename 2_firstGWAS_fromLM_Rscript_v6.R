args = commandArgs(trailingOnly=TRUE)
if(length(args)==0) args <- c(
  "prepared_BEEHIVE_LVL_full_17.RData", # prepared phenotypic data
  #"~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/Gmatrix_data/" # and where to find it--on this depends whether we run the discovery or replication GWAS 
  "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/Gmatrix_data/" # for discovery
  )
library(MASS)
library(Matrix)
library(methods)
#wd <- "~/DID/BEEHIVE_Hackathon/" #wd <- "/Users/Christophe/Dropbox (Infectious Disease)/PROJECTS/HIV/BEEHIVE_Hackathon/"
#if(grepl(pattern = "/Users/fblanqua", getwd())) setwd(paste0(wd, "Code/DevelopMethods/LMM/SolvingLMMinR/"))

source("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/function_mixed_model.R")
source("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/functions_2_draft_GWAS_v2.R")

do_prot_GWAS <- TRUE
do_kmer_GWAS <- TRUE
do_len3_GWAS <- TRUE
do_len_GWAS <- TRUE
do_freq_GWAS <- TRUE
do_rare_GWAS <- TRUE
do_drm_GWAS <- TRUE
remove_effects <- FALSE # optional; to remove the effects from the GWAS table that is saved (useful for the replication GWAS)

cat("Some options have to be set manually, TO CHECK: ", c(do_prot_GWAS, do_kmer_GWAS, do_len3_GWAS, do_len_GWAS, do_freq_GWAS, do_rare_GWAS, do_drm_GWAS, remove_effects), "\n")

# PARSE ARGUMENTS AND LOAD DATA
tmp <- process_args(args); attach(tmp) # note that this defines a default argument that's used when args is empty
print(best_filtering_idx)
#filename <- paste0(wd, "Code/DevelopMethods/LMM/SolvingLMMinR/Gmatrix_data/prepared_", pheno_name, "_", aln_name, "_", best_filtering_idx, ".RData")
cat("loading data ", RDatafile_name, " ...\n")
load(paste0(Gmatrixdata_folder, RDatafile_name)) # load the prepared data
Xrand <- lapply(Xrand, as.matrix)
X <- as.matrix(X)

# LOAD OPTIMAL LMM RANDOM EFFECT PARAMETERS
optim_results <- load_optim(args = args, myXrand = Xrand, mysuba = suba, myX = X, mynind = nind) # load opimal parameters and check them 
myattach(optim_results) # this defines variables opt_pars, npar, fixed_effects, log.det.Sigma, inverse.Sigma
n_fixed_effects <- length(fixed_effects)

# filename_prot <- paste0(wd, "Code/DevelopMethods/LMM/SolvingLMMinR/Gmatrix_data/prepared_", pheno_name, "_prot.RData")
# load(filename_prot)

# WHILE WE ARE AT IT, COMPUTE THE BLUPS FOR EACH INDIVIDUAL
inverse.Sigma <- as(inverse.Sigma, "dgeMatrix")
blups <- opt_pars[1] * Xrand$Kmat %*% inverse.Sigma %*% (suba[, pheno_name] - X %*% fixed_effects) # from Lippert thesis paragraph 2.2.2
pred <- X %*% fixed_effects + blups
as.vector(pred) -> pred
R2 <- 1 - sum((suba[, pheno_name] - pred)^2) / sum((suba[, pheno_name] - mean(suba[, pheno_name], na.rm = T))^2) # R2 for full model with fixed effects

######                                     MAIN GWAS WITH PROTEINS                                                ######

if(do_prot_GWAS){
  load(paste0(Gmatrixdata_folder, "prot.RData"))
  Lvec_prot <- cbind(window.starts_prot, window.ends_prot)
  max_df <- max(table(paste(window.starts_prot, window.ends_prot, sep = "_"))) # maximum df (maximum number of alternative alleles)
  
  dim(G_prot)
  stopifnot(all(idsG %in% idsG_prot)) # check that all beehive ids of individuals in G matrix are also in prot matrix
  stopifnot(all(idsG[sub_ind] == suba$PATIENT))
  reorder_G_prot <- match(idsG, idsG_prot) # reorder the G_prot matrix
  
  GWAS_prot <- do_GWAS(my_G = G_prot[reorder_G_prot, ], my_Lvec = Lvec_prot, my_reference = reference_prot, my_alternative = alternative_prot,
                       my_epidata = suba, my_sub_ind = sub_ind, my_X = X, my_inverse.Sigma = inverse.Sigma, my_log.det.Sigma = log.det.Sigma, my_pheno_name = pheno_name, my_n_fixed_effects = n_fixed_effects, max_df = max_df, return_design_matrix = T)
  myattach(GWAS_prot)
  GWAS_prot <- all_lik_beta
  full_design_matrix_prot <- full_design_matrix
  full_ind_tokeep_prot <- full_ind_tokeep
  
  rm(G_prot, full_design_matrix, full_ind_tokeep)
}

######                                     MAIN GWAS WITH K-MERS, K IN [1;6]                                                ######
if(do_kmer_GWAS){
  
  load(paste0(Gmatrixdata_folder, "kmers.RData"))
  Lvec_kmer <- cbind(window.starts_kmer, window.ends_kmer)
  max_df <- max(table(paste(window.starts_kmer, window.ends_kmer, sep = "_"))) # maximum df (maximum number of alternative alleles)
  
  dim(G_kmer)
  stopifnot(all(idsG %in% idsG_kmer)) # check that all beehive ids of individuals in G matrix are also in kmer matrix
  stopifnot(all(idsG[sub_ind] == suba$PATIENT))
  reorder_G_kmer <- match(idsG, idsG_kmer) # reorder the G_kmer matrix
  
  GWAS_kmer <- do_GWAS(my_G = G_kmer[reorder_G_kmer, ], my_Lvec = Lvec_kmer, my_reference = reference_kmer, my_alternative = alternative_kmer,
                       my_epidata = suba, my_sub_ind = sub_ind, my_X = X, my_inverse.Sigma = inverse.Sigma, my_log.det.Sigma = log.det.Sigma, my_pheno_name = pheno_name, my_n_fixed_effects = n_fixed_effects, max_df = max_df, return_design_matrix = T)
  myattach(GWAS_kmer)
  GWAS_kmer <- all_lik_beta
  full_design_matrix_kmer <- full_design_matrix
  full_ind_tokeep_kmer <- full_ind_tokeep
  
  rm(G_kmer, full_design_matrix, full_ind_tokeep)
  
}

######                                     MAIN GWAS WITH LENGTH MODULO 3 VARIANTS                                                ######

if(do_len3_GWAS){
  
  load(paste0(Gmatrixdata_folder, "lengths3.RData"))
  Lvec_len3 <- cbind(window.starts_len3, window.ends_len3)
  max_df <- max(table(paste(window.starts_len3, window.ends_len3, sep = "_"))) # maximum df (maximum number of alternative alleles)
  
  dim(G_len3)
  stopifnot(all(idsG %in% idsG_len3)) # check that all beehive ids of individuals in G matrix are also in kmer matrix
  stopifnot(all(idsG[sub_ind] == suba$PATIENT))
  reorder_G_len3 <- match(idsG, idsG_len3) # reorder the G_len3 matrix
  stopifnot(reorder_G_len3 == reorder_G_kmer)
  
  GWAS_len3 <- do_GWAS(my_G = G_len3[reorder_G_len3, ], my_Lvec = Lvec_len3, my_reference = reference_len3, my_alternative = alternative_len3,
                       my_epidata = suba, my_sub_ind = sub_ind, my_X = X, my_inverse.Sigma = inverse.Sigma, my_log.det.Sigma = log.det.Sigma, my_pheno_name = pheno_name, my_n_fixed_effects = n_fixed_effects, max_df = max_df, return_design_matrix = T)
  myattach(GWAS_len3)
  GWAS_len3 <- all_lik_beta
  full_design_matrix_len3 <- full_design_matrix
  full_ind_tokeep_len3 <- full_ind_tokeep
  
  rm(G_len3, full_design_matrix, full_ind_tokeep)
}

######                                     MAIN GWAS WITH LENGTH VARIANTS                                                             ######
if(do_len_GWAS){
    
  load(paste0(Gmatrixdata_folder, "lengths.RData"))
  Lvec_len <- cbind(window.starts_len, window.ends_len)
  max_df <- max(table(paste(window.starts_len, window.ends_len, sep = "_"))) # maximum df (maximum number of alternative alleles)
  
  dim(G_len)
  stopifnot(all(idsG %in% idsG_len)) # check that all beehive ids of individuals in G matrix are also in kmer matrix
  stopifnot(all(idsG[sub_ind] == suba$PATIENT))
  reorder_G_len <- match(idsG, idsG_len) # reorder the G_len matrix
  stopifnot(reorder_G_len == reorder_G_kmer)
  
  GWAS_len <- do_GWAS(my_G = G_len[reorder_G_len, ], my_Lvec = Lvec_len, my_reference = reference_len, my_alternative = alternative_len,
                      my_epidata = suba, my_sub_ind = sub_ind, my_X = X, my_inverse.Sigma = inverse.Sigma, my_log.det.Sigma = log.det.Sigma, my_pheno_name = pheno_name, my_n_fixed_effects = n_fixed_effects, max_df = max_df, return_design_matrix = T)
  myattach(GWAS_len)
  GWAS_len <- all_lik_beta
  full_design_matrix_len <- full_design_matrix
  full_ind_tokeep_len <- full_ind_tokeep
  
  rm(G_len, full_design_matrix, full_ind_tokeep)
  
}

######                                     MAIN GWAS WITH BASE FREQUENCIES VARIANTS                                                    ######
if(do_freq_GWAS){
  load(paste0(Gmatrixdata_folder, "base_frequencies.RData"))
  Lvec_freq <- cbind(window.starts_freq, window.ends_freq)
  max_df <- max(table(paste(window.starts_freq, window.ends_freq, sep = "_"))) # maximum df (maximum number of alternative alleles)
   
  stopifnot(all(idsG %in% idsG_freq)) # check that all beehive ids of individuals in G matrix are also in kmer matrix
  stopifnot(all(idsG[sub_ind] == suba$PATIENT))
  reorder_G_freq <- match(idsG, idsG_freq) # reorder the G matrix
  #stopifnot(all(reorder_G_freq==reorder_G_kmer))
   
  GWAS_freq <- do_GWAS(my_G = G_freq[reorder_G_freq, ], my_Lvec = Lvec_freq, my_reference = reference_freq, my_alternative = alternative_freq, my_epidata = suba, my_sub_ind = sub_ind, my_X = X,
                       my_inverse.Sigma = inverse.Sigma, my_log.det.Sigma = log.det.Sigma, my_pheno_name = pheno_name, my_n_fixed_effects = n_fixed_effects, max_df = max_df, return_design_matrix = T)
  myattach(GWAS_freq)
  GWAS_freq <- all_lik_beta
  full_design_matrix_freq <- full_design_matrix
  full_ind_tokeep_freq <- full_ind_tokeep
  
  rm(G_freq, full_design_matrix, full_ind_tokeep)
}

######                                        NOW RUN THE RARE VARIANTS ANALYSIS                                             ######
if(do_rare_GWAS){
  GWAS_rare <- do_rare_variants(my_freq_rare_variants = c(0.001, 0.005, 0.01, 0.05, 0.1), my_G = G, my_epidata = suba, my_sub_ind = sub_ind, my_X = X,
                                my_inverse.Sigma = inverse.Sigma, my_log.det.Sigma = log.det.Sigma, my_pheno_name = pheno_name, my_n_fixed_effects = n_fixed_effects)
}

######                                        NOW RUN THE DRM ANALYSIS                                             ######

if(do_drm_GWAS){
  GWAS_drm <- do_drm(my_epidata = suba, my_sub_ind = sub_ind, my_X = X,
                     my_inverse.Sigma = inverse.Sigma, my_log.det.Sigma = log.det.Sigma, my_pheno_name = pheno_name, my_n_fixed_effects = n_fixed_effects)
  
}

######                                        TABLE OF RESULTS POST-PROCESSING                                               ######

# re-format
if(exists("GWAS_prot")) GWAS_prot <-  clean_table(GWAS_prot, remove_effects = remove_effects)
if(exists("GWAS_kmer")) GWAS_kmer <-  clean_table(GWAS_kmer, remove_effects = remove_effects)
if(exists("GWAS_len")) GWAS_len <-    clean_table(GWAS_len, remove_effects = remove_effects)
if(exists("GWAS_len3")) GWAS_len3 <-  clean_table(GWAS_len3, remove_effects = remove_effects)
if(exists("GWAS_freq")) GWAS_freq <-  clean_table(GWAS_freq, remove_effects = remove_effects)
if(exists("GWAS_rare")) GWAS_rare <-  clean_table(GWAS_rare, remove_effects = remove_effects)
if(exists("GWAS_drm")) GWAS_drm <-    clean_table(GWAS_drm, remove_effects = remove_effects)

all_design_matrix_names <- c("full_design_matrix_prot", "full_design_matrix_kmer", "full_design_matrix_len", "full_design_matrix_len3", "full_design_matrix_freq",
                             "full_ind_tokeep_prot", "full_ind_tokeep_kmer", "full_ind_tokeep_len", "full_ind_tokeep_len3", "full_ind_tokeep_freq")
for(n in all_design_matrix_names){
  if(exists(n)) assign(n, remove_empty_rows_cols(get(n)))
}

if(draw_plots <- F){
  draw_plot(GWAS_kmer)
}

# significant hits
# if(exists("GWAS_prot")) print_bon_significant(GWAS_prot)
# if(exists("GWAS_kmer")) print_bon_significant(GWAS_kmer)
# if(exists("GWAS_len")) print_bon_significant(GWAS_len)
# if(exists("GWAS_len3")) print_bon_significant(GWAS_len3)
# if(exists("GWAS_freq")) print_bon_significant(GWAS_freq)
# if(exists("GWAS_rare")) print_bon_significant(GWAS_rare)
# if(exists("GWAS_drm")) print_bon_significant(GWAS_drm)

# SAVE RESULTS OF THE GWAS
to_save <- c("opt_pars", "fixed_effects", "n_fixed_effects", "blups", "pred", "nb_alternative", "nb_ind_per_gpos", "nb_gpos_per_ind", "unique_pos_to_test",
             "GWAS_prot", "GWAS_kmer", "GWAS_len", "GWAS_len3", "GWAS_freq", "GWAS_rare", "GWAS_drm",
             all_design_matrix_names
             )
if(!all(to_save %in% ls())){
  print("warning not all variables exist - saving those that exist only")
  to_save <- to_save[to_save %in% ls()]
}

results_folder <- gsub(pattern = "Gmatrix_data/", replacement = "GWAS_results/", x = Gmatrixdata_folder)
save(list = to_save, file = paste0(results_folder, "firstGWAS_", pheno_name, "_", aln_name, "_", best_filtering_idx, ".RData"))

quit()




