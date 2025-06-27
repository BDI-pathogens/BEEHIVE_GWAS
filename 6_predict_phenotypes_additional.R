args = commandArgs(trailingOnly=TRUE)
library(plyr)
library(MASS)

find_reference_alternative <- function(allele_label){
  # Function to find the reference and alternative alleles in discovery from the allele label in the form "type_start_end_n"
  # Note that this information is normally in files of the form kmers.RData, lengths.RData, etc.
  # in objects of the form reference_kmer, alternative_kmer, etc.
  # but can be retrieved more simply and more precisely from the GWAS table (some alleles with missing data are pooled at the GWAS step)
  # allele_label <- "len_9249_9255_1"
  # returns a vector with type of allele, start, end, idx in GWAS_table, ref, alternative, number
  
  tmp <- strsplit(x = allele_label, split = "_")[[1]]; stopifnot(length(tmp) == 4)
  type <- tmp[1]
  st <- as.numeric(tmp[2])
  en <- as.numeric(tmp[3])
  no <- as.numeric(tmp[4]) # number of the allele
  GWAS_table <- get(paste0("GWAS_", type))
  idx <- which(GWAS_table$start_position_alignment == st & GWAS_table$end_position_alignment == en); stopifnot(length(idx) == 1)
  ref <- as.character(GWAS_table[idx, c("reference")])
  alt <- as.character(GWAS_table[idx, paste0(c("alternative"), no)])
  return(c(type, st, en, idx, ref, alt, no))
}

if(length(args)==0){
  print("warning length of args is 0, replacing by default")
  args <- c("spvl_adjusted", 0.2)
  print(args)
  setwd("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/")
}

list_args <- list(
  c("CD4_slope", 0.05),
  c("CD4_slope", 0.2),
  
  c("BEEHIVE_LVL", 0.05),
  c("BEEHIVE_LVL", 0.2),
  
  c("spvl_normalised_adjusted", 0.05),
  c("spvl_normalised_adjusted", 0.2),
  
  c("spvl_adjusted", 0.05),
  c("spvl_adjusted", 0.2)
)

# introduced 19/03/2025; loop on phenotype/p-value
all_res <- c()
for(ii in 1:length(list_args)){
args <- list_args[[ii]]
print(args)

stopifnot(length(args)==2)
pheno_name <- args[1] # "spvl_normalised_adjusted"
p_t <- as.numeric(as.character(args[2])) # the p_value_threshold

####################                                        1. load results of re-optimised model to get the relevant positions to include         ####################

load(paste0("3_reoptimisation/re_optimisation_", pheno_name, "_", p_t, ".RData"))
if(use_kmer_only) my_h2 <- h2_full1 else my_h2 <- h2_full2 # choose the model to re-optimise (QUESTION: only kmer/prot or both kmer/prot and len variants ?)

# the relevant effects to consider are those:
nvariants_nofreq <- sum(grepl("kmer_", rownames(my_h2$fixed_effects))) # number of variants considered when removing frequency effects (and rare)
my_h2$fixed_effects[n_fixed_effects+1:nvariants_nofreq] # more precisely

# get list of reference and alternative alleles in discovery:
list_ref_alt_disc <- sapply(rownames(my_h2$fixed_effects)[n_fixed_effects+1:nvariants_nofreq], find_reference_alternative)
list_ref_alt_disc <- data.frame(t(list_ref_alt_disc))
names(list_ref_alt_disc) <- c("type", "start", "end", "idx_GWAS", "reference", "alternative", "allele_no")
list_ref_alt_disc$variants <- paste(list_ref_alt_disc$start, list_ref_alt_disc$end, list_ref_alt_disc$reference, list_ref_alt_disc$alternative, sep = "_")

####################                 2. load the prepared replication data & GWAS results to get design matrices                                   ####################

# remove objects of the discovery data (this will avoid mistakes later on as we will now load the same objects for replication)
objects_to_remove <- c("aln", "alternative", "blups", "fixed_covariate_names",
                       "fixed_covariate_type", "fixed_effects", "full_design_matrix_freq", "full_design_matrix_kmer",
                       "full_design_matrix_len", "full_design_matrix_len3", "full_design_matrix_prot", "full_ind_tokeep_freq", "full_design_matrix_rare",
                       "full_ind_tokeep_kmer", "full_ind_tokeep_len", "full_ind_tokeep_len3", "full_ind_tokeep_prot",  
                       "G", "Gn", "GWAS_drm", "GWAS_freq", 
                       "GWAS_kmer", "GWAS_len",  "GWAS_len3", "GWAS_prot", 
                       "GWAS_rare", "ids_tokeep_epi", "idsG", "Lvec",     
                       "n_fixed_covariates", "n_fixed_effects", "nb_alternative", "nb_gpos_per_ind",   
                       "nb_ind_per_gpos", "ngpos", "nind", "normalise_factors",
                       "opt_pars", "PC01_AE", "PC02_AG", "PCA1",     
                       "PCB", "PCC", "PCD", "PCother",  
                       "PCsd01_AE", "PCsd02_AG", "PCsdA1", "PCsdB",    
                       "PCsdC", "PCsdD", "PCsdother",  
                       "pred", "reference", "Sigma_sub", "sub_ind",  
                       "sub_pos", "suba", "subtype_01_AE", "subtype_02_AG",     
                       "subtype_A1", "subtype_B", "subtype_C", "subtype_D", 
                       "subtype_other", "unique_pos_to_test", "X", "Xrand")

rm(list = objects_to_remove[objects_to_remove %in% ls()])

# load prepared replication data and GWAS results
load(paste0("~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/Gmatrix_data/prepared_", pheno_name, "_full_17.RData"))
load(paste0("~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/GWAS_results/firstGWAS_", pheno_name, "_full_17.RData"))

# check objects are indeed from the replication data now:
stopifnot(nrow(full_design_matrix_kmer)==nind)
stopifnot(nrow(full_ind_tokeep_kmer)==nind)
stopifnot(nrow(X)==nind)

get_index_in_dm <- function(st = NULL, en = NULL, idx_in_GWAS = NULL, GWAS_tab, full_design_matrix, full_ind_tokeep){
  
  # this function get column position in design matrix of variant starting at "st" and ending at "en" position in GWAS_tab
  # and the column position in "individuals to keep" matrix of this variant
  # alternatively, the variant can be specified by the row number, idx_in_GWAS in GWAS_tab
  
  if((is.null(st) | is.null(en)) & is.null(idx_in_GWAS)) stop("st and en, or idx_in_GWAS have to be specified")
  if(is.null(idx_in_GWAS)){
    idx_in_GWAS <- which(GWAS_tab$start_position_alignment==st & GWAS_tab$end_position_alignment == en)
  }
  if(length(idx_in_GWAS) == 0){
    print("variant not found in GWAS table")
    return(c(NA, NA, NA))
  }
  stopifnot(length(idx_in_GWAS)==1)
  cum_df <- cumsum(GWAS_tab$df)
  if(idx_in_GWAS > 1){
    idx_in_dm_start <- 1 + cum_df[idx_in_GWAS-1]
    idx_in_dm_end <- 1 + cum_df[idx_in_GWAS] - 1 # position in list of design matrices
  } else {
    idx_in_dm_start <- 1
    idx_in_dm_end <- GWAS_tab$df[idx_in_GWAS]
  }
  
  idx_in_ind <- idx_in_GWAS # idx in ind is the same as index in GWAS
  
  # check that the indices are correct by checking this is the same total number of individuals considered and N1, N2 etc for variants
  stopifnot(sum(full_ind_tokeep[, idx_in_ind])==GWAS_tab$N[idx_in_GWAS])
  stopifnot(
    all(colSums(full_design_matrix[, idx_in_dm_start:idx_in_dm_end, drop = F] >= 0.01, na.rm = T) == GWAS_tab[idx_in_GWAS, paste0("N", 1:GWAS_tab$df[idx_in_GWAS])])
    # the ">= 0.01" is for frequency variants
  )
  vec <- c(idx_in_ind, idx_in_dm_start, idx_in_dm_end)
  names(vec) <- c("idx_in_ind", "idx_in_dm_start", "idx_in_dm_end")
  return(vec)
  #apply(, 2, sum, na.rm = T)
} # need to re-enter this function (because the saved version in .RData file was uncorrected for some small mistake)

# generate a list containing, for all the nucleotidic positions retained in discovery,
# the reference/alternative alleles tested in replication dataset
list_ref_alt_rep <- c()
tmp_df <- list_ref_alt_disc[list_ref_alt_disc$type == "kmer", c("start", "end")]; tmp_df <- tmp_df[!duplicated(tmp_df), ]
for(i in 1:nrow(tmp_df)){ # for each index in GWAS table (of positions retained in discovery)
  # find position in GWAS_kmer of replication data:
  tmp_idx <- which(GWAS_kmer$start_position_alignment==tmp_df$start[i] & GWAS_kmer$end_position_alignment==tmp_df$end[i])
  if(length(tmp_idx) == 0) next else if(length(tmp_idx) > 1) stop()
  n_tested_alleles <- GWAS_kmer$df[tmp_idx]
  for(k in 1:n_tested_alleles){
    tested_allele <- GWAS_kmer[tmp_idx, paste0("alternative", k)]
    if(is.na(tested_allele)) stop("error tested allele is NA")
    list_ref_alt_rep <- rbind(list_ref_alt_rep, c(unlist(unname(GWAS_kmer[tmp_idx, c("start_position_alignment", "end_position_alignment", "reference", paste0("alternative", k))])), k))
  }
  rm(tmp_idx)
}
rm(tmp_df)
list_ref_alt_rep <- data.frame(list_ref_alt_rep) 

# now get these alleles in replication dataset and test correlation between predicted and true phenotype
if(nrow(list_ref_alt_rep) > 0){
  
  # clean list_ref_alt_rep table:
  names(list_ref_alt_rep) <- c("start", "end", "reference", "alternative", "allele_no")
  list_ref_alt_rep$allele_no <- as.numeric(list_ref_alt_rep$allele_no)
  list_ref_alt_rep$variants <- paste(list_ref_alt_rep$start, list_ref_alt_rep$end, list_ref_alt_rep$reference, list_ref_alt_rep$alternative, sep = "_")
  
  # now select only kmer matching those in list_ref_alt_disc
  idx_in_rep <- match(list_ref_alt_disc$variants, list_ref_alt_rep$variants); idx_in_rep <- idx_in_rep[!is.na(idx_in_rep)]
  idx_in_disc <- match(list_ref_alt_rep$variants, list_ref_alt_disc$variants);  idx_in_disc <- idx_in_disc[!is.na(idx_in_disc)]
  
  stopifnot(all((matching_variants <- list_ref_alt_rep[idx_in_rep, "variants"]) == list_ref_alt_disc[idx_in_disc, "variants"]))
  
  # look at the variants retained:
  print(list_ref_alt_rep[idx_in_rep, ])
  list_ref_alt_disc[idx_in_disc, ]
  
  # column vector of fixed effects retained:
  matching_kmer_fixed_effects <- my_h2$fixed_effects[match(rownames(list_ref_alt_disc[idx_in_disc, ]), rownames(my_h2$fixed_effects)), 1, drop = F]
  
  # and construct corresponding design matrix for each of the variants tested for replication
  # first find the index in dm (design matrix)
  matching_kmer_dm_start_end <- t(
    sapply(1:length(idx_in_rep), function(i){
      get_index_in_dm(st = list_ref_alt_rep[idx_in_rep[i], "start"], en = list_ref_alt_rep[idx_in_rep[i], "end"],
                      GWAS_tab = GWAS_kmer, full_design_matrix = full_design_matrix_kmer, full_ind_tokeep = full_ind_tokeep_kmer)
    })
  )
  # second construct design matrix
  matching_kmer_dm <- sapply(
    1:length(idx_in_rep),
    function(i) full_design_matrix_kmer[ , matching_kmer_dm_start_end[i, "idx_in_dm_start"]+list_ref_alt_rep[idx_in_rep[i], "allele_no"]-1, drop = F]
  )
  # replace missing values with 0 (?)
  matching_kmer_dm[is.na(matching_kmer_dm)] <- 0
  pred_genetic_effect <- matching_kmer_dm %*% matching_kmer_fixed_effects
  
  print(
    summary(lm0 <- lm(suba[, pheno_name] ~ pred_genetic_effect))
  )
  # CI for the regression coefficient:
  ci0 <- confint(lm0)
  # R2 and CI for R2
  R2 <- summary(lm0)$r.s
  boot_R2 <- quantile(
    sapply(1:1000, function(i){
      samp <- sample(1:nrow(suba), replace = T)
      summary(lm(suba[, pheno_name][samp] ~ pred_genetic_effect[samp]))$r.s
    }), c(0.025, 0.975))
  
  # alternative: Fisher test for association between pred genetic effect and phenotype
  tab <- table(pred_genetic_effect>mean(pred_genetic_effect), suba[, pheno_name]>mean(suba[, pheno_name], na.rm = T))
  print(tab)
  if(all(dim(tab) == c(2, 2))) print(ft <- fisher.test(tab)) else { tab <- c(NA, NA, NA, NA); ft <- list(p.value = NA)}
}
all_res <- rbind(
  all_res, c(args, lm0$coefficients["pred_genetic_effect"], ci0["pred_genetic_effect",],  summary(lm0)[[4]]["pred_genetic_effect", "Pr(>|t|)"], R2, boot_R2, c(t(tab)), ft$p.value)
)
}
#TODO FORMAT THIS TABLE AND SAVE
all_res <- data.frame(all_res)
names(all_res) <- c("phenotype", "p_threshold", "pred_genetic_effect",  "pred_genetic_effect_low",  "pred_genetic_effect_up", "p_value", "correlation_coefficient", "correlation_coefficient_low", "correlation_coefficient_up", "N_neg_neg", "N_neg_pos", "N_pos_neg", "N_pos_pos", "p_value_Fisher")
for(cc in c("pred_genetic_effect",  "pred_genetic_effect_low",  "pred_genetic_effect_up", "p_value",
            "correlation_coefficient", "correlation_coefficient_low", "correlation_coefficient_up",
            "p_value_Fisher")) all_res[, cc] <- signif(as.numeric(all_res[, cc]), 2)

all_res$polygenic_score_regression <- paste0(all_res$pred_genetic_effect, " [", all_res$pred_genetic_effect_low, "; ", all_res$pred_genetic_effect_up, "]")
all_res$polygenic_score_R2 <- paste0(all_res$correlation_coefficient*100, "% [", all_res$correlation_coefficient_low*100, "; ", all_res$correlation_coefficient_up*100, "%]")

all_res2 <- all_res[, c("phenotype", "p_threshold", "polygenic_score_regression", "polygenic_score_R2", "p_value")]
all_res2

write.csv(x = all_res2, file = "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/prediction_polygenic_score.csv", row.names = F)
