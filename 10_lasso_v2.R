#rm(list= ls())
args = commandArgs(trailingOnly=TRUE)

# An attempt at a lasso to predict phenotype (e.g. BEEHIVE_LVL) from genotype
# Initially, this did not work well (0 mutations were retained a significant)
# Solved by selecting only alleles not too rare, not too frequent

if(length(args)==0){
  print("warning length of args is 0, replacing by default")
  args <- c("spvl_adjusted") # "BEEHIVE_LVL", "spvl_normalised_adjusted", "CD4_slope"
  print(args)
  setwd("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/10_lasso/")
}
stopifnot(length(args)==1)
mypheno <- args[1]

######################################################
###### 0) LOAD DISCOVERY DATA and GWAS RESULT
######################################################

load(paste0(
   "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/Gmatrix_data/prepared_", mypheno, "_full_17.RData")
)
ls()->truc
load(paste0(
   "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/GWAS_results/firstGWAS_", mypheno, "_full_17.RData")
)
ls()->truc2
setdiff(truc2, truc) # variables added by the GWAS

# rename variables that will be reloaded in replication dataset as _disc and remove:
var_to_rename <- c(
  "blups", "fixed_effects", "full_design_matrix_kmer", "full_ind_tokeep_kmer",
  "GWAS_kmer", "GWAS_rare", "n_fixed_effects",
  "nb_alternative", "nb_gpos_per_ind", "nb_ind_per_gpos", "opt_pars",
  "pred", "unique_pos_to_test"
)
for(myvar in var_to_rename) assign(x = paste0(myvar, "_disc"), value = get(myvar))
rm(list = var_to_rename) # remove old ones

######################################################
###### 1) LOAD REPLICATION GWAS RESULT
######################################################

ls() -> truc
load(paste0(
  "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/GWAS_results/firstGWAS_", mypheno, "_full_17.RData")
)
ls()->truc2
setdiff(truc2, truc) 

# rename and remove
for(myvar in var_to_rename) assign(x = paste0(myvar, "_rep"), value = get(myvar))
rm(list = var_to_rename) # remove old ones


# correlations between phenotype
# cor(suba[, c("spvl_adjusted", "spvl_normalised_adjusted", "BEEHIVE_LVL")], use="complete.obs")


# 80%-20% split and examine the R^2 in remaining 20%
# added 18/03/2025
# load replication data and select only variants of discovery that can be replicated

length(reference)
length(alternative)
dim(Gn)

############################################################################################################
##### 2) create a consolidated list of variants present at at least N=10, and in replication DONE=list_match
############################################################################################################

# code variants as st_en_ref_alt (for first 10 alternatives)
n_alt_to_try <- min(10, sum(grepl(pattern = "alternative", colnames(GWAS_kmer_rep))))
for(i in 1:n_alt_to_try) {
  GWAS_kmer_disc[, paste0("var", i)] <- paste(GWAS_kmer_disc$start_position_alignment, GWAS_kmer_disc$end_position_alignment, GWAS_kmer_disc$reference, GWAS_kmer_disc[, paste0("alternative", i)], sep = "_")
  GWAS_kmer_disc[, paste0("var", i)][grepl(pattern = "_NA", x = GWAS_kmer_disc[, paste0("var", i)])] <- NA
}
for(i in 1:n_alt_to_try) {
  GWAS_kmer_rep[, paste0("var", i)] <-  paste(GWAS_kmer_rep$start_position_alignment,  GWAS_kmer_rep$end_position_alignment,  GWAS_kmer_rep$reference,  GWAS_kmer_rep[, paste0("alternative", i)], sep = "_")
  GWAS_kmer_rep[, paste0("var", i)][grepl(pattern = "_NA", x = GWAS_kmer_rep[, paste0("var", i)])] <- NA
}

# look at the match between discovery and replication data
list_match <- c()
for(i in 1:n_alt_to_try){
  for(j in 1:n_alt_to_try){
    cat(i, j, "\n")
    tmp <- GWAS_kmer_disc[, paste0("var", i)] # variants in discovery
    tmp_noNA <- tmp[non_NA_tmp <- which(!is.na(tmp))] # remove NA variant
    idx <- match(tmp_noNA, GWAS_kmer_rep[, paste0("var", j)]) # match to variants in replication
    idx <- idx[non_NA_idx <- which(!is.na(idx))]
    rm(tmp)
    if(length(idx)>0){ # if any variant is matched
      for(kk in 1:length(idx)) {
        if(GWAS_kmer_disc[non_NA_tmp[non_NA_idx[kk]], paste0("N", i)] >= 10 & GWAS_kmer_rep[idx[kk], paste0("N", j)] >= 10){ # keep only if this variant was tested in discovery GWAS AND IS TESTABLE IN REPLICATION GWAS
          list_match <- rbind(list_match, c(i, j, non_NA_tmp[non_NA_idx[kk]], idx[kk], GWAS_kmer_disc[non_NA_tmp[non_NA_idx[kk]], paste0("var", i)],  GWAS_kmer_rep[idx[kk], paste0("var", j)]))
        }
      }
    }
    #print(length(idx[!is.na(idx)]))
  }
}
list_match <- data.frame(list_match)
names(list_match) <- c("alt_nb_disc", "alt_nb_rep", "idx_disc", "idx_rep", "var_disc", "var_rep")
stopifnot(all(list_match$var_disc==list_match$var_rep))
dim(list_match) # previously was 80006, now with filtering on N >= 10 for rep, is 74154

dim(list_match)
table(paste(list_match$alt_nb_disc, list_match$alt_nb_rep))

for(mycol in c("alt_nb_disc", "alt_nb_rep", "idx_disc", "idx_rep")) list_match[, mycol] <- as.numeric(as.character(list_match[, mycol]))

############################################################################################################
##### # 3) identify in Gn (of discovery) where (which columns) the variants in list_match are
############################################################################################################

# use the function get_index_in_dm of 6_predict_phenotypes_additional

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

idx_in_dm_disc <- sapply(1:nrow(list_match), function(i)   get_index_in_dm(idx_in_GWAS = list_match$idx_disc[i], GWAS_tab = GWAS_kmer_disc, full_design_matrix = full_design_matrix_kmer_disc, full_ind_tokeep = full_ind_tokeep_kmer_disc))
idx_in_dm_disc <- t(idx_in_dm_disc)
idx_in_dm_rep <-  t(sapply(1:nrow(list_match), function(i) get_index_in_dm(idx_in_GWAS = list_match$idx_rep[i],  GWAS_tab = GWAS_kmer_rep,  full_design_matrix = full_design_matrix_kmer_rep,  full_ind_tokeep = full_ind_tokeep_kmer_rep)))

dim(full_ind_tokeep_kmer_disc)
dim(full_design_matrix_kmer_disc)
length(reference)
length(alternative)

# now get all columns to select in G
all_col_toselect_disc <- c()
all_col_toselect_rep <- c()
all_pos_in_list_match <- c() # position of each selected column in list_match

for(ii in 1:nrow(list_match)){
  #list_match[ii,] # the variant in discovery and replication
  #idx_in_dm_disc[ii,] # the start and end in full design matrix
  
    # for discovery:
    col_toselect <- (idx_in_dm_disc[ii,"idx_in_dm_start"]:idx_in_dm_disc[ii,"idx_in_dm_end"])[list_match[ii,"alt_nb_disc"]]
    tmp <- full_design_matrix_kmer_disc[, col_toselect] # the design matrix in discovery
    stopifnot(sum(tmp, na.rm = T) ==
                GWAS_kmer_disc[list_match[ii, "idx_disc"], paste0("N", list_match[ii,"alt_nb_disc"])]) # the number of alternative
    stopifnot(sum(tmp, na.rm = T) >= 10) # check again the condition greater than or equal to 10
    all_col_toselect_disc <- c(all_col_toselect_disc, col_toselect)
    
    # for replication
    col_toselect <- (idx_in_dm_rep[ii,"idx_in_dm_start"]:idx_in_dm_rep[ii,"idx_in_dm_end"])[list_match[ii,"alt_nb_rep"]]
    tmp <- full_design_matrix_kmer_rep[, col_toselect] # the design matrix in replication
    stopifnot(sum(tmp, na.rm = T) ==
                GWAS_kmer_rep[list_match[ii, "idx_rep"], paste0("N", list_match[ii,"alt_nb_rep"])]) # the number of alternative
    stopifnot(sum(tmp, na.rm = T) >= 10) # check again the condition greater than or equal to 10
    all_col_toselect_rep <- c(all_col_toselect_rep, col_toselect)
    
    all_pos_in_list_match <- c(all_pos_in_list_match, rep(ii, length(col_toselect))) # record position in list_match
  
}


############################################################################################################
##### 4 ) Do the lasso 
############################################################################################################

# but now apply it to the full design matrix of k-mer instead of G (G is the multi-allelic alignement in matrix form)
# TODO: maybe test if works better when applying only to 1-mer, no k-mer variants.

library(glmnet)
dim(full_ind_tokeep_kmer_disc)
dim(full_design_matrix_kmer_disc)
dim(suba)

fixed_covariate_names %in% names(suba)

# the final design matrix:
dm <- full_design_matrix_kmer_disc[,all_col_toselect_disc]
stopifnot(all(colSums(dm, na.rm = T)>=10))
cm <- colMeans(dm, na.rm = T)

all_R2_cv <- c()
all_R2_cv_v2 <- c()
all_coeff <- c() # to store the coefficients
all_int <- c() # to store the intercepts
all_R2_model <- c()

# assign NA values to mean frequency:
for(i in 1:ncol(dm)) dm[,i][which(is.na(dm[, i]))] <- cm[i]

# create design matrix for covariates
mylm <- paste(mypheno, " ~ ", paste(fixed_covariate_names, collapse = " + "))
mm <- model.matrix(lm(mylm, data = suba))[, -1]

# create matrix of covariate design matrix + phenotype
mydata <- cbind(
  mm, # adding mm (covariates) slightly improves the fit
  dm, # design matrix of genetic effects
  suba[, mypheno]
  )
stopifnot(sum(is.na(mydata))==0) # check there is no NA
mydata <- data.frame(mydata)
names(mydata)[ncol(mydata)] <- "status" # last one is phenotype

fraction_training <- 0.8
linear_model = "simple"
  
for(kk in 1:(nrep<-100)){ # several replicates of cutting the dataset in several bits
    
    cat("kk=", kk, "\n")
    #set.seed(123)
    idx.training.samples <- sort(sample(1:dim(mydata)[1], dim(mydata)[1]*fraction_training, replace=FALSE))
    train.data  <- mydata[idx.training.samples, ]
    test.data <- mydata[-idx.training.samples, ]
    
    # dummy code categorical predictor variables
    
    x <- as.matrix(mydata[idx.training.samples, 1:(ncol(mydata)-1)])
    # for(j in 1:ncol(x)){
    #   mymean <- mean(x[, j], na.rm = T)
    #   x[ ,j][is.na(x[, j])] <- mymean
    # }
    
    x.test <- as.matrix(mydata[-idx.training.samples, 1:(ncol(mydata)-1)])
    
    # convert the outcome (class) to a numerical variable
    y <- train.data[, "status"]
    
    # I found out that there is variability in best value of lambda across runs
    # hence doing 10 runs below to find best value
    # to handle the fact that the function cv.glmnet is stochastic
    best_cv.lasso <- NULL
    best_R2 <- -0.001 # negative value to ensure a model is retained even if R2 is 0 
    for(ll in 1:5){ # this does not need to be repeated (?). It does not improve R2 in validation
      print(ll)
      # find the optimal value of lambda that minimizes the cross-validation error:
      cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "gaussian")
      # Final model with lambda.min
      lasso.model_min <- glmnet(x, y, alpha = 1, family = "gaussian", lambda = cv.lasso$lambda.min)
      if(lasso.model_min$dev.ratio > best_R2){
        best_cv.lasso <- cv.lasso
        best_R2 <- lasso.model_min$dev.ratio
        print(best_R2)
      }
    } # best R2 is typically around 28%-35% for spvl
    cv.lasso <- NULL
    #plot(best_cv.lasso)
    #cv.lasso$lambda.min
    #cv.lasso$lambda.1se
    
    
    # Final model with lambda.min, lambda.1se
    lasso.model_min <- glmnet(x, y, alpha = 1, family = "gaussian", lambda = best_cv.lasso$lambda.min)
    lasso.model_1se <- glmnet(x, y, alpha = 1, family = "gaussian", lambda = best_cv.lasso$lambda.1se)
    
    # Save coefficients (lambda.min is best)
    all_coeff <- cbind(all_coeff, lasso.model_min$beta)
    all_int <- c(all_int, lasso.model_min$a0) # don't forget the intercept
    all_R2_model <- c(all_R2_model, best_R2)
    # quantify R2 on test data:
    if(fraction_training < 1){
      y.test_min <- predict(object = lasso.model_min, newx = x.test)
      # summary(as.numeric(lasso.model_min$a0 + x.test %*% lasso.model_min$beta - y.test_min)) # to check the way the prediction is computed
      y.test_1se <- predict(object = lasso.model_1se, newx = x.test)
      R2_min <- 1 - sum((y.test_min - test.data$status)^2, na.rm = T) / sum((mean(test.data$status) - test.data$status)^2)
      R2_min_v2 <- summary(lm(y.test_min ~ test.data$status))$r.s # for comparison with next one
      R2_1se <- 1 - sum((y.test_1se - test.data$status)^2, na.rm = T) / sum((mean(test.data$status) - test.data$status)^2)
    } else {
      y.train_min <- predict(object = lasso.model_min, newx = x)
      y.train_1se <- predict(object = lasso.model_1se, newx = x)
      R2_min <- 1 - sum((y.train_min - train.data$status)^2, na.rm = T) / sum((mean(train.data$status) - train.data$status)^2)
      R2_1se <- 1 - sum((y.train_1se - train.data$status)^2, na.rm = T) / sum((mean(train.data$status) - train.data$status)^2)
    }
    # save R2 cross-validation:
    all_R2_cv <- rbind(all_R2_cv,
                    c(R2_min, R2_1se) # lambda_min always works better than lambda_se;
    )
    all_R2_cv_v2 <- c(all_R2_cv_v2, R2_min_v2) # R2 when correcting with lm
} # end of loop on replicates
  
# }
colMeans(all_R2_cv)
# R2 in replication dataset, on average:
# 0.06687072 0.03182619 # with no repeats of the lasso, phenotype "spvl_adjusted"
# 0.06285352 0.03425500 with 10 repeats of lasso

# UPDATE with new method 04/04:
# 0.05685921 0.01845737
# 0.04910188 0.01838096
# 0.05728540 0.02484719

# UPDATE only variants at high-enough frequency in replication:
# 0.04286529 0.01497168
# 0.04631755 0.01919717
# 0.05393201 0.02237511

# TRYING ONLY THE VARIANTS SUCH THAT idx_in_dm_disc[ii,"idx_in_dm_start"] == idx_in_dm_disc[ii,"idx_in_dm_end"]
# 0.02167858 0.00466745

if(mypheno=="BEEHIVE_LVL") idx_start <- 11
if(mypheno=="CD4_slope") idx_start <- 6
if(mypheno %in% c("spvl_adjusted", "spvl_normalised_adjusted")) idx_start <- 8

all_coeff_g <- all_coeff[idx_start:nrow(all_coeff),]; stopifnot(nrow(all_coeff_g)==ncol(dm))

# Get genomic design matrix for replication
full_design_matrix_kmer_rep_sub <- full_design_matrix_kmer_rep[, all_col_toselect_rep]
save(list = c("mypheno", "suba", "fixed_covariate_names", "n_fixed_effects_disc", "list_match", "idx_start", "dm", "all_coeff", "all_int", "all_coeff_g", "nrep",
              "all_R2_model", "all_R2_cv", "all_R2_cv_v2",
              "full_design_matrix_kmer_rep_sub"), file = paste0("lasso_", mypheno, ".RData"))


quit()
rm(list = ls())


############################################################################################################
######## 5) LASSO ON REPLICATION DATASET ########
############################################################################################################

# setwd("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/10_lasso/")

# a) Load lasso
for(mypheno in list_of_pheno <- c("BEEHIVE_LVL", "spvl_adjusted", "spvl_normalised_adjusted", "CD4_slope")){
  
load(paste0("lasso_", mypheno, ".RData"))

ls()

dim(all_coeff)
dim(full_design_matrix_kmer_rep_sub)

# b) save all existing variables with special names
# old_names <- ls()
# new_names <- paste0(old_names, "_tmp123")
# for(i in 1:length(old_names)) assign(x = new_names[i], value = get(old_names[i]))
# rm(list = old_names)

# c) load replication dataset
load(paste0(
  "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/Gmatrix_data/prepared_", mypheno, "_full_17.RData")
)
ls()

dim(suba)

# d) add covariates to design matrix
suba$hepacivirus_present <- as.numeric(suba$hepacivirus_present>0)
suba$hepacivirus_present <- as.factor(suba$hepacivirus_present)
mylm <- paste(mypheno, " ~ ", paste(fixed_covariate_names, collapse = " + "))
mm <- model.matrix(lm(mylm, data = suba))[, -1]

if(mypheno=="BEEHIVE_LVL"){ # eliminate covariates that are not in discovery
  mm <- mm[, -which(colnames(mm) == "age.catOTH/UNK")]
  mm <- mm[, -which(colnames(mm) == "assay.gsvl2SEQ_VL")]
  mm <- mm[, -which(colnames(mm) == "SAMP_TYPENK")]
}
colnames(mm)

# create matrix of covariate design matrix + phenotype
dm_rep <- cbind(
  mm, # adding mm (covariates) slightly improves the fit
  full_design_matrix_kmer_rep_sub # design matrix of genetic effects
)

# e) replace NA with the mean
cm_rep <- colMeans(dm_rep, na.rm = T)
for(i in 1:ncol(dm_rep)) dm_rep[,i][which(is.na(dm_rep[, i]))] <- cm_rep[i]

# f) manually remove some covariates not present in replication ("age.cat.20..40.", "age.cat.60..80.")
if(mypheno=="BEEHIVE_LVL") sel_of_columns <- c(1, 4, 6, 7, 8, 9, 10) else sel_of_columns <- c(1, 3, 5, 6, 7)
stopifnot(length(sel_of_columns)==ncol(mm))
rownames(all_coeff)[1:14]
print("think more carefully about epi covariates")
print(cbind( # test that the covariates are matched in mm and in all_coeff
  rownames(all_coeff)[sel_of_columns],
  colnames(mm)
))
# coefficients from lasso
if(mypheno=="BEEHIVE_LVL") {
  all_coeff <- all_coeff[c(sel_of_columns, 11:nrow(all_coeff)), ]
} else {
  all_coeff <- all_coeff[c(sel_of_columns, 8:nrow(all_coeff)), ] # selecting only the covariates that are in replication ("sel_of_columns")
}

# g) predict phenotype in test data
all_lm <- list()
all_R2_test <- c()
all_R2_test_v2 <- c()
for(k in 1:nrep){ # across replication of the lasso
  print(k)
  predicted_pheno <- all_int[k] + c(as.matrix(dm_rep %*% all_coeff[, k, drop = F])) # intercept + coefficients
  pheno_topredict <- suba[, mypheno]
  # plot(predicted_pheno, pheno_topredict)
  all_lm[[k]] <- lm(predicted_pheno ~ pheno_topredict)
  all_R2_test <- c(all_R2_test, 1 - sum((predicted_pheno-pheno_topredict)^2)/sum((pheno_topredict - mean(pheno_topredict))^2)) # definition of R2 for test
  all_R2_test_v2 <- c(all_R2_test_v2, summary(all_lm[[k]])$r.s)
  #print(summary(all_lm[[k]]))
  #print(summary(all_lm[[k]])$r.s)
}

# Define variables useful for plotting:
n_effects <- apply(all_coeff_g, 2, function(vec) sum(abs(vec)>0))
# all_R2_test <- unlist(lapply(all_lm, function(mylm) summary(mylm)$r.s)) # former definition but this is not correct
all_R2_model
all_p_test <- unlist(lapply(all_lm, function(mylm) summary(mylm)[[4]][2,4]))

# save as special coeff with phenotype name appended
assign(x = paste0("n_effects_", mypheno), value = n_effects)
assign(x = paste0("all_R2_test_", mypheno), value = all_R2_test) # Rs additional data
assign(x = paste0("all_R2_test_v2_", mypheno), value = all_R2_test_v2) # Rs additional data
assign(x = paste0("all_R2_model_", mypheno), value = all_R2_model) # Rs lasso model
assign(x = paste0("all_R2_cv_", mypheno), value = all_R2_cv[,1]) # Rs cross-validation
assign(x = paste0("all_R2_cv_v2_", mypheno), value = all_R2_cv_v2) # Rs cross-validation
assign(x = paste0("all_p_test_", mypheno), value = all_p_test)

} # end of loop on phenotype

# Result of lasso
# 1) Number of loci retained
# 2) R^2 values in cross-validation / test
# 3) p-value in test

plot_vec <- function(vec_to_plot, bin_width, xpos, width_scaling_factor = 0.2, point.cex = 1, plot_median = T, log_transform_median = F, mycol ="gray"){
  
  # function to draw a nice violin plot
  # from a vector of values vec_to_plot
  if(plot_median){
      med <- median(vec_to_plot, na.rm = T)
      print(med)
  }
    
  #segments(y0 = med, y1 = med, x0 = xpos-width_scaling_factor, x1 = xpos+width_scaling_factor, col = "black", lwd = 3)

  rr <- range(vec_to_plot, na.rm = T)
  bins <- hist(vec_to_plot, breaks = seq(rr[1]-bin_width, rr[2]+bin_width, bin_width), plot = F)
  max_d <- max(bins$density)
  
  for(i in 1:length(bins$mids)){ # for each bin
    sub <- which(vec_to_plot >= bins$mids[i]-bin_width / 2 & vec_to_plot < bins$mids[i]+bin_width/2)
    width <- bins$density[i] / max_d * width_scaling_factor
    points(runif(n = length(sub), min = xpos-width, max = xpos+width), vec_to_plot[sub], pch = 20, col = mycol, cex = point.cex)
  }
  
  if(plot_median){
    segments(y0 = med, y1 = med, x0 = xpos-width_scaling_factor, x1 = xpos+width_scaling_factor, col = "black", lwd = 3)
    if(!log_transform_median){ # whether the y scale is log10-transformed or not
      text(x = xpos+width_scaling_factor, y = med, labels = signif(med, 2), adj = c(0,0))
    } else {
      text(x = xpos+width_scaling_factor, y = med, labels = signif(10^(med), 2), adj = c(0,0))
    }
  }
  
}
list_of_titles <- c("GSVL", "SPVL adj.", "SPVL norm. adj.", "CD4 slope")

pdf("lasso_results.pdf", width= 8, height = 8)

par(mar = c(4,4,1,1), mfrow = c(2, 2), xpd = T)
plot(NULL, xlim = c(0, 5), ylim = c(0,0.5), axes = F, xlab = "", ylab = "Coefficient of determination, training")
axis(side = 2, at = seq(0, 0.5, 0.1), las = 1)
# plot R2
idx <- 1
for(mypheno in list_of_pheno){
  plot_vec(vec_to_plot = get(paste0("all_R2_model_", mypheno)), bin_width = 0.01, xpos = idx)     # R2 in lasso MODEL (not cross-validation)
  idx <- idx + 1
}
text(1:4, -0.05, list_of_titles, srt = 45)
unlist(lapply(list_of_pheno, function(mypheno) signif(median(get(paste0("all_R2_model_", mypheno))), 2)))
text(0.1, 0.5, "A", cex = 3)

# TODO need to plot all_R2_cv which is R2 in lasso cross-validation!!!

plot(NULL, xlim = c(0, 9+1.5), ylim = c(0,0.12), axes = F, xlab = "", ylab = "Coefficient of determination, CV and BEEHIVE add.")
axis(side = 2, at = seq(0, 0.12, 0.02), las = 1)
idx <- 1
for(mypheno in list_of_pheno){
  if(mypheno=="CD4_slope") myplot_median <- F else myplot_median <- T
  plot_vec(vec_to_plot = get(paste0("all_R2_cv_v2_", mypheno)), bin_width = 0.01, xpos = idx, mycol = "cornflowerblue", plot_median = myplot_median)     # R2 in lasso cross-validation
  plot_vec(vec_to_plot = get(paste0("all_R2_test_v2_", mypheno)), bin_width = 0.01, xpos = idx+1, plot_median = myplot_median)     # R2 in BEEHIVE additional
  idx <- idx + 2.5 # a larger gap
}
# manual plotting for CD4 slope (otherwise not visible)
med <- 0.00074
segments(y0 = med, y1 = med, x0 = 8.5-0.2, x1 = 8.5+0.2, col = "black", lwd = 3)
text(x = 8.5-0.3, y = med+0.003, labels = med, adj = c(0,0))
med <- 0.00037
segments(y0 = med, y1 = med, x0 = 9.5-0.2, x1 = 9.5+0.2, col = "black", lwd = 3)
text(x = 9.5+0.2, y = med, labels = med, adj = c(0,0))

text(seq(1.5, 7.5, 2), -0.015, list_of_titles, srt = 45)
legend("topright", col = c("cornflowerblue", "gray"), pch= 20, legend = c("Cross-validation (CV)", "BEEHIVE additional"), bty = "n")
unlist(lapply(list_of_pheno, function(mypheno) signif(median(get(paste0("all_R2_cv_", mypheno))), 2)))
unlist(lapply(list_of_pheno, function(mypheno) signif(median(get(paste0("all_R2_test_", mypheno))), 2)))

text(0.1, 0.12, "B", cex = 3)

# plot p-values in replication
plot(NULL, xlim = c(0, 5), ylim = c(-4.5, 0), axes = F, xlab = "", ylab = "p-value BEEHIVE additional")
axis(side = 2, at = log10(c(0.0001, 0.001, 0.01, 0.05, 0.1, 1)), labels = c(0.0001, 0.001, 0.01, 0.05, 0.1, 1), las = 1)
idx <- 1
for(mypheno in list_of_pheno){
  plot_vec(vec_to_plot = log10(get(paste0("all_p_test_", mypheno))), bin_width = 0.1, log_transform_median = T, xpos = idx)
  idx <- idx + 1
}
text(1:4, -4.5, list_of_titles, srt = 45)
text(0.1, 0, "C", cex = 3)

# plot N effects
plot(NULL, xlim = c(0, 5), ylim = c(0,600), axes = F, xlab = "", ylab = "N effects retained in lasso")
axis(side = 2, at = seq(0, 600, 100), las = 1)
idx <- 1
for(mypheno in list_of_pheno){
  plot_vec(vec_to_plot = get(paste0("n_effects_", mypheno)), bin_width = 10, xpos = idx)
  idx <- idx + 1
}
text(1:4, -100, list_of_titles, srt = 45)
text(0.1, 600, "D", cex = 3)

dev.off()

unlist(lapply(list_of_pheno, function(mypheno) signif(median(get(paste0("n_effects_", mypheno))), 3)))


# Lvec coordinates is wrt 10257-long reference--the correspondance with HXB2 is in 
# GlobalAln_to_HXB2coords.csv

############################################################################################################
######## 6) CHECK AND SAVE THE VARIANTS RETAINED IN LASSO ########
############################################################################################################

hxb2_map <- read.csv("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/hxb2_annotated_Chris.csv", header = T)

for(mypheno in list_of_pheno){
  
  load(paste0("lasso_", mypheno, ".RData"))
 
  # look at effects across the replicates
  non_zero <- apply(all_coeff, 1, function(vec) sum(abs(vec)>0.001) >= 0.5*nrep) # at least 50% of replicates have this position
  non_zero_g <- apply(all_coeff_g, 1, function(vec) sum(abs(vec)>0.001) >= 0.5*nrep) # at least 50% of replicates have this position
  effects_nonzero_g <- apply(all_coeff_g[non_zero_g, ], 1, function(vec) median(vec[abs(vec)>0.001])) # median of non-zero effects
  
  table(non_zero_g)
  list_match$st_position_alignement <- sapply(as.character(list_match$var_disc), function(x) strsplit(x, split = "_")[[1]][1])
  list_match$en_position_alignement <- sapply(as.character(list_match$var_disc), function(x) strsplit(x, split = "_")[[1]][2])

  annot <- c() # create table of annotations
  if(any(non_zero_g)){
    
    idx_in_hxb2_map_st <- match(list_match[non_zero_g,"st_position_alignement"], hxb2_map$HXB2.base.position)
    idx_in_hxb2_map_en <- match(list_match[non_zero_g,"en_position_alignement"], hxb2_map$HXB2.base.position)
    for(i in 1:sum(non_zero_g)){
      n_nuc <- idx_in_hxb2_map_en[i]-idx_in_hxb2_map_st[i]+1 # number of nucleotides in this pos
      
      # CTL annotation at these positions:
      any_CTL <- any(grepl(pattern = "CTL", x = hxb2_map[idx_in_hxb2_map_st[i]:idx_in_hxb2_map_en[i], "CTL.escape."]))
      CTL_annot <- unique(hxb2_map[idx_in_hxb2_map_st[i]:idx_in_hxb2_map_en[i], "CTL.escape."])
      CTL_annot <- CTL_annot[CTL_annot!=""]
      CTL_annot <- paste(CTL_annot, collapse = "/")
      
      # position/aa at these positions:
      hxb2_map[idx_in_hxb2_map_st[i]:idx_in_hxb2_map_en[i],  c("RF1.protein", "RF1.aa.position", "RF1.aa", "RF2.protein", "RF2.aa.position", "RF2.aa", "RF3.protein", "RF3.aa.position", "RF3.aa")] -> truc
      truc <- apply(truc, 1, function(vec) {vec <- vec[!is.na(vec)]; vec <- vec[vec!=""]; return(paste(vec, collapse = " "))})
      truc2 <- paste(unique(truc), collapse = ", ")
      annot <- rbind(annot, 
                     cbind(hxb2_map[idx_in_hxb2_map_st[i]:idx_in_hxb2_map_en[i], ],
                           matrix(rep(list_match[which(non_zero_g)[i], "var_disc"], n_nuc), ncol = 1), # variant
                           matrix(rep(effects_nonzero_g[i], n_nuc), ncol = 1), # inferred effect
                           matrix(rep(any_CTL, n_nuc, ncol = 1)), # any CTL
                           matrix(rep(CTL_annot, n_nuc, ncol = 1)), #  CTL annotation
                           matrix(truc, ncol = 1), # prot
                           matrix(rep(truc2, n_nuc, ncol = 1))
                           )
                     )
    }
  } else {
    print("no effect found in > 50% of splits")
  }
  names(annot)[17:22] <- c("variant", "median_effect", "any_CTL", "CTL", "protein/aa", "protein/aa cat")
  annot <- annot[order(annot$HXB2.base.position), ]
  annot$median_effect <- signif(annot$median_effect, 3)
  
  annot2 <- annot[, c("variant", "protein/aa cat", "median_effect", "CTL")]
  annot2 <- annot2[!duplicated(annot2),]
  
  write.csv(x = annot2, file = paste0("lasso_summary_", mypheno, ".csv"), row.names = F, col.names = T)
  assign(x = paste0("annot2_", mypheno), value = annot2)

}


table(grepl(pattern = "CTL escape", x = hxb2_map$CTL.escape.))
table(annot2$CTL=="")

table(annot$median_effect>0, annot$any_CTL)

freqs_nonzero <- cm[non_zero_g]
plot(effects_nonzero_g, log10(freqs_nonzero), pch = 20, ylab = "log10 frequency", xlab = "Effect")
abline(v = 0); abline(h=-1)
table(neg_eff = effects_nonzero_g < 0, small_freq = log10(freqs_nonzero)< -1)
