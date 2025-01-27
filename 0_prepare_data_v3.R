rm(list = ls())

# this code prepares epi and genomic data for likelihood inference
# output a prepared_data.RData file
# VERSION: ALL SOURCED FILES WERE UPDATED AS OF 18/03/2019

# this code takes as input the alignmnents, the epi data, the subtype data, the dual infection data
# for each phenotype, we select the intersection of epi data and genomic data that can be included in the GWAS
# for epi data, we remove missing covariates (normally none; those are coded as 'unknown' categorical variable), missing phenotype, duplicates, dual infections, non-beehive OK, failed beehive viral load 
# for genomic data, we remove individuals where more than 50% positions are missing.
# Indeed, we previously found that heritability of BEEHIVE_LVL was much higher (suggesting the population structure correction is better) when leaving those out.
# We save all we need for the GWAS in one .RData file per phenotype

# TODO INCORPORATE PHENOTYPIC RESISTANCE INFO

######                              INTRO                           #####

#wd <- "~/DID/BEEHIVE_Hackathon/"
wd <- "/Users/Christophe/Dropbox (Infectious Disease)/PROJECTS/HIV/BEEHIVE_Hackathon/"
if(!file.exists(wd)) wd <- "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/"
setwd(paste0(wd, "Code/DevelopMethods/LMM/SolvingLMMinR/"))
library(seqinr)
library(Matrix)
source(paste0(wd, "Code/DevelopMethods/LMM/SolvingLMMinR/0_createG.R"))

######                          SOURCES FILES TO BE UPDATED IF NEEDED                       #####

aln_file_list <- c(
  paste0(wd, "Data/Sequences/GlobalAlignment/Clean/",
         paste0("2019-11-18_", c(
           "BEEHIVE_GlobalAln_clean_A12.fasta", "BEEHIVE_GlobalAln_clean_A1U.fasta", "BEEHIVE_GlobalAln_clean_A23.fasta",
           "BEEHIVE_GlobalAln_clean_A2U.fasta", "BEEHIVE_GlobalAln_clean_A34.fasta", "BEEHIVE_GlobalAln_clean_A3U.fasta", "BEEHIVE_GlobalAln_clean_A4U.fasta",
           "BEEHIVE_GlobalAln_clean.fasta"
         ))
  )
)
aln_file_list <- as.list(aln_file_list)
names(aln_file_list) <- c("A12", "A1U", "A23", "A2U", "A34", "A3U", "A4U", "full")
#prot_file <- paste0(wd, "Data/Sequences/TypewriterProteinAlignments_HXB2/concatenated_PRO.fasta") # now dealt with separately in prepare_aminoacids.R
epi_file <- paste0(wd, "Data/MetaData/CleanMetaData/BEEHIVE_summary17.10.2024.csv") # FB updated 17.10.2024 to use more recent file
sub_file <- paste0(wd, "Data/subtype/subtype_2019-11-18_BEEHIVE_GlobalAln_clean_final.csv") # SUBTYPE IS PRODUCED BY THE COMET SOFTWARE. R function is ~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/subtyping/subtype.R
dual_file <- paste0(wd, "Data/duals/Clean/PerPatientDualClassification_2019-03-12.csv")     # all suba patients are in there   
drm_file <- paste0(wd, "Data/Final_Resistance_consensus/Final_DR_Beehive_consensus.csv")    # 280 individuals missing from DRM as of 2019-12-06       

stopifnot(all(c(unlist(lapply(aln_file_list, file.exists)), file.exists(sub_file), file.exists(dual_file), file.exists(drm_file)))) # check that all files exist


######                   LOOP ON PHENOTYPE                    #####

for(pheno_name in c(
  "spvl_normalised_adjusted",
  "spvl_adjusted",
  "BEEHIVE_LVL",
  "CD4_slope"
  )){
  print(pheno_name)
  
  ######                  VARIABLES TO DEFINE                 #####
  
  # pheno_name <- "spvl_normalised_adjusted"
  # do not include transmission mode (no effect & correlated with genotype / subtype)
  # coinfection, age of infection, gender
  # for GSVL, time from infection to sampling
  
  fixed_covariate_names <- c("GENDER", "age.cat", "pegivirus_present", "hepacivirus_present")
  fixed_covariate_type <- c("categorical", "categorical", "categorical", "categorical")
  
  create_within_error <- NA
  if(pheno_name == "BEEHIVE_LVL"){
    rand_covariate_names <- NULL
    fixed_covariate_names <- c(fixed_covariate_names, "Time.to.sample", "SAMP_TYPE")
    fixed_covariate_type <- c(fixed_covariate_type, "continuous", "categorical")
    create_within_error <- function(rand_covariate_names, mydim, epi_data) return(NULL)
  }
  if(pheno_name == "CD4_slope"){
    rand_covariate_names <- c("CD4_slope_var")
    create_within_error <- function(rand_covariate_names, mydim, epi_data){ # define function to generate the matrices corresponding to error
      return(diag(mydim) * epi_data[, rand_covariate_names])
    }
  }
  if(pheno_name == "spvl_normalised_adjusted" | pheno_name == "spvl_adjusted") {
    rand_covariate_names <- c("N.SPVL", pheno_name) # FB as of 23/10/2024, don't understand any more why I appended "pheno_name" to rand_covariate_names
    create_within_error <- function(rand_covariate_names, mydim, epi_data){ # define function to generate the variance-covariance matrices corresponding to error
      return(diag(mydim) / epi_data[, rand_covariate_names[1]])
    }
  }
  if(is.na(create_within_error)) stop() # if create_within_error has not been defined stop
  
  n_fixed_covariates <- length(fixed_covariate_names)
  
  ######                  LOAD EPI DATA                 #####
  
  a <- read.csv(epi_file, stringsAsFactors = F)
  for(colname in fixed_covariate_names[fixed_covariate_type=='categorical']){
    cat("fixed covariate", colname, " replacing NA with OTHUNK for ", sum(is.na(a[, colname])), "datapoints\n")
    a[, colname][is.na(a[, colname])] <- "OTHUNK"
  }
  
  ######                        SET HEPACIVIRUS 0.5 TO 1                     #####
  
  a[, "hepacivirus_present"][a[, "hepacivirus_present"]=="0.5"] <- "1"
  table(a$hepacivirus_present)
  
  ######                        REMOVE MISSING COVARIATES                     #####
  
  tokeep_epi <- rep(T, nrow(a))
  for(colname in fixed_covariate_names) tokeep_epi <- tokeep_epi & !is.na(a[, colname])
  for(colname in rand_covariate_names)  tokeep_epi <- tokeep_epi & !is.na(a[, colname])
  
  # normally this does not remove any of the BEEHIVE.OK2 patients for BEEHIVE_LVL
  if(pheno_name=="BEEHIVE_LVL") stopifnot(all(tokeep_epi[a$BEEHIVE.OK2]))
  
  ######                              REMOVE THE MINOR                        #####
  
  tokeep_epi <- tokeep_epi & (a$Year.pos - a$YOB != 3)
  
  
  ######    REMOVE MISSING PHENOTYPES, DUPLICATES       #####
  
  tokeep_epi <- tokeep_epi & !is.na(a[, pheno_name]) # remove patients with missing focal phenotype
  tokeep_epi <- tokeep_epi & !a$is_duplicate         # remove duplicates
  
  ######                       REMOVE DUALLY INFECTED PEOPLE                         #####
  
  dual <- read.csv(dual_file) # proportion dual  sum(dual$is.dual, na.rm = T)/(sum(dual$is.dual, na.rm = T) + sum(!dual$is.dual, na.rm = T))
  a$is.dual <- a$PATIENT %in% dual$PATIENT[dual$is.dual]
  tokeep_epi <- tokeep_epi & !a$is.dual
  
  ######                       REMOVE NON-BEEHIVE-OK PATIENTS FOR ALL PHENOTYPES                #####
  
  tokeep_epi <- tokeep_epi & a$BEEHIVE.OK2 # keep patients 0-24 months

  ######                          REMOVE THE FAILED BEEHIVE_LVL                      #####
  
  if(pheno_name == "BEEHIVE_LVL") tokeep_epi <- tokeep_epi & (a$BEEHIVE_LVL > 2 | (a$BEEHIVE_LVL == 2 & a$SPVL < 3)| is.na(a$SPVL))
  
  ids_tokeep_epi <- a$PATIENT[tokeep_epi]
  
  ######                          REMOVE THE EXTREME CD4 SLOPE VALUES                      #####
  #if(pheno_name == "CD4_slope") tokeep_epi <- tokeep_epi & a$CD4_slope
  
  ######                  LOAD SUBTYPE DATA                 #####
  
  sub_table <- read.csv(sub_file, header = T)
  sub_table$PATIENT <- sapply(sub_table$name, get.bee.id)
  a$sub <- sub_table$final_subtype[match(a$PATIENT, sub_table$PATIENT)]
  #table(a$sub)
  
  ######                  LOAD DRM DATA                 #####
  
  drm <- read.csv(drm_file, header = T)
  drm$PATIENT <- sapply(drm$sequence, get.bee.id)
  all_R <- c("INSTI", "NNRTI", "NRTI", "PI") # phenotypic resistance names
  for(n in all_R){
    a[, n] <- as.character(drm[match(a$PATIENT, drm$PATIENT), n])
  }
  for(colname in all_R){
    a[, colname][is.na(a[, colname])] <- "OTHUNK"
    a[, colname][a[, colname] == "I" | a[, colname] == "R"] <- 1
    a[, colname][a[, colname] == "S"] <- 0
  }
  # generate pairwise dual resistances:
  combine_R <- function(type1, type2){ # combine 2 types of resistance into multitypes
    new_vec <- rep(0, nrow(a))
    new_vec[a[, type1]=="OTHUNK" | a[, type2]=="OTHUNK"] <- "OTHUNK"
    new_vec[(a[, type1] == 1) & (a[, type2] == 1)] <- 1
    return(new_vec)
  }
  Rcomb <- combn(x = all_R, 2)
  for(i in 1:ncol(Rcomb)){
    type1 <- Rcomb[1,i]
    type2 <- Rcomb[2,i]
    a[, paste(type1, type2, sep = "_")] <- combine_R(type1, type2)
  }
  
  ######                  LOAD AND PREPARE GENETIC DATA                 #####
  for(idx_aln in 8:8){#1:length(aln_file_list)){
  
    aln_file <- aln_file_list[[idx_aln]]
    print(names(aln_file_list)[idx_aln])
    aln <- read.fasta(aln_file, as.string = T)
    names(aln) <- sapply(names(aln), get.bee.id)
    stopifnot(all(!duplicated(names(aln))))
    GSigma <- create_G_sigma(aln = aln, ids_tokeep = ids_tokeep_epi)
    
    # unwrap GSigma:
    G <- GSigma$Gmat
    Sigma <- GSigma$Sigma
    idsG <- GSigma$names
    Lvec <- GSigma$Lvec
    reference <- GSigma$reference
    alternative <- GSigma$alternative
    n_missing_per_ind <- GSigma$n_missing_per_ind
    nindG <- length(GSigma$names)
    ngposG <- ncol(G)               # number of polymorphic positions in the full alignment
    npos_aln <- nchar(aln[[1]])   # length of the alignment
    unique_pos_to_test <- unique(Lvec)
    
    # ######                                  CREATE THE PROTEIN MATRIX AND SAVE IT (now done in prepare_aminoacids.R)                              #####
    # 
    # if(names(aln_file_list)[idx_aln] == "full"){
    #   prot_aln <- read.fasta(prot_file, as.string= T)
    #   prot_mat_list <- create_protein_mat(prot_aln = prot_aln, ids_tokeep = ids_tokeep_epi) # TO CHECK HOW IT DEALS WITH N AND GAP; WHAT IS THE SYMBOL FOR MISSING DATA VERSUS GENUINE GAP
    #   idsprot <- prot_mat_list$names
    #   idx_to_remove <- which(!idsprot %in% idsG) # sequence presents in protein but not in genetic data
    #   
    #   prot_mat <- prot_mat_list$prot_mat[-idx_to_remove, ]
    #   idsprot <- idsprot[-idx_to_remove]
    #   Lvec_prot <- prot_mat_list$Lvec
    #   reference_prot <- prot_mat_list$reference
    #   alternative_prot <- prot_mat_list$alternative
    #   nind_prot <- length(idsprot)
    #   ngpos_prot <- ncol(G)
    #   
    #   stopifnot(all(idsprot %in% idsG))
    #   filename <- paste0(wd, "Code/DevelopMethods/LMM/SolvingLMMinR/Gmatrix_data/prepared_", pheno_name, "_prot.RData")
    #   save(list = c("prot_mat", "idsprot", "Lvec_prot", "reference_prot", "alternative_prot", "nind_prot", "ngpos_prot"), file = filename)
    # }
    
    ######                                  CREATE THE K MATRIX WITH FILTERING                       #####

    missing50pc <- apply(G, 2, function(col) sum(is.na(col))) / nindG < 0.5 # positions where less than 50% of individuals miss data
    ref_alt_gap <- reference != "-" & alternative != "-" # positions where reference and alternative is NOT a gap
    nb_alternativeG <- colSums(G, na.rm =T) # number of the alternative allele
    #missing_half <- apply(G, 1, function(row) sum(is.na(row))) < ngposG/2. # individuals with less than half positions missing
    missing_half <- n_missing_per_ind < npos_aln / 2                      # individuals with less than half positions missing, based on the alignment almost identical to previous version
    
    crit_tab <- expand.grid(list(pos_missing50pc = c(0, 1), pos_ref_alt_gap = c(0, 1), pos_alt_threshold = c(0, 10, 50, 100), ind_missing_half = c(0,1)))
    crit_tab$suffix_name <- paste0("posMissing50pc", crit_tab$pos_missing50pc, "_posGap", crit_tab$pos_ref_alt_gap, "_alt", crit_tab$pos_alt_threshold, "_indMissing", crit_tab$ind_missing_half)
    stopifnot(all(!duplicated(crit_tab$suffix_name)))
    
    #write.csv(x = crit_tab, file = paste0(wd, "Code/DevelopMethods/LMM/SolvingLMMinR/criteria_table.csv"), row.names = F)
    
    for(ll in 17:17){ # LOOP ON COMBINATIONS OF CRITERIA 1:nrow(crit_tab)
      
      # subset of positions
      sub_pos <- rep(TRUE, ngposG)
      if(crit_tab$pos_missing50pc[ll]) sub_pos <- sub_pos & missing50pc
      if(crit_tab$pos_ref_alt_gap[ll]) sub_pos <- sub_pos & ref_alt_gap
      sub_pos <- sub_pos & nb_alternativeG > crit_tab$pos_alt_threshold[ll]
      
      sub_ind <- rep(TRUE, nindG)
      if(crit_tab$ind_missing_half[ll]) sub_ind <- sub_ind & missing_half
      
      cat("criterion ", ll, " number of pos ", table(sub_pos), '\n') 
      cat("criterion ", ll, " number of ind ", table(sub_ind), '\n')
      
      # the fast option
      Gn <- normalise_colmeans_mat(G[sub_ind, sub_pos]) # normalised, selected matrix
      Sigma_sub <- as(Sigma[sub_pos, sub_pos], "dgCMatrix")
      Gn <- as(Gn, "dgeMatrix")
      ngpos <- ncol(Gn); stopifnot(ngpos==sum(sub_pos))
      
      # deal with missing data:
      Gn_nastatus <- !is.na(Gn) # record non-NA positions
      na_pos <- which(is.na(Gn), arr.ind = T)
      Gn[na_pos] <- 0 # replace NA with 0 (the population mean because we normalised)
      n_non_na <- Gn_nastatus %*% t(Gn_nastatus) # matrix giving the number of non-NA genetic values for each pair of sequence
      
      t0 <- Sys.time()
      K_firstpart <- Sigma_sub %*% t(Gn) # this term is going to be used for the BLUPs estimation
      K <- Gn %*% K_firstpart
      #rm(K_firstpart)
      Sys.time() - t0
      
      # optional correction for K based on missing data: NOTE DEFINITION OF K_FIRSTPART MUST BE CHANGED AS WELL IF WE DO THAT
      #K <- K * ncol(Gn) / n_non_na
      
      ######                  CHECK THE EPI AND G DATA AND REMOVE DUPLICATES                #####
      
      stopifnot(all(idsG %in% a$PATIENT))
      #suba <- a[a$PATIENT %in% idsG[sub_ind], ]
      suba <- a[match(idsG[sub_ind], a$PATIENT), ] # first match
      if(any(duplicated(suba$PATIENT))) stop('duplicated patients in epi data')
      if(any(duplicated(idsG))) stop('duplicated patients in genetic data') # TO DO DEAL WITH DUPLICATES IN GENETIC DATA WHEN PRESENT
      stopifnot(all(suba$PATIENT==idsG[sub_ind]))
      nind <- nrow(suba)
      
      ######                  DEFINE PCs per subtype                 #####
      print("generating PCs...")
      sort(table(suba$sub), decreasing = T)
      sub_to_handle <- c("A1", "B", "02_AG", "01_AE", "C", "D")
      for(subtype in sub_to_handle){ # for each subtype
        unassigned_subtype <- paste0("unassigned", subtype)
        if(subtype == "B"){
          subtype_subset <- which(suba$sub == subtype | suba$sub == unassigned_subtype | suba$sub == "unassigned12_BF") # unassigned12_BF are probably B
          } else {
            subtype_subset <- which(suba$sub == subtype | suba$sub == unassigned_subtype)
          }
        pca <- prcomp(Gn[subtype_subset, ])
        assign(paste0("PC", subtype), pca$x) # compute PCs and assign to "PC" variable
        assign(paste0("PCsd", subtype), pca$sdev)
        assign(paste0("subtype_", subtype), subtype_subset)
      }
      other <- which((!suba$sub %in% sub_to_handle) & (!suba$sub %in% paste0("unassigned", sub_to_handle)) & suba$sub != "unassigned12_BF")
      pca <- prcomp(Gn[other, ])
      assign(paste0("PC", "other"), pca$x) # compute PCs and assign to "PC" variable
      assign(paste0("PCsd", "other"), pca$sdev)
      assign(paste0("subtype_", "other"), other)
      stopifnot(nrow(PCA1) + nrow(PCB) + nrow(PC02_AG) + nrow(PC01_AE) + nrow(PCC) + nrow(PCD) + nrow(PCother) == nrow(suba))
      rm(pca) # remove so we can see if there's a problem in calculation later on
      
      #suba_withPC <- cbind(suba[subtype_B,], PCB)
      #write.csv(x = suba_withPC, file = "~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/Gmatrix_data/epi_withPC_subtypeB.csv", row.names = F, col.names = T)
      
      
      ######                  DEFINE THE LINEAR MODEL                #####
      
      linear.model <- paste(pheno_name, "~")
      for(i in 1:n_fixed_covariates){ # AS OF 02/04/2019 CHECK IT IS CATEGORICAL
        covariate <- fixed_covariate_names[i]
        if(fixed_covariate_type[i]=="categorical"){
          tab.cov <- table(suba[, covariate])
          tab.cov <- tab.cov[tab.cov > 0]
          if(length(tab.cov) > 1 & all(tab.cov > 2)){
            if(i == 1) linear.model <- paste(linear.model, covariate) else linear.model <- paste(linear.model, "+", covariate)
          }
        }
        if(fixed_covariate_type[i]=="continuous"){# continuous, always add it
          if(i == 1) linear.model <- paste(linear.model, covariate) else linear.model <- paste(linear.model, "+", covariate)
        }
      }
      if(substr(linear.model, nc <- nchar(linear.model), nc) == "~") linear.model <- paste(linear.model, "1")
      linear.model
      lm0 <- lm(linear.model, data = suba)
      X <- model.matrix(lm0)
      X <- as(X, "dgCMatrix") # convert into the right format
      
      if(all(c("pegivirus_presentOTHUNK",  "hepacivirus_presentOTHUNK") %in% colnames(X))){
        if(all(X[, "pegivirus_presentOTHUNK"] == X[, "hepacivirus_presentOTHUNK"])){ # 'unknown' viral infection status is identical for the two virus data when we have not yet updated the data
          X <- X[, -match("hepacivirus_presentOTHUNK",colnames(X))]
        }
      }
      
      ######                  DEFINE A LIST XRAND WITH THE RANDOM EFFECT STRUCTURE    #####
      
      Kmat <- as.matrix(K)                                                                                               # genetic structure
      Imat <- diag(nind)                                                                                                 # unstructured variance
      withinmat <- create_within_error(rand_covariate_names = rand_covariate_names, mydim = nind, epi_data = suba)       # measurement error
      if(!is.list(withinmat)) Xrand <- list(Kmat = Kmat, Imat = Imat, withinmat = withinmat)
      if(!exists("Xrand")) stop('Xrand has not been defined')
      Xrand <- Xrand[!unlist(lapply(Xrand, is.null))]
      normalise_factors <- vector(length = length(Xrand))
      for(i in 1:length(Xrand)){
        tmp <- normalise_mat_diag(Xrand[[i]], return_factor = T)
        Xrand[[i]] <- tmp$mat
        normalise_factors[i] <- tmp$factor
      }
      Xrand <- lapply(Xrand, function(x) as(x, "dgeMatrix"))
      
      ######                               DEFINE OTHER VARIABLES                     #####
      
      nb_ind_per_gpos <- apply(G[sub_ind,], 2, function(x) sum(!is.na(x))) # nb of individuals with data at each generalised position
      nb_gpos_per_ind <- apply(G[sub_ind,], 1, function(x) sum(!is.na(x))) # nb of generalised positions with data for each individual
      nb_alternative <- colSums(G[sub_ind, ], na.rm =T)
      
      #Xrand <- lapply(Xrand, normalise_mat_diag)
      
      ######                  SAVE EVERYTHING WE NEED FOR THE LIKELIHOOD                #####
      filename <- paste0(wd, "Code/DevelopMethods/LMM/SolvingLMMinR/Gmatrix_data/prepared_", pheno_name, "_", names(aln_file_list)[idx_aln], "_", ll, ".RData")
      save(list = c("Xrand", "normalise_factors", "aln", "G", "Gn", "Sigma_sub", "ids_tokeep_epi", "idsG", "Lvec", "sub_ind", "sub_pos", "reference", "alternative",
                    "suba", "lm0", "X", "nind", "pheno_name", "fixed_covariate_names", "fixed_covariate_type", "n_fixed_covariates",
                    "ngpos", "nb_ind_per_gpos", "nb_gpos_per_ind", "nb_alternative", "unique_pos_to_test",
                    paste0("PC", sub_to_handle), "PCother", paste0("subtype_", sub_to_handle), "subtype_other", paste0("PCsd", sub_to_handle), "PCsdother"
                    ), file = filename)
      
    } # end of loop on criteria
  } # end of loop on alignment
} # end of loop on phenotype


