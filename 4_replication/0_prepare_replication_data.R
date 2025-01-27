rm(list = ls())

# this code prepares epi and genomic data for likelihood inference
# this code is specific to the replication dataset
# output a prepared_data.RData file
# THIS IS SIMILAR TO 0_prepare_data_v3.R EXCEPT WE DO NOT HAVE ALL UPDATED FILES (DUAL INFECTION AND DRM MUTATIONS ARE MISSING)

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
library(plyr)
source(paste0(wd, "Code/DevelopMethods/LMM/SolvingLMMinR/0_createG.R"))

######                          SOURCES FILES TO BE UPDATED IF NEEDED                       #####

aln_file_list <- c(
  paste0(wd, "Data/Sequences_replication/",
         c(
           "BEEHIVE_Oxford_GlobalAln_2019-11-06_MultiSeqsPerPat_clean.fasta"
         )
  )
)
aln_file_list <- as.list(aln_file_list)
names(aln_file_list) <- c("full")
#prot_file <- paste0(wd, "Data/Sequences/TypewriterProteinAlignments_HXB2/concatenated_PRO.fasta") # now dealt with separately in prepare_aminoacids.R
epi_file <- paste0(wd, "Data/MetaData/CleanMetaData/BEEHIVE_summary24.10.2019.csv")
sub_file <- paste0(wd, "Data/subtype/subtype_BEEHIVE_Oxford_GlobalAln_2019-11-06_MultiSeqsPerPat_clean_final.csv") # SUBTYPE IS PRODUCED BY THE COMET SOFTWARE. R function is ~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/subtyping/subtype.R
dual_file <- paste0(wd, "Data/duals/Clean/PerPatientDualClassification_2019-03-12.csv")     # all suba patients are in there   
drm_file <- paste0(wd, "Data/Final_Resistance_consensus/Final_DR_Beehive_consensus.csv")    # 280 individuals missing from DRM as of 2019-12-06       
discovery_patients_file <- paste0(wd, "Code/DevelopMethods/LMM/SolvingLMMinR/Gmatrix_data/all_individuals_discovery.csv")

stopifnot(all(c(unlist(lapply(aln_file_list, file.exists)), file.exists(sub_file), file.exists(dual_file), file.exists(drm_file), file.exists(discovery_patients_file)))) # check that all files exist


######                   LOOP ON PHENOTYPE                    #####

all_replication_sequences_tokeep <- c() # a list of all sequences we will keep

for(pheno_name in c(
  "spvl_normalised_adjusted",
  "spvl_adjusted",
  "BEEHIVE_LVL",
  "CD4_slope"
  )){
  #pheno_name <- "BEEHIVE_LVL"
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
    fixed_covariate_names <- c(fixed_covariate_names, "Time.to.sample", "SAMP_TYPE", "assay.gsvl2")
    fixed_covariate_type <- c(fixed_covariate_type, "continuous", "categorical", "categorical")
    create_within_error <- function(rand_covariate_names, mydim, epi_data) return(NULL)
  }
  if(pheno_name == "CD4_slope"){
    rand_covariate_names <- c("CD4_slope_var")
    create_within_error <- function(rand_covariate_names, mydim, epi_data){ # define function to generate the matrices corresponding to error
      return(diag(mydim) * epi_data[, rand_covariate_names])
    }
  }
  if(pheno_name == "spvl_normalised_adjusted" | pheno_name == "spvl_adjusted") {
    rand_covariate_names <- c("N.SPVL", pheno_name)
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
  disc_patients <- read.csv(discovery_patients_file, header = T)$x
  
  ######                       REMOVE NON-BEEHIVE-OK PATIENTS FOR ALL PHENOTYPES                #####
  
  a <- a[a$BEEHIVE.OK2, ] # we now remove straight away the non-OK patients
  a <- a[!is.na(a[, pheno_name]), ] # and we remove straight away missing phenotypes

  ######         NEW TO THE REPLICATION DATASET -- REMOVE PATIENTS THAT ARE ALREADY IN DISCOVERY                #####
  
  a <- a[!a$PATIENT %in% disc_patients, ]
  
  ######                  CREATE PATIENT_SEQ_DATE VARIABLE                 #####
  
  a$PATIENT_SEQ_DATE <- paste(a$PATIENT, a$SEQ_DATE, sep = "_")
  a$SAMPLE_SEQ_DIFF <- difftime(a$SEQ_DATE, a$Date.sampled, units = "days")
  a <- a[with(a, order(PATIENT, SAMPLE_SEQ_DIFF)), ] # reorder by patient then how close is sample and seq date
  
  ######                        REMOVE MISSING COVARIATES                     #####
  
  a$tokeep_epi <- T
  for(colname in fixed_covariate_names) a$tokeep_epi <- a$tokeep_epi & !is.na(a[, colname])
  for(colname in rand_covariate_names)  a$tokeep_epi <- a$tokeep_epi & !is.na(a[, colname])

  ######                       REMOVE DUALLY INFECTED PEOPLE                         #####
  
  dual <- read.csv(dual_file)
  a$is.dual <- a$PATIENT %in% dual$PATIENT[dual$is.dual]
  a$tokeep_epi <- a$tokeep_epi & !a$is.dual
  
  ######                          REMOVE THE FAILED BEEHIVE_LVL                      #####
  
  if(pheno_name == "BEEHIVE_LVL") a$tokeep_epi <- a$tokeep_epi & (a$BEEHIVE_LVL > 2 | (a$BEEHIVE_LVL == 2 & a$SPVL < 3)| is.na(a$SPVL))
  
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
  ids_tokeep_epi <- unique(a$PATIENT[a$tokeep_epi]) # all the IDs we can keep
  
  ######                  LOAD AND PREPARE GENETIC DATA                 #####
    
    idx_aln <- 1 # in the replication dataset, we have only 1 alignment
  
    aln_file <- aln_file_list[[idx_aln]]
    print(names(aln_file_list)[idx_aln])
    aln <- read.fasta(aln_file, as.string = T)
    n_missing_per_ind <- create_G_sigma(aln = aln, only_n_missing_per_ind = T)
    npos_aln <- nchar(aln[[1]])                   # length of the alignment
    aln <- aln[n_missing_per_ind < npos_aln / 2]  # keep only sequences with less than half positions missing
    
    aln_dates <- sapply(names(aln), get.date.from.aln.name)
    aln_ids <- sapply(names(aln), get.bee.id)
    aln_ids_dates <- paste(aln_ids, aln_dates, sep = "_")

    # NOW FOR EACH ID, WE NEED TO KEEP A UNIQUE EPI DATA ROW & SEQUENCE IN THE ALIGNMENT
    a$matching_id_date <- match(a$PATIENT_SEQ_DATE, names(aln))
    a$matching_id_only <- is.na(a$matching_id_date) & a$PATIENT %in% aln_ids
    a_sum <- ddply(a, .(PATIENT), summarise, matching_id_date = any(!is.na(matching_id_date)), matching_id_only = any(matching_id_only))
    table(matching_id_date = a_sum$matching_id_date, matching_id_only = a_sum$matching_id_only) # only one patient has no matching id_date but has matching id
    df_tokeep <- data.frame(id = ids_tokeep_epi, id_date_epi = NA, idx_epi = NA, id_date_genome = NA, idx_genome = NA)
    for(id in ids_tokeep_epi){ # for each patient that we can keep, get the index of epi data and the sequence that are matching
      idx_a <- which(a$PATIENT==id)
      idx_a_sum <- which(a_sum$PATIENT == id); stopifnot(length(idx_a_sum) == 1)
      idx_df <- which(df_tokeep$id == id)
      if(a_sum$matching_id_date[idx_a_sum]){ # if any epi / sequence data of this patient are matching
        first_matching <- which(!is.na(a$matching_id_date[idx_a]))[1] #; if(length(unique(which(!is.na(a$matching_id_date[idx_a])))) > 1) cat("id ", id, " had several possibilities\n")
        my_a_idx <- idx_a[first_matching] # the first id_date that is matching one in the sequence data
        my_genome_idx <- a$matching_id_date[my_a_idx]
        df_tokeep$idx_epi[idx_df] <- my_a_idx
        df_tokeep$id_date_epi[idx_df] <- a$PATIENT_SEQ_DATE[my_a_idx] # the first id_date that is matching one in the sequence data
        df_tokeep$id_date_genome[idx_df] <- names(aln)[my_genome_idx]
        df_tokeep$idx_genome[idx_df] <- my_genome_idx
        
      } else if(a_sum$matching_id_only[idx_a_sum]){ # if only the ID is matching
        print(id)
      }
    } 
    df_tokeep <- df_tokeep[!is.na(df_tokeep$idx_epi), ]
    
    # FINALLY SELECT THE RESULTING EPI FILE AND ALIGNMENT
    suba <- a[df_tokeep$idx_epi, ]
    aln <- aln[df_tokeep$idx_genome]
    stopifnot(all(suba$PATIENT_SEQ_DATE == names(aln)))
    all_replication_sequences_tokeep <- c(all_replication_sequences_tokeep, names(aln)) # save the sequences kept for this phenotype
    names(aln) <- sapply(names(aln), get.bee.id)
    stopifnot(all(!duplicated(names(a))))
    
    # we now not select any id in G:
    GSigma <- create_G_sigma(aln = aln)
    
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
    unique_pos_to_test <- unique(Lvec)
    
    ######                                  CREATE THE K MATRIX WITH FILTERING                       #####
    
    ll <- 17 # index of the criterion corresponding to all positions, but only individuals with less than half positions missing
    stopifnot(all(n_missing_per_ind < npos_aln / 2))                      # check all individuals have less than half positions missing

      
      # subset of positions: now the filtering has been done already
      sub_pos <- rep(TRUE, ngposG)
      sub_ind <- rep(TRUE, nindG)
      
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
      
      if(any(duplicated(suba$PATIENT))) stop('duplicated patients in epi data')
      if(any(duplicated(idsG))) stop('duplicated patients in genetic data')
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
          if(covariate == "age.cat"){ # for age category need to redefine
            suba$age.cat[suba$age.cat %in% names(tab.cov)[tab.cov < 10]] <- "OTH/UNK"
            tab.cov <- table(suba[, covariate])
            tab.cov <- tab.cov[tab.cov > 0]
          }
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
      filename <- paste0(wd, "Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/Gmatrix_data/prepared_", pheno_name, "_", names(aln_file_list)[idx_aln], "_", ll, ".RData")
      save(list = c("Xrand", "normalise_factors", "aln", "G", "Gn", "Sigma_sub", "ids_tokeep_epi", "idsG", "Lvec", "sub_ind", "sub_pos", "reference", "alternative",
                    "suba", "X", "nind", "pheno_name", "fixed_covariate_names", "fixed_covariate_type", "n_fixed_covariates",
                    "ngpos", "nb_ind_per_gpos", "nb_gpos_per_ind", "nb_alternative", "unique_pos_to_test",
                    paste0("PC", sub_to_handle), "PCother", paste0("subtype_", sub_to_handle), "subtype_other", paste0("PCsd", sub_to_handle), "PCsdother"
                    ), file = filename)
      

} # end of loop on phenotype

all_replication_sequences_tokeep <- sort(unique(all_replication_sequences_tokeep))
write.csv(x = all_replication_sequences_tokeep, file = paste0(wd, "Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/Gmatrix_data/all_replication_sequences_tokeep.csv"), row.names = F, col.names = F)



