# R code to prepare the consensus sequence file by removing the duplicated patients, resulting in MergedIndividualConsensuses_renamed_clean_noduplicate_2019-11-20.csv
# (incidentally also removes duplicated patients that are not in the data that can be used for the GWAS)
# similarly we place duplicate frequency files in the "ExtraSeq" folder

# Rationale: we want to eliminate the duplicated individuals
# from the genomic dataset before preparing the kmer / length / length3 / base frequency data
# For that, we need to come back to the saved epi data (files of the form "prepared_pheno_adjusted_full_17.RData") outputted by 0_prepare_replication_data.R
# Indeed the 0_prepare_replication_data.R file prepares the G matrix data, etc, and saves it under the .RData file,
# but does not remove duplicates in the consensus sequence file nor moves the duplicate individual frequency files in a separate folder

######        1. first record all IDs and dates that will be used for the GWAS         ###### 

# these patients are the intersection of patients in ~/DID/BEEHIVE_Hackathon/Data/Sequences_replication/MergedIndividualConsensuses_renamed_clean_noduplicate_2019-11-20.csv
# and patients in the epi data

all_patient_date <- c() # a vector of all patient_date combinations that are in the epi data
for(pheno in c("BEEHIVE_LVL", "spvl_adjusted", "spvl_normalised_adjusted", "CD4_slope")){
  
  # need to first run the code 4_prepare_replication_data.R to generate the prepared data files and be able to create the list of patients and dates
  load(paste0(
    "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/Gmatrix_data/prepared_", pheno, "_full_17.RData")
  )
  
  print(dim(suba))
  print(dim(G))
  print(table(sub_ind))
  stopifnot(anyDuplicated(suba$PATIENT_SEQ_DATE) == 0)
  all_patient_date <- c(all_patient_date, suba$PATIENT_SEQ_DATE)

}
all_patient_date <- sort(unique(all_patient_date)) # 323 unique individuals as of 22/10/2020
stopifnot(all(all_patient_date == read.csv("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/Gmatrix_data/all_replication_sequences_tokeep.csv")$x))
write.csv(all_patient_date, "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/Gmatrix_data/all_patients_dates.csv", row.names = F, col.names = F)

######        2. identify the merged individual consensuses and eliminate duplicates ids & those without matching date        ###### 

library(stringr)

# remove sequences from the merged consensuses in the replication dataset
csv.file <- "~/DID/BEEHIVE_Hackathon/Data/Sequences_replication/MergedIndividualConsensuses_renamed_clean_2019-11-20.csv"
stopifnot(file.exists(csv.file))
df <- read.csv(csv.file, stringsAsFactors = F, check.names = F)
npat <- ncol(df) - 2
col.name.regex <- "kmer in (.*)"
keep.col.names <- grepl(col.name.regex, colnames(df)[3:(2 + npat)])
stopifnot(any(keep.col.names))
pat.ids <- str_match(colnames(df)[3:(2 + npat)], col.name.regex)[,2]
stopifnot(all(!is.na(pat.ids)))
all_dates <- sapply(pat.ids, function(x) strsplit(x, "_")[[1]][2])
all_ids <- sapply(pat.ids, function(x) strsplit(x, "_")[[1]][1])

all_ids_dates <- paste(all_ids, all_dates, sep = "_") # list of unique BEE ID_dates
stopifnot(!any(duplicated(all_ids_dates))) # those must be unique

to_remove <- c() # create a vector of indices to remove
for(id in unique(all_ids[duplicated(all_ids)])){ # for each duplicated id
  idx <- unname(which(all_ids == id))
  if(any(all_ids_dates[idx] %in% all_patient_date)){
    to_remove <- c(to_remove, idx[!all_ids_dates[idx] %in% all_patient_date]) # if one of the sequences matches the date, remove the other sequence
  } else {
    cat(id, " not having matching patient_date in prepared GWAS data\n")
    to_remove <- c(to_remove, idx[!all_ids_dates[idx] %in% all_patient_date]) # in that case both of them will be removed
  }
  #tokeep <- idx[which()]
}

to_remove <- to_remove + 2 # since the two first column of df are for the references, correct the column number
df_noduplicate <- df[, -to_remove]
View(df_noduplicate)

write.csv(df_noduplicate, "~/DID/BEEHIVE_Hackathon/Data/Sequences_replication/MergedIndividualConsensuses_renamed_clean_noduplicate_2019-11-20.csv", row.names = F, col.names = F)

######        3. similary place duplicated frequency files in separate folder         ###### 

freq_dir <- "~/DID/BEEHIVE_Hackathon/Data/Sequences_replication/BaseFreqs_renamed_clean/"
lf <- list.files(freq_dir)
lf <- lf[!lf %in% c("ExtraSeq", "ShorterSeqsFromTheSamePatandTimePoint")]

all_dates <- sapply(lf, function(x) strsplit(x, "_")[[1]][2])
all_dates <- gsub(pattern = ".csv", replacement = "", x = all_dates)
all_ids <- sapply(lf, function(x) strsplit(x, "_")[[1]][1])

all_ids_dates <- paste(all_ids, all_dates, sep = "_")
stopifnot(!any(duplicated(all_ids_dates)))

to_move <- c()
for(id in unique(all_ids[duplicated(all_ids)])){
  idx <- unname(which(all_ids == id))
  if(any(all_ids_dates[idx] %in% all_patient_date)){
    to_move <- c(to_move, idx[!all_ids_dates[idx] %in% all_patient_date])
  } else {
    cat(id, " not having matching patient_date\n")
    to_move <- c(to_move, idx[!all_ids_dates[idx] %in% all_patient_date]) # in that case both of them will be removed
  }
  #tokeep <- idx[which()]
}

# move files to ExtraSeq directory
tmp <- file.copy(from = paste0(freq_dir, lf[to_move]), to = paste0(freq_dir, "ExtraSeq"))
if(all(tmp)) file.remove(paste0(freq_dir, lf[to_move]))



