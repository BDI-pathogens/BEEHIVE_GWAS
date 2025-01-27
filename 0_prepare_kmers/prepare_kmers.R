rm(list = ls())
source("~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/0_createG.R")

# define kmers files
kmerfilelist <- paste0("~/DID/BEEHIVE_Hackathon/Data/Sequences/KmerPartitions/MergedIndividualConsensuses_2018-11-21_genotypes_", 2:6, "mers.csv")
stopifnot(all(file.exists(kmerfilelist)))

for(kmerfile in kmerfilelist){
  print(kmerfile)
  kmers <- read.csv(kmerfile)
  names(kmers)
  n_kmer <- ncol(kmers)-1
  
  # convert ids into beehive ids
  kmers$ID <- sapply(kmers$ID, get.bee.id)
  stopifnot(all(!duplicated(kmers$ID)))
  
  # convert the column names into beginning and end of kmers
  tmp <- sapply(names(kmers)[2:(n_kmer+1)], get.begin.end.kmers)
  kmer_begin <- unname(tmp[1,])
  kmer_end <- unname(tmp[2,])
  table(kmer_end-kmer_begin)
  
  # now convert into a design matrix
  tmp <- create_G_kmers(kmer_table = kmers, ids_tokeep = NULL)
  print(object.size(tmp))
  attach(tmp)
  save(list = c("Gmat_kmer", "Lvec_kmer", "reference_kmer", "alternative_kmer", "names_kmer"), file = gsub(pattern = ".csv", replacement = ".RData", x = kmerfile))
}




