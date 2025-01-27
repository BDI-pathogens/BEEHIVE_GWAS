rm(list = ls())
library(seqinr)
wd <- "~/DID/BEEHIVE_Hackathon/"
prot_file <- paste0(wd, "/Data/Typewriter/clean/concatenated_PRO.fasta")

create_protein_mat <- function(prot_aln, ids_tokeep = NULL, verbose = F){
  # eliminate patients not in ids_tokeep
  if(!is.null(ids_tokeep)) {
    names(prot_aln) %in% ids_tokeep -> tokeep
    cat("keeping ", sum(tokeep), " individuals out of ", length(prot_aln), "\n")
    prot_aln <- prot_aln[which(tokeep)]
  }
  list.of.names <- names(prot_aln)
  # convert to matrix
  npos <- nchar(prot_aln[[1]])
  nind <- length(prot_aln)
  aln2 <- matrix(unlist(lapply(prot_aln, function(char) strsplit(char, split = "")[[1]])), nrow = npos)
  
  # deal with missing data TODO CHECK THIS IS ALRIGHT
  ambiguous.aa <- c('x')
  aln2[which(aln2 %in% ambiguous.aa, arr.ind = T)] <- NA
  list_of_aa <- c("*", "a", "b",  "c", "d", "e", "f", "g", "h", "i", "k", "l", "m", "n", "p", "q", "r", "s", "t", "v", "w", "y", "z") #tolower(a(aaa()))
  all_aa_found <- names(table(c(aln2)))
  if(!all(all_aa_found %in% c('-', list_of_aa))) stop()
  
  tab_list <- apply(aln2, 1, function(vec) table(vec))
  n_alleles_list <- unlist(lapply(tab_list, length))
  if(any(n_alleles_list==0)){
    print('there are fully unspecified columns in the alignment... removing them')
    cat('positions', which(n_alleles_list==0), "\n")
  }
  
  ncol_mat <- sum(n_alleles_list[which(n_alleles_list > 0)]) - (npos - sum(n_alleles_list == 0)) # there are nallele - 1 column for a position with nallele alleles - except the pos with 0 alleles that count 0
  prot_mat <- matrix(NA, nrow = nind, ncol = ncol_mat)
  reference <- rep(NA, ncol_mat)
  alternative <- rep(NA, ncol_mat)
  Lvec <- rep(NA, ncol_mat) # a vector containing the nucleotidic locus position of each column in the alignment
  
  idx <- 0
  for(pos in 1:npos){
    if(verbose) print(pos)
    n_alleles <- n_alleles_list[pos]
    if(n_alleles < 2) next
    tab <- sort(tab_list[[pos]], decreasing = T)
    namestab <- names(tab)
    most.common <- namestab[1]
    
    for(allele in 2:n_alleles){
      idx <- idx + 1
      reference[idx] <- most.common
      alternative[idx] <- namestab[allele]
      prot_mat[ ,  idx] <- (aln2[pos, ] == namestab[allele])
      Lvec[idx] <- pos # indicate that this column idx concerns nucleotidic position pos in the alignment
    }
    #print(pos)
  }
  return(list(prot_mat = prot_mat, Lvec = Lvec, reference = reference, alternative = alternative, names = list.of.names))
}
######                                  CREATE THE PROTEIN MATRIX AND SAVE IT                               #####

prot_aln <- read.fasta(prot_file, as.string= T)
prot_mat_list <- create_protein_mat(prot_aln = prot_aln) # TO CHECK HOW IT DEALS WITH N AND GAP; WHAT IS THE SYMBOL FOR MISSING DATA VERSUS GENUINE GAP

G_prot <- prot_mat_list$prot_mat
reference_prot <- prot_mat_list$reference
alternative_prot <- prot_mat_list$alternative
window.starts_prot <- prot_mat_list$Lvec
window.ends_prot <- prot_mat_list$Lvec
idsG_prot <- prot_mat_list$names

save(list = c("G_prot", "reference_prot", "alternative_prot", "window.ends_prot", "window.starts_prot", "idsG_prot"), 
     file = "~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/Gmatrix_data/prot.RData")


