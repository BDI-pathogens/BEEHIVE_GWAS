# 15 March 2018, Francois and Chris
# code to create a G matrix from an alignment
create_blocks <- function(n_alleles){
  mat <- diag(n_alleles)
  mat[mat == 0] <- 1/2
  return(mat)
}
get.date.from.aln.name <- function(x){
  x <- as.character(x)
  if(nchar(x) < 3 | is.na(x)) return(NA)
  if(!substr(x,1,3) == "BEE") return(NA)
  x <- strsplit(x, "_")[[1]][2]
  return(x)
}
get.bee.id <- function(x){
  x <- as.character(x)
  if(nchar(x) < 3 | is.na(x)) return(x)
  if(!substr(x,1,3) == "BEE"){
    return(x)
  }else {
    return(substr(x,1,7)) # the beehive ID
  }
}
get.begin.end.kmers <- function(name){
  name <- gsub(pattern = "X", replacement = "", x = name)
  name <- strsplit(name, split = "\\.")[[1]]
  return(c(as.numeric(name[1]), as.numeric(name[2])))
}

create_G_sigma <- function(aln, ids_tokeep = NULL, verbose = F, only_n_missing_per_ind = F){
  # eliminate patients not in ids_tokeep
  if(!is.null(ids_tokeep)) {
    names(aln) %in% ids_tokeep -> tokeep
    cat("keeping ", sum(tokeep), " individuals out of ", length(aln), "\n")
    aln <- aln[which(tokeep)]
  }
  list.of.names <- names(aln)
  # convert to matrix
  npos <- nchar(aln[[1]])
  nind <- length(aln)
  aln2 <- matrix(unlist(lapply(aln, function(char) strsplit(char, split = "")[[1]])), nrow = npos)
  
  # deal with ambiguous nuc TO DO deal with ambiguous nucleotides better
  ambiguous.nuc <- c('b', 'h', 'k', 'm', 'n', 'r', 's', 'v', 'w', 'y', '?')
  aln2[which(aln2 %in% ambiguous.nuc, arr.ind = T)] <- NA
  all_nuc <- unique(c(aln2)); all_nuc <- all_nuc[!is.na(all_nuc)]
  if(! all(all_nuc %in% c('-', 'a', 'c', 'g', 't'))) stop('unknown nucleotides are present')
  
  n_missing_per_ind <- colSums(is.na(aln2)) # number of nucleotides missing per individual
  if(only_n_missing_per_ind) return(n_missing_per_ind)
  #table(c(aln2))
  
  tab_list <- apply(aln2, 1, function(vec) table(vec))
  n_alleles_list <- unlist(lapply(tab_list, length))
  if(any(n_alleles_list==0)) stop('there might be fully unspecified columns in the alignment')
  
  #create the blocks we will need for the sigma matrix (correction for multi-allelic sites)
  list_of_blocks <- list()
  for(i in 1:max(n_alleles_list)){
    list_of_blocks[[i]] <- create_blocks(i) # computes each block we need
  }
  
  ncolG <- sum(n_alleles_list) - npos # there are nallele - 1 column for a position with nallele alleles
  Gmat <- matrix(NA, nrow = nind, ncol = ncolG)
  Sigma <- matrix(0, nrow = ncolG, ncol = ncolG)
  reference <- rep(NA, ncolG)
  alternative <- rep(NA, ncolG)
  Lvec <- rep(NA, ncolG) # a vector containing the nucleotidic locus position of each column in the alignment
  
  idx <- 0
  for(pos in 1:npos){
    if(verbose) print(pos)
    n_alleles <- n_alleles_list[pos]
    if(n_alleles < 2) next
    tab <- sort(tab_list[[pos]], decreasing = T)
    namestab <- names(tab)
    most.common <- namestab[1]
    
    # here add the block to the sigma matrix (correction for multi-allelic sites)
    Sigma[(idx + 1):(idx + n_alleles - 1), (idx + 1):(idx + n_alleles - 1)] <- list_of_blocks[[n_alleles - 1]]
    
    for(allele in 2:n_alleles){
      idx <- idx + 1
      reference[idx] <- most.common
      alternative[idx] <- namestab[allele]
      Gmat[ ,  idx] <- (aln2[pos, ] == namestab[allele])
      Lvec[idx] <- pos # indicate that this column idx concerns nucleotidic position pos in the alignment
    }
    #print(pos)
  }
  return(list(Gmat = Gmat, Sigma = Sigma, Lvec = Lvec, reference = reference, alternative = alternative, names = list.of.names, n_missing_per_ind = n_missing_per_ind))
}
create_G_kmers <- function(kmer_table, ids_tokeep = NULL, verbose = F){ # create G matrix and others for kmer table
  # eliminate patients not in ids_tokeep
  if(!is.null(ids_tokeep)) {
    kmer_table$ID %in% ids_tokeep -> tokeep
    cat("keeping ", sum(tokeep), " individuals out of ", nrow(kmer_table), "\n")
    kmer_table <- kmer_table[which(tokeep),]
  }
  list.of.names <- kmer_table$ID
  
  npos <- ncol(kmer_table)-1
  kmer_table <- kmer_table[, 2:(npos+1)]
  nind <- nrow(kmer_table)
  
  # deal with NA
  missing_data <- "?"
  kmer_table[which(kmer_table == missing_data, arr.ind = T)] <- NA
  
  n_missing_per_ind <- colSums(is.na(kmer_table)) # number of nucleotides missing per individual
  #table(c(aln2))
  
  tab_list <- apply(kmer_table, 2, function(vec) table(vec))
  n_alleles_list <- unlist(lapply(tab_list, length))
  if(any(n_alleles_list==0)) stop('there might be fully unspecified columns in the alignment')
  
  
  ncolG <- sum(n_alleles_list) - npos # there are nallele - 1 column for a position with nallele alleles
  Gmat <- matrix(NA, nrow = nind, ncol = ncolG)
  reference <- rep(NA, ncolG)
  alternative <- rep(NA, ncolG)
  Lvec <- rep(NA, ncolG) # a vector containing the nucleotidic locus position of each column in the alignment
  
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
      Gmat[ ,  idx] <- (kmer_table[, pos] == namestab[allele])
      Lvec[idx] <- pos # indicate that this column idx concerns nucleotidic position pos in the alignment
    }
    #print(pos)
  }
  return(list(Gmat_kmer = Gmat,Lvec_kmer = Lvec, reference_kmer = reference, alternative_kmer = alternative, names_kmer = list.of.names, n_missing_per_ind_kmer = n_missing_per_ind))
}
normalise_mat_diag <- function(mat, return_factor = F){
  x <- mean(diag(mat)) - mean(c(mat))
  if(x!=0) mat <- mat / x else stop()
  if(return_factor) return(list(mat = mat, factor = x)) else return(mat)
}

normalise_colmeans_mat <- function(mat){
  column_means <- colMeans(mat, na.rm = T) # mean across individuals
  mat_normalised <- sweep(x = mat, MARGIN = 2, STATS = column_means, FUN = "-") # a normalised matrix
  return(mat_normalised)
}
