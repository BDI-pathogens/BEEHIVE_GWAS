args = commandArgs(trailingOnly=TRUE)
library(plyr)
library(MASS)
library(seqinr)
library(caret)
library(igraph)

# also needs igraph, caret

if(length(args)==0){
  print("warning length of args is 0, replacing by default")
  args <- c("BEEHIVE_LVL", 0.05)
  print(args)
  setwd("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/")
}
stopifnot(length(args)==2)
pheno_name <- args[1] # "spvl_normalised_adjusted"
p_t <- as.numeric(as.character(args[2])) # the p_value_threshold

#####                           DEFINE AND LOAD DATA                        ####

bonferroni_folder <- "bonferroni/"
bonferroni_files <- list.files(bonferroni_folder)
prot_to_hxb2_file <- "concatenated_PRO_coords.csv"
hxb2_reference_file <- "hxb2_annotated_Chris.csv"

stopifnot(file.exists(prot_to_hxb2_file))
stopifnot(file.exists(hxb2_reference_file))

# load the prepared data and the GWAS results
load(
  paste0("Gmatrix_data/prepared_", pheno_name, "_full_17.RData")
)
load(
  paste0("GWAS_results/firstGWAS_", pheno_name, "_full_17.RData")
)

prot_to_hxb2 <- read.csv(prot_to_hxb2_file, header = T) # file containing for each position in the protein alignment, the corresponding HIV amino-acid
hxb2_map <- read.csv(hxb2_reference_file, header = T)

#TODO PRINT WHY SOME POSITIONS ARE REMOVED MAYBE

#####                      SHOW THE VARIANTS IN THE TABLE 1 and SUP TABLES 2, 3, 4 OF PAPER                             ####

# BEEHIVE_LVL
GWAS_kmer[which(GWAS_kmer$start_position_alignment==1124 & GWAS_kmer$end_position_alignment==1126), c("reference", paste0("alternative", 1:9))]
GWAS_kmer[which(GWAS_kmer$start_position_alignment==1413 & GWAS_kmer$end_position_alignment==1416), c("reference", paste0("alternative", 1:3))]
GWAS_kmer[which(GWAS_kmer$start_position_alignment==1514 & GWAS_kmer$end_position_alignment==1514), c("reference", paste0("alternative", 1:2))]
GWAS_kmer[which(GWAS_kmer$start_position_alignment==7676 & GWAS_kmer$end_position_alignment==7678), c("reference", paste0("alternative", 1:3))]
GWAS_kmer[which(GWAS_kmer$start_position_alignment==9008 & GWAS_kmer$end_position_alignment==9008), c("reference", paste0("alternative", 1:2))]

# SPVL adjusted
GWAS_kmer[which(GWAS_kmer$start_position_alignment==660 & GWAS_kmer$end_position_alignment==664), ]
GWAS_kmer[which(GWAS_kmer$start_position_alignment==1413 & GWAS_kmer$end_position_alignment==1416), ]
GWAS_kmer[which(GWAS_kmer$start_position_alignment==1514 & GWAS_kmer$end_position_alignment==1514), ]
GWAS_kmer[which(GWAS_kmer$start_position_alignment==9008 & GWAS_kmer$end_position_alignment==9008), ]


# SPVL adjusted normalised
GWAS_kmer[which(GWAS_kmer$start_position_alignment==1413 & GWAS_kmer$end_position_alignment==1416), ]
GWAS_kmer[which(GWAS_kmer$start_position_alignment==1514 & GWAS_kmer$end_position_alignment==1514), ]
GWAS_kmer[which(GWAS_kmer$start_position_alignment==9008 & GWAS_kmer$end_position_alignment==9008), ]


# CD4 slope
GWAS_kmer[which(GWAS_kmer$start_position_alignment==2172 & GWAS_kmer$end_position_alignment==2173), ]
GWAS_kmer[which(GWAS_kmer$start_position_alignment==5553 & GWAS_kmer$end_position_alignment==5558), ]


#####                                 DEFINE FUNCTIONS                              ####
source(paste0(bonferroni_folder, "get_bonferroni_dfs.R"))

is_overlapping <- function(idx1, idx2) any(idx1 %in% idx2)

clean_table <- function(focal_tab){
  # remove empty rows and columns
  rows_allna <- apply(focal_tab, 1, function(vec) all(is.na(vec)))
  cols_allna <- apply(focal_tab, 2, function(vec) all(is.na(vec)))
  focal_tab <- focal_tab[!rows_allna,]
  focal_tab <- focal_tab[, !cols_allna]
  
  return(focal_tab)
}
get_index_in_dm <- function(st = NULL, en = NULL, idx_in_GWAS = NULL, GWAS_tab, full_design_matrix, full_ind_tokeep){
  # this function get column position in design matrix of variant starting at "st" and ending at "en" position in GWAS_tab
  # and the column position in "individuals to keep" matrix of this variant
  # alternatively, the variant can be specified by row idx_in_GWAS in GWAS_tab
  
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
  
  idx_in_ind <- idx_in_GWAS
  
  # check that the indices are correct
  stopifnot(sum(full_ind_tokeep[, idx_in_ind])==GWAS_tab$N[idx_in_GWAS])
  
  stopifnot(
    all(colSums(full_design_matrix[, idx_in_dm_start:idx_in_dm_end, drop = F] >= 0.01, na.rm = T) == GWAS_tab[idx_in_GWAS, paste0("N", 1:GWAS_tab$df[idx_in_GWAS])])
    # the ">= 0.01" is for frequency variants
  )
  vec <- c(idx_in_ind, idx_in_dm_start, idx_in_dm_end)
  names(vec) <- c("idx_in_ind", "idx_in_dm_start", "idx_in_dm_end")
  return(vec)
  #apply(, 2, sum, na.rm = T)
}
select_significant <- function(GWAS_tab, reduce_redundant = F, p_threshold){
  # this function takes as input the table of GWAS result
  # and selects bonferroni-corrected significant positions
  # and optionally select only one position among overlapping positions ("reduce_redundant" option)
  
  stopifnot("p_b" %in% colnames(GWAS_tab))
  
  # significant hits (0.05 level)
  signif_GWAS_tab <- clean_table(GWAS_tab[GWAS_tab$p_b < p_threshold,])
  
  # reduce redundant hits
  if(reduce_redundant & nrow(signif_GWAS_tab) > 0){
    
    stopifnot("k" %in% colnames(GWAS_tab)) # k = k-mer length is used in choice
    stopifnot("start_position_alignment" %in% colnames(GWAS_tab))
    stopifnot("end_position_alignment" %in% colnames(GWAS_tab))
    
    # order by start position
    signif_GWAS_tab <- signif_GWAS_tab[order(signif_GWAS_tab$start_position_alignment), ]
    # generate a matrix equal to 1 when positions are overlapping, 0 otherwise
    n_s <- nrow(signif_GWAS_tab)
    mat <- matrix(NA, nrow = n_s, ncol = n_s)
    cat("building adjacency matrix for", n_s, "positions...\n")
    for(i in 1:n_s){
      if(i/10==floor(i/10)) print(i)
      for(j in i:n_s){
        mat[i, j] <- is_overlapping(
          signif_GWAS_tab$start_position_alignment[i]:signif_GWAS_tab$end_position_alignment[i],
          signif_GWAS_tab$start_position_alignment[j]:signif_GWAS_tab$end_position_alignment[j]
        )
      }
    }
    mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
    stopifnot(isSymmetric(mat))
    
    print("building graph based on adjacency matrix...")
    mygraph <- graph_from_adjacency_matrix(adjmatrix = mat, mode = "undirected") # build a graph based on adjacency matrix
    signif_GWAS_tab$cluster <- clusters(graph = mygraph, mode = "strong")$membership    # membership to connected components
    n_clusters <- length(unique(signif_GWAS_tab$cluster))
    
    # select only one hit per cluster
    reduced_signif_GWAS_tab <- data.frame()
    for(i in 1:n_clusters){
      tmp <- signif_GWAS_tab[signif_GWAS_tab$cluster==i,]
      reduced_signif_GWAS_tab <- rbind(reduced_signif_GWAS_tab,
                                       tmp[with(tmp, order(k, p_b))[1],] # prioritise small k-mers, small p-values
      )
    }
    signif_GWAS_tab <- reduced_signif_GWAS_tab
  }
  return(signif_GWAS_tab)
}

get_dm <- function(reduced_signif, GWAS, full_design_matrix, full_ind_tokeep, type, return_full = F){
  
  # function to get the indices in design matrix and list of individual to keep
  # of the significant positions in reduced_signif
  # reduced_signif is a reduced GWAS result table whose rows are position
  # and whose columns must contains names "start_position_alignment" and "end_position_alignment"
  # but contains only (non-redundant) significant positions
  # GWAS is the full GWAS results table
  
  if(nrow(reduced_signif) == 0) return(list(idx = NULL, idxnames = NULL))
  # test that arguments are compatible with this function:
  stopifnot("start_position_alignment" %in% colnames(reduced_signif))
  stopifnot("end_position_alignment" %in% colnames(reduced_signif))
  stopifnot(ncol(full_ind_tokeep)==nrow(GWAS))
  stopifnot(sum(GWAS$df) == ncol(full_design_matrix))
  npos <- nrow(reduced_signif)
  idx <- lapply(1:npos, function(i){
    
    get_index_in_dm(st = reduced_signif$start_position_alignment[i],
                        en = reduced_signif$end_position_alignment[i],
                        GWAS_tab = GWAS,
                        full_design_matrix = full_design_matrix,
                        full_ind_tokeep = full_ind_tokeep)
    
  })
  if(return_full){
    return(idx) 
  } else {
    idx <- lapply(idx, function(vec) vec["idx_in_dm_start"]:vec["idx_in_dm_end"])
    posnames <- lapply(1:npos, function(i) {
      paste(type, reduced_signif$start_position_alignment[i], reduced_signif$end_position_alignment[i], sep = "_")
    })
    idxnames <- lapply(1:npos, function(i) {
      paste(type, reduced_signif$start_position_alignment[i], reduced_signif$end_position_alignment[i], 1:length(idx[[i]]), sep = "_")
    })
    idx <- unlist(idx)
    return(list(
      idx = unlist(idx),
      posnames = unlist(posnames),
      idxnames = unlist(idxnames)
    ))
  }
}

#####                           NOW DEFINE VARIABLES                             ####

all_phenos <- c("BEEHIVE_LVL", "spvl_normalised_adjusted", "spvl_adjusted", "CD4_slope")
all_types <- c("prot", "kmer", "len", "len3", "freq")
all_types2 <- c("drm", "rare")
names_tab <- paste0("GWAS_", all_types) # names of tables for each type of test

#####                           1. ATTACH BONFERRONI-CORRECTED P-VALUE TO EACH GWAS MATRIX                      ####

# k-mer length
GWAS_prot$k <- GWAS_prot$end_position_alignment - GWAS_prot$start_position_alignment + 1
GWAS_kmer$k <- GWAS_kmer$end_position_alignment - GWAS_kmer$start_position_alignment + 1
GWAS_len$k <- GWAS_len$end_position_alignment - GWAS_len$start_position_alignment + 1
GWAS_len3$k <- GWAS_len3$end_position_alignment - GWAS_len3$start_position_alignment + 1
GWAS_freq$k <- GWAS_freq$end_position_alignment - GWAS_freq$start_position_alignment + 1

# associate effective dfs to each GWAS
eff_dfs <- sapply(all_types, function(ty) {
  if(ty=="kmer"){
    get_bonferroni_dfs(pheno = pheno_name, type = ty)$all_dfs
  } else {
    get_bonferroni_dfs(pheno = pheno_name, type = ty)$total_dfs
  }
})
eff_dfs$drm <- nrow(GWAS_drm)
eff_dfs$rare <- nrow(GWAS_rare)

GWAS_drm$p_b <-  GWAS_drm$p  * eff_dfs[["drm"]]
GWAS_freq$p_b <- GWAS_freq$p * eff_dfs[["freq"]]
GWAS_len$p_b  <- GWAS_len$p  * eff_dfs[["len"]]
GWAS_len3$p_b  <-GWAS_len3$p  * eff_dfs[["len3"]]
GWAS_prot$p_b  <-GWAS_prot$p  * eff_dfs[["prot"]]
GWAS_rare$p_b  <-GWAS_rare$p  * eff_dfs[["rare"]]

# kmers, define effective dfs for each kmer length
GWAS_kmer$p_b <- NA
for(i in 1:6){
  subset <- GWAS_kmer$k == i
  GWAS_kmer$p_b[subset] <- GWAS_kmer$p[subset] * eff_dfs[["kmer"]][i]
}

#####                                     2. SELECT AND REDUCE SIGNIFICANT POSITIONS                                           ####

reduced_signif_GWAS_prot <- select_significant(GWAS_tab = GWAS_prot, reduce_redundant = F, p_threshold = p_t)
reduced_signif_GWAS_kmer <- select_significant(GWAS_tab = GWAS_kmer, reduce_redundant = T, p_threshold = p_t) # note with reduce kmers
reduced_signif_GWAS_len <- select_significant(GWAS_tab = GWAS_len, reduce_redundant = F, p_threshold = p_t)
reduced_signif_GWAS_len3 <- select_significant(GWAS_tab = GWAS_len3, reduce_redundant = F, p_threshold = p_t)
reduced_signif_GWAS_freq <- select_significant(GWAS_tab = GWAS_freq, reduce_redundant = F, p_threshold = p_t)
reduced_signif_GWAS_drm  <- select_significant(GWAS_tab = GWAS_drm, reduce_redundant = F, p_threshold = p_t)
reduced_signif_GWAS_rare  <- select_significant(GWAS_tab = GWAS_rare, reduce_redundant = F, p_threshold = p_t)


# remove protein hits that are already in kmer hits
if(nrow(reduced_signif_GWAS_prot) > 0){
  
  ################        MODIFICATIONS TO CHRIS' HXB2 REFERENCE        ################
  
  # 1) in hxb2_map, transform gp120 and gp41 in "env" and renumber
  hxb2_map$RF3.protein <- as.character(hxb2_map$RF3.protein)
  stopifnot(hxb2_map$RF3.aa.position[hxb2_map$RF3.protein=="gp120"]==rep(1:511, each = 3))
  stopifnot(hxb2_map$RF3.aa.position[hxb2_map$RF3.protein=="gp41"]==rep(1:345, each = 3))
  hxb2_map$RF3.aa.position[hxb2_map$RF3.protein=="gp41"] <- 511 + rep(1:345, each = 3) # re-number gp41 positions, starting with 512
  
  hxb2_map$RF3.protein[hxb2_map$RF3.protein == "gp120" | hxb2_map$RF3.protein == "gp41"] <- "env"
  stopifnot(hxb2_map$RF3.aa.position[hxb2_map$RF3.protein=="env"] == rep(1:(511+345), each = 3))
  
  # 2) MAKE POL START FROM POSITION 2085 INSTEAD OF 2253
  pol_beginning <- which(hxb2_map$HXB2.base.position == 2085):which(hxb2_map$HXB2.base.position == 2252)
  initial_pol_pos <- which(hxb2_map$HXB2.base.position == 2253):which(hxb2_map$HXB2.base.position == 5093)

  hxb2_map[pol_beginning, "RF3.protein"] <- "pol"
  hxb2_map[pol_beginning, "RF3.aa.position"] <- c(sapply(1:(length(pol_beginning)/3), rep, 3))
  hxb2_map[initial_pol_pos, "RF3.aa.position"] <- hxb2_map[initial_pol_pos, "RF3.aa.position"] + length(pol_beginning)/3
  hxb2_map[pol_beginning, "RF3.aa"] <- c(sapply(translate(seq = hxb2_map[pol_beginning, "HXB2.base"]), rep, 3))
  
  # 3) SHIFT INDICES IN prot_to_hxb2 TO ACCOUNT FOR STOP CODONS IN NEF AND TAT
  add1 <- function(cc){ # function to add 1 to the character string index in prot_to_hxb2$x
    cc <- strsplit(cc, "_")[[1]]
    cc[2] <- as.character(as.numeric(cc[2]) + 1)
    return(paste(cc, collapse = "_"))
  }
  # for nef:
  prot_to_hxb2[1480:1561, "x"] <- sapply(prot_to_hxb2[1480:1561, "x"], add1)
  # for tat:
  prot_to_hxb2[2768:2780, "x"] <- sapply(prot_to_hxb2[2768:2780, "x"], add1)

  # 4) CHECK THE GENES DEFINED IN TWO FILES HAVE SAME LENGTH
  
  # pol OK
  stopifnot(max(as.numeric(gsub(pattern = "pol_", replacement = "", x = prot_to_hxb2$x[grepl("pol", prot_to_hxb2$x)]))) ==
              max(hxb2_map[which(hxb2_map$RF3.protein=="pol"), "RF3.aa.position"]))
  # env OK
  stopifnot(max(as.numeric(gsub(pattern = "env_", replacement = "", x = prot_to_hxb2$x[grepl("env", prot_to_hxb2$x)]))) ==
              max(hxb2_map[which(hxb2_map$RF3.protein=="env"), "RF3.aa.position"]))
  # gag OK
  stopifnot(max(hxb2_map[which(hxb2_map$RF1.protein=="gag"), "RF1.aa.position"]) == 
              max(as.numeric(gsub(pattern = "gag_", replacement = "", x = prot_to_hxb2$x[grepl("gag", prot_to_hxb2$x)]))))
  
  # nef not quite the same (205 aa in my amino-acid alignment, 206 in Chris reference) (BECAUSE OF STOP CODON; CORRECTED)
  range(hxb2_map[which(hxb2_map$RF1.protein=="nef"), "RF1.aa.position"])
  stopifnot(range(hxb2_map[which(hxb2_map$RF3.protein=="nef"), "RF3.aa.position"])[2] == 
    max(as.numeric(gsub(pattern = "nef_", replacement = "", x = prot_to_hxb2$x[grepl("nef", prot_to_hxb2$x)]))))
  
  # vif OK
  stopifnot(max(hxb2_map[which(hxb2_map$RF1.protein=="vif"), "RF1.aa.position"]) == 
    max(as.numeric(gsub(pattern = "vif_", replacement = "", x = prot_to_hxb2$x[grepl("vif", prot_to_hxb2$x)]))))
  
  # tat not quite the same (100 aa in my amino-acid alignment, 101 in Chris) (BECAUSE OF STOP CODON; CORRECTED)
  range(hxb2_map[which(hxb2_map$RF2.protein=="tat"), "RF2.aa.position"])
  stopifnot(range(hxb2_map[which(hxb2_map$RF1.protein=="tat"), "RF1.aa.position"])[2] ==
    max(as.numeric(gsub(pattern = "tat_", replacement = "", x = prot_to_hxb2$x[grepl("tat", prot_to_hxb2$x)]))))
  
  # vpr OK
  stopifnot(range(hxb2_map[which(hxb2_map$RF1.protein=="vpr"), "RF1.aa.position"])[2] == 
              max(as.numeric(gsub(pattern = "vpr_", replacement = "", x = prot_to_hxb2$x[grepl("vpr", prot_to_hxb2$x)]))))
            
  # rev OK
  stopifnot(
    range(hxb2_map[which(hxb2_map$RF2.protein=="rev"), "RF2.aa.position"])[2] ==
   max(as.numeric(gsub(pattern = "rev_", replacement = "", x = prot_to_hxb2$x[grepl("rev", prot_to_hxb2$x)]))))
  
  # vpu OK
  stopifnot(range(hxb2_map[which(hxb2_map$RF2.protein=="vpu"), "RF2.aa.position"])[2] == 
    max(as.numeric(gsub(pattern = "vpu_", replacement = "", x = prot_to_hxb2$x[grepl("vpu", prot_to_hxb2$x)]))))
  
  ################        NOW LOOK AT POSITIONS OF SIGNIFICANT AA IN HXB2        ################
  
  # look for gene and amino-acid position of this variant in prot_to_hxb2 table
  idx_in_prot <- match(reduced_signif_GWAS_prot$start_position_alignment, prot_to_hxb2$X)
  stopifnot(all(!is.na(idx_in_prot)))
  genes_pos <- prot_to_hxb2$x[idx_in_prot]
  genes_pos <- strsplit(as.character(genes_pos), split = "_")
  genes <- unlist(lapply(genes_pos, "[", 1))
  pos <- as.numeric(unlist(lapply(genes_pos, "[", 2))) # amino-acid number in that gene
  
  # look for this position in all three reading frames of hxb2
  idx_in_hxb2 <- sapply(1:length(pos), function(i){
    c(
    which(hxb2_map$RF1.protein == genes[i] & hxb2_map$RF1.aa.position == pos[i]),
    which(hxb2_map$RF2.protein == genes[i] & hxb2_map$RF2.aa.position == pos[i]),
    which(hxb2_map$RF3.protein == genes[i] & hxb2_map$RF3.aa.position == pos[i])
  )})
  
  # cut to 3 indices only (for nef gene is present at two places in the genome -> can be found twice)
  if(is.list(idx_in_hxb2)){
    idx_in_hxb2 <- matrix(unlist(lapply(idx_in_hxb2, function(vec) vec[1:3])), ncol = 3, byrow = T)
  } else {
    idx_in_hxb2 <- t(idx_in_hxb2)
  }
  stopifnot(all(dim(idx_in_hxb2) == c(length(pos), 3)))
  
  # remove protein position corresponding to nucleotides already in kmer
  any_in_kmer <- apply(idx_in_hxb2, 1, function(vec_idx) any(vec_idx %in% sapply(
      1:nrow(reduced_signif_GWAS_kmer), function(i) reduced_signif_GWAS_kmer$start_position_alignment[i]:reduced_signif_GWAS_kmer$end_position_alignment[i]
    )
  ))
  reduced_signif_GWAS_prot <- reduced_signif_GWAS_prot[!any_in_kmer, ]
}

# Error in apply(idx_in_hxb2, 2, function(vec_idx) any(vec_idx %in% sapply(1:nrow(reduced_signif_GWAS_kmer),  : 
#                                                                            dim(X) doit avoir un longueur positive
#                                                                          Exécution arrêtée

# table(Xvariants[, "kmer_9008_9008_1"], Xvariants[, "prot_1427_1427_1"]) # these things were perfectly correlated when not removing the amino-acids corresponding to already-found kmers
# table(Xvariants[, "kmer_1514_1514_1"], Xvariants[, "prot_1098_1098_1"])

#####                           3. GET DESIGN MATRIX AND SELECTED INDIVIDUALS FOR EACH POSITIONS                               ####

pos_dm_kmer <- get_dm(reduced_signif = reduced_signif_GWAS_kmer, GWAS = GWAS_kmer,
       full_design_matrix = full_design_matrix_kmer, full_ind_tokeep = full_ind_tokeep_kmer, type = "kmer")
pos_dm_prot <- get_dm(reduced_signif = reduced_signif_GWAS_prot, GWAS = GWAS_prot,
                      full_design_matrix = full_design_matrix_prot, full_ind_tokeep = full_ind_tokeep_prot, type = "prot")
pos_dm_len <- get_dm(reduced_signif = reduced_signif_GWAS_len, GWAS = GWAS_len,
                      full_design_matrix = full_design_matrix_len, full_ind_tokeep = full_ind_tokeep_len, type = "len")
pos_dm_len3 <- get_dm(reduced_signif = reduced_signif_GWAS_len3, GWAS = GWAS_len3,
                     full_design_matrix = full_design_matrix_len3, full_ind_tokeep = full_ind_tokeep_len3, type = "len3")
pos_dm_freq <- get_dm(reduced_signif = reduced_signif_GWAS_freq, GWAS = GWAS_freq,
                      full_design_matrix = full_design_matrix_freq, full_ind_tokeep = full_ind_tokeep_freq, type = "freq")

# construct design matrix for variant effects:
Xvariants <- cbind(
  full_design_matrix_kmer[, pos_dm_kmer$idx],
  full_design_matrix_prot[, pos_dm_prot$idx],
  full_design_matrix_len[, pos_dm_len$idx],
  full_design_matrix_len3[, pos_dm_len3$idx],
  full_design_matrix_freq[, pos_dm_freq$idx]
)
colnames(Xvariants) <- c(
  pos_dm_kmer$idxnames,
  pos_dm_prot$idxnames,
  pos_dm_len$idxnames,
  pos_dm_len3$idxnames,
  pos_dm_freq$idxnames
)
  
# design matrix for the rare burden if present
if(nrow(reduced_signif_GWAS_rare) > 0){
  variant_frequency <- nb_alternative / nb_ind_per_gpos
  stopifnot(length(variant_frequency) == ncol(G))
  stopifnot(length(nb_gpos_per_ind) == sum(sub_ind))
  myfreq <- reduced_signif_GWAS_rare$frequency[which.min(reduced_signif_GWAS_rare$p_b)] # best frequency
  nb_rare_perind <- rowSums(G[sub_ind, which(variant_frequency < myfreq)], na.rm = T) / nb_gpos_per_ind # number of rare positions by individual divided by total number of positions
  full_design_matrix_rare <- as.matrix(nb_rare_perind, nrow = nind)
  colnames(full_design_matrix_rare) <- "rare"
  Xvariants <- cbind(Xvariants, full_design_matrix_rare) # cbind full_design_matrix_rare to Xvariants
}

#####                                               4. NOW RE-OPTIMISE FULL MODEL                                               #####

source("function_mixed_model.R")
# epi covariates design matrix (contained in the file "Gmatrix_data/prepared_PHENONAME_full_17.RData)
X <- as.matrix(X)

# misc. definitions for optimisation
npar1 <- length(Xrand)
lower1 <- c(rep(0, npar1))
upper1 <- c(rep(1, npar1))
init_fun1 <- function(){
  return(c(runif(npar1, 0, 1)))
}

# RE-OPTIMISE ON THE subset of individuals with non-NA at all these positions
# non_missing <- apply(Xvariants, 1, function(vec) all(!is.na(vec)))


# REMOVE GENETIC VARIANTS THAT ARE LINEAR COMBINATIONS OF OTHERS
# Use imputation when too many variants
if(ncol(Xvariants) > 1000) use_imputation <- TRUE else use_imputation <- FALSE
if(ncol(Xvariants) > nrow(Xvariants)){
  print("Too many variants... considering only kmer variants")
  use_kmer_only <- TRUE
  Xvariants <- Xvariants[, grepl("kmer", colnames(Xvariants))]
  if(ncol(Xvariants) > nrow(Xvariants)) stop("Too many variants...")
} else {
  use_kmer_only <- FALSE
}
if(use_imputation){
  Xvariants_imputed <- Xvariants
  Xvariants_imputed[is.na(Xvariants)] <- 0 # impute missing to 0
  tmp <- findLinearCombos(Xvariants_imputed)
  if(!is.null(tmp$remove)){
    to_remove <- tmp$remove
    cat("removing variants ", colnames(Xvariants)[to_remove], "because colinear with other variants\n")
    Xvariants <- Xvariants_imputed[, -to_remove]
  } else {
    Xvariants <- Xvariants_imputed
  }
} else {
  non_missing <- apply(Xvariants, 1, function(vec) all(!is.na(vec)))
  tmp <- findLinearCombos(Xvariants[non_missing, ])
  if(!is.null(tmp$remove)){
    to_remove <- tmp$remove
    cat("removing variants ", colnames(Xvariants)[to_remove], "because colinear with other variants\n")
    Xvariants <- Xvariants[, -to_remove]
  }
}
non_missing <- apply(Xvariants, 1, function(vec) all(!is.na(vec)))
non_missing_Xrand <- lapply(Xrand, function(mm) as.matrix(mm[non_missing, non_missing])) # conversion to matrix is important here
non_missing_nind <- sum(non_missing)
nvariants <- ncol(Xvariants)

if(nvariants > non_missing_nind - 100) stop("Too many variants included in the model")

# 0) model without variants
all_opt_full0 <- optim.fun.repeated(n.repeats = 3, lik.fun = mynegativeloglik1,
                               init.fun = init_fun1, verbose = T, lower = lower1, upper = upper1, control = list(tol = 0.001),
                               data = list(phenotype = suba[non_missing, pheno_name]),
                               X = X[non_missing, ],
                               n = non_missing_nind,
                               Xrand = non_missing_Xrand, npar = npar1)

# 1) model WITH all variants
all_opt_full1 <- optim.fun.repeated(n.repeats = 3, lik.fun = mynegativeloglik1,
                                    init.fun = init_fun1, verbose = T, lower = lower1, upper = upper1, control = list(tol = 0.001),
                                    data = list(phenotype = suba[non_missing, pheno_name]),
                                    X = cbind(X[non_missing, ], Xvariants[non_missing, ]),
                                    n = non_missing_nind,
                                    Xrand = non_missing_Xrand, npar = npar1)

# 2) model WITH variants but WITHOUT frequency variants
not_freq <- !grepl(pattern = "freq", x = colnames(Xvariants))
if(!use_kmer_only) all_opt_full2 <- optim.fun.repeated(n.repeats = 3, lik.fun = mynegativeloglik1,
                                    init.fun = init_fun1, verbose = T, lower = lower1, upper = upper1, control = list(tol = 0.001),
                                    data = list(phenotype = suba[non_missing, pheno_name]),
                                    X = cbind(X[non_missing, ], Xvariants[non_missing, not_freq]),
                                    n = non_missing_nind,
                                    Xrand = non_missing_Xrand, npar = npar1)

# 3) model WITH ONLY kmer / prot variants
only_kmer_prot <- grepl(pattern = "prot", x = colnames(Xvariants)) | grepl(pattern = "kmer", x = colnames(Xvariants))
if(!use_kmer_only) all_opt_full3 <- optim.fun.repeated(n.repeats = 3, lik.fun = mynegativeloglik1,
                                    init.fun = init_fun1, verbose = T, lower = lower1, upper = upper1, control = list(tol = 0.001),
                                    data = list(phenotype = suba[non_missing, pheno_name]),
                                    X = cbind(X[non_missing, ], Xvariants[non_missing, only_kmer_prot]),
                                    n = non_missing_nind,
                                    Xrand = non_missing_Xrand, npar = npar1)

system("say Just finished!")

opt_full0 <- clean_lik(all_opt_full0, npar1)
opt_full1 <- clean_lik(all_opt_full1, npar1)
if(!use_kmer_only) opt_full2 <- clean_lik(all_opt_full2, npar1)
if(!use_kmer_only) opt_full3 <- clean_lik(all_opt_full3, npar1)

# the below get_h2 function may produce an error for (spvl_adjusted, 0.2)
# Error in X %*% full_lik$beta : 
#   nécessite des arguments numériques/complexes matrice/vecteur
# Calls: get_h2
# caused by system being singular
# caused by colinearities in X matrix probably

h2_full0 <- get_h2(par = opt_full0$pars,
                  data = list(phenotype = suba[non_missing, pheno_name]),
                  X = X[non_missing, ],
                  n = non_missing_nind,
                  Xrand = non_missing_Xrand, npar = npar1,
                  lik_fun_name = "mynegativeloglik1", return_full = T)

h2_full1 <- get_h2(par = opt_full1$pars,
                   data = list(phenotype = suba[non_missing, pheno_name]),
                   X = cbind(X[non_missing, ], Xvariants[non_missing, ]),
                   n = non_missing_nind,
                   Xrand = non_missing_Xrand, npar = npar1,
                   lik_fun_name = "mynegativeloglik1", return_full = T)

if(!use_kmer_only) h2_full2 <- get_h2(par = opt_full2$pars,
                   data = list(phenotype = suba[non_missing, pheno_name]),
                   X = cbind(X[non_missing, ], Xvariants[non_missing, not_freq]),
                   n = non_missing_nind,
                   Xrand = non_missing_Xrand, npar = npar1,
                   lik_fun_name = "mynegativeloglik1", return_full = T)

if(!use_kmer_only) h2_full3 <- get_h2(par = opt_full3$pars,
                   data = list(phenotype = suba[non_missing, pheno_name]),
                   X = cbind(X[non_missing, ], Xvariants[non_missing, only_kmer_prot]),
                   n = non_missing_nind,
                   Xrand = non_missing_Xrand, npar = npar1,
                   lik_fun_name = "mynegativeloglik1", return_full = T)

#####                                      5. NOW GENERATE CONFIDENCE INTERVALS                                     #####

# In this part of the code, we generate confidence intervals on heritability
# This part comes from "patch_CI_on_heritability.R"
# We generate confidence intervals assuming the likelihood has approx the shape of a multivariate normal distribution
# We do so for each phenotype
# We do so for the model "0" without variants (defined in "3_reoptimise_model_v2.R")
# (Note we also have model "1" with all variants; model "2" with all variants but frequency variants)

# Hessian at the optimum for model 0:
h0 <- rootSolve::hessian(f = mynegativeloglik1, x = opt_full0$pars, pert = min(opt_full0$pars/10), # relatively larger perturbation (the effect of the 'pert' parameter on CI was systematically tested)
                        data = list(phenotype = suba[non_missing, pheno_name]),
                        X = X[non_missing, ],
                        n = non_missing_nind,
                        Xrand = non_missing_Xrand,
                        npar = npar1
)
h3 <- rootSolve::hessian(f = mynegativeloglik1, x = opt_full3$pars, pert = min(opt_full3$pars/10), # relatively larger perturbation (the effect of the 'pert' parameter on CI was systematically tested)
                         data = list(phenotype = suba[non_missing, pheno_name]),
                         X = cbind(X[non_missing, ], Xvariants[non_missing, only_kmer_prot]),
                         n = non_missing_nind,
                         Xrand = non_missing_Xrand,
                         npar = npar1
)

original_cov_matrix0 <- solve(h0) # covariance matrix is inverse of Hessian https://stats.stackexchange.com/questions/68080/basic-question-about-fisher-information-matrix-and-relationship-to-hessian-and-s?rq=1
original_cov_matrix3 <- solve(h3)

force_symmetric <- function(i, j, mymat){
  # force matrix to be symmetric:
  if(i==j) stop("i and j must be different")
  #if(abs(mymat[i, j]/mymat[j, i]) > 0.99 & abs(mymat[i, j]/mymat[j, i]) < 1.01){
    my_mean <- 0.5 * (mymat[i, j] + mymat[j, i])
    mymat[i,j] <- mymat[j, i] <- my_mean
  #} 
  mymat
}
cov_matrix0 <- original_cov_matrix0
cov_matrix3 <- original_cov_matrix3

for(i in 1:(ncol(cov_matrix0)-1)){
  for(j in ((i+1):ncol(cov_matrix0))){
    cov_matrix0 <- force_symmetric(i, j, mymat = cov_matrix0)
    cov_matrix3 <- force_symmetric(i, j, mymat = cov_matrix3)
  }
}

sds0 <- sqrt(diag(cov_matrix0))
ci_hessian0 <- rbind(opt_full0$pars - 1.96 * sds0, opt_full0$pars, opt_full0$pars + 1.96 * sds0)
rownames(ci_hessian0) <- c("2.5%", "ML", "97.5%")

sds3 <- sqrt(diag(cov_matrix3))
ci_hessian3 <- rbind(opt_full3$pars - 1.96 * sds3, opt_full3$pars, opt_full3$pars + 1.96 * sds3)
rownames(ci_hessian3) <- c("2.5%", "ML", "97.5%")

nrep <- 200

if(with_long_bootstrap <- TRUE){ # draw random parameters for CI (long version, takes 10 hours)
  
  rand_pars0 <- mvtnorm::rmvnorm(n = nrep, mean = opt_full0$pars, sigma = cov_matrix0)
  rand_pars3 <- mvtnorm::rmvnorm(n = nrep, mean = opt_full3$pars, sigma = cov_matrix3)
  
  t0 <- Sys.time()
  
  # vector of h2 in model 0 (no variant) and model 3 (kmer and prot variants)
  vec_h2_0 <- numeric(length = nrep)
  vec_h2_3 <- numeric(length = nrep)
  
  mat_fixed0 <- matrix(NA, nrow = nrep, ncol = nrow(h2_full0$fixed_effects))
  mat_fixed3 <- matrix(NA, nrow = nrep, ncol = nrow(h2_full3$fixed_effects))
  
  for(i in 1:nrep){
    
    print(i)
    mypars0 <- rand_pars0[i,]
    mypars3 <- rand_pars3[i,]
    
    if(mypars0[1]>0){ # for CD4_slope (low heritability) this sometimes goes smaller than 0 (in which case h2 is 0)
      tmp0 <- get_h2(par = mypars0, # heritability for randomly chosen parameters
                    data = list(phenotype = suba[non_missing, pheno_name]),
                    X = X[non_missing, ],
                    n = non_missing_nind,
                    Xrand = non_missing_Xrand, npar = npar1,
                    lik_fun_name = "mynegativeloglik1", return_full = T)
      vec_h2_0[i] <- tmp0$h2
      mat_fixed0[i, ] <- tmp0$fixed_effects[,1]
    }
    if(mypars3[1]>0){ # for CD4_slope (low heritability) this sometimes goes smaller than 0 (in which case h2 is 0)
      # same with model 3:
      tmp3 <- get_h2(par = mypars3, # heritability for randomly chosen parameters
                     data = list(phenotype = suba[non_missing, pheno_name]),
                     X = cbind(X[non_missing, ], Xvariants[non_missing, only_kmer_prot]),
                     n = non_missing_nind,
                     Xrand = non_missing_Xrand, npar = npar1,
                     lik_fun_name = "mynegativeloglik1", return_full = T)
      vec_h2_3[i] <- tmp3$h2
      mat_fixed3[i, ] <- tmp3$fixed_effects[,1]
    }
  }
  Sys.time() - t0
}

assign(x = paste0("ci_h2_0_", pheno_name), value = quantile(vec_h2_0, c(0.025, 0.975)))
assign(x = paste0("ci_h2_3_", pheno_name), value = quantile(vec_h2_3, c(0.025, 0.975)))

ci_fixed_mat0 <- apply(mat_fixed0, 2, quantile, c(0.025, 0.975), na.rm = T)
ci_fixed_mat3 <- apply(mat_fixed3, 2, quantile, c(0.025, 0.975), na.rm = T)

colnames(ci_fixed_mat0) <- rownames(tmp0$fixed_effects)
colnames(ci_fixed_mat3) <- rownames(tmp3$fixed_effects)

assign(x = paste0("ci_fixed_0_", pheno_name), value = ci_fixed_mat0)
assign(x = paste0("ci_fixed_3_", pheno_name), value = ci_fixed_mat3)

#####                                                 6. NOW SAVE THE RESULTS                                       #####

if(!file.exists("3_reoptimisation")) dir.create("3_reoptimisation/")
save.image(paste0("3_reoptimisation/re_optimisation_", pheno_name, "_", p_t, ".RData"))

############################################ 7. save table of fixed effects ################################################

# TODO: add CI to the fixed effects
# h2_full2 <- get_h2(par = opt_full2$pars,
#                    data = list(phenotype = suba[non_missing, pheno_name]),
#                    X = cbind(X[non_missing, ], Xvariants[non_missing, not_freq]),
#                    n = non_missing_nind,
#                    Xrand = non_missing_Xrand, npar = npar1,
#                    lik_fun_name = "mynegativeloglik1", return_full = T)

# ATTEMPTS TO VISUALIZE FIXED EFFECTS AND CI
# datavec2 <- suba[non_missing, pheno_name]
# Xmat2 <- cbind(X[non_missing, ], Xvariants[non_missing, not_freq])
# lm_ci2 <- lm(datavec2 ~ 0 + Xmat2)
# plot(h2_full2$fixed_effects, lm_ci2$coefficients) # "rare" are somewhat of an outlier here

if(!use_kmer_only){
  
  tmp0 <- data.frame(cbind(h2_full0$fixed_effects, rownames(h2_full0$fixed_effects)))
  tmp1 <- data.frame(cbind(h2_full1$fixed_effects, rownames(h2_full1$fixed_effects)))
  tmp2 <- data.frame(cbind(h2_full2$fixed_effects, rownames(h2_full2$fixed_effects)))
  tmp3 <- data.frame(cbind(h2_full3$fixed_effects, rownames(h2_full3$fixed_effects)))
  
  names(tmp0) <- names(tmp1) <- names(tmp2) <- names(tmp3) <- c("effect", "name")
  tmp0$model <- "epi"; tmp1$model <- "with_freq"; tmp2$model <- "no_freq"; tmp3$model <- "only_kmer_prot"
  tmp <- rbind(tmp0, tmp1, tmp2, tmp3)
  tmp$effect <- signif(as.numeric(tmp$effect), 3)
  tmp <- tmp[c("model", "name", "effect")]
  
  #####                           ADD CONFIDENCE INTERVALS                   #####
  
  # get sd of error for the corresponding linear model (without random effects)
  X0 <- X[non_missing, ]
  X3 <- cbind(X[non_missing, ], Xvariants[non_missing, only_kmer_prot])
  
  # compute CI, only for models 0 and 3 (the others not done)
  colnames(X0) <- rownames(h2_full0$fixed_effects)
  colnames(X3) <- rownames(h2_full3$fixed_effects)
  phenotype = suba[non_missing, pheno_name]
  
  lm0 <- lm(phenotype ~ 0 + X0)
  lm3 <- lm(phenotype ~ 0 + X3)
  
  sd_lm0 <- summary(lm0)[[4]][, "Std. Error"]
  sd_lm3 <- summary(lm3)[[4]][, "Std. Error"]
  
  sd_bs0 <- apply(mat_fixed0, 2, sd, na.rm = T) # sd of ML across bs values of genetic random effects
  sd_bs3 <- apply(mat_fixed3, 2, sd, na.rm = T) # sd of ML across bs values of genetic random effects
  
  # assume total error sd is approx the sum of these two
  #sd_tot0 <- sd_lm0 + sd_bs0
  #sd_tot3 <- sd_lm3 + sd_bs3
  
  tmp$error_lm_sd <- tmp$error_bs_sd <- NA
  tmp$error_lm_sd[tmp$model=="epi"] <- sd_lm0
  tmp$error_lm_sd[tmp$model=="only_kmer_prot"] <- sd_lm3
  tmp$error_bs_sd[tmp$model=="epi"] <- sd_bs0
  tmp$error_bs_sd[tmp$model=="only_kmer_prot"] <- sd_bs3
  tmp$error_sd <- tmp$error_lm_sd + tmp$error_bs_sd
  
  tmp$effect_lower <- tmp$effect - 1.96 * tmp$error_sd
  tmp$effect_upper <- tmp$effect + 1.96 * tmp$error_sd
  tmp$ci <- paste(signif(tmp$effect, 3), " [", signif(tmp$effect_lower, 3), "; ", signif(tmp$effect_upper, 3), "]", sep = "")
  
  # save
  tmp <- tmp[, c("model", "name", "ci", "effect", "error_lm_sd", "error_bs_sd", "error_sd", "effect_lower", "effect_upper")]
  write.csv(x = tmp, file = paste0("3_reoptimisation/fixed_effects_", pheno_name, "_", p_t, ".csv"), row.names = F, sep = ",")
  
  #plot(h2_full0$fixed_effects, fixed_effects) # just checking that they correspond
  
  
  #####                           DECOMPOSITION OF VARIANCE                   #####
  
  var_decomposition <- data.frame(name = c("phenotype", "total_variance", "h2_null_model", "h2_with_hits"), value = c(pheno_name, var(suba[non_missing, pheno_name]), h2_full0$h2, h2_full1$h2))
  var_decomposition$name <- as.character(var_decomposition$name); var_decomposition$value <- as.character(var_decomposition$value)
  
  # epidemiological effects
  colnames(X)
  X <- as.matrix(X)
  fixed_effects <- as.matrix(fixed_effects)
  h2_full1$fixed_effects <- as.matrix(h2_full1$fixed_effects)
  
  # epi effects
  for(i in 1:n_fixed_covariates){
    idx_in_X <- grep(pattern = fixed_covariate_names[i], x = colnames(X))
    #print(fixed_covariate_names[i])
    var_decomposition <- rbind(var_decomposition, c(
      paste0("fraction_variance_", fixed_covariate_names[i]),
      var(X[, idx_in_X, drop=F] %*% h2_full1$fixed_effects[idx_in_X, , drop = F]) / var(suba[, pheno_name])
    ))
  }
  var_decomposition <- rbind(var_decomposition, c(
       "fraction_variance_total_epi", var(X[, , drop=F] %*% h2_full1$fixed_effects[1:n_fixed_effects, , drop = F])/ var(suba[, pheno_name])
  ))
  
  # genetic effects
  alltypes <- c("kmer", "prot", "len", "len3", "freq", "rare")
  for(type in alltypes){
    idx_inX <- grepl(pattern = type, x = colnames(Xvariants))
    #full_design_matrix <- as.matrix(get(paste0("full_design_matrix_", type)))
  
    var_decomposition <- rbind(var_decomposition, c(
        paste0("fraction_variance_", type)
       ,
       var(Xvariants[non_missing, idx_inX, drop = F] %*%
             h2_full1$fixed_effects[n_fixed_effects+which(idx_inX), , drop = F]) / var(suba[non_missing, pheno_name])
      )
    )
  }
  if(nvariants > 0){
    var_decomposition <- rbind(var_decomposition, c(
      "fraction_variance_total_genetic", var(Xvariants[non_missing, , drop=F] %*% h2_full1$fixed_effects[n_fixed_effects + 1:nvariants, , drop = F])/ var(suba[non_missing, pheno_name])
    ))
  } else {
    var_decomposition <- rbind(var_decomposition, c("fraction_variance_total_genetic", 0))
  }
  
  
  for(type in alltypes){
    if(type != "rare"){
      variant_names <- get(paste0("pos_dm_", type))$posnames
    } else {
      variant_names <- "rare"
    }
    for(variant in variant_names){ # pos_dm_kmer$posnames
      idx_variant_inX <- grepl(pattern = variant, x = colnames(Xvariants))
      #full_design_matrix <- as.matrix(get(paste0("full_design_matrix_", type)))
      var_decomposition <- rbind(var_decomposition, c(
        paste0("fraction_variance_", variant)
        ,
        var(Xvariants[non_missing, idx_variant_inX, drop = F] %*%
              h2_full1$fixed_effects[n_fixed_effects+which(idx_variant_inX), , drop = F]) / var(suba[non_missing, pheno_name])
      ))
    }
  }
  
  # add R2 measures
  sum_square <- sum((suba[non_missing, pheno_name] - mean(suba[non_missing, pheno_name], na.rm = T))^2)
  
  pred_epi <- as.vector(X[non_missing, ] %*% h2_full0$fixed_effects) # original model without genetic effects
  pred_withgen <- as.vector(cbind(X[non_missing, ], Xvariants[non_missing, ]) %*% h2_full1$fixed_effects) # with genetic effects
  pred_withgen2 <- as.vector(cbind(X[non_missing, ], Xvariants[non_missing, not_freq]) %*% h2_full2$fixed_effects) # with genetic effects but not the frequency effects
  
  R2_epi <- 1 - sum((suba[non_missing, pheno_name] - pred_epi)^2) / sum_square # R2 for full model with fixed effects
  R2_withgen <- 1 - sum((suba[non_missing, pheno_name] - pred_withgen)^2) / sum_square # R2 for full model with fixed effects
  R2_withgen2 <- 1 - sum((suba[non_missing, pheno_name] - pred_withgen2)^2) / sum_square # R2 for full model with fixed effects
  
  var_decomposition <- rbind(var_decomposition, c("R2_epi", R2_epi))
  var_decomposition <- rbind(var_decomposition, c("R2_withgen", R2_withgen))
  var_decomposition <- rbind(var_decomposition, c("R2_withgen_nofreq", R2_withgen2))
  
  
  # save
  # h2 original model
  # h2 with genetic covariates
  # for each covariate save
  # name
  # variance explained
  
  var_decomposition$value2 <- sapply(var_decomposition$value, function(x) signif(as.numeric(x), 3))
  var_decomposition$value2[is.na(var_decomposition$value2)] <- var_decomposition$value[is.na(var_decomposition$value2)]
  
  write.csv(x = var_decomposition[, c("name", "value2")], file = paste0("3_reoptimisation/var_decomposition_", pheno_name, "_", p_t, ".csv"), row.names = F)
}
quit()



