rm(list = ls())
library(fitdistrplus)

define_filenames <- function(){
  if(grepl(pattern = "blanqua", getwd())) folder <- "BEEHIVE_Hackathon" else folder <- "BEEHIVE Hackathon"
  bonferroni_folder <- paste0("~/Dropbox (Infectious Disease)/", folder, "/Code/DevelopMethods/LMM/SolvingLMMinR/bonferroni/")
  data_folder <- paste0("~/Dropbox (Infectious Disease)/", folder, "/Code/DevelopMethods/LMM/SolvingLMMinR/GWAS_results/")
  bonferroni_files <- list.files(bonferroni_folder)
  
  all_phenos <- c("BEEHIVE_LVL", "spvl_normalised_adjusted", "spvl_adjusted", "CD4_slope")
  all_types <- c("prot", "kmer", "len", "len3", "freq")
  names_tab <- paste0("GWAS_", all_types) # names of tables for each type of test
  
  return(list(folder = folder, bonferroni_folder = bonferroni_folder, data_folder = data_folder, bonferroni_files = bonferroni_files, all_phenos = all_phenos, all_types = all_types, names_tab = names_tab))
}

########################################       FIRST SAVE GWAS-ES IN MINIMAL_DATASET       ##################################

if(write_minimal_dataset <- FALSE){
  
  # Loop mode or just focus on one pheno  / one type
  for(pheno_name in all_phenos){
    #pheno_name <- "BEEHIVE_LVL" 
    
    load(
      paste0("~/Dropbox (Infectious Disease)/", folder, "/Code/DevelopMethods/LMM/SolvingLMMinR/GWAS_results/firstGWAS_", pheno_name, "_full_17.RData")
    )
    # remove what is not necessary
    rm(list = ls()[grepl(pattern = "full_", x = ls())]) # remove the large objects starting with "full_"
    rm(list = ls()[!ls() %in% c("pheno_name", "all_phenos", "GWAS_prot", "GWAS_kmer","GWAS_len","GWAS_len3", "GWAS_freq", "GWAS_drm", "GWAS_rare",  "folder")])
    
    write.csv(x = GWAS_prot, file = paste0("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/GWAS_prot_", pheno_name, ".csv"), row.names = F)
    write.csv(x = GWAS_kmer, file = paste0("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/GWAS_kmer_", pheno_name, ".csv"), row.names = F)
    write.csv(x = GWAS_len, file = paste0("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/GWAS_len_", pheno_name, ".csv"), row.names = F)
    write.csv(x = GWAS_len3, file = paste0("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/GWAS_len3_", pheno_name, ".csv"), row.names = F)
    write.csv(x = GWAS_freq, file = paste0("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/GWAS_freq_", pheno_name, ".csv"), row.names = F)
    write.csv(x = GWAS_drm, file = paste0("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/GWAS_drm_", pheno_name, ".csv"), row.names = F)
    write.csv(x = GWAS_rare, file = paste0("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/GWAS_rare_", pheno_name, ".csv"), row.names = F)
  }
}

#####         REMOVE EVERYTHING BUT define_filenames function ####

rm(list = ls()[!ls() %in% "define_filenames"])
tmp <- define_filenames()
for(i in 1:length(tmp)) assign(x = names(tmp)[i], value = tmp[[i]])

#####         DEFINE COORDINATES OF MAIN GENES ####

gene_coords <- list(
  gag_start = 790,
  gag_end = 2289,
  pol_start = 2253,
  pol_end = 5093,
  env_start = 6225,
  env_end = 8792
)

#####         DEFINE FUNCTIONS        ####

source(paste0(bonferroni_folder, "get_bonferroni_dfs.R"))

draw_qqplot <- function(tab, limit_p, label, save_plot = F){
  if(!"p" %in% colnames(tab)) stop('tab has to include a column named p')
  punif <- (1:nrow(tab)) / nrow(tab)
  tab[tab$p ==0, "p"] <- 10^-16 # some are 0s
  max_p <- 1.03 * max(-log10(punif), -log10(tab$p), na.rm = T)
  if(save_plot) pdf(paste0(data_folder, "figures/", label, ".pdf"), width = 4, height = 4)
  par(mar = c(4,4,1,1))
  
  qqplot(-log10(punif), -log10(tab$p), type = "p", pch = 20, las = 1, cex = 0.5, xaxs = "i", yaxs = "i", xlim = c(0,max_p), ylim = c(0,max_p), xlab = expression(negative~log[10]~"p-value"~uniform), ylab = expression(negative~log[10]~"p-value"), main = "", bty = "n", axes = F)
  axis(1, at = seq(0, floor(max_p), 1))
  axis(2, at = seq(0, floor(max_p), 1), las = 1)
  segments(x0 = 0, y0 = 0, x1 = max_p, y1 = max_p)
  abline(h = limit_p, lty = 2)
  
  if(save_plot) dev.off()
}
draw_manhattan <- function(tab, signif_pos, eff_dfs, label, save_plot = F){
  
  #stopifnot(length(eff_dfs) == 6)
  if(save_plot) pdf(paste0(data_folder, "figures/", label, ".pdf"), width = 6, height = 4)
  par(mar = c(4,4,1,2), xpd = 1)
  cols <- RColorBrewer::brewer.pal(n = 8, name = "Set3") # color scheme for each kmer length (TODO this will recycle colors for length variants that can be > 6 in length)
  #for(pos in signif_pos){
    #zoomed_tab <- tab[which(tab$start_position_alignment > pos - 100 & tab$end_position_alignment < pos + 100),]
    zoomed_tab <- tab # do not zoom actually
    if((length(eff_dfs) == 6)) mycols <- sapply(zoomed_tab$end_position_alignment - zoomed_tab$start_position_alignment + 1, function(i) cols[i]) else mycols <- "cornflowerblue"
    ymax <- 1.1 * max(c(-log10(zoomed_tab$p), 1.1 * -log10(0.05 / eff_dfs)))
    xmax <- 1.01 * max(GWAS_tab$start_position_alignment)
    xmin <- min(GWAS_tab$start_position_alignment)
    plot(zoomed_tab$start_position_alignment, -log10(zoomed_tab$p), type = "p", pch = 20, xaxs = "i", yaxs = "i",
         xlab = paste0(""), # paste0("Position"),
         ylab = expression(paste("-", log[10],"(p-value)")),
         xlim = c(xmin, xmax), ylim = c(0, ymax), col = mycols, las = 1, cex= 0.5, axes = F)
    axis(side = 1, at = seq(round(xmin / 500) * 500, round(xmax / 500) * 500, 1500))
    axis(side = 2, at = seq(0, 8, 2), las = 1)
    
    # add significance level for each k of k-mer:
    if((length(eff_dfs) == 6)) {
      for(i in 1:length(eff_dfs)){
        segments(x0=xmin, x1=xmax, y0 = -log10(0.05 / eff_dfs[i]), y1 = -log10(0.05 / eff_dfs[i]), col = cols[i], lty = 2)
      }
    } else {
      segments(x0=xmin, x1=xmax, y0 = -log10(0.05 / eff_dfs[1]), y1 = -log10(0.05 / eff_dfs[1]), col = mycols, lty = 2)
    }
    if((length(eff_dfs) == 6)){
      # for k-mers add legend:
      legend(x = 4000, y = ymax, legend = paste0(1:6, "-mer"), pch = 20, col = cols, cex = 0.7, bty = "n")
      shape::Arrows(x0 = c(1413, 1514, 6570, 9008), x1 = c(1413, 1514, 6570, 9008), y0 = 8, y1 = 7.5, lwd = 3, arr.type = "triangle", col = c("gray", "black", "gray", "red"))
    }
    # add annotations of gene positions
    ybot <- -1.2; yw <- 0.4
    rect(xleft = gene_coords$gag_start, xright = gene_coords$gag_end, ybottom = ybot-yw, ytop = ybot, col = "gray", border = NA)
    rect(xleft = gene_coords$pol_start, xright = gene_coords$pol_end, ybottom = ybot-2*yw, ytop = ybot-yw, col = "gray", border = NA)
    rect(xleft = gene_coords$env_start, xright = gene_coords$env_end, ybottom = ybot-2*yw, ytop = ybot-yw, col = "gray", border = NA)
    text(x = 0.5 * (gene_coords$gag_start + gene_coords$gag_end), y = ybot-0.5*yw, labels = "gag")
    text(x = 0.5 * (gene_coords$pol_start + gene_coords$pol_end), y = ybot-1.5*yw, labels = "pol")
    text(x = 0.5 * (gene_coords$env_start + gene_coords$env_end), y = ybot-1.5*yw, labels = "env")
    
  #}
  if(save_plot) dev.off()
}
draw_rare <- function(tab, label, save_plot = F){
  if(save_plot) pdf(paste0(data_folder, "figures/", label, ".pdf"), width = 6, height = 4)
  if(save_plot) dev.off()
}
clean_table <- function(focal_tab){
  # remove empty rows and columns
  rows_allna <- apply(focal_tab, 1, function(vec) all(is.na(vec)))
  cols_allna <- apply(focal_tab, 2, function(vec) all(is.na(vec)))
  focal_tab <- focal_tab[!rows_allna,]
  focal_tab <- focal_tab[, !cols_allna]
  
  return(focal_tab)
}
df_inflation <- data.frame(expand.grid(all_phenos, all_types, 1:6))
names(df_inflation) <- c("phenos", "type", "k")
df_inflation$inflation <- NA
  
for(pheno_name in all_phenos){

  GWAS_prot <- read.csv(file = paste0("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/GWAS_prot_", pheno_name, ".csv"))
  GWAS_kmer <- read.csv(file = paste0("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/GWAS_kmer_", pheno_name, ".csv"))
  GWAS_len <-  read.csv(file = paste0("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/GWAS_len_", pheno_name, ".csv"))
  GWAS_len3 <- read.csv(file = paste0("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/GWAS_len3_", pheno_name, ".csv"))
  GWAS_freq <- read.csv(file = paste0("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/GWAS_freq_", pheno_name, ".csv"))
  GWAS_drm <- read.csv(file = paste0("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/GWAS_drm_", pheno_name, ".csv"))
  GWAS_rare <- read.csv(file = paste0("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/GWAS_rare_", pheno_name, ".csv"))

  # type <- "kmer"
  for(type in all_types){
    
    # get table of p values for focal phenotype and type of GWAS
    mytab <- paste0("GWAS_", type)
    GWAS_tab <- get(mytab)
    GWAS_tab$p[GWAS_tab$p == 0] <- 1e-16 # set 0 to 1e-16 (the limit)
    
    # check there is no 'pool' allele (only 'pooled')
    stopifnot(!"pool" %in% GWAS_tab$reference)
    for(colname in colnames(GWAS_tab)[grepl(pattern = "alternative", x = colnames(GWAS_tab))]){
      stopifnot(!"pool" %in% GWAS_tab[, colname])
    }
    stopifnot(!anyDuplicated(GWAS_tab))
    
    na_lik <- is.na(GWAS_tab$lik0)
    if(any(na_lik)){
      print("eliminating NA likelihood, rows:")
      print(GWAS_tab[which(na_lik),])
      GWAS_tab <- GWAS_tab[!na_lik,]
    }
    
    if(type == "kmer") n_mers <- 6 else n_mers <- 1
    all_GWAS_tab_focalkmer_significant <- c() # a table to gather all kmer significant positions
    all_eff_dfs <- c()
    
    for(k in 1:n_mers){ # loop on "k" in k-mers, if kmer GWAS
      
      mylabel <- paste0(pheno_name, "_", mytab, "_", k)
      print(mylabel)
      
      if(mytab=="GWAS_kmer") {
        GWAS_tab_focalkmer <- GWAS_tab[GWAS_tab$end_position_alignment - GWAS_tab$start_position_alignment == k - 1, ]
        eff_df <- get_bonferroni_dfs(pheno = pheno_name, type = type)$all_dfs[[k]]
      } else {
        GWAS_tab_focalkmer <- GWAS_tab
        eff_df <- get_bonferroni_dfs(pheno = pheno_name, type = type)$total_dfs
      }
      
      all_eff_dfs <- c(all_eff_dfs, eff_df)
      inflation <- qchisq(median(GWAS_tab_focalkmer$p, na.rm = T), 1, lower.tail = FALSE) / qchisq(0.5, 1)
      cat("effective dfs", eff_df, ", inflation ", inflation, "\n")
      df_inflation$inflation[df_inflation$phenos == pheno_name & df_inflation$type == type & df_inflation$k == k] <- inflation
      
      draw_qqplot(tab = GWAS_tab_focalkmer, limit_p = -log10(0.05 / eff_df), label = paste("qqplot", mylabel, sep = "_"), save_plot = T)
      
      significant_rows <- which(GWAS_tab_focalkmer$p < 0.05 / eff_df)
      
      GWAS_tab_focalkmer$p_bonferroni <- GWAS_tab_focalkmer$p * eff_df
      GWAS_tab_focalkmer$effective_tests <- eff_df
      GWAS_tab_focalkmer$n_tests <- nrow(GWAS_tab_focalkmer)
      GWAS_tab_focalkmer$inflation <- inflation
      
      if(length(significant_rows) > 0){
        
        GWAS_tab_focalkmer_significant <- GWAS_tab_focalkmer[significant_rows, ]
        all_GWAS_tab_focalkmer_significant <- rbind(all_GWAS_tab_focalkmer_significant, GWAS_tab_focalkmer_significant) # append the significant positions to a table that will contain significant positions for all kmers
        GWAS_tab_focalkmer_significant <- clean_table(GWAS_tab_focalkmer_significant)
        tmp <- GWAS_tab_focalkmer_significant[order(GWAS_tab_focalkmer_significant$p),]
        print(tmp)
        write.csv(x = tmp, file = paste0(data_folder, mylabel, ".csv"))
        
        # plot(GWAS_tab_focalkmer$start_position_alignment, -log10(GWAS_tab_focalkmer$p), type = "l", main = mylabel)
        # abline(h = -log10(0.05 / eff_df))

      }
      
      
      #cols_to_show <- c('lik0', 'lik1', paste0('effect',1:5), 'df', 'start_position_alignment', 'end_position_alignment',
      #                  'reference', paste0('alternative', 1:5), paste0('N', 1:5), "N", "gamma", "p", "p_bonferroni", "effective_tests", "n_tests", "inflation")
      #cols_to_show <- cols_to_show[cols_to_show %in% colnames(GWAS_tab_focalkmer)]
      #GWAS_tab_focalkmer_significant <- GWAS_tab_focalkmer[significant_rows, ]
      
      # if(nrow(GWAS_tab_focalkmer_significant) > 0){
      #   tmp <- GWAS_tab_focalkmer_significant[order(GWAS_tab_focalkmer_significant$p),]
      # } else {
      #   tmp <- GWAS_tab_focalkmer_significant
      # }
      
      #if(pheno_name == "BEEHIVE_LVL" & mytab == "GWAS_freq")
      
      
    } # end of loop on number of kmers
    
      # clean the table containing all significant positions:
      if(!is.null(all_GWAS_tab_focalkmer_significant)){
        all_GWAS_tab_focalkmer_significant <- clean_table(all_GWAS_tab_focalkmer_significant)
      }
        
      # draw manatthan plot of all positions
      mylabel <- paste0(pheno_name, "_", mytab)
      draw_manhattan(tab = GWAS_tab, # the overall table
                     signif_pos = sort(unique(all_GWAS_tab_focalkmer_significant$start_position_alignment)), # significant positions for kmer length 1 to 6
                     eff_dfs  = all_eff_dfs, # effective dfs for each kmer (1 to 6)
                     label = paste("manhattan", mylabel, sep = "_"), save_plot = T)
  } # end of loop on type of GWAS
  # deal with GWAS_drm and GWAS_rare
  for(GWAS_type in c("GWAS_drm", "GWAS_rare")){
    mylabel <- paste0(pheno_name, "_", GWAS_type)
    GWAS_tab_focalkmer <- get(GWAS_type)
    significant_rows <- which(GWAS_tab_focalkmer$p < 0.05)
    GWAS_tab_focalkmer_significant <- GWAS_tab_focalkmer[significant_rows,]
    write.csv(x = GWAS_tab_focalkmer_significant, file = paste0(data_folder, mylabel, ".csv"))
  }
  
  
  # fit distribution to frequency for further simulations
  singleSNP <- GWAS_kmer[which(GWAS_kmer$end_position_alignment==GWAS_kmer$start_position_alignment), ]
  singleSNP$freq <- singleSNP$N1 / singleSNP$N
  print(min(singleSNP$freq))
  myfreq <- singleSNP$freq;  myfreq <- myfreq-(0.004444444 - 1e-3); myfreq <- myfreq * 2 # transform frequency
  f <- fitdist(data = myfreq, distr = "beta") # 0.5160741 3.1197464 
  hist(myfreq, breaks = 100, freq = F)
  
} # end of loop on phenotype

df_inflation <- df_inflation[!is.na(df_inflation$inflation), ]
write.csv(df_inflation, file = paste0(data_folder, "inflation_table.csv"), row.names = F)

df_inflation -> inf

if(compute_h2 <- F){
  
  dim(full_design_matrix_freq)
  dim(full_ind_tokeep_freq)
  stopifnot(sum(GWAS_freq$df)==ncol(full_design_matrix_freq))
  wd <- "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/"
  setwd(paste0(wd, "Code/DevelopMethods/LMM/SolvingLMMinR/"))
  RDatafile_name <- "prepared_BEEHIVE_LVL_full_17.RData"
  load(paste0("Gmatrix_data/", RDatafile_name))
  
  # CAREFUL THE FOLLOWING IS WORKING ONLY WHEN 1 MINORITY FREQUENCY VARIANT IS FITTED FOR THE FOCAL POSITION
  for(my_pos in c(4551, 7513, 9085)){
    
    idx_in_tab <- which(GWAS_freq$start_position_alignment==my_pos)
    if(idx_in_tab==1) stop("not working for idx_in_tab=1")
    my_df <- GWAS_freq$df[idx_in_tab]
    idx_in_design <- (sum(GWAS_freq$df[1:(idx_in_tab-1)]) + 1):(sum(GWAS_freq$df[1:(idx_in_tab-1)]) + my_df)
    my_design <- full_design_matrix_freq[, idx_in_design, drop = F] # design matrix
    
    # individuals to keep only
    my_tokeep <- full_ind_tokeep_freq[, idx_in_tab, drop = F]
    stopifnot(all(!is.na(my_design[my_tokeep,])))
    stopifnot(all(colSums(my_design[my_tokeep, , drop = F] > 0) == GWAS_freq[idx_in_tab, paste0("N", 1:my_df)])) # check the numners are matching
    
    # heritability of this position
    print(
      var(my_design[my_tokeep, , drop=F] %*% t(GWAS_freq[idx_in_tab, paste0("effect", 1:my_df)])) / var(suba[my_tokeep, "BEEHIVE_LVL"])
    )
  }
  #my_pos <- 7513
}

# # "negative control" with non-significant positions
# for(my_pos in sample(GWAS_freq$start_position_alignment, size = 10)){
#   
#   idx_in_tab <- which(GWAS_freq$start_position_alignment==my_pos)
#   my_df <- GWAS_freq$df[idx_in_tab]
#   idx_in_design <- (sum(GWAS_freq$df[1:(idx_in_tab-1)]) + 1):(sum(GWAS_freq$df[1:(idx_in_tab-1)]) + my_df)
#   my_design <- full_design_matrix_freq[, idx_in_design, drop = F] # design matrix
#   
#   # individuals to keep only
#   my_tokeep <- full_ind_tokeep_freq[, idx_in_tab, drop = F]
#   stopifnot(all(!is.na(my_design[my_tokeep,])))
#   stopifnot(all(colSums(my_design[my_tokeep, , drop = F] > 0) == GWAS_freq[idx_in_tab, paste0("N", 1:my_df)])) # check the numners are matching
#   
#   # heritability of this position
#   print(
#     var(my_design[my_tokeep, , drop=F] %*% t(GWAS_freq[idx_in_tab, paste0("effect", 1:my_df)])) / var(suba[my_tokeep, "BEEHIVE_LVL"])
#   )
# }

