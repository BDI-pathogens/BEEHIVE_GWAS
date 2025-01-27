get_bonferroni_dfs <- function(pheno, type, lf = bonferroni_files, bonferroni_folder = "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/bonferroni/"){
  
  if(type == "len"){
    sel_neffective <- grepl(pattern = "n_effective", x = lf) & grepl(pattern = pheno, x = lf) & grepl(pattern = type, x = lf) &
      !grepl(pattern = "len3", x = lf)  # also eliminate the len3 which might be selected
  } else {
    sel_neffective <- grepl(pattern = "n_effective", x = lf) & grepl(pattern = pheno, x = lf) & grepl(pattern = type, x = lf)
  }
  n_effective_files <- lf[sel_neffective]
  
  if(length(n_effective_files) == 0){ # effective df file not present -> SET TO DEFAULTS_DFS (WHICH ARE NOT DEFINED NOW)
    warning(paste("effective df set to default of ", default_dfs, ", for phenotype", pheno, "GWAS on", type))
    all_dfs <- NA
    total_dfs <- default_dfs
  } else {
    list_of_dfs <- lapply(n_effective_files, function(filename) read.csv(paste0(bonferroni_folder, filename), header = T))
    if(type == "freq"){
      all_dfs <- unlist(lapply(list_of_dfs, function(ll) ll$n_tests[1])) # USE TOTAL NUMBER OF TESTS FOR FREQ VARIANTS
    } else {
      all_dfs <- unlist(lapply(list_of_dfs, function(ll) ll$n_effective_tests[1]))
    }
    
    total_dfs <- sum(all_dfs)
  }
  return(list(all_dfs = all_dfs, total_dfs = total_dfs))
}