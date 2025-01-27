rm(list = ls())

get_power <- function(mydisc, myphenovar, nrep = 1000){
  nvar <- nrow(mydisc)
  powers <- c()
  powers1 <- c()
  alternatives <- c()
  for(i in 1:nvar){
    
    # check the alternative alleles are the same... at least the first one
    truc1 <- as.character(unlist(disc_withrep[i, paste0("alternative", 1:nallelemax)]))
    truc2 <- as.character(unlist(disc_withrep[i , paste0("rep_alternative", 1:nallelemax)]))
    is_common <- intersect(truc1, truc2)
    is_common <- is_common[is_common!="pooled" & !is.na(is_common)]
    #is_identical <- which(truc1 == truc2)
    
    if(length(is_common) == 0) {
      powers <- c(powers, NA)
      powers1 <- c(powers1, NA)
      alternatives <- c(alternatives, NA)
      next
    }
    nalternative <- length(is_common)
    #if(length(is_identical) == 1) if(is_identical == 1) nalternative <- 1
    #if(length(is_identical) == 2) if(all(is_identical == c(1, 2))) nalternative <- 2
    #if(length(is_identical) == 3) if(all(is_identical == c(1, 2, 3))) nalternative <- 3
    #if(is.na(nalternative)) stop("pattern of allele identity between discovery and replication is weird")
    idx_in_disc <- match(is_common, truc1)
    idx_in_rep <- match(is_common, truc2)
    if(length(idx_in_disc) == 0) stop()
    if(length(idx_in_rep) == 0) stop()
    
    Ns <- mydisc[i, paste0("rep_N", c("", idx_in_rep))]
    #Ns <- Ns[!is.na(Ns)]
    #nalternative <- length(Ns)
    #stopifnot(nalternative >= 2)
    alleles <-  unlist(mydisc[i, c("rep_reference", paste0("rep_alternative", idx_in_rep))])
    effects <- unlist(mydisc[i, paste0("effect", idx_in_disc)]) # effect of these alleles from discovery
    
    # create fake data; allele status and mean is always the same...
    fakedata <- rep(x = alleles[1], Ns[1] - sum(Ns[2:(nalternative+1)]))
    for(j in 2:(nalternative+1)) fakedata <- c(fakedata, rep(x = alleles[j], Ns[j]))
    fakepheno <- rep(0, Ns[1] - sum(Ns[2:(nalternative+1)]))
    for(j in 2:(nalternative+1)) fakepheno <- c(fakepheno, rep(x = effects[j-1], Ns[j]))
    stopifnot(length(fakedata) == Ns[1])
    stopifnot(length(fakepheno) == Ns[1])
    
    n_signif <- 0
    n_signif1 <- 0
    ps <- c()
    ps1 <- c()
    for(k in 1:nrep){
      # ... but some noise is added
      fakepheno_withnoise <- fakepheno + rnorm(n = unlist(Ns[1]), mean = 0, sd = sqrt(myphenovar))
      ll1 <- as.numeric(logLik(lm1 <- lm(fakepheno_withnoise ~ as.factor(fakedata))))
      ll0 <- as.numeric(logLik(lm0 <- lm(fakepheno_withnoise ~ 1)))
      gamma <- 2*(ll1 - ll0) # twice the diff in log-lik
      p <- 1- pchisq(gamma, df = nalternative)
      if(p <= 0.05) n_signif <- n_signif + 1
      ps <- c(ps, p)
      
      if(nalternative == 1){ # if only one alternative allele apply a 1-sided t-test
        tvalue <- summary(lm1)[[4]][2, "t value"]
        p1 <- 2 * pt(abs(tvalue), df = unlist(Ns[1])-nalternative-1, lower.tail = F) # the two-sided t-test
        p1 <- pt(abs(tvalue), df = unlist(Ns[1])-nalternative-1, lower.tail = F) # the one-sided t-test -> simply half the p-value
        ps1 <- c(ps1, p1)
        if(p1 <= 0.05) n_signif1 <- n_signif1 + 1
      }
    }
    #hist(ps, breaks = 100)
    power <- n_signif / nrep
    power1 <- n_signif1 / nrep
    powers <- c(powers, power)
    powers1 <- c(powers1, power1)
    alternatives <- c(alternatives, nalternative)
    
  }
  return(
    cbind(powers, powers1, alternatives)
  )
}

for(pheno_name in c("BEEHIVE_LVL", "spvl_normalised_adjusted", "spvl_adjusted", "CD4_slope")){
  
  myphenovar <- as.numeric(as.character(read.csv(paste0("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/3_reoptimisation/var_decomposition_", pheno_name, "_0.05.csv"))[2, 2])) # variance of this phenotype
  
  # load replication dataset
  load(paste0("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/GWAS_results/firstGWAS_", pheno_name, "_full_17.RData"))
  
  types <- c("kmer_1", "kmer_2", "kmer_3", "kmer_4", "kmer_5", "kmer_6", "len_1", "len3_1") # forget about prot for now
  
  for(GWAS_type in types){
    # NEED TO LOOP ON KMER LENGTH AS WELL
    print(pheno_name)
    print(GWAS_type)
    # load discovery results:
    tryCatch(
      {
        myfile <- paste0("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/GWAS_results/", pheno_name, "_GWAS_", GWAS_type, ".csv")
        if(file.exists(myfile)) disc <- read.csv(myfile) else stop("file does not exist")
        
        # get the appropriate table of the GWAS 
        if(grepl(pattern = "kmer", GWAS_type)){
          my_GWAS <- get("GWAS_kmer")
          kmer_length <- as.numeric(gsub(pattern = "kmer_", replacement = "", GWAS_type))
        } else {
          kmer_length <- 1
        }
        if(GWAS_type=="len_1") my_GWAS <- get("GWAS_len")
        if(GWAS_type=="len3_1") my_GWAS <- get("GWAS_len3")
  
        my_GWAS$start_end <- paste(my_GWAS$start_position_alignment, my_GWAS$end_position_alignment, sep = "_")
        disc$start_end <- paste(disc$start_position_alignment, disc$end_position_alignment, sep = "_")
        
        # modify column names to avoid duplicate column names when we bind the two tables
        colnames(my_GWAS) <- gsub(pattern = "reference", replacement = "rep_reference", x = colnames(my_GWAS))
        colnames(my_GWAS) <- gsub(pattern = "alternative", replacement = "rep_alternative", x = colnames(my_GWAS))
        colnames(my_GWAS) <- gsub(pattern = "N", replacement = "rep_N", x = colnames(my_GWAS))
        
        # HERE NEED TO ADJUST TO NUMBER OF ALLELES
        nallelemax <- sum(grepl(pattern = "N[0-9]", x = names(disc)))
        
        index_in_GWAS <- match(disc$start_end, my_GWAS$start_end)
        
        if(!is.na(index_in_GWAS)){
          disc_withrep <- cbind(disc, my_GWAS[index_in_GWAS, c("rep_reference", paste0("rep_alternative", 1:nallelemax), "rep_N", paste0("rep_N", 1:nallelemax))])
          print(disc_withrep)
          if(any(disc_withrep[, "reference"] != disc_withrep[, "rep_reference"])) stop("not same references")
          
          tmp_power <- get_power(mydisc = disc_withrep, myphenovar = myphenovar)
          stopifnot(nrow(tmp_power) == nrow(disc_withrep))
          disc_withrep$power <- tmp_power[,1]
          disc_withrep$power1 <- tmp_power[,2]
          disc_withrep$nalternative <- tmp_power[,3]
          
          print("writing down the file in GWAS_results with replication power:")
          stop("pause here")
          write.csv(x = disc_withrep, file = gsub(x = myfile, pattern = ".csv", replacement = "_withreplicationpower.csv"), row.names = F)
        } else {
          print("index of these positions not found in the GWAS... not returning any result regarding power")
        }

      }, error = function(e){
        print(e)
      }
    )
  }
  
} # end of loop on phenotype

# NOW COMPUTE POWER IN REPLICATION DATASET

# loop on replicates
# (1) generate a fake dataset
# (2) re-do inference on this fake dataset
# (3) record p-value and effect sizes 

# disc_withrep <- disc_withrep[, c("lik0", "lik1", "effect1", "effect2", "effect3", "df", "start_position_alignment", "end_position_alignment", "reference", "alternative1", "alternative2", "alternative3", "N1", "N2", "N3", "N",
#                                  "rep_reference", "rep_alternative1", "rep_alternative2", "rep_alternative3", "rep_N1", "rep_N2", "rep_N3", "rep_N", "p", "p_bonferroni", "effective_tests", "n_tests", "inflation", "gamma")]

#write.csv(x = disc_withrep, file = paste0("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/power_replication_", pheno_name, ".csv"), row.names = F, col.names = T)

