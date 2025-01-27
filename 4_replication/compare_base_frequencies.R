rm(list = ls())
# load the replication data
load("~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/Gmatrix_data/base_frequencies.RData")

# rename replication objects with suffix "_rep"
for(oo in ls()){
  assign(x = paste0(oo, "_rep"), value =  get(oo)) # rename
  rm(list = c(oo))
}

#... and load the discovery
load("~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/Gmatrix_data/base_frequencies.RData")

#... and load the positions of interest: those that are significant at the 0.05 level in the within-host frequency GWAS
pheno_name <- "BEEHIVE_LVL"
GWAS_type <- "freq"
disc <- read.csv(paste0("~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/GWAS_results/", pheno_name, "_GWAS_", GWAS_type, "_1.csv"))

# compare freq GWAS significant positions between replication and discovery dataset
# for a focal position 2398, for example, this is:
idx <- which(window.ends_freq_rep == 2398 & window.starts_freq_rep == 2398) # index of this position in the replication dataset
apply(G_freq_rep[, idx], 2, summary) # what's at this position in replication dataset
apply(G_freq_rep[, idx] > 0.001, 2, sum, na.rm = T) # how many individuals have greater than 0.001 frequency at that position  
reference_freq_rep[idx] # what is the allelic diversity  
alternative_freq_rep[idx]
# ... now do it systematically

########################## SYSTEMATICALLY COMPARE ALLELE FREQUENCIES AT HIGHLY SIGNIFICANT POSITIONS, IN GENOMES THAT ARE BOTH IN DISCOVERY AND REPLICATION DATASET ##########################

idsG_freq_both <- intersect(idsG_freq, idsG_freq_rep) # BEEHIVE IDs which are in both datasets (134 individuals)
# index of these people in G:
i   <-   match(idsG_freq_both, idsG_freq)
i_rep <- match(idsG_freq_both, idsG_freq_rep)

stopifnot(all(idsG_freq[i]==idsG_freq_rep[i_rep]))

pdf("~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/within_host_frequencies_significantpositions.pdf", paper = "a4", width = 0, height = 0)
par(mfrow = c(4, 3))

# now look at specific positions
for(my_pos in disc$start_position_alignment){ # for each of the 12 significant positions
  
  idx_rep <- which(window.ends_freq_rep == my_pos & window.starts_freq_rep == my_pos)
  idx <- which(window.ends_freq == my_pos & window.starts_freq == my_pos)
  stopifnot(all(idx_rep == idx)) # positions should be the same
  
  # select subsets of G_freq matrix corresponding to these people, this position
  G <- G_freq[i, idx]
  G_rep <- G_freq_rep[i_rep, idx]
  
  
  # select alternative alleles present in both datasets
  alt_both <- intersect(alternative_freq[idx], alternative_freq_rep[idx]) 
  n_alt_alleles <- length(alt_both) # number of alternative alleles common to the two datasets
  j <- match(alt_both, alternative_freq[idx]) # position of these common alternative alleles
  j_rep <- match(alt_both, alternative_freq_rep[idx])
  
  for(k in 1:n_alt_alleles){
    plot(G[, j[k]], G_rep[, j_rep[k]], pch = 20, main = paste("position", my_pos, "- allele", alt_both[k], "(ref", reference_freq[idx][1], ")"), xlab = "within-host frequency in old genomes", ylab = "within-host frequency in new genomes")
    abline(0, 1)
    lm0 <- lm(G_rep[, j_rep[k]] ~ G[, j[k]])
    if(!is.na(lm0$coefficients[1]) & !is.na(lm0$coefficients[2]) & lm0$coefficients[2] != 0){
      abline(lm0$coefficients[1], lm0$coefficients[2], col = "red", lwd = 1)
    }
  }
}
dev.off()

########################## SYSTEMATICALLY COMPARE ALL ALLELE FREQUENCIES IN GENOMES THAT ARE BOTH IN DISCOVERY AND REPLICATION DATASET ##########################

# create a dataframe to store results
corr <- data.frame(pos = rep(NA, length(window.starts_freq_rep) * 5)) # each row is a position x allele (common to both sequencing runs) combination
corr$var_1 <- NA
corr$var_2 <- NA
corr$N_1 <- NA
corr$N_2 <- NA
corr$int <- NA
corr$beta <- NA
corr$R2 <- NA
corr$p <- NA
corr$allele <- NA
l <- 1
# now look at all positions
for(my_pos in unique(window.starts_freq_rep)){
  # for each position
  
  idx <- which(window.ends_freq == my_pos & window.starts_freq == my_pos)
  idx_rep <- which(window.ends_freq_rep == my_pos & window.starts_freq_rep == my_pos)
  
  stopifnot(all(idx_rep == idx)) # should be the same
  
  # select subsets of G_freq matrix corresponding to these people, this position
  G <- G_freq[i, idx]
  G_rep <- G_freq_rep[i_rep, idx]
  
  
  # select alternative alleles in both
  alt_both <- intersect(alternative_freq[idx], alternative_freq_rep[idx])
  n_alt_alleles <- length(alt_both) # number of alternative alleles common to the two datasets
  j <- match(alt_both, alternative_freq[idx]) # position of these alternative alleles that are common to both sequencing runs
  j_rep <- match(alt_both, alternative_freq_rep[idx])
  
  for(k in 1:n_alt_alleles){
    #plot(G[, j[k]], G_rep[, j_rep[k]], pch = 20, main = paste("position", my_pos, "- allele", alt_both[k], "(ref", reference_freq[idx][1], ")"), xlab = "within-host frequency in old genomes", ylab = "within-host frequency in new genomes")
    #abline(0, 1)
    
    var_1 <- var(G[, j[k]], na.rm = T)
    var_2 <- var(G_rep[, j_rep[k]], na.rm = T)
    
    mean_1 <- mean(G[, j[k]], na.rm = T)
    mean_2 <- mean(G_rep[, j_rep[k]], na.rm = T)
    
    N_1 <- sum(G[, j[k]] > 0, na.rm = T) # number of individuals with non-zero frequency of this allele
    N_2 <- sum(G_rep[, j_rep[k]] > 0, na.rm = T) 
    
    N_1_low <- sum(G[, j[k]] < 0.1, na.rm = T) # number of individuals with low frequency of this allele
    N_2_low <- sum(G_rep[, j_rep[k]] < 0.1, na.rm = T) 

    N_1_high <- sum(G[, j[k]] > 0.9, na.rm = T) # number of individuals with high frequency of this allele
    N_2_high <- sum(G_rep[, j_rep[k]] > 0.9, na.rm = T)
   
    
    corr[l, c("pos", "allele", "mean_1", "mean_2", "N_1", "N_2", "N_1_low", "N_2_low", "N_1_high", "N_2_high", "var_1", "var_2")] <- c(my_pos, alt_both[k], mean_1, mean_2, N_1, N_2, N_1_low, N_2_low, N_1_high, N_2_high, var_1, var_2)
    
    if(var_1 > 0 & var_2 > 0 &  N_1 > 5 & N_2 > 5){ # when at least 5 individuals are non 0 in both sequencing runs
      lm0 <- lm(G_rep[, j_rep[k]] ~ G[, j[k]]) # correlation between the two frequencies
      if(all(!is.na(lm0$coefficients))){
        corr[l, c("int", "beta", "R2", "p")] <- c(lm0$coefficients[1], lm0$coefficients[2], summary(lm0)$r.squared, summary(lm0)[[4]]["G[, j[k]]", "Pr(>|t|)"])
        l <- l + 1
      } else {
        corr[l, c("int", "beta", "R2", "p")] <- c(lm0$coefficients[1], lm0$coefficients[2], summary(lm0)$r.squared, NA)
        l <- l + 1
      }
      
    } else {
      l <- l + 1
    }
  }
}
for(col in c("var_1",  "var_2",  "int",    "beta",   "R2",     "p", "mean_1", "mean_2")){
  # var_1 and var_2 are the variance in allele frequency
  corr[, col] <- as.numeric(corr[, col])
}
which_na <- apply(corr, 1, function(x) all(is.na(x)))
write.csv(corr[!which_na,], "~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/all_correlations_withinhostfrequencies.csv", row.names = F, col.names = T)
corr <- corr[!which_na,]
plot(corr$R2, corr$var_1, pch = ".", ylab= "variance of within-host frequency across individuals", xlab = "R2 correlation between frequencies")
plot(corr$R2, corr$beta, pch = ".", xlim = c(0, 1), ylim = c(0,1), xlab = "R2 correlation between frequencies", ylab= "regression coefficient between frequencies")

hist(corr$R2[which(corr$N_1_high == 0 & corr$N_2_high == 0)], breaks = 100, main = "distribution R2 for those without high frequency") # without the sites that have frequency 1, very few high correlations
View(corr[which(corr$R2 > 0.8 & corr$N_1_high == 0 & corr$N_2_high == 0),]) # positions with high correlation that are not caused by having just 0s and 1s within-host frequency

#### zoom on one of the positions where correlation exists (for alternative allele A that is the 1st position in the alternative alleles) #### 
my_pos <- 8679
idx <- which(window.ends_freq == my_pos & window.starts_freq == my_pos)
idx_rep <- which(window.ends_freq_rep == my_pos & window.starts_freq_rep == my_pos)
stopifnot(all(idx_rep == idx)) # should be the same
# select subsets of G_freq matrix corresponding to these people, this position
G <- G_freq[i, idx]
G_rep <- G_freq_rep[i_rep, idx]
plot(G[, 1], G_rep[, 1])
abline(0,1)

