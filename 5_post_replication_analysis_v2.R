rm(list = ls())

if(write_minimal_dataset <- TRUE){
  
  # this code compares effect sizes of variants inferred in discovery and replication datasets
  # and draws several figures
  
  tab_res <- c() # a table storing results of main lm correlating effects
  tab_res2 <- c() # same but filtering by small p-value in main analysis
  
  hits <- c() # a table storing the estimated effects for hits
  phenos <-  c("BEEHIVE_LVL", "spvl_adjusted", "spvl_normalised_adjusted", "CD4_slope")
  
  for(mytrait in phenos){
  
    load(paste0("~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/GWAS_results/firstGWAS_", mytrait, "_full_17.RData"))
    GWAS_disc <- GWAS_kmer
    load(paste0("~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/GWAS_results/firstGWAS_", mytrait, "_full_17.RData"))
    GWAS_rep <- GWAS_kmer
    
    objects_to_keep <- c("GWAS_disc", "GWAS_rep", "mytrait", "tab_res", "tab_res2", "hits", "phenos")
    objects_to_remove <- ls()
    objects_to_remove <- objects_to_remove[!objects_to_remove %in% objects_to_keep]
    rm(list = objects_to_remove) # clear objects to avoid mistakes
    
    #start-end:
    GWAS_disc$se <- paste(GWAS_disc$start_position_alignment, GWAS_disc$end_position_alignment, sep = "_")
    GWAS_rep$se <- paste(GWAS_rep$start_position_alignment, GWAS_rep$end_position_alignment, sep = "_")
    
    # merge discovery and replication
    GWAS <- merge(x = GWAS_disc, y = GWAS_rep, by = "se") # columns labelled .x are discovery, labelled .y are replication
    dim(GWAS)
    names(GWAS)
    
    ############################################################ LOOK AT THE HITS ###############################################
    
    if(mytrait %in% c("BEEHIVE_LVL", "spvl_adjusted", "spvl_normalised_adjusted")){
      
      hit1 <- which(GWAS$start_position_alignment.x == 1413 & GWAS$end_position_alignment.x == 1416)
      hit2 <- which(GWAS$start_position_alignment.x == 1514 & GWAS$end_position_alignment.x == 1514)
      hit3 <- which(GWAS$start_position_alignment.x == 6570 & GWAS$end_position_alignment.x == 6574)
      hit4 <- which(GWAS$start_position_alignment.x == 9008 & GWAS$end_position_alignment.x == 9008)
      # table to look at hits:
      hit_tab <- GWAS[c(hit1, hit2, hit3, hit4), c("se",
                                        "reference.x", "alternative1.x", "alternative2.x", "N1.x", "N2.x", "effect1.x", "effect2.x", "lik0.x", "lik1.x", "gamma.x", "p.x",
                                        "reference.y", "alternative1.y", "alternative2.y", "N1.y", "N2.y",  "effect1.y", "effect2.y", "lik0.y", "lik1.y", "gamma.y", "p.y"
      )]
      hit_tab$phenotype <- mytrait
      hits <- rbind(hits, hit_tab)
      
      GWAS$se[order(GWAS$p.x)[1:27]] # these 4 hits correspond to the first 27 lowest p-values
      
      GWAS[c(hit1, hit2, hit3, hit4), c("p.x")]
    }
    
    
    # keep only 1-mer:
    GWAS_tmp <- GWAS[which(GWAS$end_position_alignment.x==GWAS$start_position_alignment.x & GWAS$alternative1.x==GWAS$alternative1.y), ]
    GWAS_tmp2 <- GWAS[which(GWAS$end_position_alignment.x==GWAS$start_position_alignment.x & (GWAS$alternative1.x==GWAS$alternative1.y | GWAS$alternative1.x==GWAS$alternative2.y)), ]
    
    GWAS <- GWAS_tmp
    
    # SAVE IN MINIMAL DATASET
    write.csv(GWAS, paste0("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/GWAS_main_additional_SNPonly_combined_", mytrait, ".csv"))
    
    GWAS$effect <- GWAS$effect1.x # name this column "effect"
    lm1 <- lm(effect1.y ~ effect, data = GWAS) # replication vs. discovery
    print(
      summary(lm1)
    )
    ci1 <- confint(lm1)
    tab_res <- rbind(tab_res,
                     c(lm1$coefficients, ci1["(Intercept)",], ci1["effect",], summary(lm1)[[4]]["effect", "Pr(>|t|)"], summary(lm1)$r.s)
    )
    lm2 <- lm(effect1.y ~ effect, data = GWAS[GWAS$p.x<0.05,]) # replication vs. discovery
    print(
      summary(lm2)
    )
    ci2 <- confint(lm2)
    tab_res2 <- rbind(tab_res2,
                     c(lm2$coefficients, ci2["(Intercept)",], ci2["effect",], summary(lm2)[[4]]["effect", "Pr(>|t|)"], summary(lm2)$r.s)
    )
  }
  
  ################################################################################################################
  ####################################### save table of all correlations  ########################################
  ################################################################################################################
  
  tab_res <- data.frame(tab_res)
  colname <- c("int", "effect", "int_low", "int_high", "effect_low", "effect_high", "p", "R2")
  names(tab_res) <- colname
  for(mycol in colname) tab_res[mycol] <- signif(tab_res[mycol], 3)
  tab_res$pheno <- phenos
  
  tab_res2 <- data.frame(tab_res2)
  names(tab_res2) <- colname
  for(mycol in colname) tab_res2[mycol] <- signif(tab_res2[mycol], 3)
  tab_res2$pheno <- phenos
  
  write.csv(tab_res, file = "~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/discovery_replication_correlations.csv", row.names = F)
  
  ################################################################################################################
  #######################################     save table of all hits      ########################################
  ################################################################################################################
  
  # finally write down a table of the hits in discovery and replication dataset:
  hits2 <- hits[c("se", "reference.x", "alternative1.x", "alternative2.x", "effect1.x", "effect2.x", "p.x", "alternative1.y", "effect1.y", "p.y", "phenotype")]
  for(mycol in c("effect1.x", "effect2.x", "p.x", "effect1.y", "p.y")) hits2[,mycol] <- signif(hits2[,mycol], 2)
  hits2$se <- sapply(hits2$se, function(x) strsplit(x, split = "_")[[1]][1])
  names(hits2) <- c("start", "reference", "allele 1", "allele 2", "effect 1", "effect 2", "p", "rep. allele", "rep. effect", "rep. p", "phenotype")
  hits2$phenotype[hits2$phenotype == "BEEHIVE_LVL"] <- "GSVL"
  hits2$phenotype[hits2$phenotype == "spvl_adjusted"] <- "SPVL adjusted"
  hits2$phenotype[hits2$phenotype == "spvl_normalised_adjusted"] <- "SPVL adjusted normalised"
  
  write.csv(hits2, file = "~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/hits_discovery_replication.csv", row.names = F)
  
}

################################################################################################################
################################################ DRAW THE FIGURES NOW   ########################################
################################################################################################################
rm(list = ls())

phenos <-  c("BEEHIVE_LVL", "spvl_adjusted", "spvl_normalised_adjusted", "CD4_slope")

for(mytrait in phenos){
  
  GWAS <- read.csv(paste0("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/GWAS_main_additional_SNPonly_combined_", mytrait, ".csv"))
  
  View(GWAS[c("se", "effect1.x", "lik0.x", "lik1.x", "gamma.x", "p.x", "effect1.y", "lik0.y", "lik1.y", "gamma.y", "p.y")])
  addmargins(
    table(GWAS$effect1.x > 0, GWAS$effect1.y > 0)
  )
  # overall correlation (weak but existing for BEEHIVE_LVL)
  GWAS$effect <- GWAS$effect1.x # name this column "effect"
  
  lm1 <- lm(effect1.y ~ effect, data = GWAS) # replication vs. discovery
  print(
    summary(lm1)
  )
  ci1 <- confint(lm1)
  
  pdf(paste0("~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/replication_correlation_effect_sizes_", mytrait, ".pdf"), width = 4*1.5, height = 4*1.5)
  
  par(mfrow = c(1,1))
  if(mytrait == "BEEHIVE_LVL" | mytrait == "spvl_adjusted" | mytrait == "spvl_normalised_adjusted") myrange <- c(-0.4, 0.4) else if(mytrait == "CD4_slope") myrange <- c(-0.2, 0.2)
  
  plot(NULL, xlim = range(GWAS$effect), ylim = range(GWAS$effect1.y), main = mytrait, las = 1, xlab = "Effect Main", ylab = "Effect Additional")
  
  abline(h=0)
  abline(v=0)
  #abline(0,1)
  
  points(GWAS$effect, GWAS$effect1.y, pch = 20)
  
  # abline(lm1)
  plot_lm <- function(mylm , mydata, line_color = "black", b_color = "gray"){
    stopifnot("effect" %in% names(mydata))
    my_effect <- seq(quantile(mydata$effect, 0.01), quantile(mydata$effect, 0.99), length.out = 100)
    pp1 <- predict.lm(object = mylm,
                      newdata = data.frame(effect=my_effect),
                      interval = "confidence", level = 0.95)
    polygon(x = c(my_effect, rev(my_effect)), y = c(pp1[,"lwr"], rev(pp1[, 'upr'])), col = b_color, border = NA)
    points(my_effect, pp1[, "fit"], type = "l", lwd = 3, col = line_color)
  }
  
  plot_lm(lm1, GWAS) # plot general model
 
  if(mytrait %in% c("BEEHIVE_LVL", "spvl_adjusted", "spvl_normalised_adjusted")){
    hits_idx <- which(GWAS$start_position_alignment.x %in% c(1514, 9008))
    points((GWAS$effect1.x[hits_idx]), (GWAS$effect1.y[hits_idx]), pch = 20, cex = 2, col = "black")
    text((GWAS$effect1.x[hits_idx]), (GWAS$effect1.y[hits_idx]) + 0.03, c(1514, 9008), col = "black")
  }
  
  # check these are all 1-mer and position is same
  stopifnot(all(GWAS$start_position_alignment.x-GWAS$end_position_alignment.x == 0))
  stopifnot(all(GWAS$start_position_alignment.x==GWAS$start_position_alignment.y))
  oo <- order(GWAS$start_position_alignment.x)
  GWAS$negative_log_distance <- -log10(((GWAS$effect1.y-GWAS$effect)^2)) # negative log-distance
  plot(GWAS$start_position_alignment.x[oo], GWAS$negative_log_distance[oo], type = "l", xlab = "Position in alignment", ylab = "Negative log-distance")
  loess_fit <- loess(negative_log_distance ~ start_position_alignment.x, data = GWAS[oo, ], span = 0.05)
  points(GWAS$start_position_alignment.x[oo], loess_fit$fitted, type = "l", col = "red", lwd = 3)
  
  #TODO maybe could look into Syn vs. NonSyn mutations, and the correlation in effect in these two classes
  
  #abline(h=0, lty = 2); abline(v = 0, lty = 2)
  

  dev.off()

}




