rm(list = ls())

library(MetBrewer)

a <- read.table("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/9_Gabrielaite/column_ID_BETA_P_analysis_with_plink_plink_analysis_plink2.VL.glm.linear", header = T)
#b <- read.table("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/9_Gabrielaite/matched_marc_analysis_plink_bonf_p_beta.tsv", header = T)

head(a)
#head(b)
dim(a)
#dim(b)
hist(a$BETA)
plot(a$BETA, -log10(a$P), pch = 20)

a[which(a$P*nrow(a) < 0.05), ]
abline(h = -log10(0.05/nrow(a)), col = "red")

#par(mar = c(4,4,1,1))
dev.off()
original_par <- par()

tab_res <- c()


# hits in Gabrielaite
idx_hit_Gab <- which(a$ID %in% c("nef_071_K", "gag_242_N"))
idx_hit_main_all <- c()

phenos <-  c("BEEHIVE_LVL", "spvl_adjusted", "spvl_normalised_adjusted", "CD4_slope")
for(pheno_name in phenos){
  
  load(paste0("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/GWAS_results/firstGWAS_", pheno_name, "_full_17.RData"))
  
  head(GWAS_prot)
  
  range(GWAS_prot$start_position_alignment)
  range(GWAS_prot$end_position_alignment)
  table(GWAS_prot$end_position_alignment-GWAS_prot$start_position_alignment)
  
  # the 3150-bp long protein alignment comes from concatenated_PRO.fasta
  # the coordinate file is concatenated_PRO_coords.csv
  coords <- read.csv("/Users/francois.blanquart/Infectious Disease Dropbox/Francois Blanquart/BEEHIVE_Hackathon/Data/Typewriter/clean/concatenated_PRO_coords.csv")
  
  # COMPARE THE TWO:
  
  # coordinates of protein alignment in BEEHIVE:
  coords$gene <- sapply(coords$x, function(x) {
    if(grepl("tat", x) | grepl("rev", x)){ # tat and rev
      mygene <- strsplit(x, "_")[[1]][1]
      myexon <- as.numeric(gsub("exon", "", strsplit(x, "_")[[1]][2]))
      return(paste0(mygene, myexon))
    } else {
      return(strsplit(x, "_")[[1]][1])
    }
  })
  coords$pos <- as.numeric(sapply(coords$x, function(x) {
    if(grepl("tat", x) | grepl("rev", x)){ # tat and rev
      return(strsplit(x, "_")[[1]][3])
    } else {
      return(strsplit(x, "_")[[1]][2])
    }
  }))
  # the two hits of the main GWAS have p-values < 0.0001:
  print("The two hits in main GWAS:")
  print(GWAS_prot[idx_hit_main <- which(GWAS_prot$p < 0.0001),])
  if(length(idx_hit_main)<2) idx_hit_main <- c(idx_hit_main, rep(NA, 2-length(idx_hit_main)))
  idx_hit_main_all <- rbind(idx_hit_main_all, c(pheno_name, idx_hit_main))
  print(coords[match(GWAS_prot$start_position_alignment[which(GWAS_prot$p < 0.001)], coords$X), ])
  
  # coordinates of protein alignment in Gab:
  a$gene <- sapply(a$ID, function(x) strsplit(x, "_")[[1]][1])
  a$pos <- as.numeric(sapply(a$ID, function(x) strsplit(x, "_")[[1]][2]))
  a$var <- (sapply(a$ID, function(x) strsplit(x, "_")[[1]][3]))
  
  unique(a$gene)
  for(mygene in c("env", "gag", "nef", "pol", "vpu")){
    print(mygene)
    print(max(coords$pos[coords$gene==mygene]))
    print(max(a$pos[a$gene==mygene]))
    print(range(a$pos[a$gene==mygene]))
  }
  
  coords[coords$gene=="tat1",]
  a[a$gene=="tat1",]
  
  coords[coords$gene=="tat2",]
  a[a$gene=="tat2",]
  
  coords[coords$gene=="rev1",]
  a[a$gene=="rev1",]
  
  coords[coords$gene=="rev2",]
  a[a$gene=="rev2",]
  
  tab <- c()
  not_same_ref <- c()
  for(i in 1:nrow(GWAS_prot)){
    # for each amino-acid tested in main GWAS
    
    #i <- 1202
    #i <- 740
    global_pos <- GWAS_prot[i,"start_position_alignment"]
    mygene <- coords[which(coords$X==global_pos), "gene"]
    mypos <- coords[which(coords$X==global_pos), "pos"]
    
    idx0 <- which(a$gene==mygene & a$pos==mypos) # coordinate in Gabrielaite GWAS of this variant
    
    if(length(idx0) > 0){
      
      if(toupper(GWAS_prot[i, "reference"]) %in% a[idx0, "var"]){ # check the reference of main is not in alternative
        not_same_ref <- c(not_same_ref, i) # TODO improve this (check reference is same as HXB2, not necesarily alternative1)
      } else { # record only if (likely) same reference used for Gab and BEEHIVE
        for(k in 1:18){ 
          myvar <- GWAS_prot[i, paste0("alternative", k)]
          idx <- which(a$gene==mygene & a$pos==mypos & a$var == toupper(myvar))
          if(length(idx==1)){
            tab <- rbind(tab,
                         c(i, global_pos, mygene, mypos, myvar, a$BETA[idx], a$P[idx], GWAS_prot[i, paste0("effect", k)], GWAS_prot[i, paste0("N", k)],  GWAS_prot[i, "p"])
            )
          }
        }
      }
    }
  }
  tab <- data.frame(tab)
  head(tab)
  names(tab) <- c("GWAS_row", "global_pos", "gene", "pos", "var", "effect_Gab", "p_Gab", "effect", "N", "p")
  head(tab)
  for(colname in c("effect", "effect_Gab", "N", "global_pos", "pos", "p", "p_Gab")){
    tab[ , colname] <- as.numeric(as.character(tab[, colname]))
  }
  tab$gene[grepl("rev", tab$gene)] <- "rev"
  tab$gene[grepl("tat", tab$gene)] <- "tat"
  
  #gene_list <- c("pol", "env", "gag", "nef", "vpu", "rev", "tat")
  gene_list <- c("gag", "pol", "env", "tat", "vpu", "rev", "nef")
  
  colset <- met.brewer(name = "Nizami", n = 7)
  colset <- col2rgb(colset); colset <- colset/256
  colnames(colset) <- gene_list
  
  # general lm:
  assign(paste0("lm1_", pheno_name), lm(effect_Gab ~ effect, data = tab))
  
  # lm for each gene
  for(gg in gene_list){
    assign(
      paste0("lm_", pheno_name, "_", gg), lm(effect_Gab ~ effect, data = tab[tab$gene==gg,])
    )
  }
  
  # summary(lm_pol <- lm(effect_Gab ~ effect, data = tab[tab$gene=="pol",]))
  # summary(lm_env <- lm(effect_Gab ~ effect, data = tab[tab$gene=="env",]))
  # summary(lm_gag <- lm(effect_Gab ~ effect, data = tab[tab$gene=="gag",]))
  # summary(lm_nef <- lm(effect_Gab ~ effect, data = tab[tab$gene=="nef",]))
  # summary(lm_vpu <- lm(effect_Gab ~ effect, data = tab[tab$gene=="vpu",]))
  # summary(lm_rev <- lm(effect_Gab ~ effect, data = tab[grepl("rev", tab$gene),]))
  # summary(lm_tat <- lm(effect_Gab ~ effect, data = tab[grepl("tat", tab$gene),]))
  # 
  
  # META-ANALYSIS OF GABRIELAITE ET AL
  
  # Fisher's method https://en.wikipedia.org/wiki/Fisher%27s_method
  
  tab$Fisher_chi <- -2 * (log(tab$p_Gab) + log(tab$p))
  tab$Fisher_p <- 1-pchisq(tab$Fisher_chi, df = 2*2)
  tab[which(tab$Fisher_p < 0.05/nrow(tab)), ]
  tab[which(tab$Fisher_chi > 15), ]
  
  # Save table of correlations in minimal dataset
  write.csv(tab, paste0("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/GWAS_main_INSIGHT_combined_", pheno_name, ".csv"))
  
  # overall correlation (weak but existing for BEEHIVE_LVL)

  ci1 <- confint(get(paste0("lm1_", pheno_name)))
  tab_res <- rbind(tab_res,
                   c(nrow(tab), get(paste0("lm1_", pheno_name))$coefficients, ci1["(Intercept)",], ci1["effect",], summary(get(paste0("lm1_", pheno_name)))[[4]]["effect", "Pr(>|t|)"], summary(get(paste0("lm1_", pheno_name)))$r.s)
  )

}

# clean hit tables:
hit_table <- rbind(a[which(a$ID=="nef_071_K"),], a[which(a$ID=="gag_242_N"),])
hit_table <- data.frame(hit_table)
hit_table$BETA <- signif(as.numeric(hit_table$BETA), 2)
hit_table$P <- signif(as.numeric(hit_table$P), 2)
hit_table[hit_table$ID == "nef_071_K", "ID"] <- "Nef71K"
hit_table[hit_table$ID == "gag_242_N", "ID"] <- "Gag242N"
hit_table <- apply(hit_table, 2, as.character)
write.table(hit_table, file = "~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/9_Gabrielaite/hit_table_Gabrielaite.csv", row.names = F, sep = "\t")


# clean correlation table:
tab_res <- data.frame(tab_res)
colname <- c("N", "int", "effect", "int_low", "int_high", "effect_low", "effect_high", "p", "R2")
names(tab_res) <- colname
for(mycol in colname) tab_res[mycol] <- signif(tab_res[mycol], 3)
tab_res$pheno <- phenos
tab_res
#tab_res$R2 <- tab_res$R2*100

tab_res$clean <- paste(tab_res$effect, " [", tab_res$effect_low, "; ", tab_res$effect_high, "]", sep = "")
write.table(tab_res, file = "~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/9_Gabrielaite/discovery_Gabrielaite_correlations.csv", row.names = F, sep = "\t")


# write index of hits in minimal dataset
write.table(x = idx_hit_main_all, file = "~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/idx_prot_hits_main.csv")

#################################################################### FIGURE 3, right panel ####################################################################

rm(list = ls())

idx_hit_main_all <- read.table(file = "~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/idx_prot_hits_main.csv")

pdf("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/9_Gabrielaite/figure_discovery_INSIGHT.pdf", width = 4*1.5, height = 4*1.5)

phenos <-  c("BEEHIVE_LVL", "spvl_adjusted", "spvl_normalised_adjusted", "CD4_slope")
cute_names <- c("GSVL", "SPVL adjusted", "SPVL normalised adjusted", "CD4 slope"); names(cute_names) <- phenos

for(pheno_name in phenos){
  
  tab <-   read.csv(paste0("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/GWAS_main_INSIGHT_combined_", pheno_name, ".csv"))
  
  #par(original_par)
  
  plot(NULL, xlim = range(tab$effect)*1., ylim = range(tab$effect_Gab)*1., main = cute_names[pheno_name], las = 1, xlab = "Effect Main", ylab = "Effect INSIGHT START")
  
  abline(h=0)
  abline(v=0)
  #abline(0,1)
  
  points(tab$effect, tab$effect_Gab, pch = 20)
  
  # general lm:
  assign(paste0("lm1_", pheno_name), lm(effect_Gab ~ effect, data = tab))
  
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
  
  plot_lm(get(paste0("lm1_", pheno_name)), tab) # plot general model
  
  # write down actual effects
  if(pheno_name != "CD4_slope"){
    idx_hit_main <- as.numeric(idx_hit_main_all[idx_hit_main_all$V1==pheno_name,c("V2", "V3")])
    points(tab[which(tab$GWAS_row==idx_hit_main[1]), c("effect", "effect_Gab")], cex = 2, pch = 20)
    text(tab[which(tab$GWAS_row==idx_hit_main[1]), "effect"], tab[which(tab$GWAS_row==idx_hit_main[1]), "effect_Gab"]+0.005, "GagT242N")
    
    points(tab[which(tab$GWAS_row==idx_hit_main[2]), c("effect", "effect_Gab")], cex = 2, pch = 20)
    text(tab[which(tab$GWAS_row==idx_hit_main[2]), "effect"], tab[which(tab$GWAS_row==idx_hit_main[2]), "effect_Gab"]+0.005, "Nef71K")
  }
  
  #legend("bottomright", legend = gene_list, col = sapply(gene_list, function(gg) rgb(colset[1, gg], colset[2, gg], colset[3, gg])), lwd = 2, lty = 1, pch = 20, bty = "n")
  
  #TODO: add the hit
  # hits_idx <- which(GWAS$start_position_alignment.x %in% c(1514, 9008))
  # points((GWAS$effect1.x[hits_idx]), (GWAS$effect1.y[hits_idx]), pch = 20, cex = 2, col = cols[1])
  # text((GWAS$effect1.x[hits_idx]), (GWAS$effect1.y[hits_idx]) + 0.03, c(1514, 9008), col = cols[1])
  
  gene_list <- c("gag", "pol", "env", "tat", "vpu", "rev", "nef")
  
  colset <- met.brewer(name = "Nizami", n = 7)
  colset <- col2rgb(colset); colset <- colset/256
  colnames(colset) <- gene_list
  
  # plot by-gene analysis as insets
  # these three values are tuned for the pdf:
  if(plot_bygene <- FALSE){ # plot by gene or not
    xleft <- 0.15
    height <- 0.1
    ybottom <- 0.27
    for(ii in 1:length(gene_list)){
      gg <- gene_list[ii]
      par(fig = c(xleft, ybottom, xleft + height*(ii-1), ybottom + height*(ii-1)), new = T, mar = c(0,0,1,0))  
      plot(tab$effect[tab$gene==gg], tab$effect_Gab[tab$gene==gg], pch = 20, cex = 0.5, main = gg, col.main = rgb(colset[1, gg], colset[2, gg], colset[3, gg]), col = rgb(colset[1, gg], colset[2, gg], colset[3, gg], 0.5), axes=  F, xlab = "", ylab = "")
      plot_lm(mylm = get(paste0("lm_", pheno_name, "_", gg)), mydata = tab[tab$gene==gg, ],
              line_color = rgb(colset[1, gg], colset[2, gg], colset[3, gg]), b_color = rgb(colset[1, gg], colset[2, gg], colset[3, gg] ,0.5))
    }
  }
  if(pheno_name=="spvl_adjusted") dev.off() # finish the pdf here
}


