rm(list = ls())
library(plyr)
library(RColorBrewer)

if(write_minimal_dataset <- TRUE){
  
  load("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/Gmatrix_data/prepared_BEEHIVE_LVL_full_17.RData")
  load("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/GWAS_results/firstGWAS_BEEHIVE_LVL_full_17.RData")
  
  ls()
  
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
    
    idx_in_ind <- idx_in_GWAS # idx in ind is the same as index in GWAS
    
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
  } # need to re-enter this (because saved version was uncorrected for some small mistake)
  
  ################################################################################################################################################ 
  ################################################    1) TEMPORAL TRAJECTORIES OF THE HITS        ################################################
  ################################################################################################################################################ 
  
  # position 1:
  st1 <- 1514
  en1 <- 1514
  idx1 <- get_index_in_dm(st = st1, en = en1, GWAS_tab = GWAS_kmer, full_design_matrix = full_design_matrix_kmer, full_ind_tokeep = full_ind_tokeep_kmer)
  dm_pos1 <- full_design_matrix_kmer[, idx1["idx_in_dm_start"]:idx1["idx_in_dm_end"]]
  colSums(dm_pos1, na.rm = T)
  
  # position 2:
  st2 <- 9008
  en2 <- 9008
  idx2 <- get_index_in_dm(st = st2, en = en2, GWAS_tab = GWAS_kmer, full_design_matrix = full_design_matrix_kmer, full_ind_tokeep = full_ind_tokeep_kmer)
  dm_pos2 <- full_design_matrix_kmer[, idx2["idx_in_dm_start"]:idx2["idx_in_dm_end"]]
  colSums(dm_pos2, na.rm = T)
  
  dim(dm_pos1)
  dim(dm_pos2)
  dim(suba)
  names(suba)
  
  as.numeric.date <- function(date) as.numeric(as.Date(date, tryFormats = c("%d/%m/%Y"), origin = "01/01/1970")) - as.numeric(as.Date("1985-01-01", origin = "1970-01-01")) + 1
  
  process_tab <- function(dm_pos){
    
    tab <- table(suba$Year.pos, dm_pos[, 1])
    tab <- as.matrix(tab, ncol = 2)
    rn <- rownames(tab)
    tab <- unname(tab)
    tab <- as.data.frame.matrix(tab)
    colnames(tab) <- c("WT", "MUT")
    
    tab$year <- as.numeric(as.character(rn))
    
    tab$WT <- as.numeric(tab$WT)
    tab$MUT <- as.numeric(tab$MUT)
    tab$N <- tab$WT + tab$MUT
    tab$year[tab$year %in% 1985:1989] <- 1990 # put 1985-1989 into 1990 for the sake of having enough data points
    
    tab$year_cat <- floor(tab$year / 5) * 5
    tab <- ddply(tab, .(year_cat), summarise, MUT = sum(MUT), WT = sum(WT), N = sum(N))
    
    # remove beginning and end
    # tab0 <- tab[tab$year %in% 1985:1995, ]
    # tab0_2 <- tab[tab$year %in% 2012:2014, ]
    # tab <- tab[!tab$year %in% 1985:1995 & !tab$year %in% 2012:2014, ] 
    # 
    # 
    # truc0 <- matrix(c(1990, sum(tab0$WT), sum(tab0$MUT), sum(tab0$N)), nrow = 1)
    # colnames(truc0) <- c("year", "WT", "MUT", "N")
    # rm(tab0)
    # 
    # truc0_2 <- matrix(c(2013, sum(tab0_2$WT), sum(tab0_2$MUT), sum(tab0_2$N)), nrow = 1)
    # colnames(truc0_2) <- c("year", "WT", "MUT", "N")
    # rm(tab0_2)
    # 
    # 
    # tab <- rbind(truc0, tab, truc0_2)
    
    # compute frequencies
    tab$f <- tab$MUT / tab$N
    
    tab$f_low <- tab$f  - 1.96 * sqrt(tab$f * (1 - tab$f) / tab$N)
    tab$f_up <- tab$f  + 1.96 * sqrt(tab$f * (1 - tab$f) / tab$N)
    
    tab
  }

  frequency_through_time_hit1 <- process_tab(dm_pos1)
  frequency_through_time_hit2 <- process_tab(dm_pos2)
  
  write.csv(x = frequency_through_time_hit1, file = "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/frequency_through_time_hit1.csv", row.names = F)
  write.csv(x = frequency_through_time_hit2, file = "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/frequency_through_time_hit2.csv", row.names = F)
  
  # NO linkage disequilibrium between the two positions
  dim(dm_pos1)
  dim(dm_pos2)
  colSums(dm_pos1, na.rm = T)
  colSums(dm_pos2, na.rm = T)
  
  table(dm_pos1[,1], dm_pos2[,1])
  
  
  dm_pos1[
    which(suba$Year.pos==1996), 1
  ]
  dm_pos2[
    which(suba$Year.pos==1996), 1
  ]
  
  ################################################################################################################################################
  ############################################################### ENRICHMENT ANALYSES ############################################################
  ################################################################################################################################################
  
  ############################################### 1. COMPLETE ANNOTATION OF HXB2 REF ##########################################################
  
  # We look at whether positions have more signifant effects when they are CTL escape mutations
  
  # import reference for enrichment analysis
  hxb2 <- read.csv(
    file = "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/ReferenceData/hxb2_annotated_Chris.csv"
  )
  
  all(hxb2$HXB2.base.position==1:9719) # position
  hxb2[which(hxb2$HXB2.base.position==1514), "CTL.escape."]
  hxb2[which(hxb2$HXB2.base.position==1514), "Other.features"]
  hxb2_CTL <- hxb2[grepl(pattern = "CTL", x = hxb2$CTL.escape.), ]
  
  table(hxb2$RF1.protein)
  table(hxb2$RF2.protein)
  table(hxb2$RF3.protein)
  
  # try to sort out CTL escape positions a bit more
  table(table(hxb2_CTL$CTL.escape.))
  
  # add protein_position columns
  hxb2$prot_pos_1 <- paste(hxb2$RF1.protein, hxb2$RF1.aa.position, sep = "_")
  hxb2$prot_pos_2 <- paste(hxb2$RF2.protein, hxb2$RF2.aa.position, sep = "_")
  hxb2$prot_pos_3 <- paste(hxb2$RF3.protein, hxb2$RF3.aa.position, sep = "_")
  
  which(table(hxb2$prot_pos_1) == 1)
  which(table(hxb2$prot_pos_2) == 2)
  which(table(hxb2$prot_pos_3) == 1) # four nucleotides correspond to vpr_72 (slippage at multiple Ts?)
  
  # add nucleotide position for each line
  # gag = RF1
  # vif = RF1
  # gp120 = RF3
  # gp41 = RF3
  # pol = RF3
  # nef = RF1 (la partie sur RF3 est redondante)
  # tat = RF2 puis RF1
  # vif = RF1 
  # vpr = RF3 puis RF1
  # rev = RF3 puis RF2
  # vpu = RF2
  
  # view genes that have complicated structures:
  View(hxb2[which(hxb2$RF1.protein == "nef" | hxb2$RF3.protein == "nef"), ])
  View(hxb2[which(hxb2$RF2.protein == "tat" | hxb2$RF1.protein == "tat"), ])
  View(hxb2[which(hxb2$RF1.protein == "vpr" | hxb2$RF3.protein == "vpr"), ])
  View(hxb2[which(hxb2$RF2.protein == "rev" | hxb2$RF3.protein == "rev"), ])
  
  # https://www.hiv.lanl.gov/content/sequence/HIV/MAP/landmark.html
  # position of each nucleotide in codon, for each reading frame
  hxb2$pos1 <- NA
  hxb2$pos2 <- NA
  hxb2$pos3 <- NA
  
  for(mygene in c("gag", "vif", "gp120", "gp41", "pol", "nef", "tat", "vpr", "rev", "vpu")){
    
    aa_nb <- 0
    flag <- T
    while(flag){
      
      aa_nb <- aa_nb + 1
      idx1 <- which(hxb2$RF1.protein==mygene & hxb2$RF1.aa.position == aa_nb)
      idx2 <- which(hxb2$RF2.protein==mygene & hxb2$RF2.aa.position == aa_nb)
      idx3 <- which(hxb2$RF3.protein==mygene & hxb2$RF3.aa.position == aa_nb)
      all_idx <- c(idx1, idx2, idx3)
      n_nuc <- length(idx1)+length(idx2)+length(idx3)
      
      if(n_nuc == 0) flag <- F
      
      if(n_nuc == 3){
        if(length(idx1)==3) {if(any(!is.na(hxb2$pos1[idx1]))) stop("pos 1 already written"); hxb2$pos1[idx1] <- 1:3}
        if(length(idx2)==3) {if(any(!is.na(hxb2$pos2[idx2]))) stop("pos 2 already written"); hxb2$pos2[idx2] <- 1:3}
        if(length(idx3)==3) {if(any(!is.na(hxb2$pos3[idx3]))) stop("pos 3 already written"); hxb2$pos3[idx3] <- 1:3}
        if(length(idx1) != 3 & length(idx2) != 3 & length(idx3) != 3){
          #cat("overlapping two RF", mygene, aa_nb, "\n")
          if(mygene == "tat" & aa_nb == 72){
            stopifnot(length(idx2) == 2 & length(idx1) == 1)
            hxb2$pos2[idx2] <- 1:2
            hxb2$pos1[idx1] <- 3
          }
          if(mygene == "rev" & aa_nb == 26){
            stopifnot(length(idx2) == 2 & length(idx3) == 1)
            hxb2$pos3[idx3] <- 1
            hxb2$pos2[idx2] <- 2:3
          }
        }
      }
      if(n_nuc != 3 & n_nuc != 0){
        if(n_nuc == 6 & mygene == "nef"){
          hxb2$pos1[idx1] <- 1:3 # for nef just consider RF 1
        } else {
          if(n_nuc == 4 & mygene == "vpr" & aa_nb == 72){ 
            # vpr gene has four nucleotides, 1 in RF3 then 3 in RF1
            hxb2$pos3[idx3] <- 1; stopifnot(length(idx3)==1)
            hxb2$pos1[idx1] <- 2:4; stopifnot(length(idx1)==3)
          } else {
            cat("number of nucleotides different from 3:", mygene, aa_nb, "\n")
          }
        }
      }
    }
  }
  
  hxb2$prot_pos_CTL <- gsub(pattern = "CTL escape ", replacement = "", x = hxb2$CTL.escape.)
  hxb2$RF_CTL <- NA # which is the relevant reading frame for CTL escape
  hxb2$RF_CTL[hxb2$prot_pos_CTL == hxb2$prot_pos_1] <- 1
  hxb2$RF_CTL[hxb2$prot_pos_CTL == hxb2$prot_pos_2] <- 2
  hxb2$RF_CTL[hxb2$prot_pos_CTL == hxb2$prot_pos_3] <- 3
  hxb2$pos_CTL <- NA # nucleotidic position (1, 2 or 3) for each CTL escape mutation
  hxb2$pos_CTL[which(hxb2$RF_CTL == 1)] <- hxb2$pos1[which(hxb2$RF_CTL == 1)]
  hxb2$pos_CTL[which(hxb2$RF_CTL == 2)] <- hxb2$pos2[which(hxb2$RF_CTL == 2)]
  hxb2$pos_CTL[which(hxb2$RF_CTL == 3)] <- hxb2$pos3[which(hxb2$RF_CTL == 3)]
  
  hxb2[which(is.na(hxb2$RF_CTL) & grepl(pattern = "CTL escape ", x = hxb2$CTL.escape.)), ] # those for which we don't know in which gene is the epitope
  
  write.csv(x = hxb2, file = "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/ReferenceData/hxb2_annotated_Chris_FB.csv")
  # add the multiplicity of the CTL escape noted (investigate further the 1, 2, 4, 6)
  
  
  ############################################### 2. LOOK AT ENRICHMENT ##########################################################
  
  # focuses on SNPs
  GWAS1 <- GWAS_kmer[GWAS_kmer$start_position_alignment==GWAS_kmer$end_position_alignment, ] 
  GWAS1$annotation <- hxb2$CTL.escape.[GWAS1$start_position_alignment]
  GWAS1$pos_CTL <- hxb2$pos_CTL[GWAS1$start_position_alignment]
  GWAS1$drm_ias <- (hxb2$DRM.IAS.=="DRM_IAS")[GWAS1$start_position_alignment]
  GWAS1$drm_la <- (hxb2$DRM.Los.Alamos.=="DRM_LosAlamos")[GWAS1$start_position_alignment]
  
  GWAS1$is_CTL <- grepl(pattern = "CTL escape", x = GWAS1$annotation)
  table(GWAS1$is_CTL, GWAS1$pos_CTL)
  
  hist(log10(GWAS1$p[GWAS1$is_CTL]), breaks = 50, freq = F, las = 1, col = rgb(1,0,0,0.5))
  hist(log10(GWAS1$p[!GWAS1$is_CTL]), breaks = 50, freq = F, las = 1, col = rgb(0,0,1,0.5), add = T)
  
  enrichment_df <- c() # dataframe of enrichment results
  
  # look at association for different p-value thresholds
  
  for(pthreshold in ptlist <- c(0.0001, 0.01, 0.05)){
    cat("pvalue threshold = ", pthreshold, "\n")
    tmp_tab <- table(low_p = GWAS1$p < pthreshold, CTL = GWAS1$is_CTL)
    print(tmp_tab)
    enrichment_df <- rbind(enrichment_df, 
                           c(pthreshold, "CTL", - log10(fisher.test(tmp_tab)$p))
    )
  }
  
  # for CTL-decomposed by 1st, 2nd, 3rd nucleotide
  
  for(codon_position in 1:3){
    cat("\nposition in codon: ", codon_position, "\n")
    for(pthreshold in ptlist){
      cat("pvalue threshold = ", pthreshold, "\n")
      tmp_tab1 <- table(low_p = GWAS1$p < pthreshold, CTL1 = GWAS1$pos_CTL==codon_position)
      print(tmp_tab1)
      enrichment_df <- rbind(enrichment_df, 
                             c(pthreshold, paste0("CTL_", codon_position), - log10(fisher.test(tmp_tab1)$p))
      )
    }
  }
  
  # for DRM - IAS
  for(pthreshold in ptlist <- c(0.0001, 0.01, 0.05)){
    cat("pvalue threshold = ", pthreshold, "\n")
    tmp_tab <- table(low_p = GWAS1$p < pthreshold, drm = GWAS1$drm_ias)
    print(tmp_tab)
    enrichment_df <- rbind(enrichment_df, 
                           c(pthreshold, "DRM_IAS", - log10(fisher.test(tmp_tab)$p))
    )
  }
  
  # for DRM - LA
  for(pthreshold in ptlist <- c(0.0001, 0.01, 0.05)){
    cat("pvalue threshold = ", pthreshold, "\n")
    tmp_tab <- table(low_p = GWAS1$p < pthreshold, drm = GWAS1$drm_la)
    print(tmp_tab)
    enrichment_df <- rbind(enrichment_df, 
                           c(pthreshold, "DRM_LA", - log10(fisher.test(tmp_tab)$p))
    )
  }
  
  enrichment_df <- data.frame(enrichment_df)
  names(enrichment_df) <- c("p_threshold", "type", "logp")
  enrichment_df$logp <- as.numeric(enrichment_df$logp)
  
  write.csv(x = enrichment_df, file = "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/enrichment_df.csv", row.names = F)
  
  ################################################################################################################################################
  ####################################### 3. LOOK AT POSITIONS NOTED IN BARTHA ET AL ELIFE PAPER #################################################
  ################################################################################################################################################
  
  # we look at the positions noted in figure 3 (significantly associated with human HLA variation) to see what effect they have on VL
  
  # test a few of the Bartha et al. positions:
  GWAS1[which(GWAS1$start_position_alignment %in% 
                hxb2[which(hxb2$RF1.protein=="gag" & hxb2$RF1.aa.position==147), "HXB2.base.position"]
  ), ]
  
  #RT to be found
  
  GWAS1[which(GWAS1$start_position_alignment %in% 
                hxb2[which(hxb2$RF1.protein=="vif" & hxb2$RF1.aa.position==51), "HXB2.base.position"]
  ), ]
  
  GWAS1[which(GWAS1$start_position_alignment %in% 
                hxb2[which(hxb2$RF1.protein=="nef" & hxb2$RF1.aa.position==126), "HXB2.base.position"]
  ), ]
  
  GWAS1[which(GWAS1$start_position_alignment %in% 
                hxb2[which(hxb2$RF1.protein=="nef" & hxb2$RF1.aa.position==135), "HXB2.base.position"]
  ), ]
  
}

rm(list = ls())

################################################################################################################################################
############################################################## GENERATE FIGURES     ############################################################
################################################################################################################################################

# retrieve tables:
tab1 <- read.csv(file = "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/frequency_through_time_hit1.csv")
tab2 <- read.csv(file = "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/frequency_through_time_hit2.csv")
enrichment_df <- read.csv(file = "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/enrichment_df.csv")


pdf("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/temporal_trajectories_hits.pdf", width = 4, height = 4)
par(mar = c(6,4,1,1), xpd = T)
xoffset <- 0.3
# temporal trajectories of the two alleles:
plot(tab1$year, tab1$f, type = "o", pch = 20, las = 1, xlim = c(1989, 2011), ylim = c(0,0.3), axes = F, ylab = "Frequency", xlab = "", lwd = 2)
points(tab2$year+xoffset, tab2$f, type = "o", pch = 20, col= "red", lwd = 2)
segments(x0 = tab1$year, x1 = tab1$year, y0 = tab1$f_low, y1 = tab1$f_up, col = "black", lwd = 2)
segments(x0 = tab2$year+xoffset, x1 = tab2$year+xoffset, y0 = tab2$f_low, y1 = tab2$f_up, col = "red", lwd = 2)

axis(side = 1, at = seq(1990, 2010, 5), label = rep("", 5), cex = 1)
text(x  = seq(1990, 2010, 5), y = -0.1 , labels = c("1985-1994", "1995-1999", "2000-2004", "2005-2009", "2010-2014"), srt = 45)
axis(side = 2, at = seq(0, 0.3, 0.1), las = 1)
legend("topright", legend = c("C1514A", "G9008A"), col = c("black", "red"), lty = 1, pch = 20, bty = "n")
dev.off()


pdf("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/enrichment_analysis.pdf", width = 4 * 1.2, height = 4 * 1.2)

cols <- c(rep(brewer.pal(n = 5, name = "Set3")[4], 12), rep(brewer.pal(n = 5, name = "Set3")[5], 12))
par(mar = c(4.5,4,3,1))
ymax <- 1.6
plot(1:18, y = enrichment_df$logp, las = 1, col = cols, pch = 20, xlim = c(1,18), ylim = c(0, ymax), xlab = "", ylab = "p-value", axes = F) #ylab =  expression("-"~log[10]~"(p-value)"), axes = F)
segments(x0 = 0, x1 = 22, y0 = -log10(0.05), y1 = -log10(0.05), lty = 2)
segments(x0 = 1:18, x1 = 1:18, y0 = rep(0,18), y1 = enrichment_df$logp, col = cols, lwd = 3)
axis(side = 2, at = -log10(c(0.025, 0.05, 0.1, 0.5, 1)), labels = c(0.025, 0.05, 0.1, 0.5, 1), las = 1)

par(xpd = TRUE) # Draw outside plot area is possible
segments(x0 = 3.5, y0 = 0, x1 = 3.5, y1 = ymax, col = "gray")
segments(x0 = 6.5, y0 = 0, x1 = 6.5, y1 = ymax, col = "gray")
segments(x0 = 9.5, y0 = 0, x1 = 9.5, y1 = ymax, col = "gray")
segments(x0 = 12.5, y0 = 0, x1 = 12.5, y1 = ymax, col = "gray")
segments(x0 = 15.5, y0 = 0, x1 = 15.5, y1 = ymax, col = "gray")

enrichment_df$pth <- rep(c("0.0001", "0.01", "0.05"), 6)
text(x = -0.3 + (1:18), y = -0.15, labels = enrichment_df$pth, srt = 45)
text(x = c(2, 5, 8, 11, 14, 17), y = -0.32, c("all", "pos 1", "pos 2", "pos 3", "IAS", "LANL"))

enrichment_df$label <- c("CTL - 0.0001", "CTL - 0.01", "CTL - 0.05",
                         "CTL 1st - 0.0001", "CTL 1st - 0.01", "CTL 1st - 0.05",
                         "CTL 2nd - 0.0001", "CTL 2nd - 0.01", "CTL 2nd - 0.05",
                         "CTL 3rd - 0.0001", "CTL 3rd - 0.01", "CTL 3rd - 0.05",
                         "DRM - IAS - 0.0001", "DRM - IAS - 0.01", "DRM - IAS - 0.05",
                         "DRM - LANL - 0.0001", "DRM - LANL - 0.01", "DRM - LANL - 0.05"
                         )
rect(xleft = 1, xright = 12.3, ybottom = ymax, ytop = ymax + 0.02, border = NA, col = cols[1])
text(x = 6.75, y = ymax + 0.1, labels = "CTL escape", col = cols[1], cex = 2)

rect(xleft = 12.7, xright = 18, ybottom = ymax, ytop = ymax + 0.02, border = NA, col = cols[13])
text(x = 15.5, y = ymax + 0.1, labels = "DRM", col = cols[13], cex = 2)

dev.off()






