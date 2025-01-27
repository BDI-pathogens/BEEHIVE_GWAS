rm(list = ls())
library(ape)
library(phytools)

# THE BEGINNING IS JUST TO WRITE THE RELEVANT OBJECTS IN MINIMAL_DATASET FOLDER
# THIS IS NOT NEEDED ONCE THESE OBJECTS HAVE BEEN WRITTEN
if(write_minimal_dataset <- FALSE){
  
  get.bee.id <- function(x){
    x <- as.character(x)
    if(nchar(x) < 3 | is.na(x)) return(x)
    if(!substr(x,1,3) == "BEE"){
      return(x)
    }else {
      return(substr(x,1,7)) # the beehive ID
    }
  }
  wd <- "/Users/Christophe/Dropbox (Infectious Disease)/PROJECTS/HIV/BEEHIVE_Hackathon/"
  if(!file.exists(wd)) wd <- "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/"
  
  ################################## LOAD RELEVANT FILES ##################################
  
  tr <- read.tree("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Data/trees/2019-11-18_BEEHIVE_GlobalAln_clean_GSVL_GWAS_subset.fasta.treefile")
  sub_file <- paste0(wd, "Data/subtype/subtype_2019-11-18_BEEHIVE_GlobalAln_clean_final.csv") # SUBTYPE IS PRODUCED BY THE COMET SOFTWARE. R function is ~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/subtyping/subtype.R
  sub_table <- read.csv(sub_file, header = T)
  
  tmp <- ls()
  load("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/Gmatrix_data/prepared_BEEHIVE_LVL_full_17.RData")
  ls()[!ls() %in% tmp]
  
  # numeric date
  as.numeric.date <- function(date) as.numeric(as.Date(date, origin = "1970-01-01")) - as.numeric(as.Date("1985-01-01", origin = "1970-01-01")) + 1
  suba$Date.sampled.numeric <- sapply(suba$Date.sampled, as.numeric.date)
  
  # update 24/01/2025--BEE2926 was excluded from "suba", probably the minor excluded at a later stage
  tr <- drop.tip(tr, "BEE2926")
  
  ################################## ROOT THE TREE PROPERLY ##################################
  
  sub_table$PATIENT <- sapply(sub_table$name, get.bee.id)
  
  idx <- match(tr$tip.label, sub_table$PATIENT)
  stopifnot(length(idx)==length(tr$tip.label))
  stopifnot(all(idx) %in% 1:nrow(sub_table))
  
  subtypes <- sub_table$final_subtype[idx] # use final subtype (as defined in subtype.R code)
  sort(table(subtypes), decreasing = T)[1:20]
  # subtypes of tip labels
  
  ntips <- Ntip(tr)
  tr$edge
  tr$Nnode
  length(tr$tip.label)
  
  # find MRCA of set of tips
  
  findMRCA(tree = tr, tips = tr$tip.label[which(subtypes == "B")])
  findMRCA(tree = tr, tips = tr$tip.label[which(subtypes == "01_AE")])
  findMRCA(tree = tr, tips = tr$tip.label[which(subtypes == "02_AG")])
  findMRCA(tree = tr, tips = tr$tip.label[which(subtypes == "C")])
  findMRCA(tree = tr, tips = tr$tip.label[which(subtypes == "A1")])
  findMRCA(tree = tr, tips = tr$tip.label[which(subtypes == "D")])
  
  all_clade_sizes <- c()
  for(nn in nodes_to_test <- ((ntips + 1):(2 * ntips - 2))[255:748]){ # test all roots (here directly subset of all roots that are good for subtype B)
    
    tr_tmp <- root(phy = tr, node = nn) # re-root the tree as such
    mrca_B_node <- (
      findMRCA(tree = tr_tmp, tips = tr_tmp$tip.label[which(subtypes == "B")])
    )
    mrca_01_AE_node <- (
      findMRCA(tree = tr_tmp, tips = tr_tmp$tip.label[which(subtypes == "01_AE")])
    )
    mrca_02_AG_node <- (
      findMRCA(tree = tr_tmp, tips = tr_tmp$tip.label[which(subtypes == "02_AG")])
    )
    mrca_C_node <- (
      findMRCA(tree = tr_tmp, tips = tr_tmp$tip.label[which(subtypes == "C")])
    )
    mrca_A1_node <- (
      findMRCA(tree = tr_tmp, tips = tr_tmp$tip.label[which(subtypes == "A1")])
    )
    mrca_D_node <- (
      findMRCA(tree = tr_tmp, tips = tr_tmp$tip.label[which(subtypes == "D")])
    )
    mrca_G_node <- (
      findMRCA(tree = tr_tmp, tips = tr_tmp$tip.label[which(subtypes == "G")])
    )
    mrca_F1_node <- (
      findMRCA(tree = tr_tmp, tips = tr_tmp$tip.label[which(subtypes == "G")])
    )
    # number of tips+nodes in each clade
    n_subtype_tmp <- c(
      nn, 
      length(getDescendants(tree = tr_tmp, node = mrca_B_node)),
      length(getDescendants(tree = tr_tmp, node = mrca_01_AE_node)),
      length(getDescendants(tree = tr_tmp, node = mrca_02_AG_node)),
      length(getDescendants(tree = tr_tmp, node = mrca_C_node)),
      length(getDescendants(tree = tr_tmp, node = mrca_A1_node)),
      length(getDescendants(tree = tr_tmp, node = mrca_D_node)),
      length(getDescendants(tree = tr_tmp, node = mrca_G_node)),
      length(getDescendants(tree = tr_tmp, node = mrca_F1_node))
    )
    all_clade_sizes <- rbind(all_clade_sizes, n_subtype_tmp) # number of nodes+tips descending from the MRCA node of sets of tips from each subtypes
  }
  rm(tr_tmp)
  all_clade_sizes <- data.frame(all_clade_sizes)
  names(all_clade_sizes) <- c("node", "NB", "N01_AE", "N02_AG", "NC", "NA1", "ND", "NG", "NF1")
  
  # find optimal rooting to minimise these things:
  table(rs <-rowSums(all_clade_sizes[ , 2:9]))
  good_nodes <- nodes_to_test[which(rs == min(rs))]
  
  # AND REROOT:
  tr2 <- root(phy = tr, node = good_nodes[1], resolve.root = T)
  stopifnot(is.rooted(tr2))
  
  ##################################IMPROVE ROOT BASED ON MOLECULAR CLOCK ##################################
  
  # compute R2 of RTT-time plot for each possible node for rooting:
  rs <- c()
  for(myroot in good_nodes){
    tr2 <- root(phy = tr, node = myroot, resolve.root = T) # re-root the tree
    stopifnot(is.rooted(tr2))
    
    # get dates and rtt for this tree:
    idx2_suba <- match(tr2$tip.label, suba$PATIENT)
    stopifnot(all(!is.na(idx2_suba)))
    dates2 <- suba$Date.sampled.numeric[idx2_suba]
    rtt2 <- diag(vcv.phylo(tr2))
    # plot(dates2, rtt2, pch = 20, xlab = "Date (days since 1985-01-01)", ylab = "Rtt distance", las = 1, main = paste("root: ", myroot))
    # print(myroot)
    rs <- c(rs, (summary(lm(rtt2 ~ dates2))$r.s))
  }
  rm(myroot)
  # select final best root, root the tree and plot
  final_root <- good_nodes[which.max(rs)]
  tr2 <- root(phy = tr, node = final_root, resolve.root = T) # re-root the tree
  stopifnot(is.rooted(tr2))
  
  # get dates and rtt for this tree:
  idx2_suba <- match(tr2$tip.label, suba$PATIENT); stopifnot(all(!is.na(idx2_suba)))
  dates2 <- suba$Date.sampled.numeric[idx2_suba]
  rtt2 <- diag(vcv.phylo(tr2))
  plot(dates2, rtt2, pch = 20, xlab = "Date (days since 1985-01-01)", ylab = "Rtt distance", las = 1, main = paste("root: ", final_root))
  summary(lm(rtt2 ~ dates2))$r.s
  
  # save this rooted tree
  write.tree(phy = tr2, file = "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Data/trees/2019-11-18_BEEHIVE_GlobalAln_clean_GSVL_GWAS_subset.fasta_rooted.treefile")
  
  
  ################################## REMOVE OUTLIERS BASED ON ROOT TO TIP DISTANCE WITHIN CLADE B ##################################
  
  # get list of subtypes ordered as tree:
  idx2 <- match(tr2$tip.label, sub_table$PATIENT)
  subtypes2 <- sub_table$final_subtype[idx2]
  
  # step what subtypes are assigned to the descent of MRCA of each subtype:
  mrca_B_node <- ( # a mess, mixture of B and recombinant forms and F1
    findMRCA(tree = tr2, tips = tr2$tip.label[which(subtypes2 == "B")])
  )
  sort(table(subtypes2[getDescendants(tree = tr2, node = mrca_B_node)]))
  
  mrca_01_AE_node <- ( # perfect, only 01_AE
    findMRCA(tree = tr2, tips = tr2$tip.label[which(subtypes2 == "01_AE")])
  )
  table(subtypes2[getDescendants(tree = tr2, node = mrca_01_AE_node)])
  
  mrca_02_AG_node <- ( # mixture of 02_AG / G
    findMRCA(tree = tr2, tips = tr2$tip.label[which(subtypes2 == "02_AG")])
  )
  table(subtypes2[getDescendants(tree = tr2, node = mrca_02_AG_node)])
  
  mrca_C_node <- ( # good (mostly C)
    findMRCA(tree = tr2, tips = tr2$tip.label[which(subtypes2 == "C")])
  )
  table(subtypes2[getDescendants(tree = tr2, node = mrca_C_node)])
  
  mrca_A1_node <- ( # this one is a mess (MRCA is too deep, includes 02_AG)
    findMRCA(tree = tr2, tips = tr2$tip.label[which(subtypes2 == "A1")])
  )
  table(subtypes2[getDescendants(tree = tr2, node = mrca_A1_node)])
  
  # exclude B subtypes that present quite long branches
  tr_cladeB <- extract.clade(phy = tr2, node = mrca_B_node)
  hist(
    rtt_B <- diag(vcv.phylo(tr_cladeB)), breaks = 100 # root to tip distances
  )
  B_to_exclude <- names(which(rtt_B > 0.2))
  
  tr2
  tr2 <- drop.tip(phy = tr2, tip = B_to_exclude, trim.internal = T, rooted = T)
  
  # update ntips and subtypes:
  ntips <- length(tr2$tip.label)
  idx2 <- match(tr2$tip.label, sub_table$PATIENT)
  subtypes2 <- sub_table$final_subtype[idx2]
  
  # update meta-data (dates and VL)
  idx2_suba <- match(tr2$tip.label, suba$PATIENT); stopifnot(all(!is.na(idx2_suba)))
  dates2 <- suba$Date.sampled.numeric[idx2_suba]
  vl2 <- suba$BEEHIVE_LVL[idx2_suba]
  rtt2 <- diag(vcv.phylo(tr2))
  summary(lm(rtt2 ~ dates2))$r.s
  
  ######################################################## SAVING MINIMAL DATA ############################################################
  
  write.csv(x = subtypes2, file = "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/subtypes2.csv", row.names = F)
  write.csv(x = vl2, file = "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/vl2.csv", row.names = F)
  write.tree(phy = tr2, file = "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/tr2.tr")

} # End of the condition "if(write_minimal_dataset <- FALSE){" at line 6

######################################################## DRAW THE TREE ############################################################

# remove everything
rm(list = ls())

# need the following data
subtypes2 <- read.csv(file = "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/subtypes2.csv")[,1]
vl2 <- read.csv(file = "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/vl2.csv")[,1]
tr2 <- read.tree(file = "~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/minimal_dataset/tr2.tr")
ntips <- length(tr2$tip.label)

source("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/function_plot.phylo.R")

# define colors for tips:
# 1) colors by subtypes:

colset <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
colset <-c(colset[c(2, 6, 4, 8, 10, 1)], "gray")
names(colset) <- c("B", "01_AE", "02_AG", "C", "A1", "D", "Others/unassigned")
cols2 <- sapply(subtypes2, function(ss){ # colors for each tip depending on subtype
  if(ss %in% names(colset)[names(colset) != "Others"]) return(colset[ss])
  return(colset["Others/unassigned"])
})

# 2) color by VL
vl2_cat <- cut(vl2, breaks = quantile(vl2, seq(0, 1, 0.1)))
vl2_cat[which.min(vl2)] <- "(2.1,3.74]"
stopifnot(all(!is.na(vl2_cat)))
colset_vl <- viridisLite::viridis(length(unique(vl2_cat))); names(colset_vl) <- names(table(vl2_cat))
cols2_vl <- colset_vl[vl2_cat]

# r_theta_tocoord <- function(r, theta){ # transform a radius and angle to x, y coords (useful for plotting)
#   stopifnot(theta >= -pi & theta <= pi)
#   stopifnot(r > 0)
#   return(list(x = r * cos(theta), y = r * sin(theta)))
# }


tip.labels2 <- tr2$tip.label
tr2$tip.label<- rep("", ntips)

pdf("~/Dropbox (Infectious Disease)//BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/tree_VL_Figure1.pdf", width = 10, height = 10)

par(xpd = TRUE, mar = c(1,1,2,1))
tree_plot <- plot.phylo.FB(tr2, type = "fan", show.node.label = F, show.tip.label = T, x.lim = c(-0.4, 0.5), y.lim = c(-0.41, 0.5))  # need to show tip label to get angles
#tr2$tip.label <- tip.labels2

# extract tip coords and draw colors for subtypes:
tip_xx <- tree_plot$xx.tips # x of all tips
tip_yy <- tree_plot$yy.tips # y of all tips
tip_r2 <- tip_xx^2 + tip_yy^2 # r2 of all tips

tip_angles <- atan2(tip_yy, tip_xx) # angles of all tips between -pi and pi

radius2 <- max(tip_r2)
#factor2 <- radius2 / tip_r2 # factor to get coordinate on a circle
circle_xx <- tip_xx * sqrt(radius2) / sqrt(tip_r2) # * sqrt(factor2)
circle_yy <- tip_yy * sqrt(radius2) / sqrt(tip_r2) # * sqrt(factor2)

#points(x = circle_xx, y = circle_yy, col = cols2, pch = 20, cex = 0.3) # 
#points(x = circle_xx * 1.03, y = circle_yy * 1.03, col = cols2_vl, pch = 20, cex = 0.3)

width <- 0.08
segments(x0 = circle_xx * 1.0, y0 = circle_yy * 1.0, x1 = circle_xx * (1 + width), y1 = circle_yy * (1 + width), col = cols2, lwd = 1)
segments(x0 = circle_xx * (1+width+0.01), y0 = circle_yy * (1+width+0.01), x1 = circle_xx * (1+2*width+0.01), y1 = circle_yy * (1+2*width+0.01), col = cols2_vl, lwd = 1)

legend(x = -0.44, 0.54, pch = 20, col = colset, legend = names(colset), bty = "n")
legend(x = 0.34, 0.54, pch = 20, col = colset_vl, legend = names(colset_vl), bty = "n")
text(x = -0.44, y = 0.55, "Subtype", adj = 0, cex = 1.5)
text(x = 0.2, y = 0.55, "GSVL decile (log10 copies/mL)", adj = 0, cex = 1.5)

dev.off()



