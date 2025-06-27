rm(list = ls())
library(plyr)
library(ggplot2)

draw_qqplot <- function(tab, limit_p, label, save_plot = F){
  if(!"p" %in% colnames(tab)) stop('tab has to include a column named p')
  punif <- (1:nrow(tab)) / nrow(tab)
  tab[tab$p ==0, "p"] <- 10^-16 # some are 0s
  max_p <- 1.03 * max(-log10(punif), -log10(tab$p), na.rm = T)
  
  if(save_plot) pdf(paste0(data_folder, "figures/", label, ".pdf"), width = 4, height = 4)
  
  par(mar = c(4,4,1,1), mfrow = c(1,1))
  if(save_plot) main <- "" else main <- label
  qqplot(-log10(punif), -log10(tab$p), type = "p", pch = 20, las = 1, cex = 0.5, xaxs = "i", yaxs = "i", xlim = c(0,max_p), ylim = c(0,max_p), xlab = expression(log[10]~"p-value"~uniform), ylab = expression(log[10]~"p-value"), bty = "n", axes = F, main = main)
  axis(1, at = seq(0, floor(max_p), 1))
  axis(2, at = seq(0, floor(max_p), 1), las = 1)
  segments(x0 = 0, y0 = 0, x1 = max_p, y1 = max_p)
  abline(h = limit_p, lty = 2)
  
  if(save_plot) dev.off()
}
55
nind <- 2250; nrep <- 322 # number of individuals in discovery and replication GWAS
nlocus_total <- 5530 # 5000 loci tested
n_effective_tests <- 1291
sde <- 0.7 # this is tuned to reprodcue the right heritabinlity of hits; 0.4 = 25% h2; 0.7 =  ; 0.8 = 8.5%; 0.85 = 8%; 1 = 6% h2

# we want to reproduce the effect of the ~30 "top hits" that are beyond? or the larger 
# - heritability explained by hits (around 10%, + 20% being explained by small effects controlled by the random effects for a total of 30%) 
# - shape of qqplot ? -> difficult because somehow the qqplot controls for many effects through the regression
# - overall correlation between effect size discovery and replication (regression coefficient, p-value)


res <- c()
# all_inferred_p_discovery <- c()
for(rep in 1:50){
  
  print(rep)
  for(nlocus_effect in locuslist <- c(10, 20, 50, 100, 200, 500)){
    
    print(nlocus_effect)
    effect_size <-  1  / sqrt(nlocus_effect) # 0.1 this is to ensure constant across n loci, as V[g] = n effectsize^2 p(1-p); coefficient inferred by trial and error for around 10% heritability caused by these locis 
    
    # effects of alleles:
    effects <- as.matrix(c(rnorm(n = nlocus_effect, mean = 0, sd = effect_size), rep(0, nlocus_total - nlocus_effect)), ncol = 1)
    
    # population frequency of each allele
    # 0.5160741 3.1197464  estimate from GWAS_tab (done in inspect_GWAS_result.R)
    f <- rbeta(n = nlocus_total, shape1 = 0.5160741, shape2 = 3.1197464)/2 + (0.004444444 - 1e-3) # mean frequency is around 8%
    
    flag <- T
    i  <- 1
    while(flag){
      #print(i)
      design_mat_discovery <- sapply(1:nlocus_total, function(i) rbinom(n = nind, size = 1, prob = f[i]))
      design_mat_replication <- sapply(1:nlocus_total, function(i) rbinom(n = nrep, size = 1, prob = f[i]))
      if(all(colSums(design_mat_discovery) > 0)) flag <- FALSE
      i<-i+1
    }
    
    g <- design_mat_discovery %*% effects # breeding values
    e <- rnorm(nind, mean = 0, sd = sde)
    z <- g + e
    h2 <- var(g) / var(z)
    h2
    
    # discovery GWAS:
    inferred_effect_discovery <- c()
    inferred_p_discovery <- c()
    for(i in 1:nlocus_total){
      lm_tmp <- lm(z ~ design_mat_discovery[, i])
      inferred_effect_discovery <- c(inferred_effect_discovery,
                                     lm_tmp$coefficients[2]
      )
      inferred_p_discovery <- c(inferred_p_discovery,
                                     summary(lm_tmp)[[4]]["design_mat_discovery[, i]", "Pr(>|t|)"]
      )
    }
    # draq qqplot
    if(rep == 1) draw_qqplot(tab = data.frame(p = inferred_p_discovery), limit_p = 0.05, label = paste("N loci = ", nlocus_effect))
    #all_inferred_p_discovery <- rbind(all_inferred_p_discovery, inferred_p_discovery)
    
    # replication GWAS:
    g_rep <- design_mat_replication %*% effects
    e_rep <- rnorm(nrep, mean = 0, sd = sde)
    z_rep <- g_rep + e_rep
    h2_rep  <- var(g_rep) / var(z_rep)
    
    # predicted g in rep GWAS (STRICT = 0.05)
    signif_pos <- which(inferred_p_discovery < 0.05/n_effective_tests) # bonferonni-selected significance
    if(length(signif_pos)>0){
      pred_g_rep <- design_mat_replication[ , signif_pos, drop = F] %*% as.matrix(inferred_effect_discovery, ncol = 1)[signif_pos, , drop = F]
      lmz <- lm(z_rep ~ pred_g_rep)
      slmz <- summary(lmz)
      pred_g_corr <- c(slmz$r.s, slmz[[4]]["pred_g_rep", "Estimate"],  slmz[[4]]["pred_g_rep", "Pr(>|t|)"])
    } else {
      pred_g_corr <- c(NA, NA, NA)
    }
    
    # predicted g in rep GWAS with ALL effects
    pred_g_rep2 <- design_mat_replication[ , , drop = F] %*% as.matrix(inferred_effect_discovery, ncol = 1)[, , drop = F]
    lmz2 <- lm(z_rep ~ pred_g_rep2)
    slmz2 <- summary(lmz2)
    pred_g_corr2 <- c(slmz2$r.s, slmz2[[4]]["pred_g_rep2", "Estimate"],  slmz2[[4]]["pred_g_rep2", "Pr(>|t|)"])
    
    
    inferred_effect_replication <- c()
    for(i in 1:nlocus_total){
      inferred_effect_replication <- c(inferred_effect_replication,
                                       lm(z_rep ~ design_mat_replication[, i])$coefficients[2]
      )
    }
    
    if(length(signif_pos)>0){
      # fraction significant effects going in same direction
      frac_same_direction <- 1-sum(inferred_effect_replication[signif_pos]*inferred_effect_discovery[signif_pos]<0)/length(signif_pos)
    } else {
      frac_same_direction <- NA
    }

    # correlation effect sizes discovery and replication
    # plot(inferred_effect_discovery, inferred_effect_replication)
    lm0 <- lm(inferred_effect_replication ~ inferred_effect_discovery)
    slm0 <- summary(lm0)
    
    res <- rbind(res, 
                 c(
                   nlocus_effect, length(signif_pos), frac_same_direction, h2, h2_rep,
                   slm0$r.s, slm0[[4]]["inferred_effect_discovery", "Estimate"],  slm0[[4]]["inferred_effect_discovery", "Pr(>|t|)"],
                   pred_g_corr, # STRICT
                   pred_g_corr2 # ALL EFFECTS
                   )
    )
  }
  print(res)
}
res <- data.frame(res)
names(res) <- c(
                "nlocus_effect", "n_hits", "frac_same_direction", "h2", "h2_rep",
                "R2_eff", "beta_eff", "p_eff",
                "R2_pred", "beta_pred", "p_pred",
                "R2_pred2", "beta_pred2", "p_pred2"
                )


save.image("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/tmp_GWAS_toy_simulations.RData")


# can re-start from here:
load("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/tmp_GWAS_toy_simulations.RData")

library(plyr)

res2 <- ddply(.data = res, .variables = .(nlocus_effect), summarise,
             n_hits = mean(n_hits),
             h2 = mean(h2), h2_rep = mean(h2_rep),
             
             R2_eff = mean(R2_eff), beta_eff = mean(beta_eff), p_eff = mean(p_eff),
             R2_eff_low = quantile(R2_eff, 0.025, na.rm = T), R2_eff_high = quantile(R2_eff, 0.975, na.rm = T),
             beta_eff_low = quantile(beta_eff, 0.025, na.rm = T), beta_eff_high = quantile(beta_eff, 0.975, na.rm = T),
             p_eff_low = quantile(p_eff, 0.025, na.rm = T), p_eff_high = quantile(p_eff, 0.975, na.rm = T),
             
             R2_pred = mean(R2_pred, na.rm = T), beta_pred = mean(beta_pred, na.rm = T), p_pred = mean(p_pred, na.rm = T),
             R2_pred_low = quantile(R2_pred, 0.025, na.rm = T), R2_pred_high = quantile(R2_pred, 0.975, na.rm = T),
             beta_pred_low = quantile(beta_pred, 0.025, na.rm = T), beta_pred_high = quantile(beta_pred, 0.975, na.rm = T),
             p_pred_low = quantile(p_pred, 0.025, na.rm = T), p_pred_high = quantile(p_pred, 0.975, na.rm = T),
             
             R2_pred2 = mean(R2_pred2, na.rm = T), beta_pred2 = mean(beta_pred2, na.rm = T), p_pred2 = mean(p_pred2, na.rm = T),
             R2_pred2_low = quantile(R2_pred2, 0.025, na.rm = T), R2_pred2_high = quantile(R2_pred2, 0.975, na.rm = T),
             beta_pred2_low = quantile(beta_pred2, 0.025, na.rm = T), beta_pred2_high = quantile(beta_pred2, 0.975, na.rm = T),
             p_pred2_low = quantile(p_pred2, 0.025, na.rm = T), p_pred2_high = quantile(p_pred2, 0.975, na.rm = T)
          )

plot(log10(res2$nlocus_effect), res2$h2, type = "l", ylim = c(0,1), ylab = "Heritability", xlab = "Number locus")
points(log10(res$nlocus_effect), res$h2, type = "p", pch = 20)
summary(res2$h2)

plot(log10(res2$nlocus_effect), res2$R2_eff, type = "o", pch = 20, las = 1, ylim = c(0, 1), ylab =  expression(paste(R^2, " correlation")), xlab = "Number of loci")
points(log10(res$nlocus_effect), res$R2_eff, type = "p", pch = 20)
abline(h=0)

# correlation coefficient between effect sizes
plot(log10(res2$nlocus_effect), res2$beta_eff, type = "o", pch = 20, las = 1, ylim = c(-0.1, 0.2), ylab =  c("Regression coefficient effects"), xlab = "Number of loci")
points(log10(res$nlocus_effect), res$beta_eff, pch = 20)
abline(h = 0); abline(h = 0.1, lty=  2, col = "blue", lwd = 3) # obtained for our data

# p-value correlation between effect sizes
plot(log10(res2$nlocus_effect), log10(res2$p_eff), type = "o", pch = 20, las = 1, ylim = c(-6, 0), ylab =  expression(paste(log[10], " p-value - effect correlation")), xlab = "Number of loci")
points(log10(res$nlocus_effect), log10(res$p_eff), pch = 20)
abline(h = log10(0.018), lty=  2, col = "blue", lwd = 3) # obtained for our data

# correlation polygenic score STRICT
plot(log10(res2$nlocus_effect), log10(res2$p_pred), type = "o", pch = 20, las = 1, ylim = c(-6,0), ylab =  expression(paste(log[10], " p-value - polygenic STRICT")), xlab = "Number of loci")
points(log10(res$nlocus_effect), log10(res$p_pred), pch = 20)
abline(h = log10(0.05), lty=  2, col = "black", lwd = 1)

# correlation polygenic score with ALL effects
plot(log10(res2$nlocus_effect), log10(res2$p_pred2), type = "o", pch = 20, las = 1, ylim = c(-6,0), ylab =  expression(paste(log[10], " p-value - polygenic ALL")), xlab = "Number of loci")
points(log10(res$nlocus_effect), log10(res$p_pred2), pch = 20)
abline(h = log10(0.05), lty=  2, col = "black", lwd = 1)

#save.image("~/Downloads/tmp_GWAS_toy_simulations.RData")

# p-value correlation polygenic score (STRICT) / effects for polygenic 
pdf("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/toy_model_simulations.pdf", width = 6, height = 6)
par(mar = c(4,4,1,1))
plot(NULL, xlim = c(-6, 0), ylim = c(-6, 0), xlab = "Significance of correlation polygenic score", ylab = "Significance of correlation effects", axes = F)
pal1 <- colorRampPalette(c('blue', 'red'))(7)
i <- 1
for(nn in c(10, 20, 50, 100, 200, 500)){
  sub <- which(res$nlocus_effect==nn)
  if(nn>=200){
    points(log10(res$p_pred)[sub], log10(res$p_eff)[sub], col =pal1[i], pch = 17, cex = 1)
  } else {
    points(log10(res$p_pred)[sub], log10(res$p_eff)[sub], col =pal1[i], pch = 20)
  }
  i <- i + 1
}
axis(side = 1, at = c(-6, -4, -2, log10(0.05), 0), labels = c(expression(10^-6), expression(10^-4), expression(10^-2), 0.05, 1))
axis(side = 2, at = c(-6, -4, -2, log10(0.05), 0), labels = c(expression(10^-6), expression(10^-4), expression(10^-2), 0.05, 1), las = 1)
abline(v = log10(0.05), lty = 2)
abline(h = log10(0.05), lty = 2)
rect(xleft = log10(0.05)+0.02, xright = 0, ybottom = -6, ytop = log10(0.05)-0.02, lwd = 2, border = "gray")
legend("bottomleft", legend = c(10, 50, 100, 200, 500), pch = c(20, 20, 20, 17, 17), col = pal1)
dev.off()
# h2 is about 10-11%
summary(res2$h2)

# the "zone of interest" correponding to what we observe is when the p-value of polygenic score is not significant but that of correlation of effects can be

res[zone_of_interest<-which(res$p_eff<0.05 & res$p_pred>0.05),]
dim(res)
table(res$nlocus_effect[zone_of_interest])
table(res$nlocus_effect)

table(signif_effect_correlation = res$p_eff<0.05, signif_polygenic = res$p_pred<0.05)

table(signif_effect_correlation = res$p_eff[res$nlocus_effect==200]<0.05, signif_polygenic = res$p_pred[res$nlocus_effect==200]<0.05)
table(signif_effect_correlation = res$p_eff[res$nlocus_effect==500]<0.05, signif_polygenic = res$p_pred[res$nlocus_effect==500]<0.05)
table(signif_effect_correlation = res$p_eff[res$nlocus_effect==1000]<0.05, signif_polygenic = res$p_pred[res$nlocus_effect==1000]<0.05)
