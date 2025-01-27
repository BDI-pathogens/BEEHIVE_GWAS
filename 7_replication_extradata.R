rm(list =ls())
setwd("~/Dropbox (Infectious Disease)/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/")
library(seqinr)
gag_aln <- read.fasta("7_extradata_replication/nucleotides_wHXB2/gag.fasta")
nef_aln <- read.fasta("7_extradata_replication/nucleotides_wHXB2/nef.fasta")
a <- read.csv("7_extradata_replication/G2G_SHCS/G2G_shcs.spvl.txt", sep = "\t")

idx_in_gag_aln <- match(
  paste0("SHCS$", a$ID),
  names(gag_aln)
)
idx_in_nef_aln <- match(
  paste0("SHCS$", a$ID),
  names(nef_aln)
)

gag_aln <- gag_aln[idx_in_gag_aln]
nef_aln <- nef_aln[idx_in_nef_aln]

stopifnot(all(names(gag_aln) == paste0("SHCS$", a$ID)))
table(names(nef_aln) == paste0("SHCS$", a$ID), useNA = "ifany")

a$gag_nuc <- unlist(lapply(gag_aln, function(x)x[725]))
a$gag_nuc <- factor(a$gag_nuc, levels = c("c", "a", "g", "n"))

a$nef_nuc <- unlist()


nef_tmp <- lapply(nef_aln, function(x)x[212])
a$nef_nuc <- unlist(lapply(nef_tmp, function(x) if(is.null(x)) return(NA) else return(x)))
a$nef_nuc <- factor(a$nef_nuc, levels = c("g", "a", "c", "r", "n"))

# added 06 Oct 2022: eliminate the ids already in BEEHIVE discovery
ids_in_beehive_discovery <- read.csv("~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/7_extradata_replication/overlapBEEHIVEID_in_discovery_SHCS.csv")$ID
dim(a)
print(length(ids_in_beehive_discovery))
a <- a[!a$ID %in% ids_in_beehive_discovery, ]
dim(a)

s1 <- summary(
  lm1 <- lm(a$spVL ~ a$gag_nuc)
)

s2 <-summary(
  lm2 <- lm(a$spVL ~ a$nef_nuc)
)

ci1 <- confint(lm1)["a$gag_nuca", ]
ci2 <- confint(lm2)["a$nef_nuca", ]

t.test(a$spVL[a$gag_nuc=="c"], a$spVL[a$gag_nuc=="a"], alternative = "less")
t.test(a$spVL[a$nef_nuc=="g"], a$spVL[a$nef_nuc=="a"], alternative = "greater")


hist(a$spVL[a$gag_nuc=="c"])
hist(a$spVL[a$gag_nuc=="a"], add = T, col = "red")

# summarise results in a data frame
df <- data.frame(position = c(1514, 9008), reference = c("c", "g"), alternative = c("a", "a"),
                 N_reference = c(sum(a$gag_nuc=="c"), sum(a$nef_nuc=="g", na.rm = T)),
                 N_alternative = c(sum(a$gag_nuc=="a"), sum(a$nef_nuc=="a", na.rm = T)),
                 effect = c(s1[[4]]["a$gag_nuca", "Estimate"], s2[[4]]["a$nef_nuca", "Estimate"]),
                 lower = c(ci1["2.5 %"], ci2["2.5 %"]), upper = c(ci1["97.5 %"], ci2["97.5 %"]),
                 p = c(s1[[4]]["a$gag_nuca", "Pr(>|t|)"], s2[[4]]["a$nef_nuca", "Pr(>|t|)"])
                 )
for(mycol in c("effect", "lower", "upper", "p")) df[mycol] <- signif(df[mycol], 3)
df
write.csv(x = df, file = "7_extradata_replication/summary_replication_effects.csv", row.names = F)


