# this script uses 2_firstGWAS_fromLM_Rscript_v6.R
# first argument is the name of the RData file
# second argument is directory where data files are stored
# note: for now, we need to manually remove the protein GWAS in the R code because replication alignment has not been translated yet
# note: a number of options need to be tuned directly in the 2_firstGWAS_fromLM_Rscript_v6.R file, importantly whether information on effects / pvalues are saved or not
cd ~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/
Rscript 2_firstGWAS_fromLM_Rscript_v6.R prepared_BEEHIVE_LVL_full_17.RData "~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/Gmatrix_data/"
Rscript 2_firstGWAS_fromLM_Rscript_v6.R prepared_spvl_normalised_adjusted_full_17.RData "~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/Gmatrix_data/"
Rscript 2_firstGWAS_fromLM_Rscript_v6.R prepared_spvl_adjusted_full_17.RData "~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/Gmatrix_data/"
Rscript 2_firstGWAS_fromLM_Rscript_v6.R prepared_CD4_slope_full_17.RData "~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/Gmatrix_data/"
