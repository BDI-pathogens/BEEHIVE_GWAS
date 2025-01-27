# this script uses 1_do_inference_Rscript_v2.R
# first argument is the name of the RData file
# second argument is directory where data files are stored
cd ~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/
Rscript 1_do_inference_Rscript_v2.R prepared_BEEHIVE_LVL_full_17.RData "~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/Gmatrix_data/"
Rscript 1_do_inference_Rscript_v2.R prepared_spvl_normalised_adjusted_full_17.RData "~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/Gmatrix_data/"
Rscript 1_do_inference_Rscript_v2.R prepared_spvl_adjusted_full_17.RData "~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/Gmatrix_data/"
Rscript 1_do_inference_Rscript_v2.R prepared_CD4_slope_full_17.RData "~/DID/BEEHIVE_Hackathon/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/Gmatrix_data/"
