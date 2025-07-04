# BEEHIVE GWAS

Code and minimal dataset for the BEEHIVE GWAS presented in "The genetic architecture of HIV-1 virulence". The minimal dataset (in the `minimal_dataset` folder) enables to reproduce the figures.

## Minimal dataset

The folder `minimal_dataset` contains the data necessary to reproduce the figures in the paper.

### Figure 1:

This figure is drawn by `draw_tree_Figure1.R`.

`minimal_dataset/tr2.tr` is a phylogenetic tree of the 2,249 sequences used for the GSVL GWAS.

`minimal_dataset/subtypes2.csv` contains the subtypes in the same order as the tip labels.

`minimal_dataset/vl2.csv` contains the viral load values in the same order as the tip labels.

### Figure 2A

This figure is drawn by `2.5_inspect_GWAS_results.R`.

The tables of the form `minimal_dataset/GWAS_type_phenotype.csv` contain the results of GWASes. "type" can be "prot", "kmer", "len", "len3", "freq". Phenotype can be "BEEHIVE_LVL", "spvl_normalised_adjusted", "spvl_adjusted", "CD4_slope".


### Figure 2B, C

The figure is drawn by `8_temporal_trends_enrichment.R`. 

`minimal_dataset/frequency_through_time_hit1.csv` and `minimal_dataset/frequency_through_time_hit2.csv`

Tables containing the frequency of mutations through time with confidence intervals (Figure 2B).

`minimal_dataset/enrichment_df.csv`

Table containing the results of the enrichment analyses (Figure 2C).

### Figure 3

This figure is drawn by `5_post_replication_analysis_v2.R` (BEEHIVE additional, left panels) and `9_Gabrielaite/analyse_Gabrielaite.R` (INSIGHT, right panels).

Tables of the form `minimal_dataset/GWAS_main_additional_SNPonly_combined_phenotype.csv`

Tables containing, for each phenotype ("BEEHIVE_LVL", "spvl_normalised_adjusted", "spvl_adjusted", "CD4_slope"), the combined results of the main GWAS and the GWAS in additional BEEHIVE dataset. These are used for the left panels of figure 3.

Tables of the form `minimal_dataset/GWAS_main_INSIGHT_combined_CD4_phenotype.csv`

Tables containing, for each phenotype ("BEEHIVE_LVL", "spvl_normalised_adjusted", "spvl_adjusted", "CD4_slope"), the combined results of the main GWAS and the GWAS on INSIGHT dataset. These are used for the right panels of figure 3.

`minimal_dataset/idx_prot_hits_main.csv`

Small table containing the index of rows corresponding to the 2 hits in the protein GWAS of the BEEHIVE main dataset.


## 0. COMPUTING THE G MATRICES AND PREPARING THE DATA

`0_prepare_kmers/prepare_kmer_Chris.R`	

A code that takes as input the merged individual consensuses, and prepares the `lengths.RData`, `lengths3.RData`, `kmers.RData` files containing all information on length variants, length modulo 3 variants, and kmer variants at all positions in the genome

`0_prepare_aminoacids/prepare_aminoacids.R`		

A code that takes as input the protein alignment and outputs the `prot.RData` file containing all information on protein variants at all positions in the genome. Protein alignment is generated by Typewriter.

`0_prepare_base_frequencies/prepare_base_frequencies.R`

A code that takes as input the base frequency files and outputs the `base_frequencies.RData` file containing all information on base frequencies at all positions in the genome

`0_prepare_data_v3.R`

A code that prepares the genomic and epidemiological data and generates the `prepared_phenotypename_alnname_condition_RData` files in the `Gmatrix_data` folder. These files contain, most importantly, the matrices defining the variance-covariance structure of the phenotype data. This depends on the genetic structure and the (unstructured) error variance. This code outputs various .RData files corresponding to different phenotypes and different filtering conditions on the genomic data

`Gmatrix_data`

A folder containing:

* a set of .RData file containing prepared genomic (G, K matrices) and epi data to run the LMM on, for various phenotypes and filtering conditions on the genomic data. This is produced by the code `0_prepare_data_v3.R`
									
* kmers.RData, a RData file containing kmer info in the form (mainly) of a G matrix. This is produced by `prepare_kmers/prepare_kmer_Chris.R`. This contains "G_kmer", "reference_kmer", "alternative_kmer", "window.ends_kmer", "window.starts_kmer", "idsG_kmer".
									
* lengths.RData, a RData file containing length info in the form (mainly) of a G matrix. This is produced by `prepare_kmers/prepare_kmer_Chris.R`. This contains "G_len", "reference_len", "alternative_len", "window.ends_len", "window.starts_len", "idsG_len"
									
* lengths3.RData, a RData file containing length info modulo 3 in the form (mainly) of a G matrix. This is produced by `prepare_kmers/prepare_kmer_Chris.R`. This contains "G_len3", "reference_len3", "alternative_len3", "window.ends_len3", "window.starts_len3", "idsG_len3".

* base_frequencies.RData, a RData file containing base frequency info in the form (mainly of a G matrix). This is produced by `prepare_base_frequencies/prepare_base_frequencies.R`. This contains "G_freq", "reference_freq", "alternative_freq", "window.ends_freq", "window.starts_freq", "idsG_freq"

`0_createG.R`						

Functions used to create the G matrix - the G matrix is a matrix of 1, 0 and NA describing the genomic alignment.

## 1. PERFORMING LINEAR MIXED MODELS INFERENCE FOR VARIOUS G MATRICES

`1_do_inference_Rscript_v2.R`		

Performs inference of the random effect parameters, taking as input the filename of the RData file, that gives the criterion index (corresponding to the different filtering conditions in the "criteria_table.csv" file), the phenotype, the alignment. It then loads the various RData files (in folder `Gmatrix_data`) containing the variant information among others. This code can be used on a cluster. It outputs a simple csv file with explicit name, with the parameters, likelihood, heritability for these criteria.

`results_Gmatrix_comparison`	

Folder containing the results of LMM inference for each of the filtering conditions, in the form of csv files with explicit names. `merge_results.R` produces a table "all_results.csv" compiling all these csv files.

`function_mixed_model.R`

A set of functions used for the linear mixed model inference.

## 2. PERFORMING THE GWAS

`2_firstGWAS_fromLM_Rscript_v6.R`					

Performs a GWAS. Takes as argument the name of the "prepared_phenotype_aln_$criterion.RData" file. Uses the maximum likelihood estimation of random effects for all positions tested, which considerably speeds up the calculation. Perform GWAS of k-mer variants, length variants, length%3 variants, base frequency variants.


`functions_2_draft_GWAS_v2.R`

Functions used in the GWAS code			

`GWAS_results/firstGWAS_phenotype_aln_idx.RData`

Files containing the results of the GWAS, for different filtering used. phenotype is "BEEHIVE_LVL", "spvl_adjusted", "spvl_normalised_adjusted", "CD4_slope"; aln is the alignment used ('full', 'A12', etc.).

`2.5_inspect_GWAS_results.R`								

Loads the `firstGWAS_phenotype_aln_idx.RData` files and draw qqplots and tables of the significant positions within R. It considers all phenotypes and all variants in turn within that one code. It also applies a Bonferroni correction to the p-value. The Bonferroni correction with an EFFECTIVE number of tests based on the variance-covariance structure of the data is used. It also computes inflation.



## 3. REOPTIMISING THE MODEL

`3_optimise_model_v2.R`								

Takes as argument a phenotype and a p-value threshold. Loads the prepared data (full genome, criterion #17). Selects all hits that are significant at this p-value threshold (bonferroni corrected). Re-optimise the full linear model that includes all these hits. Decomposes the variance of the phenotype explained by epidemiological and genomic factors. Results are stored in the folder `3_reoptimisation`.

## 4. REPLICATE THE FINDINGS

`4_replication/0_prepare_replication_data.R` 				

Prepares the phenotype data for the additional BEEHIVE dataset. It takes as argument the consensus alignment `/Data/Sequences_replication/BEEHIVE_Oxford_GlobalAln_2019-11-06_MultiSeqsPerPat_clean.fasta` and the epi data in `/Data/MetaData/CleanMetaData/BEEHIVE_summary24.10.2019.csv` and the subtype file `Data/subtype/subtype_BEEHIVE_Oxford_GlobalAln_2019-11-06_MultiSeqsPerPat_clean_final.csv`(+ miscellaneous files: dual infection, DRM, which have not been updated for the additional dataset). It outputs 4 files of the form `/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/Gmatrix_data/prepared_pheno_adjusted_full_17.RData` with the 4 phenotypes "BEEHIVE_LVL", "spvl_adjusted", "spvl_normalised_adjusted", "CD4_slope". This code also exports a list of all individuals-dates kept in the alignment for the additional analysis, in `/Code/DevelopMethods/LMM/SolvingLMMinR/4_replication/Gmatrix_data/all_replication_sequences_tokeep.csv`

`4_replication/0_prepare_kmers/deduplicate_merged_consensuses_replication.R`

This is the code that outputs the file `MergedIndividualConsensuses_renamed_clean_noduplicate_2019-11-20.csv` and copies the within-host frequency files that are duplicated ids in a separated "ExtraSeq" folder. It needs to be ran before `script_prepare_kmers.sh` and `script_prepare_base_frequencies.sh`.

`4_replication/0_prepare_kmers/script_prepare_kmers.sh`

Prepare the kmer data and length data for the additional BEEHIVE dataset. It uses the same script as for the discovery dataset, called `/0_prepare_kmers/prepare_kmers_Chris.R`, and the genomic dataset in `/Data/Sequences_replication/MergedIndividualConsensuses_renamed_clean_noduplicate_2019-11-20.csv`. The result is in 3 files in `/4_replication/Gmatrix_data/`. The 3 files are `kmers.RData`, `lengths.RData`, `lengths3.RData`

`4_replication/0_prepare_base_frequencies/script_prepare_base_frequencies.sh`

Prepare the base frequency data for the additional BEEHIVE dataset. It uses the same script as for the discovery dataset, called `/0_prepare_base_frequencies/prepare_base_frequencies.R`, and the dataset in `/Data/Sequences_replication/BaseFreqs_renamed_clean/`. The results are in `/Data/Sequences_replication/BaseFreqs_renamed_clean/` and `/4_replication/Gmatrix_data/base_frequencies.RData`

`4_replication/compare_base_frequencies.R`				

Compare the within-host base frequencies in the discovery and additional BEEHIVE datasets, for the patients that overlap in discovery and additional BEEHIVE datasets, and at the significant positions(those with p < 0.05 with bonferroni correction). This uses the `prepared base_frequencies.RData` files (for discovery and additional BEEHIVE datasets). This results in the figure `within_host_frequencies_significantpositions.pdf`. Also systematically compares the allele frequencies at all positions and store results in the table `all_correlations_withinhostfrequencies.csv`.

`4_replication/1_script_do_inference.sh` 					

Does the inference of random effet of the generalised linear model for each phenotype. Stores the results in `/4_replication/Gmatrix_data/`, in 4 files of the form `prepared_pheno_full_17.csv` (for the 4 phenotypes)

`4_replication/2_script_do_GWAS.sh` 					

Does the GWAS for each phenotype by calling the script `2_firstGWAS_fromLM_Rscript_v6.R`. Need to check the options of `2_firstGWAS_fromLM_Rscript_v6.R` at the beginning of the file before running. The results are in `/4_replication/GWAS_results/` in 4 files of the form `firstGWAS_pheno_full_17.RData` (for the 4 phenotypes)

`4_replication/3_power_calculation.R` 					

Calculates the power of the additional dataset to find associations identified in the discovery dataset. Details of each variants number is found in the results of the GWAS, under `GWAS_results/firstGWAS_pheno_full_17.RData`. Outputs the result of the power calculation in file of the form, e.g. `GWAS_results/spvl_normalised_adjusted_GWAS_kmer_1_withreplicationpower.csv`

## 5. EXAMINE THE RESULTS OF ADDITIONAL ANALYSIS

`5_post_replication_analysis.R`

Looks systematically at effects found in additional BEEHIVE analysis, correlates them with effects found in discovery, and produces figures of the form `4_replication/replication_correlation_effect_sizes_phenotype.pdf`

## 6. PREDICT PHENOTYPE IN BEEHIVE ADDITIONAL DATA

`6_predict_phenotypes_additional.R` 

Uses the results from discovery GWAS to predict the phenotype in the BEEHIVE additional dataset using the polygenic score. The arguments are the name of the phenotype and the p-value threshold used for this prediction. This returns the table `prediction_polygenic_score.csv`

## 7. ADDITIONAL ANALYSIS OF SHCS DATA

`7_replication_extradata.R`

Uses the gag and nef alignments of the extra additional data (the subset of Swiss sample from Bartha et al. eLife 2016 paper, not already in BEEHIVE). These alignments are in `7_extradata_replication/nucleotides_wHXB2/gag.fasta` and `7_extradata_replication/nucleotides_wHXB2/nef.fasta`. And the corresponding SPVL data in `7_extradata_replication/G2G_SHCS/G2G_shcs.spvl.txt`. From these, identify positions 1514 and 9008 and looks at their impact on viral load. Exports the data in `7_extradata_replication/summary_replication_effects.csv`.


`7_extradata_replication`					

Folder that contains the data from the SHCS with HIV sequences and human HLA genotypes. Used as an extra additional dataset.

## 8. ANALYSES OF TEMPORAL TRENDS IN THE TWO MAIN VARIANTS, AND ENRICHMENT IN CTL ESCAPE / DRM

`8_temporal_trends_enrichment.R`

Looks at temporal trends in main alternative allele at positions 1514 and 9008. Results in the figure `temporal_trajectories_hits.pdf`. Second, looks at enrichment in CTL escape mutations or drug resistance mutations among the set of SNPs variants significant at some threshold. Results in the figure `enrichment_analysis.pdf`.

## 9. CORRELATION OF DISCOVERY EFFECT SIZES WITH GABRIELAITE ET AL DATA

`9_Gabrielaite`

Folder containing the code and data necessary to correlate effect sizes in our discovery dataset and the additional Gabrielaite dataset.

`column_ID_BETA_P_marc_analysis_with_plink_plink_analysis_plink2.VL.glm.linear`

Additional data from Gabrielaite et al., showing effect size for all tested variants.

`analyse_Gabrielaite.R`

Code matching discovery effect sizes with Gabrielaite effect sizes, and producing the figure `figure_discovery_INSIGHT` and the table `discovery_Gabrielaite_correlations.csv`.

## 10. ALTERNATIVE GWAS BASED ON LASSO REGRESSION

`10_lasso_v2.R`

Does the alternative GWAS for each phenotype. The argument is the phenotype (either "spvl_adjusted", "BEEHIVE_LVL", "spvl_normalised_adjusted", "CD4_slope"). Loads the data and the first GWAS results. Loads the additional data. Then creates a list of variants present in at least N=10 individuals, and also present in replication. Does the lasso regression on the main dataset, then tests the lasso algorithm on the replication dataset. Draws the figure 4 in the paper.

## 11. TOY MODEL FOR THE GWAS

`11_toy_simulation_GWAS.R`

Performs a GWAS on simulated data to understand the difficulties to replicate hits and the failure of the polygenic score, while the inferred effects of variants are weakly correlated between main and additional datasets. This is done for a number of replicates, and testing various number of loci explaining about 10% of the genetic variance (on top of the variance explained by population structure). Draws the Supplementary figure 2 in the paper.

## MISCELLANEOUS

`bonferroni`

A folder containing tables and code to compute the effective number of tests for each type of variant tested, as well as the resulting effective number of tests. The tables in starting in `n_effective_tests_` contain the number of effective tests. The distance tables that are used to compute this number are not saved because they are enormous. The code `get_bonferroni_dfs.R` uses the `n_effective_tests_` tables to compute the total number of effective tests.

`function_mixed_model.R`

A set of function used to resolve linear mixed models.

`draw_tree_Figure1.R`

A code to draw the phylogenetic tree in figure 1 of the paper. This calls the set of function in `function_plot.phylo.R`.

`function_plot.phylo.R`

A set of functions to draw a nice phylogenetic tree with metadata.

`GlobalAln_to_HXB2coords.csv`

A table with the coordinates of the global alignment and the corresponding coordinates of the HXB2 reference.

`hxb2_annotated_Chris.csv`

An HXB2 reference with custom annotation.

`concatenated_PRO_coords.csv` a list giving the amino-acid coordinates of each amino-acid in the protein alignment `concatenated_PRO.fasta`.







