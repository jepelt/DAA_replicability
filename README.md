# DAA_replicability

In this this repository are the R codes that can be used to perform all the analyses and to create the figures included in the manuscript of a research article entitled "Elementary methods provide more replicable results in microbial differential abundance analysis". A preprint of the article is available in arXiv (https://arxiv.org/abs/2404.02691).

Data_curation folder includes the codes used in the curation of data, e.g. to filter and split datasets (needed in the split-data analyses), and to combine the datasets in one file. Running the script for shotgun datasets (data_curation_shotgun.R) automatically downloads the raw datasets but the raw 16S datasets need to be separately downloaded from https://zenodo.org/records/840333 before running the script data_curation_16s.R. The curated data used in the analyses are, however, also included in the Analysis folder (data_171023.rds, data_meta_041023.rds, data_prevl_061023.rds).

The scripts run_daa_methods.R and run_daa_methods_covariates.R in the Analysis folder can be used to run the differential abundance analysis methods and to produce files res_main.rds and res_covariates.rds, respectively. (Running these scripts may take rather long time.) The res_... files are needed to calculate the evaluation metrics and to create the figures appearing in the article. To create each figure and to calculate the metrics needed in it, the script named after the figure should be run.

Note that the method LogR is now commented in the scripts as uncommenting the logistf::logistf function may cause the base paste function to work improperly. That would further affect the working of some functions. LogR may thus need to be run in a separate R session.

Scripts needed for the additional analyses (related to absolute abundances, Figure A6) are indicated with the suffix “abs” or “absolute”.  

 
