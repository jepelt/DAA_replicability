# DAA_replicability

In this this repository are the R codes that can be used to perform all the analyses and to create the figures included in the manuscript of a research article entitled "Elementary methods provide more replicable results in microbial differential abundance analysis". A preprint of the article is available in arXiv (https://arxiv.org/abs/2404.02691).

Data_curation folder includes the codes used in the curation of data, e.g. to filter and split datasets (needed in the split-data analyses), and to combine the datasets in one file. Running the script for shotgun datasets (data_curation_shotgun.R) automatically downloads the raw datasets but the raw 16S datasets need to be separately downloaded from https://zenodo.org/records/840333 before running the script data_curation_16s.R. The curated data used in the analyses are, however, also included in the Analysis folder (data_171023.rds, data_meta_041023.rds, data_prevl_061023.rds).

The main results of the paper can be replicated by using the following versions of the R packages:
ALDEx2 1.34.0, ANCOMBC 2.8.1 (for ANCOM-BC2), corncob 0.3.2, DESeq2 1.40.2, edgeR 3.42.0, fastANCOM 0.0.4, LDM 6.0.1, limma 3.56.2 (for limma-voom), LinDA 0.1.0, logistf 1.26.0 (for Firth logistic regression), MaAsLin2 1.14.1 (for t-test/linear regression), metagenomeSeq 1.42.0, rms 6.5-0 (for ORM/Wilcoxon), GUniFrac 1.8 (for ZicoSeq)


