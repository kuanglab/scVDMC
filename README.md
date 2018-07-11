scVCMD (Variance-Driven Multitask Clustering of scRNAseq Data) is a multi-task learning method for clustering multiple scRNAseq data from a patient cohort.

* Citation: 
Zhang, Huanan, Catherine AA Lee, Zhuliu Li, John R. Garbe, Cindy R. Eide, Raphael Petegrosso, Rui Kuang, and Jakub Tolar. "A multitask clustering approach for single-cell RNA-seq analysis in Recessive Dystrophic Epidermolysis Bullosa." PLoS computational biology 14, no. 4 (2018): e1006053.

** Lung Data:

nature13173-s4.txt: This file contains single cell RNA-seq expression data (log2(FPKM) values) for all 80 lung epithelial cells at E18.5 together with the   putative cell type of each cell in a .txt file. The last two are bulk values.

** mESC data:
- cell_states_condition.txt: the label files for mESC data. (the data matrix is needed to be downloaed at below)
- nCountGeneBatchAdjusted.csv: Download from http://compbio.cs.umn.edu/wp-content/uploads/2018/07/nCountGenesBatchAdjusted.csv

** Code:

*** To test the lung data, run the following two scripts

- process_Lung_data.m: Process Lung data with Matlab to generate the processed Lung data in Lung_data.mat.

- Lung_data_test.m: Run scVDMC algorithm on Lung_data.mat.

*** To test the mESC data; run the following two scripts

- process_mESC_data.m: Process mESC data with Matlab to generate the processed mESC data in mESC_data.mat.

- mESC_data_test.m: Run scVDMC algorithm on mESC_data.mat.

- scVDMC.m: scVDMC function.
