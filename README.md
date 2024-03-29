# Intro

This is a repo with the (preprocessed) data and code to produce the figures in (add elif citation).

## R version and packages

The analysis was done using R version 4.1.2 with the following packeges: 

- ComplexHeatmap_2.10.0 
- ggalluvial_0.12.5 
- ggpubr_0.4.0    
- RColorBrewer_1.1-2
- testthat_3.1.4 
- tidyverse_1.3.1 
- valr_0.6.6  
- vdiffr_1.0.5 

A singularity image exists and can be pulled from (to be added)

## Data in data folder

The subfolder data contains a set of processed datasets (for documentation on how to generate those data from the raw data see methods in paper). Those sets are plotted in the paper. In addition the gff file used throughout the analysis is include.

### **seedling_26**
Data from the 26 states model of WT in seedlings. 
1. State assignments across the genome
2. Emission matrix for the model

### **seedling_leaf_15**
Data from the 15 states concatinated model of WT in seedlings and in leaves.   
1. State assignments for leaves and seedlings respectally 
2. Emission matrix for the model

### **leaf_16**
Data from the 16 states concatinated model of WT and ddm1 mutant in leaves.   
1. State assignments for WT and ddm respectally 
2. Emission matrix for the model


### **expression**
Expression data (mean TPM):
1. WT and ddm1 in leaves
2. WT in seedlings 

### **methyl_dnase_mnase**
Data from seedlings, bed files generated by bedtools 
1. CG/CHG/CHH 
2. DNase
3. MNase


### **EvalSubset**
Output from EvalSubset runs from ChromHMM.
1. allMarks.txt: 
2. noH2A.txt 
3. noH2B.txt 
4. noH3.txt
5. noH3mod.txt
6. noVariants.txt
7. onlySeq_Mendes.txt:

See paper for details.

## Reproduce Figures 

The figures (except the deeptools generated ones) used in the paper (before illustrator polishing) can be generated either by running the separete scripts (Fig2_1.R, Fig3_1,R, Fig3_2.R, Fig4.R ) or by running the script all_and_test.R. The second option will also check that the resulting figures are indeed 100% the same as the one used in the paper. Note that the tests are very sensitive and can fail due to small variations introduced by eg. package versions. 

