# Xpress

## Basic Expression Gene List Goals

### Yardstick Genes

Most common problem is determining baselines and comparing experiments.
 * Were these experiments measured in a different unit?  TPM, RPKM, FPKM, log2 or log10 transformed?
 
 Use the large GTEx data set to find genes that express at consistent levels across the tissue types.
 Some desired features:
  * Common Genes - ones that are often measured.
  * Consistent expression across all tissue types and samples:
  ** Poor information for tissue type classification.
  * Avoid very low expression genes.
  * Have a range of expressions levels - for internal, within sample normalization.
  * Have good correlation between mRNA and protein expression levels.
  ** this is not always that common.

## Random Point notes

 1. Expression can be complicated:
     * Cite GSC UoA - decrease mRNA resulting more expression - Overall ratios of all transcripts.
 1. Toxic Expression Levels - mutually exclusive genes - overactive pathway.
 1. Pathway level expressions and categories?
 
 Other effects
 1. Low expression levels are irrelvant for very active proteins.
     * Activating mutations.
     * Activating fusions
     * Can pathway activation / deactivation be used to identify?
 1. Micro-environment effects - cell vs tumour micro env
 1. PDL1 vs mutation burden - patterns for tumour immune attack.
