# AS-reprogramming-KD-OE
Analysis of the AS changes in MEF reprogramming regulated by the knockdown of CPSF3/hnRNP UL1 (KDs dataset) or by the overexpression of TIA1 (OE dataset).
Published in [Vivori et al., 2020](https://www.biorxiv.org/content/10.1101/2020.09.17.299867v1).
Datasets available in GEO Database: Accession number GSE158633. 

### 00. Run vast-tools and edgeR and apply initial filtering
- Alternative splicing analysis ([00_Run_VTS.sh](00_Run_VTS.sh)):
  - `Vast-tools` v2.2.2 align (mm10), combine, compare. See [vast-tools webpage](https://github.com/vastgroup/vast-tools) for details.
  - Filtering of vast-tools output for reads coverage with [VTS_INCL_filtering.R](https://github.com/cvivori/useful-cluster-scripts/blob/master/README.md#vts_incl_filteringr), to consider only AS events with a minimum of 10 actual reads per sample (0 “N” values allowed for each dataset).
  - Extension of the filtered INCLUSION tables, including all the AS events differentially spliced in at least one comparison of each dataset (and their dPSI in all the comparisons). See [VTS_add_dPSI_toINCL.R](https://github.com/cvivori/useful-cluster-scripts/blob/master/README.md#vts_add_dpsi_toinclr).

- Gene Expression analysis ([00_Run_edgeR.R](00_Run_edgeR.R)):
  - Import of gene counts (from STAR mapping) 
  - Filtering for minimum 1 cpm per sample and 5 cpm in at least in 33% of samples (8 for KDs dataset and 2 for OE dataset).
  - Calculation of cpm values and differential expression analysis.

### 01. Define sets of CPSF3/UL1/TIA1-dependent and -independent AS events
- Calculate efficiency of CPSF3/hnRNP UL1 knockdown or TIA1 overexpression ([01_EdgeR_KDOE_efficiency.R](01_EdgeR_KDOE_efficiency.R)).
- Define sets of CPSF3/UL1/TIA1-dependent and -independent AS events and their overlaps ([01_KDOE_VTS_Analysis.R](01_KDOE_VTS_Analysis.R)).

### 02. Study CPSF3/UL1/TIA1-dependent and -independent AS events
- Analysis of sequence features and RNA-binding protein motif enrichment with [Matt](http://matt.crg.eu/) ([02_Run_Matt_VTS.sh](02_Run_Matt_VTS.sh)).
- Enrichment of senescence/SASP-associated genes in GO terms enriched in TIA1-dependent AS exons ([02_GO_TIA1dep_SenescSASPgenes.R](02_GO_TIA1dep_SenescSASPgenes.R)).

