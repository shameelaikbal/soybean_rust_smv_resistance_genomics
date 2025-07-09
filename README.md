# Soybean Rust and SMV Resistance Genomics

## Overview
Comprehensive genomics analysis investigating resistance to soybean rust (*Phakopsora pachyrhizi*) and soybean mosaic virus (SMV) using bioinformatics approaches.

**Phenotype data:**
- Soybean rust data from USDA-GRIN for 2815 soybean accessions
- Soybean mosaic virus data from He et al., 2025, originally retrieved from USDA-GRIN

**Genotype data:**
- **Reference panel for genotype imputation**: Valliyodan et al. (2021) - 399 soybean PIs
- **Target VCF**: Song et al. dataset

---

## Analysis Pipeline

### **1: Imputation**
- Beagle 5.4 phasing and imputation
- Quality filtering (DR2 > 0.8, MAF > 0.03)
- **Scripts**: `01_imputation.sh`

### **2: GWAS**  
- Rust resistance association analysis using rMVP
- Comparison of original vs imputed results
- **Scripts**: `02_gwas.R`

### **3: Linkage disequilibrium analysis**  
- LD analysis around genomic regions of interest
- Selection of haplotype windows
- **Scripts**: `03_LD.sh`

### **4: Local Haplotyping for soybean rust**
- Fine-mapping around significant SNPs
- Separate analysis for rust resistant and susceptible populations
- **Scripts**: `03_local_haplotyping.R`

### **4: Cross-Population Analysis**
- Modified CrossHap for haplotype transferability
- Resistant marker group transfer → susceptible 
- **Scripts**: `04_cross_population_analysis.R`

### **5: Comparative Genomics**
- Synteny analysis with common bean
- **Scripts**: `05_synteny_analysis.R`


## Code Attribution

- **crosshap**: Modified for cross-population marker group transfer analysis. Original software: [Marsh et al., 2023]. Licensed under MIT License.


