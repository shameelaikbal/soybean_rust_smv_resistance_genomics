# Soybean Rust and SMV Resistance Genomics

## Overview
Comprehensive genomics analysis investigating resistance to soybean rust (*Phakopsora pachyrhizi*) and soybean mosaic virus (SMV) using bioinformatics approaches.

**Phenotype data:**
- Soybean rust data from USDA-GRIN for 2815 soybean accessions [Miles et al. 2006](https://doi.org/10.1094/PHP-2006-0104-01-RS)
- Soybean mosaic virus data from [He et al., 2025](https://doi.org/10.3390/ijms26052106) originally retrieved from USDA-GRIN

**Genotype data:**
- **Reference panel for genotype imputation**: [Valliyodan et al. (2021)](https://www.nature.com/articles/s41597-021-00834-w) - 399 soybean PIs
- **Target VCF**: [Song et al. (2015)](https://www.soybase.org/tools/snp50k/) dataset

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

### **4: LDBlockShow** 
- **Scripts**: `04_LDBlockShow.sh`

### **5: Local Haplotyping for soybean rust and SMV**
- Fine-mapping around significant SNPs
- Separate analysis for rust resistant and susceptible populations and SMV
- **Scripts**: `05_local_haplotyping.R`

### **6: Cross-Population Analysis**
- Modified crosshap for haplotype transferability
- Resistant marker group transfer → susceptible 
- **Scripts**: `06_cross_population_analysis.R`

### **7: Comparative Genomics**
- Synteny analysis with common bean
- **Scripts**: `07_synteny_analysis.R`


## Code Attribution

- **crosshap**: Modified for cross-population marker group transfer analysis. Original software: [Marsh et al., 2023](https://doi.org/10.1093/bioinformatics/btad518). Licensed under MIT License.


