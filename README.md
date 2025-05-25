# Soybean Rust and SMV Resistance Genomics

## Overview
Comprehensive genomics analysis investigating resistance to soybean rust (*Phakopsora pachyrhizi*) and soybean mosaic virus (SMV) using advanced population genomics approaches.

**Phenotype data:**
- Soybean rust data from USDA-GRIN for 2815 soybean accessions
- Soybean mosaic virus data from He et al., 2025, originally retrieved from USDA-GRIN

**Genotype data:**
- **Reference panel for genotype imputation**: Valliyodan et al. (2021) - 399 soybean PIs
- **Target VCF**: Song et al. dataset

---

## Analysis Pipeline

### 📊 **Phase 1: Imputation**
- Beagle 5.4 phasing and imputation
- Quality filtering (DR2 > 0.8, MAF > 0.05)
- **Scripts**: `scripts/01_imputation/`

### 🔬 **Phase 2: GWAS Analysis**  
- Rust resistance association mapping
- Comparison of original vs imputed results
- **Scripts**: `scripts/02_gwas/`

### 🧬 **Phase 3: Local Haplotyping**
- Fine-mapping around significant SNPs
- Separate analysis for rust resistant and susceptible populations
- **Scripts**: `scripts/03_local_haplotyping/`

### 🔄 **Phase 4: Cross-Population Analysis**
- Modified CrossHap for haplotype transferability
- Resistant marker group transfer → susceptible 
- **Scripts**: `scripts/04_cross_population_analysis/`

### 🌐 **Phase 5: Comparative Genomics**
- Synteny analysis with common bean
- **Scripts**: `scripts/05_synteny_analysis/`


## Code Attribution

- **CrossHap**: Modified for cross-population marker group transfer analysis. Original software: [citation]. Licensed under MIT License.


