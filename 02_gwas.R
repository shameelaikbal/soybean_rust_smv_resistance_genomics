#GWAS using target VCF and imputed VCF using rMVP R-package
#soybean rust association mapping
#phenotype data for TAN transforming to a 5-9 scale from 1-5
pheno <-  read.table("phenotype_tan_crosshap.txt", header = T)
summary(pheno$Pheno)

#transform susceptible phenotypes 1-5 to 5-9
#load original phenotype data for TAN, pheno
pheno$Pheno <- pheno$Pheno + 4
write.table(pheno, "pehnotype_tan_transformed.txt", row.names = F, quote = F, sep = "\t")

tan <- read.table("pehnotype_tan_transformed.txt", header = T)
rb <- read.table("RB_crosshap_phenotype.txt", header = T)
combined <- rbind (rb, tan)

#data preparation
library (rMVP)

MVP.Data (fileVCF="target_imputed_filtered_for_DR0.8_maf0.03_merged.vcf",
          filePhe="rust_phenotype.txt",
          fileKin=FALSE,
          filePC=FALSE,
          out="song_rb_tan")

library(rMVP)
genotype <- attach.big.matrix("song_rb_tan.desc")
phenotype <- read.table("song_rb_tan.phe",head=TRUE)
map <- read.table("song_rb_tan.geno.map", head = TRUE)


fMVP <- MVP(
  phe=phenotype,
  geno=genotype,
  map=map,
  #K=Kinship,
  #CV.GLM=Covariates,     ##if you have additional covariates, please keep there open.
  #CV.MLM=Covariates,
  #CV.FarmCPU=Covariates,
  nPC.GLM=5,      ##if you have added PC into covariates, please keep there closed.
  nPC.MLM=5,
  nPC.FarmCPU=4,
  priority="memory",       ##for Kinship construction
  #ncpus=10,
  vc.method="BRENT",      ##only works for MLM
  maxLoop=10,
  method.bin="static",      ## "FaST-LMM", "static" (#only works for FarmCPU)
  #permutation.threshold=TRUE,
  #permutation.rep=100,
  threshold=0.05,
  method=c("GLM", "MLM", "FarmCPU"),
  out="combined_song_GWAS"
)
