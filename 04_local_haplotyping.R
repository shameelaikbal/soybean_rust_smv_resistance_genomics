# local haplotyping around regions of interest associated with soybean rust done using crosshap R-package Marsh et al. 2023
# documentation for crosshap here  https://jacobimarsh.github.io/crosshap/articles/Getting_started.html

#crosshap libraries
library(crosshap)
library (ggplot2)
library (tidyverse)
library (clustree)
library (dplyr)
library (tidyr)

#data preparation
# input the delimited vcf of local genomic regions of interest around the GWAS-SNP on chromosome 18, Gm18_56100116, within only the rust-resistant RB population
#removed all non bi allelic snps using bcftools
#bcftools view -v snps -m2 -M2 -o output.vcf filtered_data_20_11056210_nld_200kb.vcf 

vcf_18 <- read_vcf ("region_200kb_RB_chr18_56100116.vcf")
head (vcf_18)

#Update the ID column
vcf_18 <- vcf_18 %>%
  mutate (ID = paste0("Gm18_", POS))
head(vcf_18)

#input LD matrix
#./plink --allow-extra-chr --r2 square --vcf region_200kb_chr20_46854737.vcf 
#replace nans in plink.ld with zeroes on CL sed 's/nan/0/g' plink.ld > chr15_plink.ld
LD_18 <-  read_LD ("chr18_56100116.ld", vcf = vcf_18)

pheno_18 <- read_pheno("RB_crosshap_phenotype.txt")
head(pheno_18)
metadata_18 <- read_metadata("RB_metadata.txt")

#Add minimum marker group for SNP count
MGmin_gm <- 8 # depending on the SNP density in the region
str(MGmin_gm)

##Add list of epsilon values to run haplotyping on
epsilon_gm <- c(0.2,0.4,0.6,0.8,1.0,1.2)
str(epsilon_gm)

#Run the haplotyping at all provided epsilon value
run_haplo_18 <- run_haplotyping (vcf = vcf_18,
                                 LD = LD_18,
                                 pheno = pheno_18,
                                 metadata = metadata_18,
                                 epsilon = epsilon_gm ,
                                 MGmin = MGmin_gm)

run_haplo_18

#Provide phenotype data and parameters used to create haplotype objects
#Add type = 'MG' to ensure it summarizes Marker Groups rather than haplotypes and type hap to visualise haplotypes
hap_clustree_chr18 <- clustree_viz(HapObject = run_haplo_18,
                                   type = 'hap')
hap_clustree_chr18

mg_clustree_chr18 <- clustree_viz(HapObject = run_haplo_18,
                                  type = 'MG')
mg_clustree_chr18


hap_viz_18 <- crosshap_viz (HapObject = run_haplo_18, epsilon = 0.2)
hap_viz_18

#choosing 0.2

ggsave(plot = hap_viz_18,"crosshap_chr18_56100116_RB_eps0.2_mgmin8.jpg", width = 10, height = 8, dpi = 300, units = "in" )

#the position plot of these mg grups are here
hap_viz_18_pos <- crosshap_viz(HapObject = run_haplo_18, epsilon = 0.2, plot_left = "pos")
hap_viz_18_pos

ggsave(plot = hap_viz_18_pos, "crosshap_chr18_RB_pos_eps0.2_mgmin8.jpg", width = 10.5, height = 8, dpi = 300, units = "in" )

#SAve these files

#Print first lines of haplotype object buckets for epsilon = 0.6 results.
#The Indfile reports the haplotype assigned to each individual
#The Hapfile reports the identified haplotypes, and their Marker Group combinations
#The Varfile reports information for each SNP, including Marker Group assignments
ind_ffile <- run_haplo_18$Haplotypes_MGmin8_E0.2$Indfile
ind_ffile
write.csv(ind_ffile, file = "ind_file_E0.2_combined_mgmin8_chr18.csv", quote = FALSE, row.names = FALSE)

hap_file <- run_haplo_18$Haplotypes_MGmin8_E0.2$Hapfile
hap_file

write.csv(hap_file, file = "hap_file_E0.2_RB_chr18.csv", quote = FALSE, row.names = FALSE)

var_file <- run_haplo_18$Haplotypes_MGmin8_E0.2$Varfile
write.csv(var_file, file = "var_file_E0.2_RB_mgmin8_chr18.csv", quote = FALSE, row.names = FALSE)
