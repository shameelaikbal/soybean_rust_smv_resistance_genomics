
#' Apply predefined marker groups from resistant individuals to susceptible individuals using modified crosshap (Marsh et al., 2023) scripts
#' original software: https://github.com/jacobimarsh/crosshap


#' run_haplotyping_predefined() applies marker groups identified in resistant individuals to susceptible individuals. 
#' It skips the clustering step in original crosshap and uses the predefined marker groups to identify haplotypes in the new population.
#' @param vcf Input VCF for region of interest.
#' @param LD Pairwise correlation matrix of SNPs in region (from PLINK).
#' @param metadata Metadata input (optional).
#' @param pheno Input numeric phenotype data for each individual.
#' @param predefined_MGs Predefined marker groups dataframe with columns ID, cluster and MGs, obtained from resistant individuals.
#' @param minHap Minimum Individuals in a haplotype combination.
#' @param prefix_MGs Whether to add "MG" prefix to marker group numbers.
#' @param hetmiss_as If hetmiss_as = "allele", heterozygous-missing SNPs './N'
#' are recoded as 'N/N', if hetmiss_as = "miss", the site is recoded as missing.
#' @param het_phenos When FALSE, phenotype associations for SNPs are calculated
#' from reference and alternate allele individuals only, when TRUE, heterozygous
#' individuals are included assuming additive allele effects.
#' @param keep_outliers When FALSE, marker group smoothing is performed to
#' remove outliers.
#'
#' @export
#'
#' @returns A comprehensive haplotyping S3 object (HapObject), needed for crosshap_viz().
#'

# load libraries
library(crosshap)
library (ggplot2)
library (tidyverse)
library (clustree)
library (dplyr)
library (tidyr)

run_haplotyping_predefined <- function(vcf, LD, pheno, predefined_MGs, metadata,
                                       minHap = 9, prefix_MGs = TRUE, hetmiss_as = 'allele',
                                       het_phenos = FALSE, keep_outliers = FALSE,
                                       verbose = TRUE){
  # Reformat VCF
  if(verbose) message("Step 1: Formatting VCF")
  
  bin_vcf <- dplyr::select(vcf, -c(1,2,4:9)) %>% tibble::column_to_rownames('ID') %>%
    dplyr::mutate_all(function(x){base::ifelse(x=='0|0',0,
                                               base::ifelse(x=='1|0'|x=='0|1',1,
                                                            base::ifelse(x=='1|1',2,
                                                                         switch(hetmiss_as, "allele" = base::ifelse(x=='1|.'|x=='.|1',1,
                                                                                                                    base::ifelse(x=='0|.'|x=='.|0',0,NA)),
                                                                                "miss" = NA))))})
  
  # Extract position information from the VCF
  if(verbose) message("Step 2: Extracting position information")
  
  vcf_info <- vcf %>%
    dplyr::select(ID, POS) %>%
    dplyr::distinct()
  
  # Join position information with predefined marker groups
  if(verbose) message("Step 3: Adding position information to marker groups")
  
  predefined_MGs$MGs <- as.character(predefined_MGs$MGs)
  
  if(!"POS" %in% colnames(predefined_MGs)) {
    predefined_MGs <- dplyr::left_join(predefined_MGs, vcf_info, by = "ID")
  }
  
  # Add MG prefix if requested
  if(prefix_MGs) {
    predefined_MGs$MGs <- ifelse(
      predefined_MGs$MGs != "0", 
      paste0("MG", predefined_MGs$MGs),
      "0"
    )
  }
  
  # Check for valid marker groups
  valid_MGs <- unique(predefined_MGs$MGs)
  valid_MGs <- valid_MGs[valid_MGs != "0"]
  
  if(length(valid_MGs) == 0) {
    message("No valid Marker Groups found in predefined_MGs")
    return(NULL)
  }
  
  if(verbose) {
    message("Found ", length(valid_MGs), " marker groups: ", paste(valid_MGs, collapse=", "))
    
    for(mg in valid_MGs) {
      snp_count <- sum(predefined_MGs$MGs == mg)
      message("  ", mg, ": ", snp_count, " SNPs")
    }
  }
  
  # Process marker groups and identify haplotypes
  if(verbose) message("Step 4: Processing marker groups and identifying haplotypes")
  
  HapObject <- identify_haplotypes_from_MGs(
    vcf = vcf,
    bin_vcf = bin_vcf, 
    predefined_MGs = predefined_MGs,
    LD = LD,
    pheno = pheno,
    metadata = metadata,
    minHap = minHap,
    het_phenos = het_phenos,
    keep_outliers = keep_outliers,
    verbose = verbose
  )
  
  if(is.null(HapObject)) {
    message("Failed to generate haplotype object. Please check input data.")
    return(NULL)
  }
  
  if(verbose) message("Haplotyping with predefined marker groups complete!")
  
  return(HapObject)
}

#' Identify haplotypes from predefined marker groups
#'
#' This function processes marker groups and identifies haplotypes
#' based on allelic patterns.
#'
#' @param vcf Original VCF data
#' @param bin_vcf Binary VCF data
#' @param predefined_MGs Predefined marker groups
#' @param LD Linkage disequilibrium matrix
#' @param pheno Phenotype data
#' @param metadata Optional metadata
#' @param minHap Minimum individuals per haplotype
#' @param het_phenos Whether to include heterozygous phenotypes
#' @param keep_outliers Whether to keep outlier SNPs
#' @param verbose Whether to print progress messages
#'
#' @return A haplotype object for use with CrossHap visualization functions
#'
identify_haplotypes_from_MGs <- function(vcf, bin_vcf, predefined_MGs, LD, pheno, metadata = NULL,
                                         minHap = 9, het_phenos = FALSE, keep_outliers = FALSE,
                                         verbose = TRUE) {
  
  # Get unique marker groups (excluding "0")
  unique_MGs <- unique(predefined_MGs$MGs)
  unique_MGs <- unique_MGs[unique_MGs != "0"]
  
  if(verbose) message("Processing ", length(unique_MGs), " marker groups")
  
  # Initialize empty pseudoSNP matrix (individuals x marker groups)
  pseudoSNP <- matrix(NA, nrow = ncol(bin_vcf), ncol = length(unique_MGs))
  rownames(pseudoSNP) <- colnames(bin_vcf)
  colnames(pseudoSNP) <- unique_MGs
  
  # For each marker group, calculate allelic state for each individual
  for(i in seq_along(unique_MGs)) {
    mg <- unique_MGs[i]
    
    if(verbose) message("  Processing marker group ", mg)
    
    # Get SNPs in this marker group
    mg_snps <- predefined_MGs$ID[predefined_MGs$MGs == mg]
    
    if(length(mg_snps) == 0) {
      if(verbose) message("    No SNPs found for marker group ", mg)
      next
    }
    
    # Filter to SNPs present in bin_vcf
    mg_snps <- mg_snps[mg_snps %in% rownames(bin_vcf)]
    
    if(length(mg_snps) == 0) {
      if(verbose) message("    No matching SNPs found in VCF for marker group ", mg)
      next
    }
    
    if(verbose) message("    Found ", length(mg_snps), " SNPs in VCF")
    
    # Extract these SNPs from binary VCF
    mg_vcf <- bin_vcf[mg_snps, , drop = FALSE]
    
    # Calculate most common allele state for each individual
    for(j in 1:ncol(mg_vcf)) {
      ind_name <- colnames(mg_vcf)[j]
      ind_values <- mg_vcf[, j]
      ind_values <- ind_values[!is.na(ind_values)]
      
      if(length(ind_values) == 0) {
        # All values are NA for this individual at these SNPs
        pseudoSNP[ind_name, i] <- NA
      } else {
        # Find most common value
        tab <- table(ind_values)
        pseudoSNP[ind_name, i] <- as.numeric(names(tab)[which.max(tab)])
      }
    }
  }
  
  if(verbose) {
    message("PseudoSNP matrix created: ", nrow(pseudoSNP), " individuals x ", ncol(pseudoSNP), " marker groups")
    # Count NA values
    na_count <- sum(is.na(pseudoSNP))
    if(na_count > 0) {
      message("Warning: ", na_count, " missing values in pseudoSNP matrix (", 
              round(100 * na_count / (nrow(pseudoSNP) * ncol(pseudoSNP)), 1), "%)")
    }
  }
  
  # Convert to data frame and add individual IDs
  pseudoSNP_df <- as.data.frame(pseudoSNP)
  pseudoSNP_df$Ind <- rownames(pseudoSNP)
  
  # Create haplotype identifier string for each individual
  hap_strings <- apply(pseudoSNP, 1, function(x) paste(x, collapse = "_"))
  
  # Count frequency of each haplotype
  hap_table <- table(hap_strings)
  
  if(verbose) {
    message("Found ", length(hap_table), " unique haplotype combinations")
    message("Top 5 haplotypes by frequency:")
    print(sort(hap_table, decreasing = TRUE)[1:min(5, length(hap_table))])
  }
  
  # Create haplotype frequency data frame
  hap_counts <- data.frame(
    haplotype = names(hap_table),
    n = as.numeric(hap_table),
    stringsAsFactors = FALSE
  )
  
  # Remove NA haplotypes (individuals with missing data)
  hap_counts <- hap_counts[!grepl("NA", hap_counts$haplotype), ]
  
  # Keep only haplotypes with frequency >= minHap
  common_haps <- hap_counts$haplotype[hap_counts$n >= minHap]
  
  if(length(common_haps) < 2) {
    if(verbose) {
      message("Warning: Only ", length(common_haps), " haplotypes with frequency >= ", minHap)
      message("Consider lowering minHap parameter (current value: ", minHap, ")")
    }
    
    if(nrow(hap_counts) >= 2) {
      # Use top 2 haplotypes regardless of frequency
      common_haps <- hap_counts$haplotype[order(hap_counts$n, decreasing = TRUE)][1:min(2, nrow(hap_counts))]
      if(verbose) message("Using top ", length(common_haps), " haplotypes regardless of frequency")
    } else {
      if(verbose) message("ERROR: Cannot proceed with fewer than 2 unique haplotypes")
      return(NULL)
    }
  }
  
  if(verbose) message("Selected ", length(common_haps), " haplotypes with frequency >= ", minHap)
  
  # Create Hapfile
  Hapfile <- data.frame(
    hap = LETTERS[1:length(common_haps)],
    n = hap_counts$n[match(common_haps, hap_counts$haplotype)],
    stringsAsFactors = FALSE
  )
  
  # Add columns for each marker group
  for(i in seq_along(unique_MGs)) {
    mg_col <- paste0("MG", i)
    Hapfile[[mg_col]] <- sapply(common_haps, function(h) {
      parts <- strsplit(as.character(h), "_")[[1]]
      if(length(parts) >= i) {
        as.numeric(parts[i])
      } else {
        NA
      }
    })
  }
  
  if(verbose) {
    message("Hapfile created with ", nrow(Hapfile), " haplotypes")
    print(Hapfile)
  }
  
  # Create mapping from haplotype string to letter code
  hap_map <- data.frame(
    haplotype = common_haps,
    hap = LETTERS[1:length(common_haps)],
    stringsAsFactors = FALSE
  )
  
  # Create nophenIndfile with haplotype assignments
  pseudoSNP_df$haplotype <- hap_strings
  
  nophenIndfile <- merge(
    pseudoSNP_df[, c("Ind", "haplotype")],
    hap_map,
    by = "haplotype",
    all.x = TRUE
  )
  
  nophenIndfile <- nophenIndfile[, c("Ind", "hap")]
  nophenIndfile$hap[is.na(nophenIndfile$hap)] <- "0"
  
  if(verbose) {
    message("nophenIndfile created with ", nrow(nophenIndfile), " individuals")
    # Count frequency of each haplotype letter
    hap_letters <- table(nophenIndfile$hap)
    message("Haplotype frequencies in individuals:")
    print(hap_letters)
  }
  
  # Calculate mean r2 for each SNP within its marker group
  r2file <- data.frame(ID = character(0), meanr2 = numeric(0), MGs = character(0), stringsAsFactors = FALSE)
  
  for(mg in unique_MGs) {
    mg_SNPs <- predefined_MGs$ID[predefined_MGs$MGs == mg]
    mg_SNPs <- mg_SNPs[mg_SNPs %in% rownames(LD)]
    
    if(length(mg_SNPs) > 1) {  # Need at least 2 SNPs to calculate correlation
      mg_ld <- LD[mg_SNPs, mg_SNPs, drop = FALSE]
      
      for(snp in mg_SNPs) {
        if(snp %in% rownames(mg_ld)) {
          r2vals <- mg_ld[snp, ]
          meanr2 <- mean(r2vals[names(r2vals) != snp], na.rm = TRUE)
          r2file <- rbind(r2file, data.frame(ID = snp, meanr2 = meanr2, MGs = mg, stringsAsFactors = FALSE))
        }
      }
    } else if(length(mg_SNPs) == 1) {
      # For single SNP marker groups, set meanr2 to 1.0
      r2file <- rbind(r2file, data.frame(ID = mg_SNPs, meanr2 = 1.0, MGs = mg, stringsAsFactors = FALSE))
    }
  }
  
  # Create MGfile with r2 values ## something's wrong here or above the R2 doesnt come in the fig
  MGfile <- merge(predefined_MGs, r2file, by = c("ID", "MGs"), all.x = TRUE)
  
  # Handle outliers if keep_outliers is FALSE
  if(!keep_outliers) {
    for(mg in unique_MGs) {
      mg_rows <- MGfile$MGs == mg
      if(sum(mg_rows) > 2) {  # Need at least 3 SNPs to detect outliers
        mg_r2 <- MGfile$meanr2[mg_rows]
        if(!all(is.na(mg_r2))) {
          median_r2 <- median(mg_r2, na.rm = TRUE)
          sd_r2 <- sd(mg_r2, na.rm = TRUE)
          
          # Identify outliers (> 2 SD from median)
          outliers <- abs(mg_r2 - median_r2) > 2 * sd_r2
          outliers[is.na(outliers)] <- FALSE
          
          if(any(outliers)) {
            if(verbose) message("Removing ", sum(outliers), " outlier SNPs from marker group ", mg)
            # Set outlier MGs to "0"
            MGfile$MGs[mg_rows][outliers] <- "0"
          }
        }
      }
    }
  }
  
  # Run tagphenos to generate Varfile
  Varfile <- tagphenos(MGfile = MGfile, bin_vcf, pheno, het_phenos = het_phenos)
  
  # Get MGmin (minimum SNPs in any marker group)
  mg_counts <- table(MGfile$MGs[MGfile$MGs != "0"])
  MGmin <- min(as.numeric(mg_counts))
  
  # Create final haplotype object
  hap_obj <- list(
    MGmin = MGmin,
    Hapfile = Hapfile,
    Indfile = merge(nophenIndfile, pheno, by = "Ind", all.x = TRUE),
    Varfile = Varfile
  )
  
  # Add metadata if provided
  if(!is.null(metadata)) {
    hap_obj$Indfile <- merge(hap_obj$Indfile, metadata, by = "Ind", all.x = TRUE)
  } else {
    hap_obj$Indfile$Metadata <- NA
  }
  
  # Return as HapObject
  return(list(Haplotypes_Predefined = hap_obj))
}

# load my inputs
# input vcf

vcf_13 <- read_vcf ("region_370kb_TAN_DR0.8_chr13_29969164.vcf")
head (vcf_13)

# update the ID column
vcf_13 <- vcf_13 %>%
  mutate (ID = paste0("Gm13_", POS))
head(vcf_13)

# input LD matrix
LD_13 <-  read_LD ("chr13_TAN_0.8.ld", vcf = vcf_13)
head(LD_13)

# load TAN phenotype data
pheno_13 <- read_pheno("phenotype_tan_crosshap.txt")
head(pheno_13)

# load TAN metadata
metadata_13 <- read_metadata("metadata_tan_crosshap.txt")

# load predefined marker groups obtained with RB haplotyping in the same genomic interval
predefined_mgs <- read.csv("MGs_from_RB.csv")

# run the Analysis
hap_obj <- run_haplotyping_predefined(
  vcf = vcf_13, 
  LD = LD_13, 
  pheno = pheno_13, 
  metadata = metadata_13,
  predefined_MGs = predefined_mgs,
  minHap = 9,
  prefix_MGs = F,
  verbose = TRUE
)

# check str of hap objected created, and create a compatible obj for crosshap_viz()
str(hap_obj)
str(hap_obj)
compatible_obj <- list()

# copy existing data with the correct naming convention
compatible_obj[["Haplotypes_MGmin9_E1"]] <- hap_obj[["Haplotypes_Predefined"]]

# add the epsilon field that crosshap_viz expects
## not sure what to do with the epsilon, kept it as 1.0 so i can run crosshap_viz
compatible_obj[["Haplotypes_MGmin9_E1"]]$epsilon <- 1.0

# try visualization with the compatible object
hap_viz_13 <- crosshap_viz(HapObject = compatible_obj, epsilon = 1.0)
hap_viz_13

# to generate final plot removing the summary table to fit
# Build individual plots
pheno_13_dot <- build_mid_dotplot(HapObject = compatible_obj, epsilon = 1.0, hide_labels = F)
right_13_pheno <- build_right_phenoplot(HapObject = compatible_obj, epsilon = 1.0, hide_labels = F)
left_13_pheno <- build_left_posplot(HapObject = compatible_obj, epsilon = 1.0, hide_labels = F)
top_13 <- build_top_metaplot(HapObject = compatible_obj, epsilon = 1.0, hide_labels = F)
bot_13 <- build_bot_halfeyeplot(HapObject = compatible_obj, epsilon = 1.0, hide_labels = T)

# Build summary tables
table_13 <- build_summary_tables(compatible_obj, epsilon = 1.0)
MGtable <- table_13[[1]]

# Required Libraries
library(ggplot2)
library(patchwork)
library(grid)
library(ggplotify)
library(gridExtra)

# Convert tables to ggplot objects
MG_ggplot <- as.ggplot(MGtable) +
  theme_void() +  # Remove titles or axis labels
  theme(plot.margin = margin(t = 0, r = 10, b = 0, l = 0))  # Add right margin

# Construct the final plot layout
final_plot10_3 <- (
  plot_spacer() +                    # Spacer (optional)
    top_13 +                           # Top 13 plot
    wrap_elements(MG_ggplot) +         # MGtable
    left_13_pheno +                    # Left pheno plot
    pheno_13_dot +                     # Dot plot
    right_13_pheno +                   # Right pheno plot
    guide_area() +                     # Guide area
    bot_13 +                           # Bottom 20 plot
    plot_layout(ncol = 3, guides = 'collect') +  # Keep the grid size
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &  # Custom tags
    theme(plot.tag = element_text(face = 'bold', size = 15))  # Define grid layout
)

ggsave(plot = final_plot10_3 ,"crosshap_chr13_TAN_basedon_RB.jpg", width = 12, height = 10, dpi = 300, units = "in" )                     
