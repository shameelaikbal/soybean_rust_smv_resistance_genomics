#syteny analysis with common bean using syntenet R-package
#documentation here:

library(syntenet)
library(stringr)

#path to directory containing the fasta files
fasta_dir <- "fasta/"

#see contents of the directory
dir(fasta_dir)

#path to directory containing the gff file 
gff_dir <- "gff/"

#see contents of  gff directory
dir(gff_dir)

# read all fasta files in fasta directory
proteomes <- fasta2AAStringSetlist(fasta_dir)
proteomes

# read all gff files in  gff dir
annotation <- gff2GRangesList(gff_dir)
annotation

# checking if data matches all their 3 criteria
check_input(proteomes, annotation, gene_field = "Name")



# so to keep only one isoform the longest for each gene
# Load required packages
library(syntenet)
library(Biostrings)
library(stringr)
library(dplyr)

# Access Pvul proteins and annotation
pvul_proteins <- proteomes$Pvul
head(pvul_proteins)
pvul_annotation <- annotation$Pvul
head(pvul_annotation)

# Get gene IDs from your annotation
annotation_gene_ids <- pvul_annotation$gene_id
length(annotation_gene_ids)
head(annotation_gene_ids)

# Extract gene IDs from protein headers
protein_gene_ids <- sapply(strsplit(names(pvul_proteins), "[.]"), function(x) x[1])
protein_gene_ids
length(protein_gene_ids)
# Find proteins that match  annotation genes
matching_indices <- which(protein_gene_ids %in% annotation_gene_ids)
length(matching_indices)

# If multiple matches per gene, keep only the first one
unique_gene_indices <- match(unique(protein_gene_ids[matching_indices]), protein_gene_ids[matching_indices])
final_indices <- matching_indices[unique_gene_indices]
length(final_indices)

# Select these proteins
filtered_proteins <- pvul_proteins[final_indices]

# Save the filtered proteins
dir.create("fasta_filtered", showWarnings = FALSE)
writeXStringSet(filtered_proteins, "fasta_filtered/Pvul.fasta")

# Access soybean data
gmax_proteins <- proteomes$Gmax

# Extract just the locus/gene ID from the protein names
#  we want to extract "Glyma.13G185100" which follows "locus="
new_names <- sapply(strsplit(names(gmax_proteins), "locus="), function(x) {
  if(length(x) > 1) {
    # Extract the locus name (gene ID)
    gene_id <- strsplit(x[2], " ")[[1]][1]
    return(gene_id)
  } else {
    # If format is different, return original
    return(x[1])
  }
})

# Apply the new names
names(gmax_proteins) <- new_names

# Save the proteins with simplified names
writeXStringSet(gmax_proteins, "fasta_filtered/Gmax.fasta")

# Check the structure of annotation object
str(annotation$Gmax, max.level=1)

# Look at available metadata columns
colnames(mcols(annotation$Gmax))

# Check a few sample entries from different fields
head(annotation$Gmax$gene_id, 3)  # Default field syntenet looks for
head(annotation$Gmax$ID, 3)       # Another common field
head(annotation$Gmax$Name, 3)     # Name field you mentioned

#now process the input
pdata <- process_input(proteomes, annotation, gene_field = "Name")

#perform all vs blast results
blast_results <- run_diamond(seq = pdata$seq)

# Infer synteny network (remember to use gene_field="Name")
synteny_network <- infer_syntenet(blast_results, pdata$annotation)

head(synteny_network)

# Get species ID mapping table
id_table <- create_species_id_table(names(proteomes))
id_table  # This shows how species names map to the short IDs used in the network

# Filter for synteny relationships involving Glycine max
gmax_synteny <- synteny_network[
  grepl("^Gma_", synteny_network$Anchor1) | grepl("^Gma_", synteny_network$Anchor2),
]

# View the first few rows
head(gmax_synteny)
gmax_synteny

