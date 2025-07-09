#Linkage disequilibrium in local genomic regions around GWAS-SNPs and plotting

# calculate LD using plink
./plink --vcf region_1Mb_from56098721_RB_chr18_56100116.vcf --allow-extra-chr --r2 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --out 1mb_ld

#then grep for this significant SNP in the above output to get its LD with all other SNPs for the region considered, can do for a 1Mb window
grep '56100116' ld.ld > gwas_ld.txt

# then plot it in R
#Load required libraries
library(ggplot2)
library(dplyr)

# Load LD data
ld <- read.table("gwas_1mb_full_ld.txt", header = FALSE)
colnames(ld) <- c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2")

# Define lead SNP
lead_pos <- 56100116

# Manually add GWAS SNP itself
gwas_point <- data.frame(
  CHR_A = "Gm18",
  BP_A = lead_pos,
  SNP_A = "ss715632290",
  CHR_B = "Gm18",
  BP_B = lead_pos,
  SNP_B = "ss715632290",
  R2 = 1
)

# Combine
ld_all <- rbind(ld, gwas_point)

# Calculate distance from lead SNP in kb
ld_all$distance_kb <- (ld_all$BP_B - lead_pos) / 1000

# Create more detailed LD categories
ld_all$LDcategory <- cut(ld_all$R2, 
                         breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
                         labels = c("Very Low", "Low", "Moderate", "High", "Very High"),
                         right = TRUE)

# Create color palette for LD categories
ld_colors <- c("Very Low" = "gray30", 
               "Low" = "gray50", 
               "Moderate" = "#6BAED6", 
               "High" = "#4292C6", 
               "Very High" = "#08519C")

# Create plot with improved aesthetics
ggplot(ld_all, aes(x = distance_kb, y = R2)) +
  # Add reference lines
  geom_hline(yintercept = c(0.2, 0.4, 0.6, 0.8), linetype = "dashed", color = "gray80", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  # Plot points with categorical colors
  geom_point(aes(color = LDcategory), size = 1.8, alpha = 0.8) +
  # Adjust colors
  scale_color_manual(values = ld_colors, name = expression(LD~(R^2))) +
  # Adjust axes
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(breaks = seq(-1000, 1000, 250)) +
  # Add labels
  labs(title = expression(paste("Linkage Disequilibrium (", R^2, ") with GWAS SNP (Chr18:56100116)")),
       x = "Distance from GWAS SNP (kb)",
       y = expression(Linkage~Disequilibrium~(R^2))) +
  # Improve theme
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    panel.border = element_rect(fill = NA, color = "gray50", linewidth = 0.5)
  )

# Save the plot at high resolution
ggsave("GWAS_SNP_LD_plot.jpg", width = 8, height = 6, dpi = 300)

# To save as PDF instead (often preferred for publications)
ggsave("GWAS_SNP_LD_plot.pdf", width = 8, height = 6)
