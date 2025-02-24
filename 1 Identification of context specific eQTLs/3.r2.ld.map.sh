#!/bin/bash

# Load necessary modules or set up environment if needed

# Step 1: Filter genotype data using PLINK
/apps/well/plink/1.90b3/plink --noweb --allow-extra-chr \
--bfile ./genotype/All_genotyping_merged_filtered_b38_refiltered_rsID \
--geno 0.02 -mind 0.02 --maf 0.01 \
--make-bed --out All_genotyping_merged_filtered_b38_refiltered_rsID_geno0.02_mind0.02

# Step 2: Get LD for SNPs in range
#/apps/well/plink/1.90b3/plink --bfile All_genotyping_merged_filtered_b38_refiltered_rsID_geno0.02_mind0.02 --chr 1 --allow-extra-chr --from-bp 11201872 --to-bp 11260293 --make-bed --out mtor.leads_for.LDheatmap
/apps/well/plink/1.90b3/plink --bfile All_genotyping_merged_filtered_b38_refiltered_rsID_geno0.02_mind0.02 --chr 1 --allow-extra-chr --from-bp 10262551 --to-bp 12262551 --make-bed --out mtor.TSS.1Mbp_for.LDheatmap

# Step 3: Calculate linkage disequilibrium (LD) (R²) for all SNPs in the region using MOTR lead as the reference
/apps/well/plink/1.90b3/plink --bfile mtor.region_snps_for.LDheatmap --ld-snp rs4845987 --r2 dprime inter-chr --ld-window-r2 0 --out ld_results

# Step 4: Generate LD matrix
/apps/well/plink/1.90b3/plink --bfile mtor.leads_for.LDheatmap \
--r2 square \
--out ld_matrix_leads

# Step 5: Generate LD heatmap using R
Rscript - <<EOF
library(LDheatmap)
library(grid)

# Read PLINK's LD matrix
ld_matrix <- as.matrix(read.table("~/MTOR.project/results/LD.heatmap/ld_matrix_leads.ld", header=FALSE))

# Read SNP positions from .bim file
bim_data <- read.table("~/MTOR.project/results/LD.heatmap/mtor.leads_for.LDheatmap.bim", header=FALSE)
snp_positions <- bim_data\$V4  # Column 4 = SNP positions
snp_names <- bim_data\$V2  # Column 2 = SNP IDs

# Assign SNP names as row and column names in the LD matrix
rownames(ld_matrix) <- snp_names
colnames(ld_matrix) <- snp_names

# Find the index of rs4845987 in the snp_names
reference_snp_index <- which(snp_names == "rs4845987")

# Identify SNPs with r² > 0.95 with rs4845987 (look at the row corresponding to rs4845987)
high_ld_indices <- which(ld_matrix[reference_snp_index, ] > 0.95)

# Extract SNP names corresponding to those with r² > 0.95
high_ld_snp_names <- snp_names[high_ld_indices]

# Create the heatmap using the custom color matrix, and specify SNP names to label
rgb.palette <- colorRampPalette(rev(c("black", "purple","yellow","green","blue","red")), space = "rgb")

LDheatmap(ld_matrix, 
          genetic.distances = snp_positions, 
          color = rgb.palette(24),
          SNP.name = high_ld_snp_names,  # Only label SNPs with r² > 0.95 with rs4845987
          name = "LDheatmap")
grid.edit("symbols", gp = gpar(cex = 1,col="cyan4"))
grid.edit(gPath("SNPnames"), gp=gpar(col="red"))

# Save the heatmap as a PNG file
png("TSS.1Mbp.png",width=10,height=10,units="in",res=400)
LDheatmap(ld_matrix, 
          genetic.distances = snp_positions, 
          color = rgb.palette(24),
          SNP.name = NULL,  # Only label SNPs with r² > 0.95 with rs4845987
          name = "LDheatmap", add.map=FALSE, add.key=FALSE)
dev.off()
EOF

