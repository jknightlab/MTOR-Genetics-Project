# Load libraries
library(DESeq2)
library(data.table)
library(reshape)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(rstatix)

#-----------------------------------------------------------------------------------#
# Step 1: Load and preprocess count data
#-----------------------------------------------------------------------------------#

# Merge count files using awk command 
# awk '{print $0 "\t" FILENAME}' *.htseq_count.txt > T.cells.ATAC.htseq_count.merge_slurm.txt

# Load merged count data
count_data <- as.data.frame(fread("./output/T.cells.ATAC.htseq_count.merge_slurm.txt.gz")) 

# Remove chrX and chrY
count_data <- count_data[!count_data$V1 %like% "chrX" & !count_data$V1 %like% "chrY",]

# Reshape data
count_data$value <- count_data$V2
count_data_reshaped <- cast(count_data, V1 ~ V3)
rownames(count_data_reshaped) <- count_data_reshaped$V1
count_data_reshaped <- count_data_reshaped[,-1]

# Clean column names
clean_colnames <- function(x, pattern, replacement) {
  colnames(x) <- gsub(pattern, replacement, colnames(x))
  return(x)
}

count_data_clean <- clean_colnames(count_data_reshaped, "_merge_ATAC.17_21.htseq_count.txt", "")
count_data_clean <- data.frame(clean_colnames(count_data_clean, "[-]", "_"))

#-----------------------------------------------------------------------------------#
# Step 2: Load and preprocess metadata
#-----------------------------------------------------------------------------------#

# Load metadata
meta_data <- fread("./output/ATACseq.meta.csv", stringsAsFactors = F)
meta_data <- meta_data[meta_data$Cell.type %in% "CD4",]

# Match count data with metadata
count_data_clean <- count_data_clean[, which(names(count_data_clean) %in% meta_data$Sample.ID)]
meta_data_matched <- meta_data[match(colnames(count_data_clean), meta_data$Sample.ID),]
all(colnames(count_data_clean) == meta_data_matched$Sample.ID)  # TRUE
colnames(count_data_clean)
#-----------------------------------------------------------------------------------#
# Step 3: DESeq2 analysis
#-----------------------------------------------------------------------------------#

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_data_clean, colData = meta_data_matched, design = ~Time.point)
dds <- DESeq(dds)

#-----------------------------------------------------------------------------------#
# Step 4: PCA visualisation
#-----------------------------------------------------------------------------------#

# Perform variance stabilizing transformation
vsd <- vst(dds, blind=TRUE)

# Plot PCA
pca_data <- plotPCA(vsd, intgroup=c("Cell.type", "Time.point", "Sample.ID"), returnData=TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

# Custom PCA plot
ggplot(pca_data, aes(PC1, PC2, col=Time.point)) +
  geom_point(size=2.5) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_bw() +
  scale_colour_manual(values = c("grey", "red")) +
  scale_shape_manual(values=c(1, 10)) +
  ggtitle("ATAC-seq") +
  theme(panel.grid.minor = element_blank()) +
  geom_text_repel(data=pca_data, aes(label = Sample.ID),
                  size = 1.5, point.padding = unit(0.2, "lines"), min.segment.length = 0)

#-----------------------------------------------------------------------------------#
# Step 5: Differential expression analysis
#-----------------------------------------------------------------------------------#

# Extract results
resultsNames(dds) # 
res <- data.frame(results(dds, name = "Time.point_Day.4_vs_Day.0"))
res$id <- rownames(res)

# Filter significant results
res_sig <- res[res$padj < 0.05 & !is.na(res$padj) & abs(res$log2FoldChange) > 1,]

#-----------------------------------------------------------------------------------#
# Step 6: Visualisation of specific peaks
#-----------------------------------------------------------------------------------#

# Example peak ID
peak_id <- "chr1:11245527-11246896_11245736_RPL39P6_13248_MTOR_16296"
padj_value <- signif(res[res$id == peak_id,]$padj, digits = 3)

# Plot counts for the specific peak
plot_data <- plotCounts(dds, gene=peak_id, intgroup=c("Time.point", "Biological.replicate", "Cell.type"), returnData=TRUE)
plot_data$gene <- "MTOR"
plot_data$value <- log2(plot_data$count + 1/2)

# Create paired plot
p <- ggpaired(plot_data, x = "Time.point", y = "value",
              color = "Time.point", line.color = "gray", line.size = 0.4,
              palette = "jco", add = "jitter", facet.by = "Cell.type",
              xlab = "", ylab = "Log2 normalised counts")

# Add p-value annotation
p + annotate("text", x=2, y=5, label=paste0("Padj = \n", padj_value)) +
  scale_color_manual(values=c("darkgrey", "darkorange"))

# Save plot
ggsave("cd4.pdf", width = 1.65, height = 2.6)