##------------------------------------------------------------------------------#
# GAinS-bulk-RNAseq-deconvolution-results
#------------------------------------------------------------------------------#
library(data.table)
library(tidyverse)
library(ggpubr) 
library(rstatix)
library(openxlsx)  
library(ggrepel)

#------------------------------------------------------------------------------#
# File path
#------------------------------------------------------------------------------#
cibersortx_results <- "./output/CIBERSORTx_Results/CIBERSORTx_Job12_Adjusted.txt"
sample_key_file <- "~/MTOR.project/GAinS/Sample_key.txt"
srs_data_file <- "~/MTOR.project/GAinS/originalSRS_2023_sorted.csv"
age_sex_file <- "~/MTOR.project/GAinS/Demo_2023-10-11_18-53-23.xlsx"

#------------------------------------------------------------------------------#
# load data
#------------------------------------------------------------------------------#
deconv <- fread(cibersortx_results)
deconv[, SangerSampleID := Mixture] 

sample_key <- fread(sample_key_file)
srs_data <- fread(srs_data_file)
meta_data <- merge(sample_key, srs_data, by.x = "GAinSID", by.y = "SampleID")

age_sex <- read.xlsx(age_sex_file, sheet = 1)
age_sex$Age = as.numeric(age_sex$Calculated.Age)

# Merge metadata
meta_data <- merge(meta_data, age_sex, by = "Subject.Id")
deconv_merged <- merge(deconv, meta_data, by = "SangerSampleID")

# clean cell subset names
deconv_3 = as.data.frame(deconv_merged)
colClean1 <- function(x){ colnames(x) <- gsub("[+]", "", colnames(x)); x } 
deconv_3 = colClean1(deconv_3)
colClean1 <- function(x){ colnames(x) <- gsub(" ", "", colnames(x)); x } 
deconv_3 = colClean1(deconv_3)
colClean1 <- function(x){ colnames(x) <- gsub("-", "_", colnames(x)); x } 
deconv_3 = colClean1(deconv_3)
colClean1 <- function(x){ colnames(x) <- gsub("/", "_", colnames(x)); x } 
deconv_3 = colClean1(deconv_3)

#------------------------------------------------------------------------------#
# Linear regression analysis
#------------------------------------------------------------------------------#
cell_types <- colnames(deconv_3)[4:26]  # Select cell type columns
results_list <- lapply(cell_types, function(cell) {
  model <- lm(as.formula(paste(cell, "~ SRS + Age + Sex")), data = deconv_3)
  coef_summary <- coef(summary(model))
  
  data.frame(
    cell_type = cell,
    beta = coef_summary[2, 1],  # SRS.x coefficient
    p_value = coef_summary[2, 4] # SRS.x p-value
  )
})

# Compile results
results_df <- rbindlist(results_list)
results_df[, padj := p.adjust(p_value, method = "bonferroni")]
results_df[, sig := fifelse(padj < 0.05 & beta > 0, "sig.up",
                            fifelse(padj < 0.05 & beta < 0, "sig.down", "non-sig"))]

#------------------------------------------------------------------------------#
# Volcano plot
#------------------------------------------------------------------------------#
ggplot(results_df, aes(x = beta, y = -log10(padj), color = sig)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "longdash", color = "grey") +
  geom_vline(xintercept = 0, linetype = "longdash", color = "grey") +
  scale_color_manual(values = c("grey", "darkorange", "cyan4")) +
  xlab("Beta") + ylab("-Log10 Adjusted P-value") +
  theme_classic() +
  xlim(-200, 250) +
  geom_label_repel(data = results_df[padj < 0.01], aes(label = cell_type),
                   size = 3, point.padding = unit(0.2, "lines"))
# ggsave("volcano_plot.png", width=8, height=6)






