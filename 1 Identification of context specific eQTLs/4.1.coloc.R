#-------------------------------------------------------------------------------
# Load libraries
#-------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(readr)
library(coloc)
library(GenomicRanges)
library(seqminer)
library(data.table)

#-------------------------------------------------------------------------------
# Load Local Tabix Paths
#-------------------------------------------------------------------------------
tabix_paths <- fread("./tabix_path_local.tsv", stringsAsFactors = FALSE) %>%
  as_tibble() 

tabix_paths$datasets <- paste0(tabix_paths$study, "_", tabix_paths$qtl_group)
tabix_paths <- filter(tabix_paths, quant_method %in% c("ge", "microarray")) ## 127 eQTL datasest from eQTL Catalogue

#-------------------------------------------------------------------------------
# Function to Import eQTL Catalogue Data
#-------------------------------------------------------------------------------
import_eQTLCatalogue <- function(ftp_path, region, selected_gene_id, column_names) {
  fetch_table <- seqminer::tabix.read.table(tabixFile = ftp_path, tabixRange = region, stringsAsFactors = FALSE) %>%
    as_tibble()
  colnames(fetch_table) <- column_names
  
  # Remove rsid duplicates and multi-allelic variants
  fetch_table %>%
    filter(gene_id == selected_gene_id) %>%
    mutate(id = paste(chromosome, position, sep = ":")) %>%
    group_by(id) %>%
    mutate(row_count = n()) %>%
    ungroup() %>%
    filter(row_count == 1) 
}

#-------------------------------------------------------------------------------
# Function to run Coloc 
#-------------------------------------------------------------------------------
run_coloc <- function(summary_stats, gwas_sumstats) {
  eQTL_dataset <- list(
    beta = summary_stats$beta,
    varbeta = summary_stats$se^2,
    N = (summary_stats$an)[1] / 2,  # Sample size is allele number (AN) divided by 2
    MAF = summary_stats$maf,
    type = "quant",
    snp = summary_stats$rsid
  )
  
  gwas_dataset <- list(
    beta = gwas_sumstats$Beta,
    varbeta = gwas_sumstats$SE^2,
    type = "quant",
    snp = gwas_sumstats$SNP,
    MAF = gwas_sumstats$MAF,
    N = gwas_sumstats$sample.size
  )
  
  coloc_res <- coloc::coloc.abf(dataset1 = eQTL_dataset, dataset2 = gwas_dataset, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5) # default setting
  as_tibble(t(as.data.frame(coloc_res$summary)))
}

#-------------------------------------------------------------------------------
# Run Coloc analysis for each dataset
#-------------------------------------------------------------------------------
results <- lapply(unique(tabix_paths$datasets), function(dataset) {
  platelet_df <- filter(tabix_paths, datasets == dataset)
  
  if (nrow(platelet_df) == 0) {
    message("Skipping dataset: ", dataset, " (not found)")
    return(NULL)
  }
  
  message("Processing dataset: ", dataset)
  
  # Extract column names from the first file
  column_names <- colnames(read_tsv(platelet_df$ftp_path, n_max = 1))
  
  # Read eQTL data
  eQTL_data <- seqminer::tabix.read.table(tabixFile = platelet_df$ftp_path,
                                          tabixRange = "1:11006407-11406407",
                                          stringsAsFactors = FALSE) %>%
    as_tibble()
  
  if (nrow(eQTL_data) == 0) {
    message("No eQTL data found for dataset: ", dataset)
    return(NULL)
  }
  
  # Import summary statistics
  summary_stats <- import_eQTLCatalogue(platelet_df$ftp_path,
                                        region = "1:11006407-11406407",
                                        selected_gene_id = "ENSG00000198793",
                                        column_names)
  
  if (nrow(summary_stats) == 0) {
    message("Skipping dataset: ", dataset, " (gene filtered out due to low expression)")
    return(NULL)
  }
  
  # Filter SNPs present in GWAS dataset
  summary_stats$rsid <- gsub("\r", "", summary_stats$rsid)
  summary_stats_filtered <- filter(summary_stats, rsid %in% Sepsis.coloc$SNP)
  
  if (nrow(summary_stats_filtered) == 0) {
    message("No matching SNPs found for dataset: ", dataset)
    return(NULL)
  }
  
  # Run coloc analysis
  run_coloc(summary_stats_filtered, Sepsis.coloc)
})

#-------------------------------------------------------------------------------
# Combine results 
#-------------------------------------------------------------------------------
results <- do.call(rbind, results)
results$dataset <- rownames(results)

write.csv(results, "coloc.res.sepsis.eQTL.csv", row.names = FALSE)

message("Results saved to: coloc.res.sepsis.eQTL.csv")
