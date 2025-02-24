#-------------------------------------------------------------------------------
# Get Command-Line Arguments
#-------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript plot.results.R <coloc.results.csv>")
}

coloc_results_csv <- args[1]

#-------------------------------------------------------------------------------
# Load required libraries
#-------------------------------------------------------------------------------
library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)
library(seqminer)

#-------------------------------------------------------------------------------
# Load coloc results
#-------------------------------------------------------------------------------
res <- read.csv(coloc_results_csv)

# Dot plot
ggplot(res, aes(x = PP.H3.abf, y = PP.H4.abf)) +
  geom_point(size = 2, shape = 21, fill = "darkgrey") +
  geom_point(data = subset(res, PP.H4.abf >= 0.95), size = 1.5, shape = 21, fill = "darkorange") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey") +
  geom_text_repel(data = subset(res, PP.H4.abf >= 0.95),
                  aes(label = dataset), size = 3, 
                  point.padding = unit(0.2, "lines"), 
                  max.overlaps = 100, min.segment.length = 0) +
  xlim(0, 1) + ylim(0, 1) +
  xlab("PP.H3 (Different Causal Variants)") + 
  ylab("PP.H4 (Shared Causal Variants)") +
  theme_bw() + theme(panel.grid.minor = element_blank())

# Filter datasets with PP.H4.abf > 0.95
res_0.9 <- filter(res, PP.H4.abf > 0.95)

#-------------------------------------------------------------------------------
# Define function to import eQTL Catalogue data
#-------------------------------------------------------------------------------
import_eQTLCatalogue <- function(ftp_path, region, selected_gene_id, column_names, verbose = TRUE) {
  if (verbose) {
    print(ftp_path)
  }
  
  # Fetch summary statistics with seqminer
  fetch_table <- seqminer::tabix.read.table(tabixFile = ftp_path, 
                                            tabixRange = region, 
                                            stringsAsFactors = FALSE) %>%
    as_tibble()
  
  colnames(fetch_table) <- column_names
  
  # Remove duplicate/multi-allelic variants
  fetch_table %>%
    filter(gene_id == selected_gene_id) %>%
    mutate(id = paste(chromosome, position, sep = ":")) %>%
    group_by(id) %>%
    mutate(row_count = n()) %>%
    ungroup() %>%
    filter(row_count == 1) 
}

#-------------------------------------------------------------------------------
# Load eQTL Catalogue tabix paths
#-------------------------------------------------------------------------------
tabix_paths <- read.delim("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv", 
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
  as_tibble()

# Fix FTP connection issue
tabix_paths$ftp_path <- gsub("ftp://ftp.ebi.ac.uk", "http://ftp.ebi.ac.uk", tabix_paths$ftp_path)

#-------------------------------------------------------------------------------
# Process datasets
#-------------------------------------------------------------------------------
tabix_paths$dataset <- paste0(tabix_paths$study_label, "_", tabix_paths$sample_group)

# Extract column names (fixed dataset for robustness)
column_names <- colnames(read_tsv("http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000010/QTD000083/QTD000083.all.tsv.gz", n_max = 1))

resul <- lapply(as.character(res_0.9$dataset), function(i) {
  platelet_df <- filter(tabix_paths, dataset == i, quant_method %in% c("ge", "microarray"))
  
  if (nrow(platelet_df) == 0) return(NULL)  # Skip if no matching dataset
  
  result <- import_eQTLCatalogue(platelet_df$ftp_path, 
                                 region = "1:11246222-11246222", 
                                 selected_gene_id = "ENSG00000198793", 
                                 column_names)
  
  result <- result %>%
    mutate(CIlow = beta - se * qnorm(0.975),
           CIup = beta + se * qnorm(0.975),
           dataset = i)
  
  return(result)
})

# Combine results into a single dataframe
resul_df <- do.call(rbind, resul) %>%
  as.data.frame() %>%
  mutate(Log10_p = -log10(pvalue)) %>%
  arrange(desc(Log10_p))

#-------------------------------------------------------------------------------
# Forest plot
#-------------------------------------------------------------------------------
ggplot(resul_df, aes(x = reorder(dataset, Log10_p), y = -beta, ymin = -CIlow, ymax = -CIup)) +
  geom_pointrange(aes(fill = Log10_p), shape = 21, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  xlab("") + ylab("Beta (95% CI) for rs4845987-G") +
  theme_bw() + theme(panel.grid.minor = element_blank()) +
  scale_fill_gradient2(low = "#FFE134", mid = "#FCCB06", high = "darkorange4", midpoint = 10, 
                       name = "-Log10 P value")
