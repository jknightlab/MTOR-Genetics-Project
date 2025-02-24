
# module load R/4.3.2-gfbf-2023a
# Rscript ~/MTOR.project/smr.analysis/1.smr.file.prep.R Schmiedel_2018 CD8_T-cell_anti-CD3-CD28

#-------------------------------------------------------------------------------
# Load Required Libraries
#-------------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(seqminer))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(readr))

#-------------------------------------------------------------------------------
# Get Command-Line Arguments
#-------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript run_smr_analysis.R <eQTL_study_label> <eQTL_sample_group> ")
}
# Assign input file and output directory from command-line arguments
study_label <- args[1]
sample_group <- args[2]

#-------------------------------------------------------------------------------
# Define Input Files and Paths
#-------------------------------------------------------------------------------
# GWAS Summary Statistics for T2D Cohort 2024
T2D_GWAS_FILE <- "~/T2D.GWAS/EUR_Metal_LDSC-CORR_Neff.v2.txt.gz"
# 1000 Genomes European Reference Panel (Allele Frequencies)
# /apps/well/plink/2.00a-20170724/plink2 -bfile /well/jknight/users/kwz374/GWAS.snps/1kg_EUR/1000genomes_EUR_chr1  --freq cols=+pos
KG_EUR_AFREQ_FILE <- "~/EUR_chr1_plink2.afreq.gz"
# Tabix file paths from eQTL Catalogue
TABIX_PATHS_FILE <- "https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv"
#TABIX_PATHS_FILE <- "./tabix_paths.csv"
# Chain file for liftover from hg38 to hg19
LIFTOVER_CHAIN_FILE <- "~/hg38ToHg19.over.chain"
# Gene annotation file (GENCODE v31, hg19)
GENE_ANNOTATION_FILE <- "~/gencode.v31lift37_gene.txt"
# Feature ID list (gene names)
OUTPUT_DIR <- paste0("esd.",args[2])
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
FEATURE_IDS_FILE <- file.path(OUTPUT_DIR, "feature_ids.txt")
# Output GWAS summary statistics (formatted for SMR)
GWAS_OUTPUT_FILE <- file.path(".", "T2D_cohort2024_subset_gwas_sum_ma")



#-------------------------------------------------------------------------------
# Read 1000 Genomes rsID Mapping
#-------------------------------------------------------------------------------
rsid <- fread(KG_EUR_AFREQ_FILE, stringsAsFactors = FALSE)
rsid <- rsid %>% filter(nchar(REF) == 1 & nchar(ALT) == 1)  # Keep only biallelic SNPs
rsid$id <- paste0(rsid$`#CHROM`, "_", rsid$POS)



#-------------------------------------------------------------------------------
#'@-1 Load GWAS Data
#-------------------------------------------------------------------------------
cat("Loading GWAS summary data...\n")
T2D_cohort2024 <- as.data.frame(fread(T2D_GWAS_FILE, stringsAsFactors = FALSE))
T2D_cohort2024$Position <- as.numeric(as.character(T2D_cohort2024$Position))

# 1:11246222 (GRCh38).  1:11306279 (GRCh37). +- 2Mbp chr1:9306281-13351841
T2D_cohort2024.subset <- T2D_cohort2024 %>%
  filter(Chromsome == "1" & Position >= 11306279 - 2000000 & Position <= 11306279 + 2000000) %>%
  mutate(id = paste0(Chromsome, "_", Position))

# Merge with GWAS data to 1KG data get rsID
T2D_cohort2024.subset_rsid <- merge(T2D_cohort2024.subset, rsid, by = "id")

#-------------------------------------------------------------------------------
# Process and Save GWAS Summary Statistics
#-------------------------------------------------------------------------------
output_GWAS <- T2D_cohort2024.subset_rsid[,c("ID","EffectAllele","NonEffectAllele","EAF","Beta","SE","Pval","Neff")]
colnames(output_GWAS) = c("SNP","A1", "A2","freq","b","se","p","n")

write.table(output_GWAS, GWAS_OUTPUT_FILE, row.names = FALSE, quote = FALSE, sep = "\t")
cat("GWAS summary saved to", GWAS_OUTPUT_FILE, "\n")

#-------------------------------------------------------------------------------
#'@-2 Load eQTL Dataset Metadata
#-------------------------------------------------------------------------------
tabix_paths <- read.delim(TABIX_PATHS_FILE, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% as_tibble()
tabix_paths$ftp_path <- gsub("ftp://ftp.ebi.ac.uk", "http://ftp.ebi.ac.uk", tabix_paths$ftp_path)
eqtl_dataset <- filter(tabix_paths, study_label == as.character(args[1]), sample_group == as.character(args[2]), quant_method == "ge")

# tabix_paths = fread(TABIX_PATHS_FILE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()
# #tabix_paths = tabix_paths[tabix_paths$quant_method %in% c("ge","microarray"),] # 
# eqtl_dataset <- filter(tabix_paths, study == as.character(args[1]), qtl_group == as.character(args[2]), quant_method == "ge")

# Read eQTL Summary Statistics for regions. # # Rs4845987. 1:11246222 (GRCh38) +- 2Mbp.  1:11306279 (GRCh37)

cat("Fetch summary statistics with seqminer...\n")
eQTL_stats <- seqminer::tabix.read.table(tabixFile = eqtl_dataset$ftp_path, tabixRange = "1:9246222-13246222", stringsAsFactors = FALSE) %>%
  as_tibble()
#
colnames(eQTL_stats) <- colnames(read_tsv(eqtl_dataset$ftp_path, n_max = 1))

##
eQTL_stats$id = paste0(eQTL_stats$chromosome, "_", eQTL_stats$position)
eQTL_stats$rsid = gsub("\r","",eQTL_stats$rsid)

print(head(eQTL_stats))
#-------------------------------------------------------------------------------
# Convert hg38 Positions to hg19 Using LiftOver
#-------------------------------------------------------------------------------
cat("Converting hg38 positions to hg19...\n")
Hg38ToHg19chain <- import.chain(LIFTOVER_CHAIN_FILE)
#
hg38 <- GRanges(seqnames = paste0("chr", eQTL_stats$chromosome),
                ranges = IRanges(start = eQTL_stats$position,
                                 end = eQTL_stats$position),
                id = eQTL_stats$id)

HG19 = liftOver(hg38, Hg38ToHg19chain)
HG19 = as.data.frame(unlist(HG19))

HG19$id.hg19 = paste0(gsub("chr", "",HG19$seqnames), "_", HG19$start)
HG19$pos.hg19 = HG19$start
HG19 = HG19[!duplicated(HG19$id.hg19),c("id","id.hg19","pos.hg19")]


# Merge hg19 positions
eQTL_stats_hg19 <- merge(eQTL_stats, HG19, by = "id")

#-------------------------------------------------------------------------------
# Merge with 1000 Genomes rsID Data
#-------------------------------------------------------------------------------
eQTL_stats_hg19_2 <- merge(eQTL_stats_hg19, rsid, by.x = "rsid", by.y = "ID")
eQTL_stats_hg19_2$effect.freq <- ifelse(eQTL_stats_hg19_2$alt == eQTL_stats_hg19_2$ALT, eQTL_stats_hg19_2$maf, 1 - eQTL_stats_hg19_2$maf)

#-------------------------------------------------------------------------------
# Load Gene Names
#-------------------------------------------------------------------------------
cat("Adding gene names...\n")
geneloc <- fread(GENE_ANNOTATION_FILE)
geneloc <- geneloc %>%
  mutate(ENSG.id = gsub("[.].*", "", V4)) %>%
  filter(ENSG.id %in% eQTL_stats_hg19_2$gene_id) %>%
  select(ENSG.id, gene.name = V7)

eQTL_stats_hg19_2_gene.name <- merge(eQTL_stats_hg19_2, geneloc, by.x = "gene_id", by.y = "ENSG.id")

#-------------------------------------------------------------------------------
# Save Processed eQTL Data
#-------------------------------------------------------------------------------

output.eQTLs <- eQTL_stats_hg19_2_gene.name[,c("gene.name","chromosome","rsid","pos.hg19","alt","ref","effect.freq","beta","se","pvalue")]
colnames(output.eQTLs) = c("gene.name", "Chr", "SNP", "Bp","A1","A2","Freq", "Beta","se","p")
#
consistent_A1_A2 <- output.eQTLs %>%
  group_by(gene.name, SNP) %>%
  filter(n_distinct(A1) == 1 & n_distinct(A2) == 1) %>%
  ungroup()

#
split(consistent_A1_A2, consistent_A1_A2$gene.name) %>%
  lapply(function(df) {
    file_name <- file.path(OUTPUT_DIR, paste0(df$gene.name[1], ".esd"))
    df <- df[, colnames(df) != "gene.name", drop = FALSE]
    write.table(df, file = file_name, sep = "\t", row.names = FALSE, quote = FALSE)
  })
#
write.table(unique(consistent_A1_A2$gene.name), file = FEATURE_IDS_FILE, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("SMR analysis complete. Outputs saved in", OUTPUT_DIR, "\n")
