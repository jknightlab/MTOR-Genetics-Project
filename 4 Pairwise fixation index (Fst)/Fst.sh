#!/bin/bash

# Load required modules
module load BCFtools/1.10.2-GCC-8.3.0

# Run PLINK to generate VCF file
/apps/well/plink/1.90b3/plink --noweb \
--bfile ~/PCA_1KG/ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.genotypes \
--recode vcf --out temp

# Filter VCF file
bcftools view --threads 24 temp.vcf --include ID=@MTOR.leads.txt -Ov > MTOR.25.leads.in.1KG.data.vcf  
rm temp*

# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped 

  # Run the R script
  Rscript -e "library(hierfstat); \
  library(vcfR); \
  library(adegenet); \
  library(pheatmap); \
  library(data.table); \
  library(dplyr); \
  library(ggplot2); \
  
  # Read VCF file
  vcf <- fread('~/MTOR.25.leads.in.1KG.data.vcf'); \
  snp = vcf[vcf$ID == '1:11246222:G:C',]; \
  
  xx1 = as.data.frame(t(snp[,c(3, 10:2557)])); \
  xx1$Individual_ID = gsub('.*[_]','',rownames(xx1)); \
  xx1 = as.data.frame(xx1[-1,]); \
  colnames(xx1) = c('genotype','Individual_ID'); \
  xx1$genotype = ifelse(xx1$genotype == '0/0', 'C/C', \
                        ifelse(xx1$genotype == '0/1', 'C/G', \
                               ifelse(xx1$genotype == '1/1', 'G/G', NA))); \
  
  # Read metadata
  PED <- fread('~/20130606_g1k.ped', stringsAsFactors = FALSE); \
  PED$Individual_ID = PED$`Individual ID`; \
  PED$IndividualID = paste0('0_', PED$Individual_ID); \
  
  population.code = fread('~/igsr_populations.tsv'); \
  population.code2 = population.code[population.code$`Population code` %in%  PED$Population,]; \
  population.code2 = population.code2[,c('Population code', 'Superpopulation code')]; \
  
  PED = merge(PED, population.code2, by.x='Population', by.y='Population code'); \
  xx2 = merge(xx1, PED, by='Individual_ID'); \
  xx2$group = xx2$`Superpopulation code`; \
  
  # Count occurrences
  df_counts <- xx2 %>% group_by(group, genotype) %>% summarise(count = n()) %>% ungroup(); \
  df_totals <- df_counts %>% group_by(group) %>% summarise(total = sum(count)); \
  
  # Generate plot
  ggplot(df_counts, aes(x = count, y = group , fill = genotype)) + \
    geom_bar(stat = 'identity', position = 'fill')  + \
    xlab('Fraction of allele carriers') + ylab('Population') + \
    theme_minimal(); \
  
  # Fst Calculation
  vcf <- read.vcfR('~/MTOR.25.leads.in.1KG.data.vcf'); \
  genind_obj <- vcfR2genind(vcf); \
  ind_names <- indNames(genind_obj); \
  PED <- PED[which(PED$IndividualID %in% ind_names), ]; \
  population_factor <- PED$`Population`[match(ind_names, PED$IndividualID)]; \
  genind_obj@pop <- as.factor(population_factor); \
  hierfstat_data <- genind2hierfstat(genind_obj); \
  pairwise_fst <- pairwise.WCfst(hierfstat_data); \
  pheatmap(pairwise_fst, display_numbers = TRUE, main = 'Pairwise Fst Heatmap');"

