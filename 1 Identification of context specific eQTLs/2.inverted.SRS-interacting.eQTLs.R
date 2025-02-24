library(data.table)
library(openxlsx)
library(seqminer)
library(readr)
library(ggplot2)
library(cowplot)

#'@1------------------------------>>>>> input data
SRS.SNPs <- read.xlsx("~/MTOR.project/Katie/Kaite_suppl.tables.xlsx", sheet = 7) # publicly available in PMID: 38897207
SRS.SNPs_2 = SRS.SNPs[SRS.SNPs$FDR <= 0.01,]
length(unique(as.factor(SRS.SNPs_2$SNP)))

## add positions
library(biomaRt)
ensembl <- useMart(host = "https://feb2023.archive.ensembl.org", "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

pos <- getBM(attributes=c(
  "refsnp_id", "chr_name", "chrom_start", "chrom_end"),
  filters="snp_filter", values=SRS.SNPs_2$SNP,
  mart=ensembl, uniqueRows=TRUE)

pos2 = pos[!pos$chr_name %like% "CHR",]
pos2$position = paste0(pos2$chr_name, ":", pos2$chrom_start, "-", pos2$chrom_end)
pos2 = pos2[,c("refsnp_id","position")]
length(unique(as.factor(pos2$refsnp_id))) # generate duplicates for non-SRS which were removed for plot
pos2 = pos2[!duplicated(pos2$refsnp_id),]
#
xxx= SRS.SNPs_2[!SRS.SNPs_2$SNP %in% pos2$refsnp_id,] # lost 3 SNPs for SRS-interacted
length(unique(as.factor(xxx$SNP)))


###
SRS.SNPs_3 = merge(SRS.SNPs_2, pos2, by.x="SNP", by.y="refsnp_id")
length(unique(as.factor(SRS.SNPs_3$SNP)))

##


#'@2------------------->>>> retrive the pairs in eQTL Catalog
#---------> function 1
import_eQTLCatalogue <- function(ftp_path, region, selected_gene_id, column_names, verbose = TRUE){
  # if(verbose){
  #   print(ftp_path)
  # }
  #Fetch summary statistics with seqminer
  fetch_table = seqminer::tabix.read.table(tabixFile = ftp_path, tabixRange = region, stringsAsFactors = FALSE) %>%
    dplyr::as_tibble()
  colnames(fetch_table) = column_names
  #Remove rsid duplicates and multi-allelic variant
  summary_stats = dplyr::filter(fetch_table, gene_id == selected_gene_id) %>%
    filter(row_number()==1)
    # dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>% 
    # dplyr::group_by(id) %>% 
    # dplyr::mutate(row_count = n()) %>% dplyr::ungroup() %>% 
    # dplyr::filter(row_count == 1) #Multialllics. # this will cause problem for microarray data with mutiple probes
}

library("dplyr")

tabix_paths = read.delim("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()
### solve the connection issue
tabix_paths$ftp_path = gsub("ftp://ftp.ebi.ac.uk", "http://ftp.ebi.ac.uk", tabix_paths$ftp_path )
tabix_paths$dataset = paste0(tabix_paths$study_label,"_" ,tabix_paths$sample_group)

# 7 T cells/neutrophil datasets
subset = tabix_paths[tabix_paths$dataset %like% "neutrophil",]
subset = subset[subset$quant_method %in% c("ge","microarray"),]

#Extract column names from first file
column_names = colnames(readr::read_tsv("http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000002/QTD000026/QTD000026.all.tsv.gz", n_max = 1))

zzz <- list()

for (a in 1:length(SRS.SNPs_3$SNP)) {
  
  input = SRS.SNPs_3[a,]

resul <- list()
j = as.character(subset$dataset)
for (i in j) {
  #
  platelet_df = dplyr::filter(tabix_paths, dataset %in% i, 
                              quant_method %in% c("ge","microarray") ) 
  #Import summary statistics
  ##--------------------------------------------##    pass the data  without SNP retrived in agiven region  in eQTL Catalog
  xxxx = seqminer::tabix.read.table(tabixFile = platelet_df$ftp_path, 
                                    tabixRange = input$position,
                                    stringsAsFactors = FALSE) %>%
    dplyr::as_tibble()
  #
  if (length(xxxx) == 0){
    next
  }
  ##--------------------------------------------## 
  
  resul[[i]] = import_eQTLCatalogue(platelet_df$ftp_path,
                                    region = input$position,
                                    selected_gene_id = input$Gene,
                                    column_names)
    resul[[i]]$dataset = i
    resul[[i]]$pair = paste0(input$SNP, "_", input$Symbol)
  
}

resul <- as.data.frame(do.call(rbind, resul))
zzz[[a]] = resul

}

zzz <- as.data.frame(do.call(rbind, zzz))
length(unique(as.factor(zzz$pair)))

# write.csv(zzz, "SRS.interacted.eQTLs_in_neutrophils.csv", row.names = F)
# write.csv(zzz, "SRS.interacted.eQTLs_in_T.cells.csv", row.names = F)


#-------------------------------------------------------------------------------
# Step 3: Analyse and plot results
#-------------------------------------------------------------------------------
df <- as.data.frame(fread("~/MTOR.project//results/SRS.interacted_Catalog.Neu.T.cell.txt.gz", stringsAsFactors = F))
df = df[!duplicated(df[c("pair","dataset")]),]
length(unique(as.factor(df$pair))) # 845 pairs

# adjust p = 0.05/4748 = 1.053075e-05
df$adjusted.p = p.adjust(df$pvalue, method = "bonferroni", n = length(df$pvalue))

# take the ones with at least one sig association with p value < 1e-05
#'@ this filtered out all associations with p > 1e05
df_2 = df[df$adjusted.p <= 0.05,]
length(unique(as.factor(df_2$pair))) # 557 pairs

#'@ this inverted phenotype is based on the beta that has p value <1e-05
df_3 <- reshape2::dcast(df_2, pair~dataset, value.var = "beta")
df_3$direction = ifelse(df_3$BLUEPRINT_neutrophil * df_3$`Schmiedel_2018_CD8_T-cell_anti-CD3-CD28` < 0 | 
                          df_3$BLUEPRINT_neutrophil * df_3$`Schmiedel_2018_CD4_T-cell_anti-CD3-CD28` < 0 |
                          df_3$CEDAR_neutrophil_CD15 * df_3$`Schmiedel_2018_CD8_T-cell_anti-CD3-CD28` < 0 | 
                          df_3$CEDAR_neutrophil_CD15 * df_3$`Schmiedel_2018_CD4_T-cell_anti-CD3-CD28` < 0 |
                          df_3$Naranbhai_2015_neutrophil_CD16 * df_3$`Schmiedel_2018_CD8_T-cell_anti-CD3-CD28` < 0 | 
                          df_3$Naranbhai_2015_neutrophil_CD16 * df_3$`Schmiedel_2018_CD4_T-cell_anti-CD3-CD28` < 0 |
                          df_3$BLUEPRINT_neutrophil * df_3$`Schmiedel_2018_CD8_T-cell_naive` < 0 | 
                          df_3$BLUEPRINT_neutrophil * df_3$`Schmiedel_2018_CD4_T-cell_naive` < 0 |
                          df_3$CEDAR_neutrophil_CD15 * df_3$`Schmiedel_2018_CD8_T-cell_naive` < 0 | 
                          df_3$CEDAR_neutrophil_CD15 * df_3$`Schmiedel_2018_CD4_T-cell_naive` < 0 |
                          df_3$Naranbhai_2015_neutrophil_CD16 * df_3$`Schmiedel_2018_CD8_T-cell_naive` < 0 | 
                          df_3$Naranbhai_2015_neutrophil_CD16 * df_3$`Schmiedel_2018_CD4_T-cell_naive` < 0, "inverted", "all-non-inverted")
table(df_3$direction) # 46 inverted
inverted = df_3[df_3$direction %in% "inverted",]
inverted$gene.symble = gsub(".*[_]","",inverted$pair)


#'@-generate-suppl.table1-sheet
SRS.SNPs_3$pair = paste0(SRS.SNPs_3$SNP, "_", SRS.SNPs_3$Symbol)
add.info = df[!duplicated(df$pair),c("variant","pair")]

sheet1 = merge(inverted, SRS.SNPs_3, by="pair")
sheet1 = merge(add.info, sheet1, by="pair")
#write.xlsx(sheet1, "45.uniq.NLR.inverted.sepsis.eQTL.xlsx")


#'@ df_2 = all the association with p > 1e-05 was filtered out
plot = df[df$pair %in% inverted$pair,]
plot$rsid = gsub("[_].*","",plot$pair)
plot$log10P = -log10(plot$adjusted.p)
plot$log10P = ifelse(plot$log10P >=50, 50, plot$log10P)
length(unique(as.factor(plot$pair)))

#
effect.size <- reshape2::dcast(plot, pair~dataset, value.var = "beta")
rownames(effect.size) = effect.size$pair
effect.size$pair = NULL
plot_size <- as.matrix(effect.size)
plot_size[is.na(plot_size)] <- 0

# make data square to calculate euclidean distance
clust <- hclust(dist(plot_size)) # hclust with distance matrix
ddgram <- as.dendrogram(clust) # create dendrogram
ggtree_plot <- ggtree::ggtree(ddgram,branch.length="none")
ggtree_plot


#glue them together with cowplot
library(dplyr)
plot.ordered <- plot %>%  
  mutate( pair = factor(pair, levels = clust$labels[clust$order]))


dotplot <- ggplot(plot.ordered, aes(x=dataset,  y=pair) ) +
  geom_point(aes(size=log10P , 
                 colour=beta)) + 
  scale_size_area(limits = c(0, 50), breaks = c(1.3,5,10,20,30,40,50)) +
  scale_colour_gradient2(  low = "darkorange4",
                           mid = "white",
                           high = "cyan4") +
  theme_bw() + 
  theme(panel.grid.major.x = element_blank() ) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) + xlab("") +
  theme(axis.text.x = element_text(size = 5 ), legend.key.size = unit(4, "mm")) +
  scale_y_discrete(position = "right")

dotplot

plot_grid(ggtree_plot,  NULL, dotplot, nrow = 1, rel_widths = c(0.5,-0.05, 2), align = 'h') # 5X8

