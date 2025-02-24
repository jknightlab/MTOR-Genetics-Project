#-------------------------------------------------------------------------------------------#
# DEseq2
#-------------------------------------------------------------------------------------------#
library(DESeq2)
library(ggplot2)
library(data.table)
# Read the data from the standard input.
countData.raw = read.table("~/1.GitHub_Code availability/output/RNAseq_featureCounts_CD4_Rapa_PRJNA532911.txt.gz", header=TRUE, sep="\t", row.names=1 )
#
colnames(countData.raw)
# clean the columns
colClean2 <- function(x){ colnames(x) <- gsub("Mapping.mapped.hisat2.", "", colnames(x)); x } 
countData.raw = colClean2(countData.raw)
colClean2 <- function(x){ colnames(x) <- gsub(".hisat2.bam", "", colnames(x)); x } 
countData.raw = colClean2(countData.raw)

#--------------------------#
######## the clinical file
#--------------------------#
####
colnames(countData.raw)
countData = countData.raw[,-c(1:6)]
colnames(countData)
#
Replicates <- as.factor(c("donor318","donor435","donor318","donor124","donor435","donor124"))
Treatment <- as.factor(c("rapa","rapa","DMSO","DMSO","DMSO","rapa"))
#
colData <- data.frame(colnames(countData), Replicates=Replicates,Treatment= Treatment, row.names=1)
colData$Treatment <- relevel(colData$Treatment, ref = "DMSO")

########
# define the contrast
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData,  design = ~Treatment + Replicates) # + 
#colData(dds)
dds <- DESeq(dds)  #

xx = as.data.frame(sizeFactors(dds))

# -------------------------------------------------------------###
# PCA PC1-PC2
# -------------------------------------------------------------###
vsd <- vst(dds, blind=TRUE)
data <- plotPCA(vsd, intgroup=c("Replicates", "Treatment"), returnData=TRUE, ntop=500)
percentVar <- round(100 * attr(data, "percentVar"))

library(ggplot2)
ggplot(data, aes(PC1, PC2, col=Treatment, shape=Replicates)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme_bw() +
  scale_colour_manual(values = c("grey", "darkblue", "red")) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())


# # -------------------------------------------------------------###
# # removeBatchEffect for PCA
# # -------------------------------------------------------------###
vsd <- vst(dds,blind=TRUE)
plotPCA(vsd, "Replicates")
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$Replicates)
data <- plotPCA(vsd, intgroup=c("Replicates", "Treatment"), returnData=TRUE, ntop=500)
percentVar <- round(100 * attr(data, "percentVar"))

library(ggplot2)
ggplot(data, aes(PC1, PC2, col=Treatment, shape=Replicates)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme_bw() +
  scale_colour_manual(values = c("grey", "red")) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())


########
resultsNames(dds)
#Pre-filtering
res <- data.frame(results(dds,name = "Treatment_rapa_vs_DMSO"))
res$gene_name = rownames(res)
res.sig= res[res$padj <0.05 & abs(res$log2FoldChange) > 1 & !is.na(res$padj),]

### 
sig.cytokines.CD4 <- fread("~/MTOR.project//results/cytokines.CD4.sig.day4.vs.day0.csv")
sig.cytokines.CD8 <- fread("~/MTOR.project/results/cytokines.CD8.sig.day4.vs.day0.csv")
sig.cytokines = merge(sig.cytokines.CD4, sig.cytokines.CD8, by="gene_name")

###
rapa.sig.cytokine = res.sig[res.sig$gene_name %in% sig.cytokines$gene_name,]

# volcano plot related to Fig.4g
res$sig.cytokines = ifelse(res$gene_name %in% sig.cytokines$gene_name, "Yes", "No")

# reduce the plot size
res = res[!is.na(res$padj),]
res$temp <- signif(res$log2FoldChange * res$padj, digits = 3)
res2 = res[!duplicated(res$temp),]

ggplot(res2, aes(x=log2FoldChange, y=-log10(padj), col=sig.cytokines)) +
  geom_point(size=0.5) + 
  geom_point(data=res[res$sig.cytokines == "Yes", ], shape=21, col="black", fill="purple") + 
  xlab("Log2 fold change") + ylab("-Log10 adjusted P value")+
  scale_color_manual(values = c("grey","purple"))+
  geom_hline(yintercept= -log10(0.05), linetype = "longdash", col='black')+
  geom_vline(xintercept= c(log2(2),-log2(2)), linetype = "longdash",col='black')+
  theme_bw() +  theme(legend.position="none") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none")  +
  ggrepel::geom_text_repel(data=res[res$gene_name %in% rapa.sig.cytokine$gene_name, ], 
                  aes(label = gene_name),
                  size = 4, colour="darkblue",point.padding = unit(0.05, "lines"),
                  max.overlaps =100, min.segment.length = 0 )


