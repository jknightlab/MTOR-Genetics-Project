#-------------------------------------------------------------------------------------------#
# DEseq2
#-------------------------------------------------------------------------------------------#
library(DESeq2)
library(ggplot2)
library(data.table)
# Read the data from the standard input.
countData.raw = read.table("./output/RNAseq_featureCounts_Tcell.day0.day4.txt.gz", header=TRUE, sep="\t", row.names=1 )
#
colnames(countData.raw)
# clean the columns
colClean2 <- function(x){ colnames(x) <- gsub("Mapping.mapped.hisat2.", "", colnames(x)); x } 
countData.raw = colClean2(countData.raw)
colClean2 <- function(x){ colnames(x) <- gsub(".hisat2.bam", "", colnames(x)); x } 
countData.raw = colClean2(countData.raw)

colClean2 <- function(x){ colnames(x) <- gsub("[.]", "_", colnames(x)); x } 
countData.raw = colClean2(countData.raw)

#--------------------------#
######## the meta file
#--------------------------#
####
colnames(countData.raw)
countData = countData.raw[,-c(1:6)]

##
meta.data <- fread("./output/RNAseqseq.meta.csv", stringsAsFactors = F) 
meta.data = meta.data[meta.data$Cell.type %in% "CD8",]

########
countData = countData[,which(names(countData) %in% meta.data$Sample.ID)]
# check if match !!!!!!!!!!!!!!!!!!!!
meta.data_2 = meta.data[match(colnames(countData), meta.data$Sample.ID),]
all(colnames(countData) == meta.data_2$Sample.ID)  # TRUE


# define the contrast
dds <- DESeqDataSetFromMatrix(countData = countData, colData = meta.data_2, design = ~Time.point) # + 
#colData(dds)
dds <- DESeq(dds)  #

# -------------------------------------------------------------###
# PCA PC1-PC2
# -------------------------------------------------------------###
vsd <- vst(dds, blind=TRUE)
data <- plotPCA(vsd, intgroup=c("Cell.type", "Time.point","Sample.ID"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

library(ggplot2)
ggplot(data, aes(PC1, PC2, col=Time.point, shape=Cell.type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme_bw() +
  scale_colour_manual(values = c("darkblue", "red")) +
  scale_shape_manual(values=c(1, 10)) + ggtitle("RNA-seq") +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
  # geom_text_repel(data=data, aes(label = Sample.ID),
  #                 size = 1.5, point.padding = unit(0.2, "lines") ,min.segment.length = 0) # 

##ggsave("xx.pdf",width = 3.5, height = 2.2)

########
resultsNames(dds)
#Pre-filtering
res <- data.frame(results(dds,name = "Time.point_Day.4_vs_Day.0"))
res$gene_name = rownames(res)

# write.csv(res, "Time.point_Chronic.Day.12_vs_Acute.Day.12_CD8.csv", row.names = F)

res.sig= res[res$padj <0.05 & abs(res$log2FoldChange) > 1 & !is.na(res$padj),] #
#write.csv(res.sig, "RNA-seq_CD4_Time.point_Day.4_vs_Day.0.sig.csv", row.names = F)

#####
#
library(ggplot2)
library(ggrepel)
plot = res[!is.na(res$padj) ,]
plot$sig = ifelse(plot$padj < 0.05 & plot$log2FoldChange > 1, "upregulated",
                  ifelse(plot$padj < 0.05 & plot$log2FoldChange < -1, "downregulated", "n.s"))


# downloaded the cytokines from https://www.uniprot.org/uniprotkb?query=%28keyword%3AKW-0202%29&facets=reviewed%3Atrue%2Cmodel_organism%3A9606
cytokines <- fread("~/MTOR.project//downloaded.data/uniprotkb_189.cytokines_human_2023_10_13.tsv")
cytokines$gene.symbol <- gsub("[ ].*", "", cytokines$`Gene Names`) 

plot$sig.cytokines = ifelse(plot$sig %in% c("upregulated", "downregulated")
                            & plot$gene_name %in% cytokines$gene.symbol, "TRUE", "FALSE")

table(plot$sig.cytokines)

# reduce the plot size
ggplot(plot, aes(x=log2FoldChange, y=-log10(padj), col=sig.cytokines)) +
  geom_point(size=0.5) + 
  geom_point(data=plot[plot$sig.cytokines == "TRUE", ], shape=21, col="black", fill="purple") + 
  xlab("Log2 fold change") + ylab("-Log10 adjusted P value")+
  scale_color_manual(values = c("grey","purple"))+
  geom_hline(yintercept= -log10(0.05), linetype = "longdash", col='black')+
  geom_vline(xintercept= c(log2(2),-log2(2)), linetype = "longdash",col='black')+
  theme_bw() +  theme(legend.position="none") + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position="none") 
# 3X3.3


# cytokines.CD8.sig= res.sig[res.sig$gene_name %in% cytokines$gene.symbol,]
# write.csv(cytokines.CD8.sig, "cytokines.CD8.sig.day4.vs.day0.csv", row.names = F)
# 
# #
# cytokines.CD4.sig= res.sig[res.sig$gene_name %in% cytokines$gene.symbol,]
# write.csv(cytokines.CD4.sig, "cytokines.CD4.sig.day4.vs.day0.csv", row.names = F)


#'@-VennDiagram

# venn diagram for eQTL-gene pairs

library(VennDiagram)
library(RColorBrewer)

cytokines.CD4.sig <- read.csv("~/MTOR.project//results/cytokines.CD4.sig.day4.vs.day0.csv")
cytokines.CD8.sig <- read.csv("~/MTOR.project/results/cytokines.CD8.sig.day4.vs.day0.csv")
###
cytokines.CD4.sig$direction = ifelse(cytokines.CD4.sig$log2FoldChange > 0, "upregulated", "downregulated")
cytokines.CD4.sig$id = paste0(cytokines.CD4.sig$gene_name, "-", cytokines.CD4.sig$direction)

cytokines.CD8.sig$direction = ifelse(cytokines.CD8.sig$log2FoldChange > 0, "upregulated", "downregulated")
cytokines.CD8.sig$id = paste0(cytokines.CD8.sig$gene_name, "-", cytokines.CD8.sig$direction)


xxx = merge(cytokines.CD4.sig, cytokines.CD8.sig, by="gene_name")

X=list(cytokines.CD4.sig=cytokines.CD4.sig$id, 
       cytokines.CD8.sig=cytokines.CD8.sig$id)
v=venn.diagram(X, filename = NULL,fill=c( "purple","cyan3"),col="transparent",print.mode="raw")
grid.draw(v)


library(eulerr)
fit1 <- euler(c("A" = 20,
                "B" = 11,
                "A&B" = 48), shape = "ellipse")

plot(fit1, quantities = TRUE,labels = c("cytokines.CD4+","cytokines.CD8+",""),
     fill = c("purple","cyan4"),alpha=0.6)


# Correlation plot

cor.test(xxx$log2FoldChange.x, xxx$log2FoldChange.y)

library(ggplot2)
library(ggpubr) 
ggplot(xxx, aes(x=log2FoldChange.x, y=log2FoldChange.y)) + 
  geom_point(fill="darkorange", shape=21)+
  geom_smooth(method=lm, colour="cyan4") + theme_classic() +
  stat_cor(method = "pearson", label.x = -1, label.y = 0.5) +
  xlab("Log2 fold change (CD4+)") + ylab("Log2 fold change (CD8+)") + theme_minimal()

ggsave("xx.pdf", width = 2.8, height = 1.8)

