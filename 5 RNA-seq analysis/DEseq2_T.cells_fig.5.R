#-------------------------------------------------------------------------------------------#
# DEseq2
#-------------------------------------------------------------------------------------------#
library(DESeq2)
library(ggplot2)
library(data.table)
# Read the data from the standard input.
countData.raw = read.table("~/1.GitHub_Code availability/output/RNAseq_featureCounts_Tcell.day0.day4.txt.gz", header=TRUE, sep="\t", row.names=1 )
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
meta.data <- fread("~/1.GitHub_Code availability//output/RNAseqseq.meta.csv", stringsAsFactors = F) 

########
countData = countData[,which(names(countData) %in% meta.data$Sample.ID)]
# check if match !!!!!!!!!!!!!!!!!!!!
meta.data_2 = meta.data[match(colnames(countData), meta.data$Sample.ID),]
all(colnames(countData) == meta.data_2$Sample.ID)  # TRUE
# define the contrast
meta.data_2$Time.point = as.factor(meta.data_2$Time.point )
meta.data_2$Time.point <- relevel(meta.data_2$Time.point, ref = "Day 0")

dds <- DESeqDataSetFromMatrix(countData = countData, colData = meta.data_2, design = ~Time.point) # + 
#colData(dds)
dds <- DESeq(dds)  #

# -------------------------------------------------------------###
# PCA PC1-PC2
# -------------------------------------------------------------###
vsd <- vst(dds, blind=TRUE)
data <- plotPCA(vsd, intgroup=c("Cell.type", "Time.point","Sample.ID", "Biological.replicate"), returnData=TRUE, ntop=500)
percentVar <- round(100 * attr(data, "percentVar"))
library(ggplot2)
ggplot(data, aes(PC1, PC2, col=Time.point, shape=Cell.type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme_bw() +
  scale_colour_manual(values = c("grey", "darkblue")) +
  scale_shape_manual(values=c(1, 10)) + ggtitle("RNA-seq") +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  geom_text_repel(data=data, aes(label = Sample.ID),
                  size = 1.5, point.padding = unit(0.2, "lines") ,min.segment.length = 0) #


########
resultsNames(dds)
res <- data.frame(results(dds,name = "Time.point_Day.4_vs_Day.0"))
res$gene_name = rownames(res)

##### box plot for MTOR expression related to fig.5
library("ggplot2")
library("ggpubr")
Padj = signif(res[res$gene_name == "MTOR",]$padj, digits = 3)

MTOR = plotCounts(dds, gene=c("MTOR"), intgroup=c("Time.point","Biological.replicate","Cell.type","Biological.replicate"),returnData=TRUE) 
MTOR$gene = "MTOR"

# log expression
MTOR$value = log2(MTOR$count+1/2)
#
MTOR$id = paste0(MTOR$Biological.replicate, MTOR$Cell.type)
MTOR$Tn.Ta = ifelse(MTOR$Time.point == "Day 0", "Rest", "Activated.T")

## Note that, the sample size should be at least 6. Otherwise, the Wilcoxon test cannot become significant.
my_comparisons <- list( c("Rest", "Activated.T"))

MTOR$Tn.Ta = factor(MTOR$Tn.Ta, levels = c("Rest", "Activated.T"))

ggplot(MTOR, aes(x = Tn.Ta, y = value)) + 
  geom_line(aes(group=id, colour = Cell.type),linetype=2, size=0.5) + 
  geom_jitter(aes(colour = Cell.type), shape=16, position=position_jitter(0.1), alpha=0.8)+
  geom_boxplot(outlier.alpha = 0.5, outlier.shape = NA, alpha=0.5) +
  theme_bw() + ylab("Log2 normalised counts") +
  stat_compare_means(comparisons=my_comparisons, method = "wilcox.test", # method.args = list(var.equal = TRUE)
                     label = "p.signif")  +
  theme(panel.grid.minor.y = element_blank(),panel.grid.major.x = element_blank())+
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) +
  scale_colour_brewer(palette="Dark2") +
  annotate("text", x=2, y= 10, label=paste0("Padj = \n", Padj))

ggsave("xx.pdf", width = 2.5, height = 3)


