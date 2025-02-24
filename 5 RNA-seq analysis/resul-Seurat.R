library(Seurat)
library(scCustomize)
#
seurat <- read_rds("rhapsody_wholeblood_sobj.rds")
Idents(seurat) <- "broad_annot"

#
SRS.obj <- seurat[,seurat$SRS_comparison %in% c("SRS1","non_SRS1","HV")]
#
n_colors <- length(unique(SRS.obj$broad_annot))
# Create a palette
palette <- as.character(pals::alphabet(n_colors))

DimPlot(object = SRS.obj, reduction = 'umap',cols = palette, raster=FALSE) #5X3


#'@_##_plot-with-Seurat-function 
SRS.obj$SRS_comparison <- factor(SRS.obj$SRS_comparison, levels = c("HV","non_SRS1","SRS1"))
FeaturePlot(object = SRS.obj, features = c("CD274"),  order = T,
            split.by = "SRS_comparison", raster=FALSE) 
#'@_##_plot-with-FeaturePlot_scCustom-from-scCustomize-package
### # activation markers: CD69 (early), CD25 (late)/IL2RA, and HLA-DR (even later).
feature_plot <- scCustomize::FeaturePlot_scCustom(seurat_object = SRS.obj, features = "CD274", 
                                  split.by = "SRS_comparison",raster=FALSE) & NoAxes() & theme(legend.position = "bottom")
ggsave(filename = "feature_plot.png", plot = feature_plot, width = 4, height = 3, dpi = 400)









#'@-1 subset the T cells and HV and sepsis
T_cells <- seurat[,seurat$broad_annot %in% c("CD4_T_cells","CD8_T_cells") & seurat$source %in% c("Sepsis")]


## drop level
T_cells@meta.data$source <- droplevels(T_cells@meta.data$source)
table(T_cells$source)

T_cells@meta.data$sample_id <- droplevels(T_cells@meta.data$sample_id)
table(T_cells$sample_id)

T_cells@meta.data$SRS_comparison <- droplevels(T_cells@meta.data$SRS_comparison)
table(T_cells$SRS_comparison)


#'@-2 vst normalisation
T_cells <- FindVariableFeatures(T_cells, selection.method = "vst", nfeatures = 2000)
T_cells <- ScaleData(T_cells)


##https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_clustering_analysis.html. # define the k 
#https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/trajectory-inference.html
#https://romanhaa.github.io/projects/scrnaseq_workflow/#create-seurat-object # nice figures



# Identify the 20 most highly variable genes
top20 <- head(x = VariableFeatures(object = T_cells), 
              n =20)
##Plot variable features with labels
plot1 <- VariableFeaturePlot(object = T_cells)
LabelPoints(plot = plot1, 
            points = top20, 
            repel = TRUE)


#'@-3  run PCA
T_cells <- RunPCA(T_cells, features = VariableFeatures(object = T_cells))

##Plot the elbow plot
ElbowPlot(object = T_cells, 
          ndims = 30)
##focus on the first 10 principal components that capture the majority of the variability in this data set. 
set.seed(1)
T_cells <- RunUMAP(T_cells, reduction = "pca", dims = 1:10)

## Plot the UMAP
DimPlot(T_cells,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

#'@-4 perform unsupervised clustering on these cells.
set.seed(1)
T_cells <- FindNeighbors(T_cells, dims = 1:10, reduction = "pca") 
T_cells <- FindClusters(T_cells, resolution = 1, algorithm = 2)

#
DimPlot(T_cells, reduction = "umap", label=T)






# CD8+ T and CD4+ TEMRA markers (Cano-Gamez et al., 2020)
CD8_T_cell_plot_1 <- FeaturePlot(T_cells, reduction = "umap", features = "GZMK", pt.size = 0.5)
CD8_T_cell_plot_2 <- FeaturePlot(T_cells, reduction = "umap", features = "CD8A", pt.size = 0.5)
plot_grid( CD8_T_cell_plot_1,CD8_T_cell_plot_2, ncol = 2, nrow = 2
)
# CD4+ Treg markers
Treg_plot_1 <- FeaturePlot(T_cells, reduction = "umap", features = "FOXP3", pt.size = 0.5)
Treg_plot_2 <- FeaturePlot(T_cells, reduction = "umap", features = "IL2RA", pt.size = 0.5)
plot_grid( Treg_plot_1,Treg_plot_2, ncol = 2, nrow = 2
)
# Naive CD4+ T
Treg_plot_1 <- FeaturePlot(T_cells, reduction = "umap", features = "CCR7", pt.size = 0.5)
Treg_plot_2 <- FeaturePlot(T_cells, reduction = "umap", features = "CD45RA", pt.size = 0.5)
plot_grid( Treg_plot_1,Treg_plot_2, ncol = 2, nrow = 2
)


# Activated T cells
Treg_plot_1 <- FeaturePlot(T_cells, reduction = "umap", features = "IL2RA", pt.size = 0.5)
Treg_plot_2 <- FeaturePlot(T_cells, reduction = "umap", features = "CD44", pt.size = 0.5)
plot_grid( Treg_plot_1,Treg_plot_2, ncol = 2, nrow = 2
)


#------------------------------------------------------------------------------------------#
#'@-Cell.type.Annotation.using.SingleR&HumanPrimaryCellAtlasData
#'# the Blueprint/ENCODE and Human Primary Cell Atlas reference datasets
# Loading reference data with Ensembl annotations.
# https://bioconductor.org/books/release/SingleRBook/using-multiple-references.html
library(celldex)



ref.data <- HumanPrimaryCellAtlasData() # BlueprintEncodeData HumanPrimaryCellAtlasData

# Performing predictions.
library(SingleR)
T_cells_exp <- Seurat::as.SingleCellExperiment(T_cells)
predictions <- SingleR(test=T_cells_exp, assay.type.test=1, 
                       ref=ref.data, labels=ref.data$label.fine)

table(predictions$labels)

# plotScoreHeatmap(predictions, show.labels = TRUE,  annotation_col=data.frame(
#   row.names=rownames(predictions)))

T_cells_exp$celltypes <- predictions$pruned.labels

## Plot  UMAP
umap_pl <- scater::plotReducedDim(T_cells_exp, dimred = "UMAP", colour_by="celltypes", text_by = "celltypes", 
                                  text_size = 3, point_size=0.5) +
  guides(fill="none")

nh_graph_pl <- scater::plotReducedDim(T_cells_exp, dimred = "UMAP",  text_by = "seurat_clusters", colour_by="seurat_clusters",
                                text_size = 3, point_size=0.5)
umap_pl + nh_graph_pl +
  patchwork::plot_layout(guides="collect")


#------------------------------------------------------------------------------------------#






gene_List <- data.frame(rownames(T_cells))

#'@_##_plot-with-Seurat-function 
FeaturePlot(object = T_cells, features = c("PDCD1"),  order = T,
            split.by = "SRS_comparison", raster=FALSE) & NoAxes()
#'@_##_plot-with-FeaturePlot_scCustom-from-scCustomize-package
### # activation markers: CD69 (early), CD25 (late)/IL2RA, and HLA-DR (even later).
scCustomize::FeaturePlot_scCustom(seurat_object = T_cells, features = "TNF", 
                                  split.by = "SRS_comparison",raster=FALSE) & NoAxes()

VlnPlot(T_cells, features = c( "IFNG"),split.by = "SRS_comparison",
        pt.size = 0.2, ncol = 4)



## heatmap for top5 Differential expression testing identifies the following genes as markers of each cluster.
cluster_markers <- FindAllMarkers(T_cells, min.pct = 0.25, logfc.threshold = 0.25)

cluster_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(T_cells, features = top10$gene) + NoLegend()


#



#---------------------------------------------------------------------------------------------#
#'@-----differential-abundance-analysis
# https://www.bioconductor.org/packages/devel/bioc/vignettes/miloR/inst/doc/milo_gastrulation.html#3_Differential_abundance_testing
# https://github.com/REdahiro/JPN_COVID-19_scRNAseq/blob/main/scripts/01.1_Milo_Case_vs_Contrl_harmony.R
#---------------------------------------------------------------------------------------------#
library(miloR)
library(SummarizedExperiment)
library(SingleCellExperiment)

T_cells_exp <- Seurat::as.SingleCellExperiment(T_cells)
table(T_cells_exp$SRS_comparison)

milo <- Milo(T_cells_exp)
milo



# Construct KNN graph & Defining representative neighbourhoods on the KNN graph
# buildFromAdjacency

reducedDim(milo, "PCA", withDimnames=TRUE) <- T_cells[['pca']]@cell.embeddings      
reducedDim(milo, "UMAP", withDimnames=TRUE) <- T_cells[['umap']]@cell.embeddings


embryo_milo <- buildGraph(milo, k = 30, d = 30, reduced.dim = "PCA")
embryo_milo <- makeNhoods(embryo_milo, prop = 0.2, k = 30, d=30, refined = TRUE, reduced_dims = "PCA")

## > 5 x N_samples
plotNhoodSizeHist(embryo_milo) + 
               geom_vline(xintercept=50, linetype="dashed", color="red") +
               geom_vline(xintercept=100, linetype="dashed", color="red") +
               geom_vline(xintercept=(5*length(unique(embryo_milo$sample_id))), 
                          linetype = "dashed", color = "grey")

# Counting cells in neighbourhoods
embryo_milo <- countCells(embryo_milo, meta.data = as.data.frame(colData(embryo_milo)), 
                          sample="sample_id")




#'@: Defining experimental design
design <- data.frame(colData(embryo_milo))[,c("sample_id", "SRS_comparison","age", "sex", "batch")] # SRS_comparison
str(design)
design$age = as.numeric(as.character(design$age))
design <- distinct(design)
rownames(design) <- design$sample_id
design

#'@: Computing neighbourhood connectivity---##
milo <- calcNhoodDistance(embryo_milo, d=30, reduced.dim = "PCA")

#'@: Testing -----##
# The last component of the formula or last column of the model matrix are by default the test variable.
library(statmod)
da_results <- testNhoods(milo, design = ~ age + sex + batch + SRS_comparison, design.df = design)

da_results %>%
  arrange(-SpatialFDR) %>%
  head() 




milo <- buildNhoodGraph(milo)
## Plot single-cell UMAP
umap_pl <- scater::plotReducedDim(milo, dimred = "UMAP", colour_by="SRS_comparison", text_by = "seurat_clusters", 
                                  text_size = 3, point_size=0.5) +
  guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(milo, da_results, layout="UMAP",alpha=0.1) 

umap_pl + nh_graph_pl +
  patchwork::plot_layout(guides="collect")


## 
da_results <- annotateNhoods(milo, da_results, coldata_col = "seurat_clusters")
head(da_results)
str(da_results)

# da_results$seurat_clusters <- ifelse(da_results$seurat_clusters_fraction < 0.7, "Mixed", da_results$seurat_clusters)
plotDAbeeswarm(da_results, group.by = "seurat_clusters")











#' test.subset <- subset(x = seurat, subset = (SRS_comparison == "HV" | SRS_comparison == "SRS1"))
#' 
#' 
#' test.subset <- subset(x = seurat, subset = (SRS_comparison == "HV"))
#' 
#' metadata <- test.subset@meta.data
#' 
#' 
#' write.csv(100 * (exp(as.matrix(GetAssayData(object = object, assay = "RNA", slot = "data"))) - 1), "em.csv", row.names = T) # Seurat 3.X
#' 
#' write.csv(Idents(object = seurat),"metadata.csv", row.names = T) # Seurat 3.X
#' 
#' 
#' # https://satijalab.org/seurat/archive/v3.1/interaction_vignette.html
#' # https://github.com/satijalab/seurat/issues/1883
#' # https://github.com/satijalab/seurat/issues/2539
#' 
#' # https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_clustering_analysis.html
#' 
#' 
#' 
#' 
#' # extract the cell annotation of interest
#' 
#' metadata <- seurat@meta.data
#' 
#' table(metadata$fine_annot)
#' table(metadata$diagnosis)
#' table(metadata$source)
#' table(metadata$SRS_comparison)
#' table(metadata$broad_annot)
#' table(metadata$fine_annot)
#' seurat@meta.data$seurat_clusters
#' 
#' 
#' Idents(seurat) <- "fine_annot"
#' 
#' DimPlot(seurat, label=TRUE, raster=TRUE, repel=TRUE) & NoLegend()
#' 
#' 
#' 
#' ####
#' metadata %>%
#'   group_by(sample_id) %>%
#'   mutate(total_cells = n()) %>%
#'   ungroup() %>%
#'   group_by(sample_id, fine_annot) %>%
#'   mutate(fine_annot_cellnum = n(),
#'          prop = (fine_annot_cellnum/total_cells) * 100) %>%
#'   distinct_at(.vars = vars(sample_id, fine_annot),.keep_all=TRUE) %>%
#'   dplyr::filter(fine_annot == "CD8_T_cells") %>%
#'   ggplot(aes(x = SRS_comparison, y = prop)) +
#'   geom_boxplot(outlier.shape = NA) +
#'   geom_jitter(aes(colour = SRS_comparison), shape=16, position=position_jitter(0.1), alpha=0.8, size=2)+
#'   theme_classic() 
#'  
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' # Set color palette
#' pal <- viridis(n = 10, option = "C", direction = -1)
#' # Create Plots
#' # https://satijalab.org/seurat/articles/visualization_vignette
#' 
#' Idents(seurat) <- 'broad_annot'
#' 
#' seurat_T.cells <- subset(seurat, idents = c("CD4_T_cells",
#'                                             "CD8_T_cells"))
#' 
#' 
#' seurat_T.cells$SRS_comparison = factor(seurat_T.cells$SRS_comparison, levels = c("SRS1","non_SRS1","Sepsis_conv","HV", "CS"))
#' 
#' 
#' #'@_##_plot-with-Seurat-function 
#' 
#' FeaturePlot(object = seurat_T.cells, features = c("CD69","PDCD1"), col=pal, order = T,
#'             split.by = "SRS_comparison", raster=FALSE) & NoAxes()
#' 
#' 
#' 
#' cell_typeA_marker_gene_list <- list(c("PDCD1", "CD69"))
#' object <- AddModuleScore(object = seurat_T.cells, features = cell_typeA_marker_gene_list, name = "cell_typeA_score")
#' FeaturePlot(object = object, features = "cell_typeA_score1", col=pal, order = T,
#'             split.by = "SRS_comparison", raster=FALSE) & NoLegend() & NoAxes()
#' 
#' #'@_##_plot-with-FeaturePlot_scCustom-from-scCustomize-package
#' ###
#' scCustomize::FeaturePlot_scCustom(seurat_object = seurat_T.cells, features = "PDCD1", 
#'                      split.by = "SRS_comparison",raster=FALSE) & NoAxes()
#' 
#' 
#' ggsave("FCGR1A.png", width = 10, height = 3)
#' 
#' #### AddModuleScore
#' cell_typeA_marker_gene_list <- list(c("TNF","PDCD1"))
#' object <- AddModuleScore(object = seurat_T.cells, features = cell_typeA_marker_gene_list, name = "CD11b_CD64")
#' scCustomize::FeaturePlot_scCustom(seurat_object = object, features = "CD11b_CD641", col=pal, order = T,
#'             split.by = "SRS_comparison",raster=FALSE, na_cutoff = 0.1) & NoAxes()
#' ggsave("1.pdf", width = 10, height = 3)
#' # 
#' # ####
#' # cell_typeA_marker_gene_list <- list(c("ITGAM", "IL3RA"))
#' # object <- AddModuleScore(object = seurat, features = cell_typeA_marker_gene_list, name = "CD11b_CD123")
#' # scCustomize::FeaturePlot_scCustom(seurat_object = object, features = "CD11b_CD1231", col=pal, order = T,
#' #                      split.by = "SRS_comparison",raster=FALSE, na_cutoff = 1) & NoAxes()
#' # 
#' # ggsave("PD-L1-CD10.pdf", width = 10, height = 3)
#' # 
#' # ####
#' # cell_typeA_marker_gene_list <- list(c("CD69", "PDCD1"))
#' # object <- AddModuleScore(object = seurat, features = cell_typeA_marker_gene_list, name = "CD69_PD1")
#' # scCustomize::FeaturePlot_scCustom(seurat_object = object, features = "CD69_PD11", col=pal, 
#' #                      split.by = "SRS_comparison",raster=FALSE, na_cutoff = 1) & NoAxes()
#' # 
#' # ggsave("3.png", width = 10, height = 3)
#' 
#' 
#' 
#' 
#' 
#' #'@_#_Create_Plots_voilin
#' library(ggpubr)
#' my_comparisons <- list( c("SRS1", "HV"))
#' 
#' cell_typeA_marker_gene_list <- list(c("FCGR1A", "IL3RA"))
#' object <- AddModuleScore(object = seurat, features = cell_typeA_marker_gene_list, name = "CD64_CD123")
#' p <- VlnPlot(object = seurat, features = c("MTOR"),  raster=FALSE, 
#'                              split.by = "SRS_comparison")
#' p
#' 
#' data = as.data.frame(p$data)
#' colnames(data)
#' 1.1# stats
#' 
#' stat.test <- data %>%
#'   group_by(ident) %>%
#'   rstatix::wilcox_test(MTOR ~ split, ref.group = "SRS1", p.adjust.method = "BH")
#' stat.test$sig = ifelse(stat.test$p.adj< 0.001, "***", ifelse(stat.test$p.adj < 0.01, "**", ifelse(stat.test$p.adj < 0.05, "*", "ns")))
#' stat.test
#' 1.1 # plot
#' stat.test <- stat.test %>% 
#'   rstatix::add_xy_position(x = "ident", dodge = 0.8)
#' bxp <- ggviolin(data, x = "ident", y = "MTOR", draw_quantiles = 0.5,
#'                 color = "split",  palette = "Dark2", size = 0.5,
#'                  x.text.angle = 30) 
#' bxp + stat_pvalue_manual(stat.test, label = "sig", tip.length = 0.01) + 
#'   ylab("Expression level CD64/CD123") + xlab("")# 10X5
#' 
#' 
#' 
#' 
#' # seurat@meta.data -> metadata
#' # umap_coord <- Embeddings(seurat, reduction = "umap") %>% as.data.frame()
#' # 
#' # metadata$UMAP_1 <- umap_coord$UMAP_1
#' # metadata$UMAP_2 <- umap_coord$UMAP_2
#' # 
#' # metadata %>%
#' #   ggplot(aes(x = UMAP_1, y = UMAP_2)) +
#' #   geom_point(aes(color = pct_counts_mito), size = 0.1) +
#' #   theme_classic() +
#' #   scale_color_viridis(option ="B")





