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
feature_plot <- scCustomize::FeaturePlot_scCustom(seurat_object = SRS.obj, features = "CD274", 
                                  split.by = "SRS_comparison",raster=FALSE) & NoAxes() & theme(legend.position = "bottom")
ggsave(filename = "feature_plot.png", plot = feature_plot, width = 4, height = 3, dpi = 400)
