
#1. 데이터 로드 및 전처리
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(zellkonverter)
library(dplyr)

# Batch Correction using Seurat

# H5AD 데이터 로드
normal_BM_file_path <- "/data/workbench/scRSEQ_AML/exdata/CancerFinder/Stemcell_DEG/DISCO_BM_downsampled.h5ad"
aml_BM_file_path <- "/data/workbench/scRSEQ_AML/exdata/CancerFinder/Stemcell_DEG/AML_GSE116256_BM_converted.h5ad"

normal_sce_obj <- readH5AD(normal_BM_file_path)
aml_sce_obj <- readH5AD(aml_BM_file_path)


# 유전자, 셀 방향 확인
dim(normal_sce_obj)          # 확인: 행 = 세포, 열 = 유전자
rownames(normal_sce_obj)     # 유전자 이름
colnames(normal_sce_obj)     # 세포 이름 


dim(aml_sce_obj)          # 확인: 행 = 세포, 열 = 유전자
rownames(aml_sce_obj)     # 유전자 이름
colnames(aml_sce_obj)     # 세포 이름


# SingleCellExperiment → Seurat 변환
# RNA 표현 데이터가 적절히 설정되어 있는지 확인 후 변환
normal_seurat_obj <- as.Seurat(normal_sce_obj, counts = "X", data = NULL)
aml_seurat_obj <- as.Seurat(aml_sce_obj, counts = "X", data = NULL)

# 병합 전에 orig.ident 설정
normal_seurat_obj$orig.ident <- "Normal"
aml_seurat_obj$orig.ident <- "AML"

# Merge datasets
combined <- merge(normal_seurat_obj, aml_seurat_obj, add.cell.ids = c("Normal", "AML"))

# 병합 후 orig.ident 확인
table(combined$orig.ident)

# 배치 정보 추가
combined$batch <- ifelse(combined$orig.ident == "Normal", "Batch1", "Batch2")

# 결과 확인
table(combined$orig.ident)
table(combined$batch)

# Normalize data
combined <- NormalizeData(combined)

# Find variable features for integration
combined <- FindVariableFeatures(combined)

# Perform batch correction
combined <- ScaleData(combined, vars.to.regress = "batch")  # Replace "batch" with the actual batch identifier column


# UMAP Clustering
# Run PCA for dimensionality reduction
combined <- RunPCA(combined, npcs = 30)

# Run UMAP for visualization
combined <- RunUMAP(combined, dims = 1:30)

# Perform clustering
combined <- FindNeighbors(combined, dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)  # Adjust resolution as needed

# Visualize UMAP
DimPlot(combined, reduction = "umap", group.by = "ident")  # Default clustering visualization


# Find markers for each cluster
markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


# Path to marker gene files
marker_dir <- "/data/workbench/scRSEQ_AML/exdata/CancerFinder/Stemcell_DEG/CellMarker2.0/stemcells"

# List all marker gene files
marker_files <- list.files(marker_dir, pattern = "*.csv", full.names = TRUE)

# Read marker files into a named list
marker_gene_lists <- lapply(marker_files, function(file) {
  df <- read.csv(file, stringsAsFactors = FALSE, sep = '\t')
  return(df$marker)
})

names(marker_gene_lists) <- basename(marker_files)

# Annotate clusters based on marker genes
annotate_clusters <- function(markers, cluster_markers) {
  annotated_clusters <- list()
  
  for (cluster in unique(cluster_markers$cluster)) {
    cluster_genes <- cluster_markers %>%
      filter(cluster == !!cluster) %>%
      pull(gene)
    
    for (marker_file in names(markers)) {
      matched_genes <- intersect(cluster_genes, markers[[marker_file]])
      
      if (length(matched_genes) > 0) {
        if (!marker_file %in% names(annotated_clusters)) {
          annotated_clusters[[marker_file]] <- list()
        }
        annotated_clusters[[marker_file]][[as.character(cluster)]] <- matched_genes
      }
    }
  }
  return(annotated_clusters)
}

# Run annotation
marker_annotations <- annotate_clusters(marker_gene_lists, markers)

# Save cell IDs and annotations to files
output_dir <- "/data/workbench/scRSEQ_AML/exdata/CancerFinder/Stemcell_DEG/Cluster_Annotations"
dir.create(output_dir, showWarnings = FALSE)


for (marker_file in names(marker_annotations)) {
  for (cluster in names(marker_annotations[[marker_file]])) {
    cluster_cells <- WhichCells(combined, idents = as.numeric(cluster))
    output_path <- file.path(output_dir, paste0(gsub(".csv", "", marker_file), "_Cluster", cluster, "_Cells.txt"))
    write.table(cluster_cells, file = output_path, row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}
print("Cluster annotation and cell ID extraction completed!")
