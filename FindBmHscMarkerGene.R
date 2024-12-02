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
