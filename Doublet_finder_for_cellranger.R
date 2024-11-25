#!/usr/bin/Rscript

library(Seurat)
library(dplyr)
library(SeuratDisk)
library(SeuratWrappers)
library(DoubletFinder)

options(Seurat.object.assay.version = "v3")

set.seed(123)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: script.R <sample_list.txt>")
}

sample_list_file <- args[1]

# sample_list_file
#/path/to/filtered_feature_bc_matrix1
#/path/to/filtered_feature_bc_matrix2
#/path/to/filtered_feature_bc_matrix3
.
.
.


# Sample 리스트 파일 읽기
sample_list <- readLines(sample_list_file) 
sample_list <- trimws(sample_list) # 공백 제거

# Loop를 통해 각 샘플 처리
for (sample in sample_list) {
  cat("Processing sample:", sample, "\n")
  tryCatch({
    # 10X 데이터 읽기
    data_path <- sample
    if (!file.exists(data_path)) stop("Directory provided does not exist")
    seurat_obj <- Read10X(data.dir = data_path)
    
    # Seurat 객체 생성
    seurat_obj <- CreateSeuratObject(counts = seurat_obj)
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
    
    ## DoubletFinder 관련 전처리
    # Normalized data 가져오기
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
    seurat_obj <- ScaleData(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj)
    
    # PCA 기준 계산
    stdv <- seurat_obj[["pca"]]@stdev
    sum.stdv <- sum(stdv)
    percent.stdv <- (stdv / sum.stdv) * 100
    cumulative <- cumsum(percent.stdv)
    co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
    co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] -
                        percent.stdv[2:length(percent.stdv)]) > 0.1),
               decreasing = TRUE)[1] + 1
    min.pc <- min(co1, co2)
    
    # UMAP 및 클러스터링
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:min.pc)
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:min.pc)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.1)
    
    ## pK Identification (DoubletFinder)
    sweep.res.list <- paramSweep_v3(seurat_obj, PCs = 1:min.pc, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric), ]
    optimal.pk <- as.numeric(levels(bcmvn.max$pK))[bcmvn.max$pK]
    
    ## Homotypic doublet proportion estimate
    annotations <- seurat_obj@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)
    nExp.poi <- round(optimal.pk * nrow(seurat_obj@meta.data))
    nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
    
    # DoubletFinder 실행
    seurat_obj <- doubletFinder_v3(seu = seurat_obj,
                                   PCs = 1:min.pc,
                                   pK = optimal.pk,
                                   nExp = nExp.poi.adj)
    
    # Meta data 업데이트
    metadata <- seurat_obj@meta.data
    colnames(metadata)[grep("DF.classifications", colnames(metadata))] <- "doublet_finder"
    seurat_obj <- AddMetaData(seurat_obj, metadata$doublet_finder, col.name = "doublet_finder")
    
    # Singlet만 추출
    seurat_obj.filtered <- subset(seurat_obj, doublet_finder == "Singlet")
    
    # 결과 저장
    saveRDS(seurat_obj.filtered, file = paste0(sample, "_filtered_singlets.rds"))
    cat("Sample processed successfully:", sample, "\n")
  }, error = function(e) {
    cat("Error processing sample", sample, ":", e$message, "\n")
  })
}

cat("Processing completed.\n")
