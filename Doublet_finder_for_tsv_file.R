#!/usr/bin/Rscript

library(Seurat)
library(DoubletFinder)
library(dplyr)

# 파일 경로 지정
args <- commandArgs(trailingOnly = TRUE)
file_path <- args[1]  # 첫 번째 인자로 파일 경로 받기

# 데이터를 읽어옴 (구분자가 탭이라 가정)
data <- read.table(file_path, header = TRUE, row.names = 1, sep = "\t")

# Seurat 객체 생성
seurat_obj <- CreateSeuratObject(counts = data)

# counts 슬롯이 정상적으로 할당되었는지 확인하기
seurat_obj[["RNA"]] <- CreateAssayObject(counts = data)  # counts 슬롯 명시적으로 추가

# Seurat 객체 정보 확인
print(seurat_obj)

# 1. 전처리: MT 비율 계산 및 추가
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# 2. 데이터 정규화 및 변수 유전자 찾기
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# 3. 데이터 스케일링 및 PCA
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)

# 4. PCA 기준 계산
stdv <- seurat_obj[["pca"]]@stdev
sum.stdv <- sum(stdv)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                    percent.stdv[2:length(percent.stdv)]) > 0.1),
             decreasing = TRUE)[1] + 1
min.pc <- min(co1, co2)

# 5. UMAP 및 클러스터링
seurat_obj <- RunUMAP(seurat_obj, dims = 1:min.pc)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:min.pc)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.1)

# 6. pK Identification (DoubletFinder)
DefaultAssay(seurat_obj) <- "RNA"  # RNA 어세이 사용
sweep.res.list <- paramSweep_v3(seurat_obj, PCs = 1:min.pc, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric), ]
optimal.pk <- as.numeric(levels(bcmvn.max$pK))[bcmvn.max$pK]

# 7. Homotypic doublet proportion estimate
annotations <- seurat_obj@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp.poi <- round(optimal.pk * nrow(seurat_obj@meta.data))
nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))

# 8. DoubletFinder 실행
seurat_obj <- doubletFinder_v3(seu = seurat_obj,
                               PCs = 1:min.pc,
                               pK = optimal.pk,
                               nExp = nExp.poi.adj)

# 9. Meta data 업데이트
metadata <- seurat_obj@meta.data
colnames(metadata)[grep("DF.classifications", colnames(metadata))] <- "doublet_finder"
seurat_obj <- AddMetaData(seurat_obj, metadata$doublet_finder, col.name = "doublet_finder")

# 10. Singlet만 필터링
seurat_obj.filtered <- subset(seurat_obj, doublet_finder == "Singlet")

# 11. 결과 저장
saveRDS(seurat_obj.filtered, file = "filtered_singlets.rds")

cat("Processing completed.\n")
