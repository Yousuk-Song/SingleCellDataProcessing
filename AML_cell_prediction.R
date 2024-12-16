# 1. 데이터 로드 및 전처리
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(zellkonverter)
library(dplyr)
library(pheatmap)
library(dittoSeq)


# 데이터 로드
#load_csv_to_seurat <- function(file_path) {
#  data <- read.csv(file_path, row.names = 1)
#  seurat_obj <- CreateSeuratObject(counts = data)
#  return(seurat_obj)
#}
#aml_seurat_obj <- lapply(aml_files, load_csv_to_seurat)
#normal_seurat_obj <- lapply(bm_files, load_csv_to_seurat)


# 파일 경로
file_path <- "/data/workbench/scRSEQ_AML/exdata/CancerFinder/Stemcell_DEG/Data/GSE130756"

# Normal과 AML 파일 목록
normal_files <- c("GSM4135770_N01_dge.txt.gz", "GSM4135771_N02_dge.txt.gz", "GSM4135772_N03_dge.txt.gz")
aml_files <- setdiff(list.files(file_path, pattern = "dge.txt.gz"), normal_files)

# Function to load data from txt.gz files
load_data <- function(file) {
  file_full_path <- file.path(file_path, file)
  data <- read.table(file_full_path, header = TRUE, row.names = 1)
  return(data)
}
# Normal 파일 로드
normal_list <- lapply(normal_files, function(f) {
  data <- load_data(f)
  seurat_obj <- CreateSeuratObject(counts = data, project = f, min.cells = 3, min.features = 200)
  seurat_obj@meta.data$batch <- f  # 배치 정보로 파일명 추가
  seurat_obj@meta.data$group <- "Normal"
  return(seurat_obj)
})

# AML 파일 로드
aml_list <- lapply(aml_files, function(f) {
  data <- load_data(f)
  seurat_obj <- CreateSeuratObject(counts = data, project = f, min.cells = 3, min.features = 200)
  seurat_obj@meta.data$batch <- f  # 배치 정보로 파일명 추가
  seurat_obj@meta.data$group <- "AML"
  return(seurat_obj)
})

aml_seurat_obj <- merge(x = aml_list[[1]], y = aml_list[-1], project = "AML")
normal_seurat_obj <- merge(x = normal_list[[1]], y = normal_list[-1], project = "Normal")


# Normal 파일 로드
#normal_list <- lapply(normal_files, function(f) {
#  data <- load_data(f)
#  CreateSeuratObject(counts = data, project = "Normal", min.cells = 3, min.features = 200)
#})
# AML 파일 로드
#aml_list <- lapply(aml_files, function(f) {
#  data <- load_data(f)
#  CreateSeuratObject(counts = data, project = "AML", min.cells = 3, min.features = 200)
#})


# 결과 출력
print(normal_seurat_obj)
print(aml_seurat_obj)


# 정상과 AML 데이터를 통합
combined <- merge(normal_seurat_obj, aml_seurat_obj)
# 원본 메타데이터 병합
meta_data_normal <- normal_seurat_obj@meta.data
meta_data_aml <- aml_seurat_obj@meta.data
# 병합된 메타데이터 재생성
combined_meta <- rbind(meta_data_normal, meta_data_aml)
combined@meta.data <- combined_meta
# 확인
head(combined@meta.data)
# 결과 확인
table(combined$batch)
table(combined$group)

# 1. Seurat 객체 리스트 생성
seurat_list <- SplitObject(combined, split.by = "batch")

# 2. 각 Seurat 객체에 대해 데이터 전처리 수행
seurat_list <- lapply(seurat_list, function(x) {
#  x <- NormalizeData(x)          # 정규화   ### 인풋 테이블이 이미 정규화 되어있어서 생략 
  x <- FindVariableFeatures(x)   # 가변 유전자 선택
#  x <- ScaleData(x)              # 데이터 스케일링
  x <- RunPCA(x)                 # PCA 실행
  return(x)
})

# 데이터 통합
combined <- IntegrateData(anchorset = integration_anchors)



#i) RPCA 기반 통합 앵커 찾기
integration_anchors <- FindIntegrationAnchors(
  object.list = seurat_list,      # Seurat 객체 리스트
  anchor.features = 3000,        # 통합에 사용할 가변 유전자 수 증가
  reduction = "rpca"             # RPCA 사용
)


# 기존 메타데이터의 batch 및 group 정보 다시 추가
combined$batch <- combined_meta$batch
combined$group <- combined_meta$group

# 메타데이터 확인
head(combined@meta.data)
table(combined$batch)
table(combined$group)

# 5. 통합 데이터 스케일링 
#combined <- ScaleData(combined)
combined <- ScaleData(combined, vars.to.regress = "batch")
# 6. PCA 및 UMAP 수행
combined <- RunPCA(combined)
combined <- RunUMAP(combined, dims = 1:30, n.neighbors = 50, min.dist = 0.1)

combined <- FindNeighbors(combined, dims = 1:30)
combined <- FindClusters(combined, resolution = 0.13)

 # 7. 결과 확인
# UMAP 시각화 
DimPlot(combined, reduction = "umap", group.by = "batch", raster=FALSE)

DimPlot(combined, reduction = "umap", group.by = "group", raster=FALSE)

# 클러스터링 결과 확인
DimPlot(combined, reduction = "umap", label = TRUE, raster=FALSE)

# 클러스터 별로 group(AML/Normal)의 비율 확인
table(combined$group, combined$seurat_clusters)

# 클러스터 마커 유전자 찾기
cluster_markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# 클러스터 별 상위 마커 확인
head(cluster_markers)




## ii) Harmony

library(harmony)

# 배치 정보 확인
head(combined@meta.data)

# 배치 정보가 있는지 확인 (batch 컬럼)
if (!"batch" %in% colnames(combined@meta.data)) {
  stop("Batch information is missing in the meta data.")
}

# PCA 수행
combined <- NormalizeData(combined)     # 정규화 (필요시)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = 30)


# Harmony 수행
combined <- RunHarmony(
  object = combined,          # Seurat 객체
  group.by.vars = "batch",    # 배치 정보 컬럼
  dims = 1:30                 # PCA 차원 수 지정
)

# Harmony 기반 UMAP 생성
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30)

# 최근접 이웃 찾기 및 클러스터링
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)

# 배치 제거 전 (원본 PCA)
DimPlot(combined, reduction = "pca", group.by = "batch", label = TRUE) +
  ggtitle("PCA by Batch (Before Harmony)")

# 배치 제거 후 (Harmony 적용)
DimPlot(combined, reduction = "umap", group.by = "batch", label = TRUE) +
  ggtitle("UMAP by Batch (After Harmony)")

# Group 기준 UMAP 시각화
DimPlot(combined, reduction = "umap", group.by = "group", label = TRUE) +
  ggtitle("UMAP by Group (After Harmony)")

DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Clusters after Harmony Batch Correction")

# RNA 레이어 활성화
DefaultAssay(combined) <- "RNA"

# DEG 분석 수행
deg <- FindMarkers(
  combined, ident.1 = "AML", ident.2 = "Normal", group.by = "group"
)

# 결과 확인
head(deg)















### 분석 준비 완료!! ###

bm_seurat <- combined

# [1] 클러스터 1,3,4 를 그룹 1 cancer로, 0을 그룹 2 Normal로 지정
bm_seurat$group <- ifelse(bm_seurat$seurat_clusters %in% c(0, 1, 4), "AML",
                          ifelse(bm_seurat$seurat_clusters %in% c(2, 3, 6), "Normal", NA))

# 그룹 확인
table(bm_seurat$group)

# DEG 계산 (Wilcoxon Rank Sum Test 적용)
deg <- FindMarkers(
  object = bm_seurat,
  ident.1 = "AML",
  ident.2 = "Normal",
  group.by = "group",      # 그룹 정보를 기반으로 비교
  min.pct = 0.25,          # 최소 발현 비율 (25% 이상의 세포에서 발현)
  logfc.threshold = 0.25,  # 최소 Log2 Fold Change 기준
  test.use = "wilcox"      # Wilcoxon Rank Sum Test 사용
)

# DEG 결과 확인
head(deg)

# 유의미한 DEG 필터링 (p-value < 0.05)
significant_deg <- deg[deg$p_val_adj < 0.05, ]

# 유의미한 DEG 확인
head(significant_deg)

# 유의미한 DEG 수 확인
cat("Number of significant DEGs:", nrow(significant_deg), "\n")


## GO 분석 
aml_specific_genes <- significant_deg %>%
  filter(avg_log2FC > 0) %>%  # AML에서 발현이 높은 유전자만 선택
  rownames()

library(clusterProfiler)
library(org.Hs.eg.db) # 인간 유전자 주석 데이터베이스

# 유전자 ID 변환 (필요 시)
gene_ids <- bitr(aml_specific_genes, 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = org.Hs.eg.db)

# GO Enrichment 분석
go_results <- enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",           # Biological Process (BP)
  pvalueCutoff = 0.05,  # 유의성 기준
  qvalueCutoff = 0.2,   # FDR 기준
  readable = TRUE       # 결과를 읽기 쉽게 반환
)

# GO 결과 확인
head(go_results)

# Hallmarks of Cancer 관련 키워드
cancer_keywords <- c("apoptosis", "proliferation", "angiogenesis", "cell cycle", 
                     "DNA repair", "immune response", "metastasis", "tumor", 
                     "oncogenesis", "differentiation")

# cancer_keywords를 포함하는 항목만 필터링
cancer_pattern <- paste(cancer_keywords, collapse = "|")  # 정규식 패턴 생성

filtered_go <- go_results@result %>%
  filter(grepl(cancer_pattern, Description, ignore.case = TRUE))

# 필터링된 GO term 결과 확인
print(paste("Number of filtered GO terms:", nrow(filtered_go)))
head(filtered_go)

# 암 관련 GO Terms의 유전자
cancer_genes <- unique(unlist(strsplit(filtered_go$geneID, "/")))

# 결과 확인
head(cancer_genes)
cat("Number of cancer-related GO terms:", nrow(cancer_go), "\n")
cat("Number of cancer-related genes:", length(cancer_genes), "\n")

write.csv(cancer_genes, "/data/workbench/scRSEQ_AML/exdata/CancerFinder/Stemcell_DEG/Cancer_related_gene_set.csv", row.names = FALSE)

# GO 항목 점수 막대 그래프
library(ggplot2)

# 상위 10개의 GO term 선택
filtered_go_top10 <- filtered_go %>%
  arrange(p.adjust) %>%
  head(10)

# Barplot 생성 (글씨 크기 조정)
ggplot(filtered_go_top10, aes(x = reorder(Description, -p.adjust), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # 가로 막대 그래프
  labs(title = "Filtered GO Terms", 
       x = "GO Term", 
       y = "Gene Count") +
  theme_minimal() +
  theme(
    text = element_text(size = 16),          # 전체 텍스트 크기 설정
    axis.text.x = element_text(size = 14),   # x축 눈금 레이블 크기
    axis.text.y = element_text(size = 14),   # y축 눈금 레이블 크기
    axis.title.x = element_text(size = 16),  # x축 제목 크기
    axis.title.y = element_text(size = 16),  # y축 제목 크기
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5) # 제목 크기와 정렬
  )






## i) cancer score 

# 암세포 점수 계산
valid_genes <- intersect(cancer_genes, rownames(GetAssayData(bm_seurat, slot = "scale.data")))
valid_logfc <- significant_deg[valid_genes, "avg_log2FC"]

# 가중치 기반 cancer score 계산
cancer_score_weighted <- rowSums(
  t(GetAssayData(bm_seurat, slot = "scale.data")[valid_genes, ]) * valid_logfc
)
bm_seurat$cancer_score <- cancer_score_weighted
#bm_seurat$cancer_score <- colMeans(GetAssayData(bm_seurat, slot = "scale.data")[valid_genes, ])

# 암세포로 판별된 세포 데이터 저장
#aml_cancer_cells <- subset(bm_seurat@meta.data, orig.ident == "AML")
aml_cancer_cells <- subset(bm_seurat, subset = group == "AML")
head(aml_cancer_cells)

# 암세포와 정상세포 구분
# Cancer Score 분포 확인
library(pROC)

# ROC Curve 및 AUC 계산
roc_obj <- roc(bm_seurat@meta.data$group, bm_seurat$cancer_score)
plot(roc_obj, print.thres = TRUE, print.auc = TRUE)

# 최적 Threshold 찾기
best_threshold <- coords(roc_obj, "best", ret = "threshold")
print(best_threshold)

threshold <- best_threshold$threshold
bm_seurat$cancer_cell <- ifelse(bm_seurat$cancer_score > threshold, "AML", "Normal")

library(ggplot2)

ggplot(bm_seurat@meta.data, aes(x = cancer_score, fill = cancer_cell)) +
  geom_density(alpha = 0.5) +
  ggtitle("Cancer Score Distribution by Classification") +
  labs(x = "Cancer Score", fill = "Cell Type") +
  geom_vline(xintercept = threshold, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = threshold, y = max(density(bm_seurat@meta.data$cancer_score, adjust = 1)$y), 
           label = paste0("Threshold = ", threshold), color = "red", vjust = -0.5) +
  theme_minimal()

# 결과 확인
table(bm_seurat$cancer_cell)
# orig.ident를 정답으로 사용한 교차표 생성
confusion_matrix <- table(Truth = bm_seurat@meta.data$group, Prediction = bm_seurat$cancer_cell)
# 교차표 출력
print(confusion_matrix)

# 성능 지표 계산
true_positive <- confusion_matrix["AML", "AML"] # AML을 AML로 맞춘 수
false_positive <- confusion_matrix["Normal", "AML"] # Normal을 AML로 잘못 분류
true_negative <- confusion_matrix["Normal", "Normal"] # Normal을 Normal로 맞춘 수
false_negative <- confusion_matrix["AML", "Normal"] # AML을 Normal로 잘못 분류

accuracy <- (true_positive + true_negative) / sum(confusion_matrix)
sensitivity <- true_positive / (true_positive + false_negative) # 민감도
specificity <- true_negative / (true_negative + false_positive) # 특이도

# 결과 출력
cat("Accuracy:", accuracy, "\n")
cat("Sensitivity:", sensitivity, "\n")
cat("Specificity:", specificity, "\n")


# 잘못 분류된 세포 필터링
misclassified_cells <- bm_seurat@meta.data %>%
  filter(group != cancer_cell)

# 분류된 세포 필터링
misclassified_normal_as_cancer <- bm_seurat@meta.data %>%
  filter(group == "Normal" & cancer_cell == "AML")

misclassified_cancer_as_normal <- bm_seurat@meta.data %>%
  filter(group == "AML" & cancer_cell == "Normal")

correctly_classified_normal <- bm_seurat@meta.data %>%
  filter(group == "Normal" & cancer_cell == "Normal")

correctly_classified_cancer <- bm_seurat@meta.data %>%
  filter(group == "AML" & cancer_cell == "AML")

# UMAP 좌표 추가
umap_coordinates <- as.data.frame(Embeddings(bm_seurat, reduction = "umap"))
misclassified_normal_as_cancer <- cbind(misclassified_normal_as_cancer, umap_coordinates[rownames(misclassified_normal_as_cancer), ])
misclassified_cancer_as_normal <- cbind(misclassified_cancer_as_normal, umap_coordinates[rownames(misclassified_cancer_as_normal), ])
correctly_classified_normal <- cbind(correctly_classified_normal, umap_coordinates[rownames(correctly_classified_normal), ])
correctly_classified_cancer <- cbind(correctly_classified_cancer, umap_coordinates[rownames(correctly_classified_cancer), ])

# UMAP 시각화: Normal -> Cancer로 잘못 분류된 세포
#p1 <- DimPlot(bm_seurat, group.by = "cancer_cell", reduction = "umap", label = TRUE) +
# geom_point(data = misclassified_normal_as_cancer, aes(x = UMAP_1, y = UMAP_2), color = "blue", size = 0.5) +
#  ggtitle("Normal Misclassified as Cancer (Blue Points)")

# UMAP 시각화: Normal -> Cancer로 잘못 분류된 세포
#p1 <- DimPlot(bm_seurat, group.by = "orig.ident", reduction = "umap", label = TRUE) +
#  geom_point(data = misclassified_normal_as_cancer, aes(x = UMAP_1, y = UMAP_2), color = "blue", size = 0.5) +
#  geom_point(data = correctly_classified_normal, aes(x = UMAP_1, y = UMAP_2), color = "orange", size = 0.5) +
#  ggtitle("Normal Misclassified as Cancer (Blue Points) and Correctly Classified Normal (Orange Points)") +
#  theme_minimal()

# UMAP 시각화: Cancer -> Normal로 잘못 분류된 세포
p2 <- DimPlot(bm_seurat, group.by = "group", reduction = "umap", label = TRUE) +
  geom_point(data = correctly_classified_cancer, aes(x = UMAP_1, y = UMAP_2), color = "green", size = 0.5) +
  ggtitle("Correctly Classified Cancer (Green Points)") +
  theme_minimal()

# 결과 출력
#print(p1)
print(p2)

# orig.ident와 cancer_cell을 각각 색상으로 표시하여 두 그룹 간의 차이를 시각화합니다 -> 두 그림을 비교하면 분류 모델과 실제 데이터의 차이를 확인할 수 있습니다.
#DimPlot(bm_seurat, group.by = "orig.ident", reduction = "umap") + ggtitle("Original Identity")
#DimPlot(bm_seurat, group.by = "cancer_cell", reduction = "umap") + ggtitle("Predicted Identity")



























# ii) random forest 
library(randomForest)
library(caret)
library(ggplot2)

# 1. 유효한 유전자 선택
#valid_genes <- intersect(cancer_genes, rownames(GetAssayData(bm_seurat, slot = "scale.data")))

# 2. Feature Matrix 생성
feature_matrix <- as.data.frame(t(GetAssayData(bm_seurat, slot = "scale.data")[valid_genes, ]))

feature_matrix$cancer_binary <- ifelse(bm_seurat$orig.ident == "AML", 1, 0) # AML: 1, Normal: 0

# 3. 학습 데이터와 테스트 데이터로 분할
set.seed(123)
train_index <- createDataPartition(feature_matrix$cancer_binary, p = 0.8, list = FALSE)
train_data <- feature_matrix[train_index, ]

common_cols <- intersect(colnames(feature_matrix), colnames(train_data))
train_data <- train_data[, common_cols, drop = FALSE]

test_data <- feature_matrix[-train_index, ]

# 반응 변수를 factor로 변환
train_data$cancer_binary <- as.factor(train_data$cancer_binary)
test_data$cancer_binary <- as.factor(test_data$cancer_binary)

# Feature Matrix 준비
X <- train_data[, colnames(train_data) != "cancer_binary"]  # 독립 변수
y <- train_data$cancer_binary  # 종속 변수

#모델 학습 RandomForest 모델에 입력되는 데이터가 열 이름을 인식하지 못하는 경우가 있습니다. 
#이를 방지하기 위해 model.matrix를 사용해 수동으로 독립 변수를 정의해줍니다.
# Random Forest 모델
#rf_model <- randomForest(
#  cancer_binary ~ ., 
#  data = train_data, 
#  ntree = 500,            # 트리 수
#  mtry = sqrt(ncol(train_data) - 1), # 분할 기준 변수 수
#  importance = TRUE       # 변수 중요도 계산
#)
rf_model <- randomForest(
  x = X,
  y = y,
  ntree = 500,
  mtry = sqrt(ncol(X)),
  sampsize = c(50, 50),  # 클래스별 샘플 크기
  importance = TRUE
)

print(rf_model)

#rf_model <- randomForest(
#  x = X,
#  y = y,
#  ntree = 500,
#  mtry = sqrt(ncol(X)),
#  importance = TRUE
#)

# 모델 성능 확인
print(rf_model)


control <- trainControl(method = "cv", number = 5)  # 5-Fold Cross Validation
rf_cv_model <- train(
  x = X,
  y = Y,
  method = "rf", 
  trControl = control,
  importance = TRUE
)
print(rf_cv_model)





# 2. 모델 성능 평가
rf_predictions <- predict(rf_model, test_data)
confusion_matrix <- confusionMatrix(as.factor(rf_predictions), as.factor(test_data$cancer_binary))
print(confusion_matrix)

  # 3. 중요 유전자 확인
importance <- as.data.frame(importance(rf_model))
importance <- importance[order(-importance$MeanDecreaseGini), ]
head(importance)

# 1. UMAP 좌표에 예측 결과 추가
bm_seurat$rf_prediction <- ifelse(predict(rf_model, feature_matrix) == 1, "AML", "Normal")

# 2. 예측 결과 시각화
DimPlot(bm_seurat, group.by = "rf_prediction", reduction = "umap", label = TRUE) +
  ggtitle("Random Forest Prediction on UMAP")


# Grid Search 설정
rf_grid <- expand.grid(
  mtry = c(sqrt(ncol(train_data) - 1), log2(ncol(train_data) - 1)), 
  ntree = c(300, 500, 700)
)

# 모델 학습 및 평가 반복
grid_results <- lapply(1:nrow(rf_grid), function(i) {
  params <- rf_grid[i, ]
  model <- randomForest(
    cancer_binary ~ ., 
    data = train_data, 
    mtry = params$mtry, 
    ntree = params$ntree
  )
  pred <- predict(model, test_data)
  acc <- mean(pred == test_data$cancer_binary)
  return(cbind(params, Accuracy = acc))
})

grid_results <- do.call(rbind, grid_results)
print(grid_results)
  

# Random Forest 모델 재훈련
rf_model <- randomForest(
  cancer_binary ~ ., 
  data = train_data, 
  ntree = 300,            # 트리 수
  mtry = 1.414214, # 분할 기준 변수 수
  importance = TRUE       # 변수 중요도 계산
)
















## RF Model을 사용해서 Prediction
library(Seurat)
library(dplyr)
library(randomForest)
library(ggplot2)
library(zellkonverter)

try_readH5AD <- function(file_path) {
  tryCatch({
    readH5AD(file_path)
  }, error = function(e) {
    cat("Error reading file:", file_path, "\n", e, "\n")
    return(NULL)
  })
}

# .h5ad 파일 로드 (aml bm, healthy bm)
aml_patient_file_path <- "/data/workbench/scRSEQ_AML/exdata/CancerFinder/Stemcell_DEG/Data/GSE116256/AML_GSE116256_BM_converted.h5ad"
healthy_BM_file_path <- "/data/workbench/scRSEQ_AML/exdata/CancerFinder/Stemcell_DEG/Data/DISCO/DISCO_BM_downsampled_20000.h5ad"

aml_patient_sce_obj<- try_readH5AD(aml_patient_file_path)
healthy_BM_sce_obj<- try_readH5AD(healthy_BM_file_path)

# SingleCellExperiment → Seurat 변환
aml_patient_seurat_obj <- as.Seurat(aml_patient_sce_obj, counts = "X", data = NULL)
healthy_BM_seurat_obj <- as.Seurat(healthy_BM_sce_obj, counts = "X", data = NULL)

aml_patient_seurat_obj$orig.ident <- "AML"
healthy_BM_seurat_obj$orig.ident <- "Normal"
combined_bm <- merge(healthy_BM_seurat_obj, aml_patient_seurat_obj, add.cell.ids = c("Normal", "AML"))


# 데이터 슬롯에서 유전자 발현 데이터 가져오기
patient_bm_gene_data <- GetAssayData(combined, assay = "originalexp", slot = "counts")[valid_genes, ]

# Transpose and convert to data.frame
patient_bm_df <- as.data.frame(t(patient_bm_gene_data))

# 랜덤 포레스트 모델로 암세포 예측
patient_bm_predictions <- predict(rf_model, patient_bm_df)

# 예측 결과를 메타데이터에 추가
combined_bm@meta.data$cancer_prediction <- patient_bm_predictions

# 결과 확인
table(combined_bm@meta.data$cancer_prediction)

# UMAP 계산 (필요한 경우)
combined_bm <- RunUMAP(combined_bm, dims = 1:30)

# 계산 후 사용 가능한 reduction 확인
Reductions(combined_bm)

# UMAP 시각화: 예측 결과
p <- DimPlot(combined_bm, reduction = "umap", group.by = "cancer_prediction") +
  ggtitle("Predicted Cancer Cell Classification") +
  theme_minimal()

#print(p)








