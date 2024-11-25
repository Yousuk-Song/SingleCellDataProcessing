#!/usr/bin/Rscript

library(Seurat)
library(dplyr)
library(stringr)
library(DoubletFinder)

set.seed(123)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: script.R <sample_list.txt>")
}

sample_list_file <- args[1]
samples <- readLines(sample_list_file)
samples <- str_trim(samples)

#setwd("/data/processed_data/scRSEQ_AML/metastasis/")

# 함수 정의
process_sample <- function(sample_name) {
  tryCatch({
    # 데이터 읽기
    file_path <- sample_name
    
    if (!file.exists(file_path)) {
      stop("File does not exist")
    }
    
    # 데이터 로드 (.txt.gz 형식)
    sample_data <- read.table(gzfile(file_path), header = TRUE, row.names = 1)
    
    # Seurat 객체 생성
    if (packageVersion("Seurat") >= "5.0.0") {
      sample <- CreateSeuratObject(counts = as.matrix(sample_data))
    } else {
      sample <- CreateSeuratObject(raw.data = sample_data)
    }
    
    sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-")
    
    # DoubletFinder 전처리
    sample <- NormalizeData(sample)
    sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)
    sample <- ScaleData(sample)
    sample <- RunPCA(sample)
    
    # PCA 분석 후 PC 수 계산
    stdv <- sample[["pca"]]@stdev
    sum.stdv <- sum(stdv)
    percent.stdv <- (stdv / sum.stdv) * 100
    cumulative <- cumsum(percent.stdv)
    co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
    co2 <- sort(which((percent.stdv[1:(length(percent.stdv) - 1)] -
                         percent.stdv[2:length(percent.stdv)]) > 0.1),
                decreasing = TRUE)[1] + 1
    min.pc <- min(co1, co2)
    
    sample <- RunUMAP(sample, dims = 1:min.pc)
    sample <- FindNeighbors(sample, dims = 1:min.pc)
    sample <- FindClusters(sample, resolution = 0.1)
    
    # DoubletFinder 실행
    sample <- doubletFinder(sample, PCs = 1:min.pc, pK = 0.05, nExp = 10)
    
    # Singlet 추출 후 저장
    singlets <- subset(sample, subset = doublet_finder == "Singlet")
    saveRDS(singlets, file = paste0(sample_name, ".filtered.singlets.rds"))
    cat("Processed:", sample_name, "\n")
    
  }, error = function(e) {
    cat("Error processing sample", sample_name, ":", e$message, "\n")
  })
}

# 모든 샘플 처리
for (sample_name in samples) {
  process_sample(sample_name)
}

cat("Processing completed.\n")

