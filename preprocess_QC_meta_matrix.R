#!/usr/bin/Rscript

library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(SeuratDisk)
##for seurat 4x version

# 파일 경로 및 매개변수 지정
args <- commandArgs(trailingOnly = TRUE)
file_path <- args[1]  # 첫 번째 인자로 파일 경로 받기

merged.obj<-readRDS(file_path)


merged.obj$sample_tag<-merged.obj$orig.ident
Idents(merged.obj)<-"sample_tag"
merged.obj[["percent.mt"]] <- PercentageFeatureSet(merged.obj, pattern = "^MT-")

# generate QC metric
merged.obj$log10GenesPerUMI <- log10(merged.obj$nFeature_RNA) / log10(merged.obj$nCount_RNA)
metadata <- merged.obj@meta.data
metadata$cells <- rownames(metadata)
metadata <- metadata %>%
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=sample_tag, x=nUMI, fill= sample_tag)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = c(200,500))


# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=sample_tag, x=nGene, fill= sample_tag)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = c(200,4000,6000))

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=percent.mt)) + 
  geom_point(cex=0.5) + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 200) +
  geom_hline(yintercept = 200) +
  facet_wrap(~sample_tag)

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample_tag, x=percent.mt, fill=sample_tag)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = c(10,30))

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample_tag, fill=sample_tag)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.80)


# Filter out low quality reads using selected thresholds - these will change with experiment
merged.obj2 <- subset(x = merged.obj, 
                          subset= (nCount_RNA >= 200) & 
                            (nFeature_RNA >= 200) & 
                            (nFeature_RNA <= 6000) & 
                            (percent.mt < 10)&
                            (log10GenesPerUMI>0.8))

# .rds를 .QC_pass.rds로 변경
out_path <- sub("\\.rds$", ".QC_pass.rds", file_path)

saveRDS(merged.obj2, out_path)


##Re-assess QC metrics
# generate QC metric
merged.obj2$log10GenesPerUMI <- log10(merged.obj2$nFeature_RNA) / log10(merged.obj2$nCount_RNA)
metadata <- merged.obj2@meta.data
metadata$cells <- rownames(metadata)
metadata <- metadata %>%
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)
# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=sample_tag, x=nUMI, fill= sample_tag)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 200)

# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=sample_tag, x=nGene, fill= sample_tag)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 200)

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=percent.mt)) + 
  geom_point(cex=0.5) + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 200) +
  geom_hline(yintercept = 200) +
  facet_wrap(~sample_tag)

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample_tag, x=percent.mt, fill=sample_tag)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 25)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample_tag, fill=sample_tag)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
