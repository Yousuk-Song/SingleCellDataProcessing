#!/usr/bin/python

import pandas as pd
import scanpy as sc

# 1. 메타데이터 파일 읽기 및 조건에 맞는 셀 ID 필터링
metadata_path = "/data/processed_data/scRSEQ_AML/TISCH/AML_GSE116256_CellMetainfo_table.tsv"
metadata = pd.read_csv(metadata_path, sep='\t')

# (a) 'HSC', 'Prog', 'GMP', 'ProMono', 'earlyEry' 포함, PredictionRefined가 'normal'인 셀 선택
stem_cells = metadata[
    (metadata['Celltype (original)'].isin(['HSC', 'Prog', 'GMP', 'ProMono', 'earlyEry'])) &
    (metadata['PredictionRefined'] == 'normal')
]['Cell']

# (b) 'HSC', 'Prog', 'GMP', 'ProMono', 'earlyEry' 제외, PredictionRefined가 'normal'인 셀 선택
non_stem_cells = metadata[
    (~metadata['Celltype (original)'].isin(['HSC', 'Prog', 'GMP', 'ProMono', 'earlyEry'])) &
    (metadata['PredictionRefined'] == 'normal')
]['Cell']


# 2. h5ad 파일에서 각각의 셀들만 필터링
adata_path = "/data/processed_data/scRSEQ_AML/train/AML_GSE116256.h5ad"
adata = sc.read_h5ad(adata_path)

# (a) 'HSC', 'Prog', 'GMP', 'ProMono', 'earlyEry' 셀 ID만 필터링
adata_stem_filtered = adata[:, adata.var_names.isin(stem_cells)]

# (b) 그 외 세포들만 필터링
adata_non_stem_filtered = adata[:, adata.var_names.isin(non_stem_cells)]


# 3. 필터링된 데이터 저장
output_path_stem = "/data/workbench/scRSEQ_AML/exdata/CancerFinder/data/train4/AML_GSE116256_StemCells.h5ad"
output_path_non_stem = "/data/workbench/scRSEQ_AML/exdata/CancerFinder/data/train5/AML_GSE116256_nonStemCell.h5ad"

adata_stem_filtered.write(output_path_stem)
adata_non_stem_filtered.write(output_path_non_stem)

print(f"필터링된 데이터가 {output_path_stem} 및 {output_path_non_stem}로 저장되었습니다.")
