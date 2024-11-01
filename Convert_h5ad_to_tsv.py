#!/usr/bin/python

import scanpy as sc
import pandas as pd

# 1. .h5ad 파일 로드
adata_path = "AML_GSE116256.h5ad"
adata = sc.read_h5ad(adata_path)

# 2. 유전자 발현 데이터 및 유전자 심볼, 세포 ID 추출
# .h5ad 파일의 세포 ID는 adata.obs_names, 유전자 ID는 adata.var_names에 저장됨
# adata.X는 세포 x 유전자 형식의 행렬 데이터로, sparse 행렬일 경우 dense로 변환
expression_data = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)

# 3. 원하는 포맷으로 재구성
# 첫 번째 열에 'SYMBOL' 열 이름을 추가하고, 나머지 열은 세포 ID로 설정
# gene names are set as row names
expression_data.insert(0, "SYMBOL", expression_data.index)

# 4. 저장할 TSV 파일 경로 설정
output_path = "AML_GSE116256_expression_data.tsv"
expression_data.to_csv(output_path, sep='\t', index=False)

print(f"데이터가 {output_path}로 성공적으로 저장되었습니다.")
