#!/usr/bin/python

import pysam
import os
import sys
import csv

def find_and_extract_reads(bam, chrom, pos, variant, reference, gene):
    # BAM 파일 열기
    bamfile = pysam.AlignmentFile(bam, "rb")

    # 출력 파일 이름 설정
    variant_bam = bam.replace('.bam', f'.{chrom}.{pos}.{variant}.{gene}.var.bam')
    reference_bam = bam.replace('.bam', f'.{chrom}.{pos}.{reference}.{gene}.ref.bam')

    variant_bamfile = pysam.AlignmentFile(variant_bam, "wb", template=bamfile)
    reference_bamfile = pysam.AlignmentFile(reference_bam, "wb", template=bamfile)

    # 특정 위치의 reads 처리
    for pileupcolumn in bamfile.pileup(chrom, pos - 1, pos):
        if pileupcolumn.pos == pos - 1:  # 0-based 위치이므로 -1 해줌
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:  # deletion과 refskip 처리
                    base = pileupread.alignment.query_sequence[pileupread.query_position]
                    if base == variant:  # 변이 염기
                        variant_bamfile.write(pileupread.alignment)
                    elif base == reference:  # 기준 염기
                        reference_bamfile.write(pileupread.alignment)

    # BAM 파일 닫기
    bamfile.close()
    variant_bamfile.close()
    reference_bamfile.close()

    # BAM 파일 인덱싱
    os.system(f'samtools index {variant_bam}')
    os.system(f'samtools index {reference_bam}')

    # 읽은 정보 추출 및 저장
    extract_read_info(variant_bam, f'{variant_bam}.cell_barcodes.umi.tsv')
    extract_read_info(reference_bam, f'{reference_bam}.cell_barcodes.umi.tsv')

def extract_read_info(bam_file, output_file):
    # BAM 파일 열기
    samfile = pysam.AlignmentFile(bam_file, "rb")

    # 결과를 텍스트 파일로 저장
    with open(output_file, 'w') as f:
        # 파일의 첫 번째 줄에 헤더 추가
        f.write("ReadID\tCellBarcode\tUMI\n")

        # 모든 read를 순회하며 정보 추출
        for read in samfile:
            read_id = read.query_name

            # Cell barcode 추출
            cell_barcode = read.get_tag('CB') if read.has_tag('CB') else 'NA'

            # UMI 추출
            if read.has_tag('im'):
                umis_raw = read.get_tag('im')
                umis_list = umis_raw.split(',')
                umi = ','.join(umis_list)  # UMIs를 쉼표로 구분하여 문자열로 변환
            else:
                umi = 'NA'

            # 정보를 파일에 작성
            f.write(f"{read_id}\t{cell_barcode}\t{umi}\n")

    # BAM 파일 닫기
    samfile.close()

if __name__ == "__main__":
    bam = sys.argv[1]  # 입력 BAM 파일
    mutation_csv = sys.argv[2]  # 변이 정보 CSV 파일

    with open(mutation_csv, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            sample = row['sample']
            gene = row['gene']
            chrom = row['chrom']
            pos = int(row['pos'])  # CSV 파일에서 pos 열 사용
            ref = row['ref']
            alt = row['alt']

            find_and_extract_reads(bam, chrom, pos, alt, ref, gene)
