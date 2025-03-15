#!/usr/bin/env sh
INPUT=$1

cp /home/dell/ncbi/feng/XQ-20240820-00151/${INPUT}/${INPUT}.R1.fq.gz .
cp /home/dell/ncbi/feng/XQ-20240820-00151/${INPUT}/${INPUT}.R2.fq.gz .
mv ${INPUT}.R2.fq.gz ${INPUT}_2.fq.gz
mv ${INPUT}.R1.fq.gz ${INPUT}_1.fq.gz
gzip -d ${INPUT}_1.fq.gz
gzip -d ${INPUT}_2.fq.gz
mv ${INPUT}_1.fq ${INPUT}_1.fastq
mv ${INPUT}_2.fq ${INPUT}_2.fastq
