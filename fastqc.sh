#!/bin/bash
#$ -cwd
#$ -j y
#$ -l m_mem_free=5G

mkdir -p ./fastqc
/sonas-hs/lippman/hpc/home/zlemmon/bin/FastQC/fastqc -o ./fastqc *.fastq.gz
