# LazyMapping

Mapping of three LAZY mutants for [Soonju Park](<sjpark75@wonkwang.ac.kr>) at Wonk Wang University South Korea. These were previously mapped to chromosomes 2, 8, and 10. Included are scripts to analyze primary bulk segregant analysis (BSA) data using Illumina sequencing. Scripts were designed for use and submission on an SGE high performance computing (HPC) cluster. Most shell scripts can be submitted via `qsub SCRIPT.sh`. Also will work from command line with `bash SCRIPT.sh`.

requisite programs:
1. samtools
2. fastqc
3. bwa
4. trimmomatic
5. SNPeff

### Pipeline outline:
1. Download reads from sequence supplier using `wget.sh`
2. Check read quality for obvious issues using `fastqc.sh`
3. Trim and align reads with trimmomatic and bwa using `bwa_lz.sh`
4. Call SNPs using samtools mpileup using `multisample_snpcall_lazy.sh`
5. Annotate SNPs for functional effects using `SNPann_lazy.sh`
6. Filter SNPs for those within genes using `genic_filt.sh`
7. Take filtered snp list and plot SNP index in R using the Rmd scripts.


