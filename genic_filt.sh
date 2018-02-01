#!/bin/bash
#$ -cwd
#$ -j y
#$ -l m_mem_free=10G

# Define commonly used files/paths
ref=$HOME/indexes/SolycGenome_2.50
GATK=$HOME/bin/GenomeAnalysisTK.jar
picard=$HOME/bin/picard-tools-1.126/bin/picard.jar
trimmomatic=$HOME/bin/Trimmomatic-0.32/trimmomatic-0.32.jar

# Set up file paths
#mutbase=lazy_mut
#mutbam=lazy_mut_dups.bam
mutbase=lazy_segsites
multisamplevcf=lazy_segsites.vcf.gz
#pimpfullvcf=/sonas-hs/lippman/hpc/data/Spimp/Heinz_Alignment/final_pimpHeinz_all.vcf.gz
#pimpvarvcf=/sonas-hs/lippman/hpc/data/Spimp/Heinz_Alignment/final_pimpHeinz_variantsites.vcf.gz
#M82fullvcf=/sonas-hs/lippman/hpc/data/Bolger2014_M82Genome/M82_genomic_refHeinz.vcf.gz
#M82varvcf=/sonas-hs/lippman/hpc/data/Bolger2014_M82Genome/M82_genomic_refHeinz_variantsites.vcf.gz
#utmfullvcf=/sonas-hs/lippman/hpc/data/MicrotomGenome/microtom_genomic.vcf.gz
#utmvarvcf=/sonas-hs/lippman/hpc/data/MicrotomGenome/microtom_genomic_variantsites.vcf.gz
#targetsfile=/sonas-hs/lippman/hpc/home/zlemmon/indexes/ITAG2.4_gene_regions.txt

# Filter novel SNPs for genic regions.
bcftools view --regions-file $HOME/indexes/ITAG2.4_gene_regions.txt lz10_mut_novel.vcf.gz
bcftools view --regions-file $HOME/indexes/ITAG2.4_gene_regions.txt lz10_mut_novel.vcf.gz | grep -v "#" | grep "SL2.50ch10" | egrep "HIGH|MODERATE"

bcftools view --regions-file $HOME/indexes/ITAG2.4_gene_regions.txt lz8_mut_novel.vcf.gz
bcftools view --regions-file $HOME/indexes/ITAG2.4_gene_regions.txt lz8_mut_novel.vcf.gz | grep -v "#" | grep "SL2.50ch08" | egrep "HIGH|MODERATE"

bcftools view --regions-file $HOME/indexes/ITAG2.4_gene_regions.txt lz5_mut_novel.vcf.gz
bcftools view --regions-file $HOME/indexes/ITAG2.4_gene_regions.txt lz5_mut_novel.vcf.gz | grep -v "#" | grep "SL2.50ch05" | egrep "HIGH|MODERATE"
