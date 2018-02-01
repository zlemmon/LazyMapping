#!/bin/bash
#$ -cwd 
#$ -j y
#$ -l m_mem_free=2G
#$ -pe threads 16

####
# This script takes as input aligned and duplicatemarked bam files to do a 
# multisample SNP calling procedure with samtools mpileup.
###

# Define commonly used files/paths
ref=$HOME/indexes/SolycGenome_2.50
GATK=$HOME/bin/GenomeAnalysisTK.jar
picard=$HOME/bin/picard-tools-1.126/bin/picard.jar
trimmomatic=$HOME/bin/Trimmomatic-0.32/trimmomatic-0.32.jar

# Set up names for files specific to this analysis
lz10_mut_bam=lz10_mut_dups.bam
lz10_wt_bam=15_2_lz10_dups.bam
lz8_mut_bam=lz8_mut_dups.bam
lz8_wt_bam=20_3_lz8_dups.bam
lz5_mut_bam=lz5_mut_dups.bam
lz5_wt_bam=lz5_wt_dups.bam

lz10_mut_base=${lz10_mut_bam/_dups.bam/}
lz10_wt_base=${lz10_wt_bam/_dups.bam/}
lz8_mut_base=${lz8_mut_bam/_dups.bam/}
lz8_wt_base=${lz8_wt_bam/_dups.bam/}
lz5_mut_base=${lz5_mut_bam/_dups.bam/}
lz5_wt_base=${lz5_mut_bam/_dups.bam/}

####################
# Set up file paths for general files such as M82, pimp, utm, etc.
####################
pimpbam=/sonas-hs/lippman/hpc/data/Spimp/Heinz_Alignment/final_pimpHeinz_recal.bam
M82bam=/sonas-hs/lippman/hpc/data/Bolger2014_M82Genome/M82_genomic_refHeinz_dups.bam
#utm_jap_bam=/sonas-hs/lippman/hpc/data/MicrotomGenome/microtom_genomic_dups.bam
#utm_fre_bam=/sonas-hs/lippman/hpc/data/MicrotomGenome/
#pimpfullvcf=/sonas-hs/lippman/hpc/data/Spimp/Heinz_Alignment/final_pimpHeinz_all.vcf.gz
#pimpvarvcf=/sonas-hs/lippman/hpc/data/Spimp/Heinz_Alignment/final_pimpHeinz_variantsites.vcf.gz
#M82fullvcf=/sonas-hs/lippman/hpc/data/Bolger2014_M82Genome/M82_genomic_refHeinz.vcf.gz
#M82varvcf=/sonas-hs/lippman/hpc/data/Bolger2014_M82Genome/M82_genomic_refHeinz_variantsites.vcf.gz
#utmfullvcf=/sonas-hs/lippman/hpc/data/MicrotomGenome/microtom_genomic.vcf.gz
#utmvarvcf=/sonas-hs/lippman/hpc/data/MicrotomGenome/microtom_genomic_variantsites.vcf.gz
#targetsfile=/sonas-hs/lippman/hpc/home/zlemmon/indexes/ITAG2.4_gene_regions.txt

############# Portion of script that calls SNPs among all bam files for MAPPING of interval
# call a parallel mpileup on all regions (chromosomes) for all samples simultaneously
echo "Running mpileup for 8 samples ( pimp, M82, lz10 (mut and wt), lz8 (mut and wt), lz5 (mut and wt) )"
samtools view -H "$M82bam" | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | xargs -I LINE -n 1 -P 14 sh -c "samtools mpileup --ignore-RG -d 1000000 -t DP,AD -Q0 -Bugf "$ref" -r 'LINE' "$M82bam" "$pimpbam" "$lz10_mut_bam" "$lz10_wt_bam" "$lz8_mut_bam" "$lz8_wt_bam" "$lz5_mut_bam" "$lz5_wt_bam" | bcftools call -mvO z -o lazy_segsites_'LINE'.vcf.gz"
# index all individual chromosome vcf files
echo "Indexing chromosome vcf.gz files" ; date
for file in lazy_segsites_SL2.50*.vcf.gz ; do tabix -fp vcf "$file" ; done
# merge SNPs
echo "Merging SNPs called on the various chromosomes" ; date
bcftools concat -O z -o lazy_segsites.vcf.gz lazy_segsites_SL2.50*.vcf.gz
# index final variants vcf.gz file
echo "Indexing final vcf.gz file" ; date
tabix -p vcf lazy_segsites.vcf.gz
echo "Done calling SNPs"
# Remove individual chromosome files.
echo -ne "Removing individual chromosome files... "
rm -f lazy_segsites_SL2.50*.vcf.gz*
echo -ne 'Done!!\n\n'

# Rename samples to something more manageable than a full filename/path
bcftools reheader -s newsamplenames lazy_segsites.vcf.gz > lazy_segsites_rehead.vcf.gz
mv -f lazy_segsites_rehead.vcf.gz lazy_segsites.vcf.gz
bcftools index lazy_segsites.vcf.gz

# Filter final VCF by a few simple metrics
bcftools view -m2 -M2 lazy_segsites.vcf.gz | bcftools filter --SnpGap 100 -i ' TYPE="SNP" & MQ>=50 ' -Oz -o lazy_segsites_filt.vcf.gz

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP[\t%SAMPLE %GT %AD]\n' lazy_segsites_filt.vcf.gz | awk 'BEGIN{FS="\t"}\
	{printf("%s\t%s\t%s\t%s\t%s",$1, $2, $3, $4, $5);\
        for(i=6;i<=NF;i++){\
		string=$i; split(string, SGD, " "); 
		sample=SGD[1]; gt=SGD[2]; ad=SGD[3];\
		split(ad,ad2,",");\
		printf("\t%s %s %s %s", sample, gt, ad2[1], ad2[2]);\
	}\
	printf("\n");\
}' > lazy_segsites_filt.txt

gzip -v lazy_segsites_filt.txt

echo -ne 'Done!\n\n' ; date
