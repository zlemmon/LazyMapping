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

########## Portion of script that calls effects for novel mutant SNPs.
# Use snpEff to annotate the variants for functional effect
echo -ne "Annotating novel mutant SNPs..." 
annvcf="$mutbase"_ann.vcf.gz
java -jar $HOME/bin/snpEff/snpEff.jar ann SL2.50 "$multisamplevcf" | bcftools view -Oz -o $annvcf
echo -ne " Done!\n\n"

# parse full file into easier to read stuff.
anntxt="$mutbase"_ann.txt
echo -ne "Parsing $annvcf into $anntxt... "
java -jar $HOME/bin/snpEff/SnpSift.jar extractFields $annvcf CHROM POS REF ALT "ANN[*].GENE" "ANN[*].IMPACT" "ANN[*].EFFECT" "ANN[*].HGVS_P" > $anntxt
echo -ne 'Done!\n\n'

# Filter mutant variants to only include those that are homozygous variant in mutant, novel (REF in M82/pimp), on specified chromosome, and .
#MappedChrom=SL2.50chXX # Chromosome where the putative mutation is located based on previous mapping.
echo "Filtering novel mutant SNPs for SNPs also in other haplotypes."
# Samples are: M82, pimp, Micro-Tom, fas_mut, fas_wt, fas_ind, fas_2ndAllele, fas_Tohru
# For lz10... filter so homoygous reference everything except the lz10 samples, which should be 1/1 and 0/1 respectively.
java -jar ~/bin/snpEff/SnpSift.jar filter " isRef(GEN[M82]) & isRef(GEN[pimp]) & ( isHom(GEN[lz10_mut]) & isVariant(GEN[lz10_mut]) ) & isHet(GEN[lz10_wt]) & isHom(GEN[lz8_mut]) & isRef(GEN[lz8_mut]) & isHom(GEN[lz5_mut]) & isRef(GEN[lz5_mut]) " $annvcf | bcftools view -Oz -o lz10_mut_novel.vcf.gz
# Likewise for lz8, but with lz5 and lz10 as ref.
java -jar ~/bin/snpEff/SnpSift.jar filter " isRef(GEN[M82]) & isRef(GEN[pimp]) & ( isHom(GEN[lz8_mut]) & isVariant(GEN[lz8_mut]) ) & isHet(GEN[lz8_wt]) & isHom(GEN[lz10_mut]) & isRef(GEN[lz10_mut]) & isHom(GEN[lz5_mut]) & isRef(GEN[lz5_mut]) " $annvcf | bcftools view -Oz -o lz8_mut_novel.vcf.gz
# Likewise for lz5, but with lz8 and lz10 as ref.
java -jar ~/bin/snpEff/SnpSift.jar filter " isRef(GEN[M82]) & isRef(GEN[pimp]) & ( isHom(GEN[lz5_mut]) & isVariant(GEN[lz5_mut]) ) & isHet(GEN[lz5_wt]) & isHom(GEN[lz8_mut]) & isRef(GEN[lz8_mut]) & isHom(GEN[lz10_mut]) & isRef(GEN[lz10_mut]) " $annvcf | bcftools view -Oz -o lz5_mut_novel.vcf.gz

# Parse annotated vcf for specific groups of effects
#java -jar ~/bin/snpEff/SnpSift.jar extractFields test.vcf CHROM POS REF ALT "ANN[0].IMPACT" "ANN[0].EFFECT" ANN[0].HGVS_P | egrep -ve "synonymous_variant|[53]_prime_UTR_variant" > novel_missense_startstop_variants.txt
