#!/bin/bash
#$ -cwd
#$ -j y
#$ -l m_mem_free=3G
#$ -pe threads 16

# set the task IDs and output dir for logs.
#$ -t 1-4
#$ -o /sonas-hs/lippman/hpc/data/LazyMutants_SJP/MutantMappingDec2016/log

echo "Reading in parameters for SGE_TASK_ID:$SGE_TASK_ID";date
# Read in the parameter line and parse
Parameters=$(sed -n -e "$SGE_TASK_ID p" Samples)
read1=$( echo "$Parameters" | awk '{print $1}' )
read2=$( echo "$Parameters" | awk '{print $2}' )
lane=$( echo "$Parameters" | awk '{print $3}' )
lanelibrary=$( echo "$Parameters" | awk '{print $4}' )
library=$( echo "$Parameters" | awk '{print $5}' )

# Define commonly used files/paths
ref=$HOME/indexes/SolycGenome_2.50
GATK=$HOME/bin/GenomeAnalysisTK.jar
picard=$HOME/bin/picard-tools-1.126/bin/picard.jar
trimmomatic=$HOME/bin/Trimmomatic-0.32/trimmomatic-0.32.jar

# Begin parallel job that will trim reads from the various libraries, align with bwa, calculate some basic stats, and finish. Output will be a final dups marked bam file ( "$lanelibrary"_dups.bam ) that has all alignments included for use in a mpileup.
echo "Processing lane:$lane from library:$library" ; date

# Trim reads
java -jar "$trimmomatic" PE -threads 16 ./"$read1" ./"$read2" "$TMPDIR"/"$lanelibrary"_P1.fastq "$TMPDIR"/"$lanelibrary"_U1.fastq "$TMPDIR"/"$lanelibrary"_P2.fastq "$TMPDIR"/"$lanelibrary"_U2.fastq ILLUMINACLIP:$HOME/bin/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:2:40:15:1:FALSE LEADING:30 TRAILING:30 MINLEN:75 TOPHRED33
# clean out unpaired reads
rm -fv "$TMPDIR"/"$lanelibrary"_U1.fastq "$TMPDIR"/"$lanelibrary"_U2.fastq
# clean out raw reads
#rm -fv ./"$read1" ./"$read2"

# Align library reads (PE) to $ref genome.
echo "Aligning lane:$lane from library:$library (PE) to $ref" ; date
bwa mem -M -t 16 -R "@RG\tID:$lanelibrary\tSM:$lane\tLB:$library\tPL:illumina" "$ref" "$TMPDIR"/"$lanelibrary"_P1.fastq "$TMPDIR"/"$lanelibrary"_P2.fastq > "$TMPDIR"/"$lanelibrary".sam
rm -fv "$TMPDIR"/"$lanelibrary"_P1.fastq "$TMPDIR"/"$lanelibrary"_P2.fastq
echo "Alignment for lane:$lane from library:$library complete" ; date

# convert to bam, sort, and index
echo "converting to bam, sorting bam, and indexing bam for lane:$lane from library:$library" ; date
samtools view -b -S -T $ref -o "$TMPDIR"/"$lanelibrary".bam "$TMPDIR"/"$lanelibrary".sam
samtools sort -O bam -o "$TMPDIR"/"$lanelibrary"_sorted.bam -T $TMPDIR/temp "$TMPDIR"/"$lanelibrary".bam
samtools index "$TMPDIR"/"$lanelibrary"_sorted.bam
# Clean out sam and unsorted bam
rm -fv "$TMPDIR"/"$lanelibrary".sam "$TMPDIR"/"$lanelibrary".bam
echo "Completed bam conversion and sorting for lane:$lane from library:$library" ; date

echo "Calculating metrics for alignment of lane:$lane from library:$library" ; date
samtools flagstat "$TMPDIR"/"$lanelibrary"_sorted.bam > "$lanelibrary".alignment_metrics
cat "$lanelibrary".alignment_metrics

# For read groups from the same library, mark duplicates with picard and merge bam files.
echo "Marking duplicate reads in lane:$lane and library:$library with picard" ; date
java -Xmx2g -jar "$picard" MarkDuplicates VALIDATION_STRINGENCY=LENIENT INPUT="$TMPDIR"/"$lanelibrary"_sorted.bam OUTPUT="$library"_dups.bam METRICS_FILE="$library"_dups.Metrics VERBOSITY=WARNING
## index the duplicates marked bam files
samtools index "$library"_dups.bam
# clean out the sorted bam
rm -fv "$TMPDIR"/"$lanelibrary"_sorted.bam

echo "Finished trimming, aligning, and marking duplicates for lane:$lane from library:$library" ; date
