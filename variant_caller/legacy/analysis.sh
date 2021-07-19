#!/bin/sh
#!/usr/bin/perl


##environment
# Number of available CPU cores
NO_CORES=20

# This stores the JAVA path
JAVA="/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java -jar"


##folder structure
# This stores the workdir path
WORK_LOCATION="/home/genomica/analysis/files/170273"

SCRIPT_LOCATION="/home/genomica/analysis/script/UMFT"
# This script will execute the variant calling for all files in the directory
#not called later. keep it still?


##files and reference
no_FQ_files=0
#paired fw files have to be merged down to 2 and placed inside a folder bearing the same of the sample

# This stores the human reference genome build hg19
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/
#url37d5="ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
#url38="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz"

#alternatively iGenomes Illumina at https://support.illumina.com/sequencing/sequencing_software/igenome.html
#or even better the GATK bundle https://gatk.broadinstitute.org/hc/en-us/articles/360035890811

#preparation steps are needed, as instructed by bwa, gatk CreateSequenceDictionary -R and samtools faidx

HG19="/home/genomica/reference/bwa_GRCh37/hs37d5.fa"


##tools
# https://gatk.broadinstitute.org/hc/en-us/articles/360041320571--How-to-Install-all-software-packages-required-to-follow-the-GATK-Best-Practices
# This variable stores the BWA 0.7.17
BWA_MEM="/home/genomica/software/bwa.kit/bwa mem"

# This variable stores the Picard path 2.18.23
# https://broadinstitute.github.io/picard/
PICARD="/home/genomica/software/picard2.18.23/picard.jar"

# This variable contains the GATK 3.8.1
# https://gatk.broadinstitute.org/hc/en-us
GATK="/home/genomica/software/GenomeAnalysisTK-3.8-1-0/GenomeAnalysisTK.jar"

# This stores BEDTOOLS path
BEDTOOLS="/home/genomica/software/bedtools2/bin/bedtools"

# This stores SAMTOOLS path
SAMTOOLS="/home/genomica/software/samtools-1.12/bin/samtools"


#Starting the loop
echo "Analysis is Active!"

while :
do
listOfFiles=$(ls -1 | wc -l)
if [ $listOfFiles -ge 1 ];
then

for file in $WORK_LOCATION*;
do
if [ -d "$file" ]
then

# Saving DIR path and extracting the DIR name (aka. ID)
dir_path=$file
dir_name="$(basename $dir_path)"

# If there is a BAM file with the directory name then it will skip this DIR
if [ -e $dir_path/$dir_name.bam ];
then
echo "INFO: "$dir_name" - is done."
continue
fi
echo $dir_path


#Verifying FASTQ files number is correct. Expected 2
#$no_FQ_files = $(find $dir_path/ -name "*.gz" | wc -l)

#echo $no_FQ_files
#if [ $(find $dir_path/ -name "*.gz" | wc -l) -eq 2 ];
#then

echo "Executing!"
#execute_analysis $dir_path $dir_name
i=1

PATH=$dir_path
echo "double check path $PATH"

ID=$dir_name
echo "double check ID $ID"

for fq_file in $PATH/*.fastq.gz;
do

if [ $i -eq 1 ];
then
fq1=$fq_file
i=2
else
fq2=$fq_file
fi
done

echo "double check fq files $fq1 $fq2"

#definitions
initial_sam_file=$PATH"/"$ID"_f01_initial_align.sam"
initial_bam_file=$PATH"/"$ID"_f02_initial_align.bam"
align_metrics=$PATH"/"$ID"_f03_align_metrics.txt"
insert_size_metrics=$PATH"/"$ID"_f03_insert_metrics.txt"
depth_file=$PATH"/"$ID"_f04_depth_out.txt"
metrics_file=$PATH"/"$ID"_f05_metrics.txt"
dedup_bam_file=$PATH"/"$ID"_f06_dedup_align.bam"
realign_targets=$PATH"/"$ID"_f07_realign_targets.list"
realign_bam_file=$PATH"/"$ID"_f08_realign.bam"
raw_variants=$PATH"/"$ID"_f09_raw_variants.vcf"
raw_snps=$PATH"/"$ID"_f10_raw_snps.vcf"
raw_indels=$PATH"/"$ID"_f11_raw_indels.vcf"
raw_filtered_snps=$PATH"/"$ID"_f12_raw_filtered_snps.vcf"
raw_filtered_indels=$PATH"/"$ID"_f13_raw_filtered_indels.vcf"
recal_data=$PATH"/"$ID"_f14_recal_data.table"
post_recal_data=$PATH"/"$ID"_f15_post_recal_data.table"
final_bam_file=$PATH"/"$ID"_final.bam"
raw_recal_variants=$PATH"/"$ID"_f16_raw_variants_recal.vcf"
raw_recal_snps=$PATH"/"$ID"_f17_raw_recal_snps.vcf"
raw_recal_indels=$PATH"/"$ID"_f18_raw_recal_indels.vcf"
final_indels=$PATH"/"$ID"_final_indels.vcf"
final_variants=$PATH"/"$ID"_final_variants.vcf"
final_snps=$PATH"/"$ID"_final_snps.vcf"
final_filtered_snps=$PATH"/"$ID"_final_filtered_snps.vcf"
final_filtered_indels=$PATH"/"$ID"_final_filtered_indels.vcf"

#alignment of the FASTQ

RG="@RG\tID:"$ID"\tPL:illumina\tPU:run\tLB:lib1\tSM:"$ID
echo RG 

$BWA_MEM -M -t $NO_CORES -R $RG $HG19 $fq1 $fq2 > $initial_sam_file

#sort
$JAVA $PICARD SortSam \
I=$initial_sam_file \
O=$initial_bam_file \
SORT_ORDER=coordinate \
CREATE_INDEX=TRUE 
#--TMP_DIR temp

#metrics
$JAVA $PICARD CollectAlignmentSummaryMetrics \
R=$HG19 \
I=$initial_bam_file \
O=$align_metrics

$JAVA $PICARD CollectInsertSizeMetrics \
I=$initial_bam_file \
O=$insert_size_metrics \
H=insert_size_histogram.pdf \
M=0.5

$SAMTOOLS depth -a $initial_bam_file > $depth_file

#duplicates
$JAVA $PICARD MarkDuplicates \
I=$initial_bam_file \
O=$dedup_bam_file \
M=$metrics_file \
CREATE_INDEX=TRUE

#realign
$JAVA $GATK -T RealignerTargetCreator \
-nt $NO_CORES \
-R $HG19 \
-I $dedup_bam_file \
-o $realign_targets

$JAVA $GATK -T IndelRealigner \
-R $HG19 \
-I $dedup_bam_file \
-targetIntervals $realign_targets \
-o $realign_bam_file


#raw variants
$JAVA $GATK -T HaplotypeCaller \
-nct $NO_CORES \
-R $HG19 \
-I $realign_bam_file \
-o $raw_variants

#select
$JAVA $GATK -T SelectVariants \
-R $HG19 \
-V $raw_variants \
-selectType SNP \
-o $raw_snps

$JAVA $GATK -T SelectVariants \
-R $HG19 \
-V $raw_variants \
-selectType INDEL \
-o $raw_indels


#filter
$JAVA $GATK -T VariantFiltration \
-R $HG19 \
-V $raw_snps \
--filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' \
--filterName "basic_snp_filter" \
-o $raw_filtered_snps

$JAVA $GATK -T VariantFiltration \
-R $HG19 \
-V $raw_indels \
--filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' \
--filterName "basic_indel_filter" \
-o $raw_filtered_indels

# knownSites - use dbSNP?

#recalibrate
$JAVA $GATK -T BaseRecalibrator \
-nct $NO_CORES \
-R $HG19 \
-I $realign_bam_file \
-knownSites $raw_filtered_snps \
-knownSites $raw_filtered_indels \
-o $recal_data

$JAVA $GATK -T BaseRecalibrator \
-nct $NO_CORES \
-R $HG19 \
-I $realign_bam_file \
-knownSites $raw_filtered_snps \
-knownSites $raw_filtered_indels \
-BQSR $recal_data \
-o $post_recal_data

#final bam
$JAVA $GATK -T PrintReads \
-nct $NO_CORES \
-R $HG19 \
-I $realign_bam_file \
-BQSR $recal_data \
-o $final_bam_file


#reassemble reads for better indel
$JAVA $GATK -T HaplotypeCaller \
-nct $NO_CORES \
-R $HG19 \
-I $final_bam_file \
-o $raw_recal_variants

#select
$JAVA $GATK -T SelectVariants \
-R $HG19 \
-V $raw_recal_variants \
-selectType SNP \
-o $raw_recal_snps

$JAVA $GATK -T SelectVariants \
-R $HG19 \
-V $raw_recal_variants \
-selectType INDEL \
-o $raw_recal_indels

#filter
$JAVA $GATK -T VariantFiltration \
-R $HG19 -V $raw_recal_snps \
--filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' \
--filterName "basic_snp_failed" \
--filterExpression 'DP < 10' \
--filterName "Low_DP" \
-o $final_snps

$JAVA $GATK -T VariantFiltration \
-R $HG19 -V $raw_recal_indels \
--filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' \
--filterName "basic_indel_failed" \
--filterExpression 'DP < 10' \
--filterName "Low_DP" \
-o $final_indels

#select
$JAVA $GATK -T SelectVariants \
-R $HG19 \
-V $final_snps \
-o $final_filtered_snps \
--excludeFiltered

$JAVA $GATK -T SelectVariants \
-R $HG19 \
-V $final_indels \
-o $final_filtered_indels \
--excludeFiltered

#final variants
$JAVA $PICARD MergeVcfs \
I=$final_filtered_indels I=$final_filtered_snps \
O=$final_variants


##clean UP
#rm $initial_sam_file
#rm $initial_bam_file
#rm $align_metrics
#rm $insert_align_metrics
#rm $depth_file
#rm $metrics_file
#rm $dedup_bam_file
#rm $realign_targets
#rm $realign_bam_file
#rm $raw_varints
#rm $raw_snps
#rm $raw_indels
#rm $raw_filtered_snps
#rm $raw_filtered_indels
#rm $recal_data
#rm $post_recal_data
#rm $raw_recal_variants
#rm $raw_recal_indels
#rm $raw_recal_snps

echo "Finished "$ID"!"

#else
#echo "ERROR:"$dir_name" Wrong number of FASTQ files. Expected 2."
#fi

fi
done
fi
#sleep 10s
break
done

echo "Analysis is Over!"

#initial script UMFT cristian_zimbru@yahoo.com


