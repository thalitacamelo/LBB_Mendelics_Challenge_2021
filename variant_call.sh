#!/usr/bin/env bash
#
#Bash script to detect genetic variants using 'freebayes', which requires two inputs: a FASTA reference sequence, and a sample.BAM-format alignment file sorted by reference position.
#The reference sequence is indexed using 'bwa index'.
#'bwa mem' algorithm mapps the sequenced reads to the reference genome.
#'samtools' is used to get a BAM-format alignment file sorted by reference position.
#'vcffilter' filter from 'vcflib' remove low quality reads from the vcf file (quality> 20 and read depth > 3).
#This script assumes that the sequence length distribution is the same for all the sample reads and the adapters were trimmed off.
#'freebayes', 'bwa', 'samtools', 'bedtools', 'bcftools' and 'libvcflib-tools' can be installed using 'sudo apt install' on Ubuntu.
#Installing 'FastQC' is as simple as unzipping the zip file it comes in into a suitable location. But it requires a suitable Java Runtime Environment (JRE) installed.
#Running 'variant_call.sh' inside a folder holding the reference sequence and the sample reads (paired end) will output a variants.vcf.gz file as well as the intermediate ones (.sam; .bam...)
#usage: variant_call.sh reference.fasta read1.fastaq read2.fastaq

reference=$1
read1=$2
read2=$3

if [[ -z "$reference" || -z "$read1" || -z "$read2" ]]; then

    echo usage: variant_call.sh reference.fasta read1.fastaq read2.fastaq
    exit 1

fi

if [[ ! -f $reference  ]]; then
    echo file $reference not found
    exit 1
fi

#The 'BWA index' command creates an index of the reference genome in fasta format
#
bwa index $reference
#
#The 'BWA mem' algorithm performs local alignment and produce alignments for different part of the query sequence
#
alignment_output=alignment_output.sam
bwa mem $reference $read1 $read2 > $alignment_output
#
#The 'samtools fixmate' tool corrects any flaws in read-pairing that may have been introduced by the aligner
#
fixmate_output=fixmate-output.bam
samtools fixmate $alignment_output $fixmate_output
#
#The 'samtools sort' command convertes data to genome chromosome and coordinate order.
#
sorted_fixmate_output=sorted-fixmate-output.bam
samtools sort -O bam -o $sorted_fixmate_output -T /tmp/ $fixmate_output
#
#The 'santools index' command index a coordinate-sorted BGZIP-compressed SAM, BAM or CRAM file for fast random access
#
samtools index $sorted_fixmate_output
#
#The 'freebayes' tool will produce a VCF file describing all SNPs, INDELs, and haplotype variants between the reference and sample
#
variants=variants.vcf
freebayes -f $reference $sorted_fixmate_output > $variants
#
#The 'vcffilter' filter by variant quality> 20 and read depth > 3
#
filtered_variants=variants_Q20_DP3.vcf
vcffilter -f "QUAL > 20 & DP > 3" $variants > $filtered_variants
#
#bgzip
bgzip variants_Q20_DP3.vcf