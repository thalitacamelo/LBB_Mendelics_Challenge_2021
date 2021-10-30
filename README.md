# silver-octo-bassoon92

## Day1
The *BWA-MEM* algorithm was used for sequence alignment, followed by *freebayes* for variant calling.

List of softwares (version) used in the **variant_call.sh** bash script:
- samtools (1.10 Using htslib 1.10.2-3)
- freebayes (v1.3.2)
- bwa (0.7.17-r1188)
- bedtools (v2.27.1)
- bcftools ( 1.10.2 Using htslib 1.10.2-3)
- vcflib (1.0.2)

usage example:
```
variant_call.sh grch38.chr22.fasta read1.fastaq read2.fastaq
```

The output is a VCF file filtered by quality > 20 and read depth > 3: **variants_Q20_DP3.vcf.gz**.

## Day2
### Question 1
The **variants_unfiltered.vcf** was filtered by quality > 20 and read depth > 3 to **variants_Q20_DP3.vcf**.
```
vcffilter -f "QUAL > 20 & DP > 3" variants_unfiltered.vcf > variants_Q20_DP3.vcf
```

### Question 2
The parts of coverage.bam that were not coverad by sorted-fixmate-output.bam were extracted using *bedtools genomecov*.
In the 'bedtools genomecov' with -bga option, regions with zero coverage are also reported. 
This allows one to quickly extract all regions of a genome with 0 coverage.
```
bedtools genomecov -ibam sorted-fixmate-output.bam -bga /
| awk '$4==0' /
| bedtools intersect -a coverage.bed -b - > not_covered.bed 
```
### Question 3
The *samtools stats* collects statistics from BAM files and outputs in a text format.
```
samtools stats sorted-fixmate-output.bam | grep ^SN | cut -f 2- > summarystats.txt
```
The following parameters were saved in questao3.tsv:
- nreads = sequences - number of processed reads
- proper_pairs = reads properly paired - number of mapped paired reads with flag 0x2 set
- maqQ_0 = reads MQ0 - number of mapped reads with mapping quality 0
