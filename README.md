## Repository Submitted to the [LBB-Mendelics-Challenge-2021](https://github.com/mendelics/lbb-mendelics-2021/tree/main)

This was my first GitHub repository. At the time, I had just learned code flow control, functions, and command-line arguments. I had no experience with sequencing analysis pipelines but understood key concepts like alignment, variant calling, and annotation.

With zero knowledge of the recommended tools, I relied on my basic **Bash** knowledge and problem-solving skills . Over three days (8 hours/day), I used standard bioinformatics tools and scripting to complete the tasks.

I didn’t win, but the challenge led to a job at Mendelics, where I spent 2 years and 7 months learning software development good practices and large-scale data analysis.

<details><summary>Day1</summary>

### Task: To extract the variants found on chromosome 22 from the sample.
  
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
$ variant_call.sh grch38.chr22.fasta read1.fastaq read2.fastaq
```

The output is a VCF file filtered by quality > 20 and read depth > 3: **variants_Q20_DP3.vcf.gz**.

</details>

<details><summary>Day2</summary>
  
### Question 1: Which variants should be disregarded in your VCF? - Any metric from the software of choice can be used. Discuss the metric used.
The **variants_unfiltered.vcf** was filtered by quality > 20 and read depth > 3 to **variants_Q20_DP3.vcf**.
```
$ vcffilter -f "QUAL > 20 & DP > 3" variants_unfiltered.vcf > variants_Q20_DP3.vcf
```

### Question 2: Discuss the regions with low coverage and what your criteria were.
The parts of coverage.bam that were not coverad by sorted-fixmate-output.bam were extracted using *bedtools genomecov*.
In the 'bedtools genomecov' with -bga option, regions with zero coverage are also reported. 
This allows one to quickly extract all regions of a genome with 0 coverage.
```
$ bedtools genomecov -ibam sorted-fixmate-output.bam -bga /
| awk '$4==0' /
| bedtools intersect -a coverage.bed -b - > not_covered.bed 
```
### Question 3: Obtain information about your alignment. How many reads? What percentage of them were mapped correctly? Did many align to more than one location in the genome with the same quality?
The *samtools stats* collects statistics from BAM files and outputs in a text format.
```
$ samtools stats sorted-fixmate-output.bam | grep ^SN | cut -f 2- > summarystats.txt
```
The following parameters were saved in questao3.tsv:
- nreads = sequences - number of processed reads
- proper_pairs = reads properly paired - number of mapped paired reads with flag 0x2 set
- maqQ_0 = reads MQ0 - number of mapped reads with mapping quality 0

</details>

<details><summary>Day3</summary>
  
### Question 1: Obtain the Ti/Tv ratio (transitions and transversions) of the variants found on chromosome 22.
TS/TV from variants_Q20_DP3.vcf
```
$ bcftools stats variants_Q20_DP3.vcf.gz | grep "TSTV"
2.52
```

### Question 2: How many variants are found in the region from 16000000 to 20000000?
Number of variants between 16000000 and 20000000 positions from variants_Q20_DP3.vcf
```
$ grep -v '^#' variants_Q20_DP3.vcf | awk '{if ($2 >= 16000000 && $2 <= 20000000) print}' | wc -l
292
```

### Question 3: Display the content of the VCF line related to a variant **Non-synonymous** and a variant in gnomAD v3.1.1 with MAF < 0.01.
Creation of avinput file format to be used in annovar
```
$ perl ~/bioinfo/app/annovar/convert2annovar.pl -format vcf4 variants_Q20_DP3.vcf > variants_Q20_DP3.avinput
```
Annotation using ANNOVAR
```
$ perl ~/bioinfo/app/annovar/table_annovar.pl variants_Q20_DP3.avinput ~/bioinfo/app/annovar/humandb/ -buildver hg38 -out ~/silver-octo-bassoon92/ -remove -protocol refGene,exac03,clinvar_20190305,gnomad_exome -operation g,f,f,f -nastring .
```
Saving the location of a random nonsynonymous variant from annovar_hg38_multianno.txt to the variable 'nonsynonymous'
```
$ nonsynonymous=$(shuf annovar_hg38_multianno.txt | grep 'nonsynonymous' | head -n1 | awk '{print $1 "\t" $2}')
```
Printing a random 'nonsynonymous' from variants_Q20_DP3.vcf
```
$ grep "$nonsynonymous" variants_Q20_DP3.vcf
chr22   36028402        .       A       C       718.463 .       AB=0;ABP=0;AC=2;AF=1;AN=2;AO=28;CIGAR=1X;DP=28;DPB=28;DPRA=0;EPP=5.80219;EPPR=0;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=0;NS=1;NUMALT=1;ODDS=43.4214;PAIRED=1;PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=830;QR=0;RO=0;RPL=16;RPP=4.25114;RPPR=0;RPR=12;RUN=1;SAF=11;SAP=5.80219;SAR=17;SRF=0;SRP=0;SRR=0;TYPE=snp   GT:DP:AD:RO:QR:AO:QA:GL 1/1:28:0,28:0:0:28:830:-74.9858,-8.42884,0
```
Printing a random 'nonsynonymous' from annovar_hg38_multianno.txt
```
$ grep "$nonsynonymous" annovar_hg38_multianno.txt
chr22   36028402        36028402        A       C       exonic  RBFOX2  .       nonsynonymous SNV       RBFOX2:NM_001082578:exon1:c.T24G:p.H8Q,RBFOX2:NM_001082579:exon1:c.T24G:p.H8Q,RBFOX2:NM_001349999:exon1:c.T24G:p.H8Q     1       1       1       1       1       1       1       1       .       .       .       .       .       1       1       1       1       0.9986  1       1       1       1       1       1       rs9607299
```

</details>
