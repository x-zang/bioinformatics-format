# Bioinformatics-format

This is a list of (oversimplified) descriptions of different file formats used in bioinformatics.

UCSC host a (relatively ) more comprehensive [list](http://genome.ucsc.edu/FAQ/FAQformat.html).



## Sequence, index

#### .fasta/.fa

```
> Id comment
AAAAATTTTTTTCCCCCCGGGGGG
GGGGCCCCNNAAAAAAAACCCCCC
```

Anyline starting with `>` is a new sequence's name and consecutive lines are its sequences.

```
# TO remove comment
sed -e 's/^\(>[^[:space:]]*\).*/\1/' my.fasta > mymodified.fasta
```


#### .fastq/.fq

Consists of multiple 4-lines: ID, sequence, + sign, quality score

```
@HWI-D00119:50:H7AP8ADXX:1:1101:1213:2058 1:N:0:TAAGGCGA
ACTCCAGCCTGGGCAACAGAGCAAGGCTCGGTCTCCCAAAAAAAAAAAAAAAAAAAAAAAATTGGAACTCATTTAAAAACACTTATGAAGAGTTCATTTCT
+
@@@D?BD?A>CBDCED;EFGF;@B3?::8))0)8?B>B@FGCFEEBC######################################################
```

See also [unmapped bam](####.unmapped.bam/.ubam:), a format containing same information as fastq.



#### SRA

It is not exactly a file format but frequently used. SRA stands for Sequence Read Archive, an NCBI database. To retrieve .fastq from SRA, you need sratoolkit.

```
# download & dump with SRR run accession
prefetch SRR0000000
fastq-dump --gzip --split-files SRR0000000
```



#### .fai

Index file for genome. See http://www.htslib.org/doc/faidx.html

```
samtools faidx Homo_sapiens_assembly38.genome.fasta
```



#### .dict

sequence dictionary file 

```shell
samtools dict Homo_sapiens_assembly38.genome.fasta -o Homo_sapiens_assembly38.genome.dict
```



## Alignment, Index

#### .sam/.bam

 http://samtools.github.io/hts-specs/SAMv1.pdf

```shell
samtools view file.bam > file.sam
santools view -b file.sam > file.sam
```

There are 11 mandatory fields in sam/bam format

```
Col Field 	Type 			Regexp/Range 									Brief description
1 	QNAME		String 		[!-?A-~]{1,254}								Query template NAME
2 	FLAG 		Int 			[0, 216 − 1] 									bitwise FLAG
3 	RNAME 	String 		\*|[:rname:∧*=][:rname:]* 		Reference sequence NAME11 
4 	POS 		Int 			[0, 231 − 1]					 				1-based leftmost mapping POSition
5 	MAPQ 		Int 			[0, 28 − 1] 									MAPping Quality
6 	CIGAR 	String 		\*|([0-9]+[MIDNSHPX=])+ 			CIGAR string
7 	RNEXT 	String 		\*|=|[:rname:∧*=][:rname:]* 	Reference name of the mate/next read 
8 	PNEXT 	Int 			[0, 231 − 1] 									Position of the mate/next read
9 	TLEN 		Int 			[−231 + 1, 231 − 1] 					observed Template LENgth
10 	SEQ 		String 		\*|[A-Za-z=.]+ 								segment SEQuence
11 	QUAL 		String 		[!-~]+ 												ASCII of Phred-scaled base QUALity+33
```



#### .unmapped.bam/.ubam:

Unmapped bam is recommended in GATK and other tools to replace fastq. It is bam file, without alignment information.

Convert fastq to ubam using GATK `FastqToSam`.

```shell
gatk FastqToSam \
-F1 reads_R1.fastq \
-F2 reads_R2.fastq \
-O reads.unmapped.bam \
--SAMPLE_NAME sample001
```

#### .bai

```shell
samtools index file.bam
```



## Graph representation

#### .gfa

#### .gam (graph-version bam)

#### .vg, .xg



## Genome features

#### .gtf/.gff
Both `.gtf` and `.gvf` are two specific types of `.gff`.

```

```



#### .bed 

```

```



#### .wig

```

```



## Variant calling

#### .vcf/.bcf

[.vcf](http://samtools.github.io/hts-specs/VCFv4.2.pdf) is Variant Calling Format. bcf is binary .vcf.



#### .vcf.idx

The index file of .vcf. Use IGV or [igvtools](https://software.broadinstitute.org/software/igv/igvtools_commandline) to generate.

```
# make file.vcf.idx with igvtools
igvtools index file.vcf 
```

#### Specifications explained from [10X genome](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/vcf)

#### `GT `field

The `GT` (genotype) field encodes allele values separated by either of `/` unphased  or `|` phased.  For phased genotypes, the allele to the left of the bar is haplotype 1, and the allele to the right of the bar is haplotype 2.

#### `PS` field

Variants with the same PS value are in the same phase block. Variants with different PS values are not phased with respect to one another. e.g. `pos 1000, GT:PS 1|0:PS1  ` and `pos 2000GT:PS 1|2:PS2` are not phased together, but `pos 1000, GT:PS 1|0:PS1` and `pos 1050, GT:PS 0|1:PS1` are phased.

#### .tbi

.tbi is index file of compressed vcf file (.vcf.gz). Use [bgzip](http://www.htslib.org/doc/bgzip.html) to compress .vcf file and use [tabix](http://www.htslib.org/doc/tabix.html) to generate .tbi file.

```
# make file.vcf.gz
bgzip file.vcf

# make file.vcf.gz.tbi
tabix -p vcf file.vcf.gz
```



## Others

#### Compressed files

```
unzip file.zip
gunzip file.gz
# zcat outputs to stdout
zcat file.gz > file
# .tar.gz and .tgz are the same thing
tar -zcvf file.tar.gz
tar -zcvf file.tgz
```
