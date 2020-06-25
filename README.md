# Bioinformatics-format

This is a list of (oversimplified) description of different file formats used in bioinformatics.



## Sequence, alignment, index

#### .fasta/.fa

```
> Id comment
AAAAATTTTTTTCCCCCCGGGGGG
```

#### .fastq/.fq

```
@Id comment
AAAAATTTTTTTCCCCCCGGGGGG
+
..dsadfsf
```

#### SRA

It is not exactly a file format but frequently used. SRA stands for Short Read Archive, a NCBI database. To retrieve .fastq from SRA, you need sratoolkit.

```
fastq-dump
```



#### .fai

```

```

#### .dict

```

```



#### .sam/.bam

 http://samtools.github.io/hts-specs/SAMv1.pdf

```

```

#### .unmapped.bam/.ubam:

unmapped bam is recommended in GATK and other tools 

```

```

```
# Convert fastq to ubam using GATK FastqToSam

```



## Genome features

.gtf/.gff

```

```



.bed 

```

```



.wig

```

```



## Variant calling

#### .vcf/.bcf

[.vcf](http://samtools.github.io/hts-specs/VCFv4.2.pdf) is Variant Calling Format. bcf is binary .vcf.

```

```

#### .vcf.idx

The index file of .vcf. Use IGV or [igvtools](https://software.broadinstitute.org/software/igv/igvtools_commandline) to generate.

```
# make file.vcf.idx with igvtools
igvtools index file.vcf 
```

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