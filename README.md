<img src="icon.png" align="right" />

# AStruct

AStruct is a powerful tool for detecting allele-specific RNA secondary structures (RiboSNitches) within one sample from structomics sequencing data.

## Introduction

Here, we present AStruct, a Java-based software for identifying the local structure difference of SNPs from RT-Stop and RT-Mut structure sequencing data. It will automatically target the comparison window according to the maximum spanning Ref and Alt reads length around SNPs from the inputted sorted BAM file(s). SNPs that do not have enough reads support will be filtered out. Replicate samples are highly recommended by removing many false positive results. The output will provide you a AStruct score, evaluating the overall allele’s structure difference. The Astruct score and other basic SNP annotation are written in the first file suffix with riboSnotch.txt. More detailed information, including the intermediate values for calculating AStruct score, the sequence, the region position, base counts, Ref/Alt each base structure score, and the significant P Value of each base structure difference, are also listed in a second file suffix with riboSnotchDetail.txt.

## Table of Contents
* [Requirements](#Requirements)
* [Installation](#Installation)
* [Quick Start Guide](#QuickStart)
* [Usage](#Usage)
* [Output Headers](#OutputHeaders)
* [License](#License)

## Requirements

* JDK 8

## Installation

* This tool can be installed by instructions as follows:

```
git clone https://github.com/canceromics/AStruct.git
cd AStruct/tool/src
sh make.sh
cd ../..
```
The tool is generated as AStruct.jar in this directory.

## QuickStart

* Start from bam file of and input sample for example.

```
java -Xmx16g -jar AStruct.jar ribosnitch -tf example/data/T1.bam example/data/T2.bam -cf example/data/C1.bam example/data/C2.bam -o example/result/Example -gf example/reference/GRCh38.p12_NC_000022.11.fa -ef example/reference/GRCh38.p12_NC_000022.11.gff 
-mutf example/reference/SNP_NC_000022.11.vcf -repli
```
Running this instruction will result in getting a file named Example_riboSNitch.txt. `example/result/Example` means output_dir/file_prefix

More commonly used

```
java [-Xmx16g] -jar AStruct.jar <method> -tf <treat.bam(s)> -gf <genome.fa> -ef <gencode.gtf/gff> -mutf <snp.cvf> -o <path/out_prefix> [-cf control.bam(s)] [options]
```

## Usage

* More details of this tool can be found with -h parameter

```
java [-Xmx16g] -jar AStruct.jar <method> [options]
Methods:
    ribosnitch
    find different structure around SNPs. (all bams need to be sorted by coordinate)
    options:
        <-tf treat.bam ...> treat bam files
        <-cf control.bam ...> control bam files
        <-ef gencode.gtf> annotation gtf/gff file
        <-gf genome.fa> genome fasta file
        <-mutf snp.vcf> snp vcf/bed file
        <-o out_prefix> output file prefix
        [-permu count] set permutation with count times. (default is 1000)
        [-repli] enable replication filter
        [-seed seed] set random seed
        [-mut] switch between RTstop/RTmut
    
    sim:
    simulate fastq seq data with circRNAs, m6A peaks, structures and SNPs.
    options:
        <-ef gencode.gtf> annotation gtf/gff file
        <-gf genome.fa> genome fasta file
        <-o out_prefix> output file prefix
        [-mutf snp.vcf] snp vcf/bed file
        [-sf struct.bed] bed file with structure in 4th column.
        [-alen length] alignment length. (default is 150, must bigger than 0)
        [-rlen length] read fragment length. (default is 300, must bigger than alignment length) (Tips: single end simulation enabled when alignment length is the same as read fragment length)
        [-minalen length] minimum alignment length which sequence lower than this length would be dropped. (default is 30)
        [-faf frequency] enable and set fixed allele frequency. (default is disabled)
        [-peak] enable simulation with m6Apeaks.
        [-genef gene.bed] only the genes in this file can be used in simulation.
        [-bases bases] allow structures show in some bases such as "AC". (default is all bases)
        [-mut] change structure simulation from RTstop to RTmutation.
```

## OutputHeaders

* Here are definitions of headers in output file named `(output_dir/file_prefix)_riboSNitch.txt`

| Field       | Description                           |
| ---------- | ------------------------------------ |
| Chr | Chromosome Name|
| Position | SNP position in the chromosome (1 based) |
| Ref | reference snp |
| Alt | variant snp |
| Name | dbSNP ID |
| Score | score of riboSNitch (-1.0 is shown not meaningful) |
| Strand | strand of this gene |
| Gene | the longest gene covers this snp |
| Transcript | the longest transcript covers this snp |
| Annotation | feature of transcript |
| Region_in_Genome/Transcriptome | judging region in chromosome/transcript |

* Here are definitions of headers in output file named `(output_dir/file_prefix)_riboSNitchDetail.txt`

| Field       | Description                           |
| ---------- | ------------------------------------ |
| Name | dbSNP ID (key with riboSNitch file) |
| eSDC | eSDC score of riboSNitch |
| eSDC_pValue | permutation pValue of eSDC score |
| SNP_Reads_Ratio | proportion of reads that can be divided into specific alleles |
| Replicate_eSDC | eSDC score of riboSNitch in replicates |
| Treat_Counts(A\|T\|C\|G) | separate read counts cover snp in treat file(s) |
| Control_Counts(A\|T\|C\|G) | separate read counts cover snp in control file(s) |
| Ref_Seq | base sequence of judging region |
| Ref_Score/Alt_Score | base scores of judging region |
| StructDiff_pValue | permutation pValues of bases in judging region |

## License
Licensed GPLv3 for open source use or contact zuoLab (zuozhx@sysucc.org.cn) for commercial use.
