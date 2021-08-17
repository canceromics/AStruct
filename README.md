<img src="icon.png" align="right" />

# diffRSS

diffRSS is a powerful tool for detection  different riboSNitch structure in alleles (SNPs).



## Table of Contents
* [Requirements](#Requirements)
* [Installation](#Installation)
* [Quick Start Guide](#QuickStart)
* [Usage](#Usage)
* [Example](#Example)
* [Output Headers](#OutputHeaders)
* [License](#License)

## Requirements

* JDK 8

## Installation

* This tool can be installed by instructions as follows:

```
git clone https://github.com/canceromics/diffRSS.git
cd diffRSS/tool/src
sh make.sh
cd ../..
```
The tool is generated as diffRSS.jar in this directory.

## QuickStart

* Start from bam file of and input sample for example.

```
java -Xmx16g -jar diffRSS.jar -tf example/data/T1.bam example/data/T2.bam -cf example/data/C1.bam example/data/C2.bam -o example/result/Example -gf example/reference/GRCh38.p12_NC_000022.11.fa -ef example/reference/GRCh38.p12_NC_000022.11.gff 
-mutf example/reference/SNP_NC_000022.11.vcf -repli
```
Running this instruction will result in getting a file named Example_riboSNitch.txt. `example/result/Example` means output_dir/file_prefix

More commonly used

```
java [-Xmx16g] -jar diffRSS.jar -tf <treat.bam(s)> -gf <genome.fa> -ef <gencode.gtf/gff> -mutf <snp.cvf> -o <path/out_prefix> [-cf control.bam(s)] [options]
```

## Usage

* More details of this tool can be found with -h parameter

```
java [-Xmx16g] -jar diffRSS.jar <method> [options]
Methods:
    rtmut
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
        
    rtmut
    find different structure between 2 samples. (all bams need to be sorted by coordinate)
    options:
        <-tf treat.bam ...> treat bam files
        <-cf control.bam ...> control bam files
        <-tf2 treat.bam ...> homo treat bam files
        <-cf2 control.bam ...> homo control bam files
        <-ef gencode.gtf> annotation gtf/gff file
        <-gf genome.fa> genome fasta file
        <-o out_prefix> output file prefix
        [-mutf snp.vcf] snp vcf/bed file (use gene region without this, overwrites -peakf parameter)
        [-peakf region.bed] region bed file
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
        [-peakf m6Apeak.bed] bed3 file that include m6Apeak regions.
        [-circf circRNA.bed] bed3 file that include circRNA regions.
        [-sf struct.bed] bed file with structure in 4th column.
        [-alen length] alignment length. (default is 150, must bigger than 0)
        [-rlen length] read fragment length. (default is 300, must bigger than alignment length) (Tips: single end simulation enabled when alignment length is the same as read fragment length)
        [-minalen length] minimum alignment length which sequence lower than this length would be dropped. (default is 30)
        [-faf frequency] enable and set fixed allele frequency. (default is disabled)
        [-peak] enable simulation with m6Apeaks.
        [-genef gene.bed] only the genes in this file can be used in simulation, and IP enrichment, IP read count, input read count can also be set in this file.
        [-lnp proportion] set linear reads proportion against circle reads. (default is 50)
        [-share] enable more than 1 circle m6Apeak in the same gene.
        [-nolinear] disable m6Apeak simulation when m6Apeaks do not overlap with circRNAs.
        [-npback proportion] set background proportion in non peak regions in IP against to input. (default is 0.01)
        [-minbr count] minimum background reads count in IP. (default is 0)
        [-enrich proportion] set m6Apeaks enrichment proportion. (default is calculated to ensure the same size between IP and input files.
        [-bases bases] allow structures show in some bases such as "AC". (default is all bases)
        [-mut] change structure simulation from RTstop to RTmutation.
```

## Example

```
java -Xmx16g -jar diffRSS.jar -tf example/data/T1.bam example/data/T2.bam -cf example/data/C1.bam example/data/C2.bam -o example/result/Example -gf example/reference/GRCh38.p12_NC_000022.11.fa -ef example/reference/GRCh38.p12_NC_000022.11.gff 
-mutf example/reference/SNP_NC_000022.11.vcf -repli
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

* Here are definitions of headers in output file named `(output_dir/file_prefix)riboSNitchDetail.txt`

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
