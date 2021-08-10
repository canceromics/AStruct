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
java -Xmx16g -jar circm6A.jar -h
Usage:
	java [-Xmx24g] -jar circm6a.jar -input <input.bam> -g <genome.fa> -o <path/out_prefix> [-ip ip.bam] [-r gencode.gtf] [options]
	
	<input.bam>	a bam/sam file of a sample mapping by bwa.
	<genome.fa>	a fasta file of the genome. The same as the file used by bwa is recommended.
	<path/out_prefix>	the path should exist and out_prefix is the first part of the output file name.
```

* Tips: -r enables exon boundary filter of circRNA detecting. -ip enables circle m6A peak detecting.

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
| Treat_Counts(A|T|C|G) | separate read counts cover snp in treat file(s) |
| Control_Counts(A|T|C|G) | separate read counts cover snp in control file(s) |
| Ref_Seq | base sequence of judging region |
| Ref_Score/Alt_Score | base scores of judging region |
| StructDiff_pValue | permutation pValues of bases in judging region |

## License
Licensed GPLv3 for open source use or contact zuoLab (zuozhx@sysucc.org.cn) for commercial use.
