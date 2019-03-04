# Varpipe

Set of scripts for analyzing NGS data

## Software dependencies:

BEDtools Version 2.17.0

Bcftools Version 1.2 

BWA Version 0.7.12

GATK Version 3.8.0

Picard Version 1.134

Trimmomatic Version 0.36

Pigz Version 2.3.3

Qualimap Version 2.1.1

Samtools Version 1.2

SnpEff Version 4.3

Vcftools Version 0.1.12b

## Command line argument
usage: Varpipeline -q STRING -r STRING -n STRING [-q2 STRING] [-o STRING]
                   [--keepfiles] [--bwa] [--all] [--gatk] [--samtools] [-a]
                   [-v] [-h] [--version]

options: -q input gzipped fastq file
         -r reference genome (H37Rv NC_000962.3) file
         -n sample name
         -q2 paired end input gzipped fastq file
         -a annotate
         -v verbose
         -h help
         -- version tool version

         default mapping tool bwa 
         default variant calling tool GATK 
