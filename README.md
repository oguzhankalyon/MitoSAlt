# MitoSAlt

This repository contains a modified version of the MitoSAlt Perl script, specifically adapted to start the analysis directly from **BAM files** instead of the original paired-end FASTQ input.

## Source & Documentation
This is a modified script for the **MitoSAlt** tool. For the full toolset, installation requirements, and detailed methodology, please visit the official SourceForge page:
👉 [MitoSAlt on SourceForge](https://sourceforge.net/projects/mitosalt/)

## Key Modifications
The internal logic of the Perl script has been updated to bypass the FASTQ alignment step and accept a single BAM file. 
* The script now treats the **BAM file** as the second argument.
* The study name (**$tag**) has been shifted to the third argument.

## Usage Comparison

### Original (Old) Usage
Originally, the tool required a configuration file and two separate FASTQ files:
```bash
perl script.pl <config_file> <fastq_1> <fastq_2> <tag>
