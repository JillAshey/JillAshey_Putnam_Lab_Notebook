---
layout: post
title: E5 deep-dive into long non-coding RNAs
date: '2023-04-05'
categories: Analysis
tags: [Bioinformatics]
projects: E5 deep-dive 
---

This post details analysis and bioinformatic steps to generate a long non-coding RNA count matrix from RNAseq data for the E5 Deep Dive Project. More information on this project can be found [here](https://github.com/urol-e5/deep-dive). 

I'll be following a workflow from [Kang & Liu (2019)](https://link.springer.com/protocol/10.1007/978-1-4939-9045-0_13). They developed a workflow to identify long non-coding RNAs in diploid strawberries. 

### Set-up

AH & I already analyzed the RNASeq data to generate a gene count matrix (i.e., all QC and trimming steps have been done - see trimming results in this [post](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/E5-Deep-Dive-RNAseq-Count-Matrix-Analysis/)). I'm going to make a new directory for the lncRNA analysis and symbolically link the files that I'll need. 

```
cd /data/putnamlab/ashuffmyer/e5-deepdive
mkdir lnc-rna
cd lnc-rna
mkdir data scripts mapped assembled refs 

cd data 
ln -s /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/trim/trim.SRR* .

cd ../refs/
ln -s /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/refs/Pver_genome_assembly_v1.0.fasta .
ln -s /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/refs/Pver_genome_assembly_v1_fixed.gff3 .
```


### Index the genome 

`nano build.sh`

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/ashuffmyer/e5-deepdive/lnc-rna/scripts            
#SBATCH --error="build_error" #if your job fails, the error report will be put in this file
#SBATCH --output="build_output" #once your job is completed, any final job report comments will be put in this file

module load Bowtie2/2.4.4-GCC-11.2.0

bowtie2-build /data/putnamlab/ashuffmyer/e5-deepdive/lnc-rna/refs/Pver_genome_assembly_v1.0.fasta /data/putnamlab/ashuffmyer/e5-deepdive/lnc-rna/refs/Pver_ref

echo "Reference genome indexed." $(date)
```

Submitted batch job 246359

### Align sequences

`nano align.sh`

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/ashuffmyer/e5-deepdive/lnc-rna/scripts            
#SBATCH --error="align_error" #if your job fails, the error report will be put in this file
#SBATCH --output="align_output" #once your job is completed, any final job report comments will be put in this file

module load GCCcore/9.3.0
module load Bowtie2/2.4.4-GCC-11.2.0
module load TopHat/2.1.2-foss-2018b

## Genome already indexed

echo "Start alignment" $(date)

# Alignment of clean reads to the reference genome
# Make array of sequences to trim 
array1=($(ls /data/putnamlab/ashuffmyer/e5-deepdive/lnc-rna/data/*_1.fastq.gz | xargs -n 1 basename))
 
for i in ${array1[@]}; do
	tophat2 -p 10 --report-secondary-alignments -G /data/putnamlab/ashuffmyer/e5-deepdive/lnc-rna/refs/Pver_genome_assembly_v1_fixed.gff3 /data/putnamlab/ashuffmyer/e5-deepdive/lnc-rna/refs/Pver_ref /data/putnamlab/ashuffmyer/e5-deepdive/lnc-rna/data/${i} /data/putnamlab/ashuffmyer/e5-deepdive/lnc-rna/data/$(echo ${i}|sed s/_1/_2/)
done 

echo "Alignment completed!" $(date)

```

Submitted batch job 246362

NOT WORKING. Issue with GCC versions