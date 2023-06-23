---
layout: post
title: Astrangia 2021 mRNA analysis 
date: '2023-03-06'
categories: Analysis
tags: [Bioinformatics, Astrangia, mRNA]
projects: Astrangia 2021
---

## Astrangia 2021 mRNA analysis

These data came from my Astrangia 2021 experiment, during which adult Astrangia colonies were exposed to ambient and high temperatures for ~9 months. 

Files were downloaded to this location: `/data/putnamlab/KITT/hputnam/20230614_Astrangia_mRNA`

### Set-up
#### Make new folders for this project in my directory 

```
cd /data/putnamlab/jillashey
mkdir Astrangia2021
cd Astrangia2021
mkdir mRNA
cd mRNA
mkdir data scripts fastqc
cd data
mkdir raw
cd raw
```

#### Copy the new files into raw data folder 

```
cp /data/putnamlab/KITT/hputnam/20230614_Astrangia_mRNA/* .
```

#### Check to make sure all files were transferred successfully 

```
md5sum *.fastq.gz > checkmd5.md5
md5sum -c checkmd5.md5

AST-1065_R1_001.fastq.gz: OK
AST-1065_R2_001.fastq.gz: OK
AST-1105_R1_001.fastq.gz: OK
AST-1105_R2_001.fastq.gz: OK
AST-1147_R1_001.fastq.gz: OK
AST-1147_R2_001.fastq.gz: OK
AST-1412_R1_001.fastq.gz: OK
AST-1412_R2_001.fastq.gz: OK
AST-1560_R1_001.fastq.gz: OK
AST-1560_R2_001.fastq.gz: OK
AST-1567_R1_001.fastq.gz: OK
AST-1567_R2_001.fastq.gz: OK
AST-1617_R1_001.fastq.gz: OK
AST-1617_R2_001.fastq.gz: OK
AST-1722_R1_001.fastq.gz: OK
AST-1722_R2_001.fastq.gz: OK
AST-2000_R1_001.fastq.gz: OK
AST-2000_R2_001.fastq.gz: OK
AST-2007_R1_001.fastq.gz: OK
AST-2007_R2_001.fastq.gz: OK
AST-2302_R1_001.fastq.gz: OK
AST-2302_R2_001.fastq.gz: OK
AST-2360_R1_001.fastq.gz: OK
AST-2360_R2_001.fastq.gz: OK
AST-2398_R1_001.fastq.gz: OK
AST-2398_R2_001.fastq.gz: OK
AST-2404_R1_001.fastq.gz: OK
AST-2404_R2_001.fastq.gz: OK
AST-2412_R1_001.fastq.gz: OK
AST-2412_R2_001.fastq.gz: OK
AST-2512_R1_001.fastq.gz: OK
AST-2512_R2_001.fastq.gz: OK
AST-2523_R1_001.fastq.gz: OK
AST-2523_R2_001.fastq.gz: OK
AST-2563_R1_001.fastq.gz: OK
AST-2563_R2_001.fastq.gz: OK
AST-2729_R1_001.fastq.gz: OK
AST-2729_R2_001.fastq.gz: OK
AST-2755_R1_001.fastq.gz: OK
AST-2755_R2_001.fastq.gz: OK
```


#### Count number of reads per file 

```
zgrep -c "@NGSNJ" *fastq.gz
```

### Raw QC

#### Run fastqc to quality check raw reads 

In scripts folder: `nano fastqc_raw.sh`

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/mRNA/scripts              
#SBATCH --error="fastqc_raw_error" #if your job fails, the error report will be put in this file
#SBATCH --output="fastqc_raw_output" #once your job is completed, any final job report comments will be put in this file

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

for file in /data/putnamlab/jillashey/Astrangia2021/mRNA/data/raw/*fastq.gz
do 
fastqc $file --outdir /data/putnamlab/jillashey/Astrangia2021/mRNA/fastqc/raw
done

multiqc --interactive fastqc_results
```

Submitted batch job 265922
ADD QC PLOTS

### Trim data 

Using [fastp](https://github.com/OpenGene/fastp/tree/master) to trim data. 

In scripts folder: `nano fastp_QC.sh`

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/mRNA/scripts              
#SBATCH --error="fastp_QC_error" #if your job fails, the error report will be put in this file
#SBATCH --output="fastp_QC_output" #once your job is completed, any final job report comments will be put in this file

# Load modules needed 
module load fastp/0.19.7-foss-2018b
module load FastQC/0.11.8-Java-1.8
module load MultiQC/1.9-intel-2020a-Python-3.8.2

# Make an array of sequences to trim in raw data directory 

cd /data/putnamlab/jillashey/Astrangia2021/mRNA/data/raw
array1=($(ls *R1_001.fastq.gz))

echo "Read trimming of adapters started." $(date)

# fastp and fastqc loop 
for i in ${array1[@]}; do
    fastp --in1 ${i} \
        --in2 $(echo ${i}|sed s/_R1/_R2/)\
        --out1 /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/trimmed.${i} \
        --out2 /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/trimmed.$(echo ${i}|sed s/_R1/_R2/) \
        --detect_adapter_for_pe \
        --qualified_quality_phred 30 \
        --unqualified_percent_limit 10 \
        --length_required 100 \
        --cut_right cut_right_window_size 5 cut_right_mean_quality 20
	 fastqc /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/trimmed.${i}
    fastqc /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/trimmed.$(echo ${i}|sed s/_R1/_R2/)
done

echo "Read trimming of adapters complete." $(date)

# Quality Assessment of Trimmed Reads
cd /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/ #go to output directory

# Compile MultiQC report from FastQC files 
multiqc --interactive ./  

echo "Cleaned MultiQC report generated." $(date)
```

Submitted batch job 267173


