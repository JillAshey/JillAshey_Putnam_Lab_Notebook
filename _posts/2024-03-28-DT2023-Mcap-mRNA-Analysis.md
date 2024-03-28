---
layout: post
title: Developmental 2023 Timeseries mRNA analysis (initial)
date: '2024-03-28'
categories: Analysis
tags: [Bioinformatics, mRNA, Montipora capitata]
projects: Developmental Timseries 2023
---

## Developmental 2023 Timeseries mRNA analysis (initial)

These data came from my developmental timeseries experiment in 2023 with *Montipora capitata* in Hawaii. In this experiment, *Montipora capitata* embryos and larvae were exposed to ambient and heat stress over 72 hours from embryo to swimming larvae. Samples were collected at 8 time points: 1, 4, 9, 14, 22, 28, 48 and 72 hpf. The github for this project is [here](https://github.com/JillAshey/Hawaii_Developmental_TimeSeries). 

To see which time points we want to sequence further, we ran an intial library prep and sequencing run on 8 samples (n=1 from each time point at ambient). The library prep was done by me using the [Zymo Ribofree](https://www.zymoresearch.com/pages/ngs-library-prep-kits#page-2) library prep kit (see my notebook [posts](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook) for more information on library prep). The sequencing was done through Oklahoma Medical Research Foundation NGS Core, which I found on Genohub. Here is the sequencing information: 

- Instrument: Illumina NovaSeq X Plus - 10B - PE 150 Cycle
- Read length: 2 x 150bp (Paired End)
- Pricing unit: per lane
- Number of samples: 8 libraries
- Guaranteed number of pass filter PE reads/sample: 20M (10M in each direction)
- Deliverables: FastQ files uploaded to Genohub project bucket

Files were downloaded to this location on Andromeda: `/data/putnamlab/KITT/hputnam/20240328_Mcap_RNASeq_Devo/`. Time to analyze! I'm going to write my notebook in chronological order by date and then will reorganize once the workflow is complete.

### 20240328

Sequences were received from sequencer today! Make a new directory in my own folder on Andromeda for this project: 

```
cd /data/putnamlab/jillashey
mkdir DT_Mcap_2023
cd DT_Mcap_2023
mkdir data scripts output
cd data 
mkdir raw trim 
cd raw
```

Copy the files into the raw data folder

```
cp /data/putnamlab/KITT/hputnam/20240328_Mcap_RNASeq_Devo/* .
```

Check md and initial raw QC was already done by Hollie. Raw QC is [here](https://github.com/JillAshey/Hawaii_Developmental_TimeSeries/blob/main/2023/data/Molecular/mRNA/20240328__Mcap_RNASeq_Devo_fastqc_multiqc_report.html). Overall, QC looks pretty good. The phred scores are high across the board. There is a decent amount of duplication but this is not unusual in early life stage samples. There is quite a bit of adapter content present so that definitely needs to come out. I'm going to use fastp to trim the reads. In the scripts folder: `nano fastp_qc.sh`

```
#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/DT_Mcap_2023/scripts            
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

# Load modules 
module load fastp/0.19.7-foss-2018b
module load FastQC/0.11.8-Java-1.8
module load MultiQC/1.9-intel-2020a-Python-3.8.2

cd /data/putnamlab/jillashey/DT_Mcap_2023/data/raw/

echo "Count number of reads in each file" $(date)

zgrep -c "@LH00" *.gz > raw_read_count.txt

echo "Start read trimming with fastp, followed by QC" $(date)

# Make an array of sequences to trim 
array1=($(ls *R1_001.fastq.gz))

# fastp and fastqc loop 
for i in ${array1[@]}; do
    fastp --in1 ${i} \
        --in2 $(echo ${i}|sed s/_R1/_R2/)\
        --out1 /data/putnamlab/jillashey/DT_Mcap_2023/data/trim/trim.${i} \
        --out2 /data/putnamlab/jillashey/DT_Mcap_2023/data/trim/trim.$(echo ${i}|sed s/_R1/_R2/) \
        --detect_adapter_for_pe \
        --qualified_quality_phred 30 \
        --unqualified_percent_limit 10 \
        --length_required 100 \
        --cut_right cut_right_window_size 5 cut_right_mean_quality 20 \
	 fastqc /data/putnamlab/jillashey/DT_Mcap_2023/output/fastqc/trim/trim.${i} \
    fastqc /data/putnamlab/jillashey/DT_Mcap_2023/output/fastqc/trim/trim.$(echo ${i}|sed s/_R1/_R2/)
done

echo "Read trimming of adapters and QC complete. Starting multiqc" $(date)

cd /data/putnamlab/jillashey/Pacuta_HI_2022/output/fastqc/trim
multiqc *

echo "MultiQC complete" $(date)
```

Submitted batch job 310216

Count number of reads in each trimmed file 

Align reads to Mcap genome using hisat2. I'm using Mcap genome V3. Download: http://cyanophora.rutgers.edu/montipora/ 