---
layout: post
title: Developmental 2023 Timeseries smRNA analysis 
date: '2024-10-24'
categories: Analysis
tags: [Bioinformatics, smRNA, Montipora capitata]
projects: Developmental Timseries 2023
---

## Developmental 2023 Timeseries smRNA analysis

These data came from my developmental timeseries experiment in 2023 with *Montipora capitata* in Hawaii. In this experiment, *Montipora capitata* embryos and larvae were exposed to ambient and heat stress over 72 hours from embryo to swimming larvae. Samples were collected at 8 time points: 1, 4, 9, 14, 22, 28, 48 and 72 hpf. The github for this project is [here](https://github.com/JillAshey/Hawaii_Developmental_TimeSeries). 

I initially sent 8 libraries (n=1 from each timepoint) to Oklahoma Medical Research Foundation NGS Core for sequencing and they looked good, so I sent out the rest of my sampels (24 samples, n=3 per timepoint). Data is now back! 

Files were downloaded to these location on Andromeda: XXXX and XXXX. Time to analyze! I'm going to write my notebook in chronological order by date and then will reorganize once the workflow is complete.

### 20241024

Sequences were received from sequencer today! Make a new directory in my own folder on Andromeda for this project: 

```
cd /data/putnamlab/jillashey
mkdir DT_Mcap_2023
cd DT_Mcap_2023
mkdir smRNA
cd dmRNA
mkdir data scripts output
cd data 
mkdir raw trim 
cd raw
```

Copy the files into the raw data folder

```
cp XXXXXX/* .
cp XXXXXX . 
```

nano fastqc_smallRNAseq_raw.sh

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/scripts            
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

echo "Start fastqc" $(date)

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw

for file in *fastq.gz
do 
fastqc $file
done

echo "Fastqc complete, start multiqc" $(date)

multiqc *

echo "multiqc complete!" $(date)
```

Submitted batch job 340170. Raw QC for smRNAseq data [here](https://github.com/JillAshey/DevelopmentalTimeseries/blob/main/data/Molecular/smRNA/multiqc_report_smRNA_raw.html). Raw QC looking very funky. There are definitely batch effects from the initial and last sequencing and a high level of duplication in the last sequencing run (not uncommon for bbs but also may be due to number of cycles during library prep). The adapter content plots look very different between the two batches as well. The Zymo library prep kit [protocol](https://files.zymoresearch.com/protocols/r3006_r3007-zymo-seq_mirna_library_kit.pdf) recommended single end sequencing at 75bp max. Oops I did paired end seq at 150bp. I believe that this means I will only be moving forward with R1. The protocol also recommends using [cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) to trim. 

In the scripts folder: `nano cutadapt_trim.sh`

```
#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/scripts            
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

# Load modules 
module load cutadapt/4.2-GCCcore-11.3.0
module load FastQC/0.11.8-Java-1.8

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/

#echo "Count number of reads in each file" $(date)

#zgrep -c "@LH00" *.gz > smRNA_raw_read_count.txt

echo "Start read trimming with cutadapt, followed by QC" $(date)

# Make an array of sequences to trim
array1=($(ls *R1_001.fastq.gz))

# cutadapt loop
for i in ${array1[@]}; do
    cutadapt -a TGGAATTCTCGGGTGCCAAGG \
        -u 3 \
        -m 15 \
        -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim/trim_${i} \
        $i
    
    fastqc /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim/trim_${i} \
        -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/fastqc
done

echo "Read trimming of adapters and QC complete. Starting multiqc" $(date)

module purge 
module load MultiQC/1.9-intel-2020a-Python-3.8.2

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/fastqc
multiqc* 

echo "MultiQC complete" $(date)

echo "Count number of reads in each trimmed file" $(date)

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim
zgrep -c "@LH00" *.gz > smRNA_trim_read_count.txt
```

Submitted batch job 345007. The `-a` denotes the adapter sequence that I want to trim and its the Illumina TruSeq small RNA adapter, which Zymo uses in its kit. 