---
layout: post
title: Pacuta HI 2022
date: '2024-01-18'
categories: Analysis
tags: [Bioinformatics, Pacuta, mRNA]
projects: Pacuta HI 2022
---

## Pacuta 2022 mRNA analysis

These data came from the Pacuta 2022 experiment in Hawaii, done by Federica and myself. In this experiment, larval and spat *Pocillopora acuta* were subjected to a combination of high pH and temperature treatments. The github for that project is [here](https://github.com/fscucchia/Hawaii2022_pH_Temp_Mcap_Pacu). 

Files were downloaded to this location: `/data/putnamlab/KITT/hputnam/20231127_Scucchia_HI`

Make a new directory in my own folder on Andromeda for this experiment

```
cd /data/putnamlab/jillashey
mkdir Pacuta_HI_2022
cd Pacuta_HI_2022
mkdir data scripts output
cd data 
mkdir raw trim 
cd raw
```

Copy the files into the raw data folder 

```
cp /data/putnamlab/KITT/hputnam/20231127_Scucchia_HI/* .
```

Check md has already been done by Hollie. Let's count how many reads each file has. 

```
zgrep -c "@LH00" *.gz

A_R1_001.fastq.gz:27485525
A_R2_001.fastq.gz:27485525
B_R1_001.fastq.gz:30895421
B_R2_001.fastq.gz:30895421
C_R1_001.fastq.gz:31829448
C_R2_001.fastq.gz:31829448
D_R1_001.fastq.gz:29779596
D_R2_001.fastq.gz:29779596
E_R1_001.fastq.gz:32423337
E_R2_001.fastq.gz:32423337
F_R1_001.fastq.gz:28875199
F_R2_001.fastq.gz:28875199
G_R1_001.fastq.gz:31635945
G_R2_001.fastq.gz:31635945
H_R1_001.fastq.gz:29455896
H_R2_001.fastq.gz:29455896
L_R1_001.fastq.gz:29250631
L_R2_001.fastq.gz:29250631
M_R1_001.fastq.gz:28581558
M_R2_001.fastq.gz:28581558
N_R1_001.fastq.gz:28609255
N_R2_001.fastq.gz:28609255
```

QC the raw data 

In scripts folder: nano fastqc_raw.sh

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Pacuta_HI_2022/scripts            
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

for file in /data/putnamlab/jillashey/Pacuta_HI_2022/data/raw/*fastq.gz
do 
fastqc $file --outdir /data/putnamlab/jillashey/Pacuta_HI_2022/output/fastqc/raw
done

cd /data/putnamlab/jillashey/Pacuta_HI_2022/output/fastqc/raw

multiqc *
```

Submitted batch job 293007

Check to make sure multiqc ran and look at the multiQC output. The phred scores look really good, as does the duplication rate. There is some substantial adapter content so I will trim. 

Trim reads w/ fastp. In scripts folder: `nano fastp_qc.sh`

```
#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Pacuta_HI_2022/scripts            
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

# Load modules needed 
module load fastp/0.19.7-foss-2018b
module load FastQC/0.11.8-Java-1.8
module load MultiQC/1.9-intel-2020a-Python-3.8.2

# Make an array of sequences to trim in raw data directory 

cd /data/putnamlab/jillashey/Pacuta_HI_2022/data/raw/
array1=($(ls *R1_001.fastq.gz))

echo "Read trimming of adapters started." $(date)

# fastp and fastqc loop 
for i in ${array1[@]}; do
    fastp --in1 ${i} \
        --in2 $(echo ${i}|sed s/_R1/_R2/)\
        --out1 /data/putnamlab/jillashey/Pacuta_HI_2022/data/trim/trim.${i} \
        --out2 /data/putnamlab/jillashey/Pacuta_HI_2022/data/trim/trim.$(echo ${i}|sed s/_R1/_R2/) \
        --detect_adapter_for_pe \
        --qualified_quality_phred 30 \
        --unqualified_percent_limit 10 \
        --length_required 100 \
        --cut_right cut_right_window_size 5 cut_right_mean_quality 20
	 fastqc /data/putnamlab/jillashey/Pacuta_HI_2022/output/fastqc/trim/trim.${i}
    fastqc /data/putnamlab/jillashey/Pacuta_HI_2022/output/fastqc/trim/trim.$(echo ${i}|sed s/_R1/_R2/)
done

echo "Read trimming of adapters and QC complete. Starting multiqc" $(date)

cd /data/putnamlab/jillashey/Pacuta_HI_2022/output/fastqc/trim
multiqc *

echo "MultiQC complete" $(date)
```

Submitted batch job 293017

Once trimming is done, look at the multiqc plots to make sure adapter content is gone and sample quality is high. Count the number of reads each trimmed file has: 

```
zgrep -c "@LH00" *.gz
```