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

Check to make sure multiqc ran and look at the multiQC output. The phred scores look really good, as does the duplication rate. There is some substantial adapter content so I will trim. Raw QC report can be found [here](https://github.com/fscucchia/Hawaii2022_pH_Temp_Mcap_Pacu/tree/main/output/QC). 

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

Submitted batch job 293017. Took a couple of hours, but looks like the fastqc step didn't work...Run fastqc as its own script. In scripts folder: `nano fastqc_trim.sh`

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

for file in /data/putnamlab/jillashey/Pacuta_HI_2022/data/trim/*fastq.gz
do 
fastqc $file --outdir /data/putnamlab/jillashey/Pacuta_HI_2022/output/fastqc/trim
done

cd /data/putnamlab/jillashey/Pacuta_HI_2022/output/fastqc/trim

multiqc *
```

Submitted batch job 293029


Count the number of reads each trimmed file has: 

```
zgrep -c "@LH00" *.gz

trim.A_R1_001.fastq.gz:22328993
trim.A_R2_001.fastq.gz:22328993
trim.B_R1_001.fastq.gz:24904172
trim.B_R2_001.fastq.gz:24904172
trim.C_R1_001.fastq.gz:25854719
trim.C_R2_001.fastq.gz:25854719
trim.D_R1_001.fastq.gz:24265173
trim.D_R2_001.fastq.gz:24265173
trim.E_R1_001.fastq.gz:26184558
trim.E_R2_001.fastq.gz:26184558
trim.F_R1_001.fastq.gz:23339918
trim.F_R2_001.fastq.gz:23339918
trim.G_R1_001.fastq.gz:25448242
trim.G_R2_001.fastq.gz:25448242
trim.H_R1_001.fastq.gz:23634641
trim.H_R2_001.fastq.gz:23634641
trim.L_R1_001.fastq.gz:23596961
trim.L_R2_001.fastq.gz:23596961
trim.M_R1_001.fastq.gz:23003891
trim.M_R2_001.fastq.gz:23003891
trim.N_R1_001.fastq.gz:22825941
trim.N_R2_001.fastq.gz:22825941
```

Once trimming is done, look at the multiqc plots to make sure adapter content is gone and sample quality is high. QC looks good! Trimmed QC report can be found [here](https://github.com/fscucchia/Hawaii2022_pH_Temp_Mcap_Pacu/tree/main/output/QC).

Align reads to Pacuta genome using hisat2. Download version 2 of the pocillopora genome [here](http://cyanophora.rutgers.edu/Pocillopora_acuta/) to Andromeda and unzip the assembly fasta file. In the scripts folder: `nano align.sh`

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Pacuta_HI_2022/scripts            
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

# load modules needed
module load HISAT2/2.2.1-foss-2019b #Alignment to reference genome: HISAT2
module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

echo "Building genome reference" $(date)

# index the reference genome for Pacuta output index to working directory
hisat2-build -f /data/putnamlab/jillashey/genome/Pacuta/V2/Pocillopora_acuta_HIv2.assembly.fasta Pacuta_ref
echo "Referece genome indexed. Starting alingment" $(date)

# This script exports alignments as bam files
# sorts the bam file because Stringtie takes a sorted file for input (--dta)
# removes the sam file because it is no longer needed

array=($(ls /data/putnamlab/jillashey/Pacuta_HI_2022/data/trim/*fastq.gz)) # call the clean sequences - make an array to align

for i in ${array[@]}; do
        sample_name=`echo $i| awk -F [.] '{print $2}'`
        hisat2 -p 8 --rna-strandness RF --dta -x Pacuta_ref -1 ${i} -2 $(echo ${i}|sed s/_R1/_R2/) -S ${sample_name}.sam
        samtools sort -@ 8 -o ${sample_name}.bam ${sample_name}.sam
                echo "${i} bam-ified!"
        rm ${sample_name}.sam
done

echo "Alignment complete!" $(date)
```

Submitted batch job 293038





