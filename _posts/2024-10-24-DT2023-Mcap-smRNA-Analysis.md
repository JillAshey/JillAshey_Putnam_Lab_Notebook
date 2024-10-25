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
multiqc * 

echo "MultiQC complete" $(date)

echo "Count number of reads in each trimmed file" $(date)

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim
zgrep -c "@LH00" *.gz > smRNA_trim_read_count.txt
```

Submitted batch job 345007. The `-a` denotes the adapter sequence that I want to trim and its the Illumina TruSeq small RNA adapter, which Zymo uses in its kit. 

### 20241025

Cutadapt finished running overnight. Cutadapt trim QC is [here](https://github.com/JillAshey/DevelopmentalTimeseries/blob/main/data/Molecular/smRNA/multiqc_report_smRNA_cutadapt.html). The initial samples that were sequenced look good, but the ones sequenced most recently don't seem to have gotten any adapter trimmed off. Run fastp instead (based on Sam's [trimming](https://github.com/urol-e5/deep-dive/blob/2b978135383271bf5026e40039fdff99307480be/E-Peve/code/06.2-Peve-sRNAseq-trimming-31bp-fastp-merged.md#4-create-adapters-fasta-for-use-with-fastp-trimming) for e5). 

In the scripts folder: `nano fastp_trim.sh`

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
module load bio/fastp/0.23.2-GCC-11.2.0
module load FastQC/0.11.8-Java-1.8
module load MultiQC/1.9-intel-2020a-Python-3.8.2

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/

#echo "Count number of reads in each file" $(date)

#zgrep -c "@LH00" *.gz > smRNA_raw_read_count.txt

echo "Start read trimming with fastp, followed by QC" $(date)

# Make an array of sequences to trim 
array1=($(ls *R1_001.fastq.gz))

# fastp and fastqc loop 
for i in ${array1[@]}; do
    fastp --in1 ${i} \
        --in2 $(echo ${i}|sed s/_R1/_R2/)\
        --detect_adapter_for_pe \
        --qualified_quality_phred 30 \
        --trim_poly_g \
        --overlap_len_require 17 \
        --length_limit 31 \
        --merge \
        --merged_out /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim/trim.fastp.31bp.${i}
    fastqc /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim/trim.fastp.31bp.${i}
done
        
echo "Read trimming of adapters and QC complete. Starting multiqc" $(date)

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/fastqc/trim
multiqc *

echo "MultiQC complete" $(date)

echo "Count number of reads in each file" $(date)

zgrep -c "@LH00" *.gz > smRNA_trim_read_count.txt
```

Submitted batch job 345041. Seems to have worked but the merge wasn't amazing...seems like a lot of reads were lost when I merged because I required a decently high overlap length. Looking at the e5 [trimming iterations] for smRNAseq data, it also seems like I might be good to just use R1, similar to what Sam did with [e5 deep dive](https://github.com/urol-e5/deep-dive/blob/2b978135383271bf5026e40039fdff99307480be/E-Peve/code/06.1-Peve-sRNAseq-trimming-R1-only.md) data. 

In the scripts folder: `nano fastp_trim_R1.sh`

```
#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH -D /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/scripts            
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

# Load modules 
module load fastp/0.19.7-foss-2018b
module load FastQC/0.11.8-Java-1.8
module load MultiQC/1.9-intel-2020a-Python-3.8.2

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/

#echo "Count number of reads in each file" $(date)

#zgrep -c "@LH00" *.gz > smRNA_raw_read_count.txt

echo "Start read trimming with fastp, followed by QC" $(date)

# Make an array of sequences to trim 
array1=($(ls *R1_001.fastq.gz))

# fastp and fastqc loop 
for i in ${array1[@]}; do
    fastp --in1 ${i} \
    	 --out1 /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim/trim.fastp.31bp.${i} \
        --qualified_quality_phred 30 \
        --trim_poly_g \
        --length_limit 31 \
    fastqc /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim/trim.fastp.31bp.${i}
done

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/fastqc/trim
multiqc *

echo "MultiQC complete" $(date)

echo "Count number of reads in each file" $(date)

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim

zgrep -c "@LH00" *.gz > smRNA_trim_read_count.txt
```

Submitted batch job 345061. Started, then stopped it. Looking at the slurm error file, I am getting this information:

```
Detecting adapter sequence for read1...
No adapter detected for read1

Read1 before filtering:
total reads: 19773964
total bases: 2985868564
Q20 bases: 2596847950(86.9713%)
Q30 bases: 2097923022(70.2617%)

Read1 after filtering:
total reads: 475
total bases: 13245
Q20 bases: 13005(98.188%)
Q30 bases: 12269(92.6312%)

Filtering result:
reads passed filter: 475
reads failed due to low quality: 901386
reads failed due to too many N: 76
reads failed due to too short: 103
reads failed due to too long: 18871924
reads with adapter trimmed: 0
bases trimmed due to adapters: 0
```

Reads are failing left and right because they are too long...I think i need to specify what adapters that I want to trim because [fastp](https://github.com/OpenGene/fastp?tab=readme-ov-file#adapters) is trying to detect the adapters by analyzing the tails of the first ~1M reads. I'm going to supply it with adapter sequences instead. In my raw QC, it says some of the reads have the illumina small RNA 3' adapter and some have the illumina universal adapter. 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data

nano illumina_adapters.fasta
>Illumina TruSeq small RNA 3' adapter
TGGAATTCTCGGGTGCCAAGG
>Illumina TruSeq Adapter Read 1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
>Illumina TruSeq Adapter Read 2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
``` 

In the scripts folder, edit: `nano fastp_trim_R1.sh`

```
#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH -D /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/scripts            
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

# Load modules 
#module load fastp/0.19.7-foss-2018b
module load bio/fastp/0.23.2-GCC-11.2.0
module load FastQC/0.11.8-Java-1.8
#module load MultiQC/1.9-intel-2020a-Python-3.8.2

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/

#echo "Count number of reads in each file" $(date)

#zgrep -c "@LH00" *.gz > smRNA_raw_read_count.txt

echo "Start read trimming with fastp, followed by QC" $(date)

# Make an array of sequences to trim 
array1=($(ls *R1_001.fastq.gz))

# fastp and fastqc loop 
for i in ${array1[@]}; do
    fastp --in1 ${i} \
    	 --out1 /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim/trim.fastp.31bp.${i} \
        --qualified_quality_phred 30 \
        --trim_poly_x \
        #--length_limit 31 \
         --adapter_fasta /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta
    fastqc /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim/trim.fastp.31bp.${i}
done
```

Submitted batch job 345073. Ran in 1.5 hours, and then ran multiQC on the output. Also looked at the individual QC htmls. Most samples have good phred scores and a peak where the miRNAs should be in terms of length. But there is a lot of weirdness... a lot of overrepresented sequences that might be adapters? And still a decent amount of adapter content. I'm not sure what I should be trimming and i am stressed haha. I need to look into this further and think about what I'm really trimming. also evaluate tools to use...maybe email some ppl???? Also the files say 31bp but they are not. 





