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

AST-1065_R1_001.fastq.gz:42852873
AST-1065_R2_001.fastq.gz:42852873
AST-1105_R1_001.fastq.gz:36082209

gzip: AST-1105_R2_001.fastq.gz: unexpected end of file
AST-1105_R2_001.fastq.gz:409
AST-1147_R1_001.fastq.gz:41008372
AST-1147_R2_001.fastq.gz:41008372
AST-1412_R1_001.fastq.gz:33078120
AST-1412_R2_001.fastq.gz:33078120
AST-1560_R1_001.fastq.gz:35909631
AST-1560_R2_001.fastq.gz:35909631
AST-1567_R1_001.fastq.gz:36789468
AST-1567_R2_001.fastq.gz:36789468
AST-1617_R1_001.fastq.gz:31764172
AST-1617_R2_001.fastq.gz:31764172
AST-1722_R1_001.fastq.gz:34941521
AST-1722_R2_001.fastq.gz:34941521
AST-2000_R1_001.fastq.gz:30315074
AST-2000_R2_001.fastq.gz:30315074
AST-2007_R1_001.fastq.gz:39369837
AST-2007_R2_001.fastq.gz:39369837
AST-2302_R1_001.fastq.gz:33326148
AST-2302_R2_001.fastq.gz:33326148
AST-2360_R1_001.fastq.gz:32840507
AST-2360_R2_001.fastq.gz:32840507
AST-2398_R1_001.fastq.gz:33309047
AST-2398_R2_001.fastq.gz:33309047
AST-2404_R1_001.fastq.gz:38349371
AST-2404_R2_001.fastq.gz:38349371
AST-2412_R1_001.fastq.gz:35436774
AST-2412_R2_001.fastq.gz:35436774
AST-2512_R1_001.fastq.gz:35323599
AST-2512_R2_001.fastq.gz:35323599
AST-2523_R1_001.fastq.gz:34086320
AST-2523_R2_001.fastq.gz:34086320
AST-2563_R1_001.fastq.gz:35799766
AST-2563_R2_001.fastq.gz:35799766
AST-2729_R1_001.fastq.gz:22735638
AST-2729_R2_001.fastq.gz:22735638
AST-2755_R1_001.fastq.gz:33057601
AST-2755_R2_001.fastq.gz:33057601
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

QC results: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/astrangia2021_bioinf/AST2021_mRNA_raw_fastqc_sequence_counts_plot.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/astrangia2021_bioinf/AST2021_mRNA_raw_fastqc_per_base_sequence_quality_plot.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/astrangia2021_bioinf/AST2021_mRNA_raw_fastqc_overrepresented_sequencesi_plot.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/astrangia2021_bioinf/AST2021_mRNA_raw_fastqc_per_sequence_gc_content_plot.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/astrangia2021_bioinf/AST2021_mRNA_raw_fastqc_per_sequence_quality_scores_plot.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/astrangia2021_bioinf/AST2021_mRNA_raw_fastqc_sequence_duplication_levels_plot.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/astrangia2021_bioinf/AST2021_mRNA_raw_fastqc_adapter_content_plot.png)

The quality scores all look really good, happy about that. There is quite a bit of adapter content though.

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

ADD QC PLOTS

### Align trimmed data to reference genome 

Should I use hisat2 or bowtie? Only considering using bowtie because that's used in a lot of ncRNA analyses. Let's do both!

Make new folders in output directory

```
cd /data/putnamlab/jillashey/Astrangia2021/mRNA/output/
mkdir hisat2 bowtie
cd hisat2
mkdir refs align
cd ../bowtie 
mkdir refs align
```

#### HISAT2 alignment

First, index the reference genome. 

In scripts folder: `nano hisat2_build.sh`

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=15
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/mRNA/scripts              
#SBATCH --error="hisat2_build_error" #if your job fails, the error report will be put in this file
#SBATCH --output="hisat2_build_output" #once your job is completed, any final job report comments will be put in this file

# load modules needed
module load HISAT2/2.2.1-foss-2019b #Alignment to reference genome: HISAT2

echo "Indexing reference genome" $(date)

# Index the reference genome for A. poculata 
hisat2-build -f /data/putnamlab/jillashey/Astrangia_Genome/apoculata.assembly.scaffolds_chromosome_level.fasta ./Apoc_ref

echo "Referece genome indexed!" $(date)
```

Submitted batch job 267749. After this is done running, move reference files from scripts folder to hisat2 refs folder: `mv Apoc_ref.* ../output/hisat2/refs/`

Use the indexed reference genome to align sequences 

In scripts folder: `nano hisat2_align.sh`

```
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=15
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/mRNA/scripts              
#SBATCH --error="hisat2_align_error" #if your job fails, the error report will be put in this file
#SBATCH --output="hisat2_align_output" #once your job is completed, any final job report comments will be put in this file

# load modules needed
module load HISAT2/2.2.1-foss-2019b #Alignment to reference genome: HISAT2
module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

echo "Aligning reads to reference genome" $(date)

# This script exports alignments as bam files, sorts the bam file because Stringtie takes a sorted file for input (--dta), and removes the sam file because it is no longer needed

array=($(ls /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/*fastq.gz)) # call the clean sequences - make an array to align

for i in ${array[@]}; do
    hisat2 -p 8 --rna-strandness RF --dta -q -x ../output/hisat2/refs/Apoc_ref -1 ${i} -2 $(echo ${i}|sed s/_R1/_R2/) -S ${i}.sam
    samtools sort -@ 8 -o ${i}.bam ${i}.sam
    	echo "${i} bam-ified!"
    rm ${i}.sam
done

echo "Alignment complete!" $(date)
```

Submitted batch job 267750 - this job was completed, but I think it didn't work properly. In the hisat2 error file, there was some kind of GCC conflict error. Purged modules from the environment, then ran job again. Submitted batch job 269176 on 6/30/23. Checked the error file 30 seconds after running and got the same errors: 

```
foss/2018b(23):ERROR:150: Module 'foss/2018b' conflicts with the currently loaded module(s) 'foss/2019b'
foss/2018b(23):ERROR:102: Tcl command execution failed: conflict foss

GCCcore/7.3.0(23):ERROR:150: Module 'GCCcore/7.3.0' conflicts with the currently loaded module(s) 'GCCcore/8.3.0'
GCCcore/7.3.0(23):ERROR:102: Tcl command execution failed: conflict GCCcore

GCCcore/7.3.0(23):ERROR:150: Module 'GCCcore/7.3.0' conflicts with the currently loaded module(s) 'GCCcore/8.3.0'
GCCcore/7.3.0(23):ERROR:102: Tcl command execution failed: conflict GCCcore

GCCcore/7.3.0(23):ERROR:150: Module 'GCCcore/7.3.0' conflicts with the currently loaded module(s) 'GCCcore/8.3.0'
GCCcore/7.3.0(23):ERROR:102: Tcl command execution failed: conflict GCCcore

GCCcore/7.3.0(23):ERROR:150: Module 'GCCcore/7.3.0' conflicts with the currently loaded module(s) 'GCCcore/8.3.0'
GCCcore/7.3.0(23):ERROR:102: Tcl command execution failed: conflict GCCcore
```

Stopping the job now. I don't think it likes the modules that I'm trying to load. Use these modules instead: `HISAT2/2.2.1-gompi-2021b` and `SAMtools/1.16.1-GCC-11.3.0`. Submitted batch job 269177. still getting a conflict error. Adding this line before loading modules: `source /usr/share/Modules/init/sh # load the module function (specific need to my computer)` from Emma's [script](https://github.com/emmastrand/EmmaStrand_Notebook/blob/master/_posts/2022-02-03-KBay-Bleaching-Pairs-RNASeq-Pipeline-Analysis.md). Not sure what it does but maybe it'll help??? Submitted batch job 269178. nope same error. Trying `HISAT2/2.2.1-gompi-2021b` and `SAMtools/1.15.1-GCC-11.2.0` (from Emma's script linked above). Submitted batch job 269179. This worked! 

View mapping percentages for the hisat2 alignment 

```
module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

for i in *.bam; do
    echo "${i}" >> mapped_reads_counts_Apoc
    samtools flagstat ${i} | grep "mapped (" >> mapped_reads_counts_Apoc
done
```

#### Bowtie2 alignment

First, index the reference genome. 

In scripts folder: `nano bowtie2_build.sh`

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=15
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/mRNA/scripts              
#SBATCH --error="bowtie2_build_error" #if your job fails, the error report will be put in this file
#SBATCH --output="bowtie2_build_output" #once your job is completed, any final job report comments will be put in this file

# load modules needed
module load Bowtie2/2.4.4-GCC-11.2.0 

echo "Indexing reference genome" $(date)

# Index the reference genome for A. poculata 
bowtie2-build /data/putnamlab/jillashey/Astrangia_Genome/apoculata.assembly.scaffolds_chromosome_level.fasta Apoc_ref.btindex

echo "Referece genome indexed!" $(date)
```

Submitted batch job 269279


NEED TO ASK HP: WAS THE SEQUENCING DONE RF, FF OR FR?