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

Let's look at the counts from the trimmed data 

```
zgrep -c "@NGSNJ" *fastq.gz
trimmed.AST-1065_R1_001.fastq.gz:0
trimmed.AST-1065_R2_001.fastq.gz:28514704
trimmed.AST-1105_R1_001.fastq.gz:0
trimmed.AST-1105_R2_001.fastq.gz:283
trimmed.AST-1147_R1_001.fastq.gz:0
trimmed.AST-1147_R2_001.fastq.gz:27401584
trimmed.AST-1412_R1_001.fastq.gz:0
trimmed.AST-1412_R2_001.fastq.gz:23324255
trimmed.AST-1560_R1_001.fastq.gz:0
trimmed.AST-1560_R2_001.fastq.gz:25261438
trimmed.AST-1567_R1_001.fastq.gz:0
trimmed.AST-1567_R2_001.fastq.gz:23711472
trimmed.AST-1617_R1_001.fastq.gz:0
trimmed.AST-1617_R2_001.fastq.gz:21247535
trimmed.AST-1722_R1_001.fastq.gz:0
trimmed.AST-1722_R2_001.fastq.gz:24754153
trimmed.AST-2000_R1_001.fastq.gz:0
trimmed.AST-2000_R2_001.fastq.gz:20483353
trimmed.AST-2007_R1_001.fastq.gz:0
trimmed.AST-2007_R2_001.fastq.gz:26160891
trimmed.AST-2302_R1_001.fastq.gz:0
trimmed.AST-2302_R2_001.fastq.gz:22054499
trimmed.AST-2360_R1_001.fastq.gz:0
trimmed.AST-2360_R2_001.fastq.gz:21377658
trimmed.AST-2398_R1_001.fastq.gz:0
trimmed.AST-2398_R2_001.fastq.gz:21796937
trimmed.AST-2404_R1_001.fastq.gz:0
trimmed.AST-2404_R2_001.fastq.gz:26894798
trimmed.AST-2412_R1_001.fastq.gz:0
trimmed.AST-2412_R2_001.fastq.gz:25372505
trimmed.AST-2512_R1_001.fastq.gz:0
trimmed.AST-2512_R2_001.fastq.gz:24221309
trimmed.AST-2523_R1_001.fastq.gz:0
trimmed.AST-2523_R2_001.fastq.gz:23751953
trimmed.AST-2563_R1_001.fastq.gz:0
trimmed.AST-2563_R2_001.fastq.gz:25242216
trimmed.AST-2729_R1_001.fastq.gz:0
trimmed.AST-2729_R2_001.fastq.gz:15551660
trimmed.AST-2755_R1_001.fastq.gz:0
trimmed.AST-2755_R2_001.fastq.gz:24720522
```

Why do the R1 files have 0 reads??

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

Move BAM files in trimmed folder to hisat align folder

```
mv *bam ../../output/hisat2/align/
```

Should the output should be R1 and R2 in one sam file...why is this not the case with the hisat2 samples? if not now, when do they get combined into a single file?

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

After this is done running, move reference files to bowtie folder: `mv *.bt2 ../output/bowtie/refs/`

Use the indexed reference genome to align sequences 

In scripts folder: `nano bowtie_align.sh`

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
#SBATCH --error="bowtie_align_error" #if your job fails, the error report will be put in this file
#SBATCH --output="bowtie_align_output" #once your job is completed, any final job report comments will be put in this file

module load Bowtie2/2.4.4-GCC-11.2.0 

echo "Aligning reads" $(date)

array=($(ls /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/*fastq.gz)) # call the clean sequences - make an array to align

for i in ${array[@]}; do
    bowtie2 -1 ${i} -2 $(echo ${i}|sed s/_R1/_R2/) -x /data/putnamlab/jillashey/Astrangia2021/mRNA/output/bowtie/refs/Apoc_ref.btindex -b --align-paired-reads
done

echo "Reads aligned!" $(date)
```

Submitted batch job 275228. Didn't work. Getting this error: 

```
0 reads
0.00% overall alignment rate
Warning: Same mate file "/data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/trimmed.AST-1065_R2_001.fastq.gz" appears as argument to both -1 and -2
0 reads
```

Maybe it didn't like the `--align-paired-reads` argument? Going to remove this and rerun the code. Submitted batch job 275230. Got the same error... Editing the code so that the array is `array=($(ls /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/*_R1_001.fastq.gz))` based on Emma's [code](https://github.com/emmastrand/EmmaStrand_Notebook/blob/master/_posts/2022-02-03-KBay-Bleaching-Pairs-RNASeq-Pipeline-Analysis.md). Submitted batch job 275237. Failed again, still saying 0% alignment rate. In the `-x` argument, changed `Apoc_ref.btindex` to `Apoc_ref`. Submitted batch job 275239. It didn't like that! Changing it back. 

Still getting 0% alignment. The error file says: "Error while reading BAM extra subfields
(ERR): bowtie2-align exited with value 1". Maybe it's an issue with the path? Changing `-x /data/putnamlab/jillashey/Astrangia2021/mRNA/output/bowtie/refs/Apoc_ref.btindex` to `-x ../output/bowtie/refs/Apoc_ref.btindex`. Submitted batch job 275248. Still getting the same error??? 

It's 20231004. Let's get back to aligning with bowtie2. This is the code I have now:

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
#SBATCH --error="bowtie_align_error" #if your job fails, the error report will be put in this file
#SBATCH --output="bowtie_align_output" #once your job is completed, any final job report comments will be put in this file

module load Bowtie2/2.4.4-GCC-11.2.0 

echo "Aligning reads" $(date)

array=($(ls /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/*_R1_001.fastq.gz)) # call the clean sequences - make an array to align

for i in ${array[@]}; do
    bowtie2 -1 ${i} -2 $(echo ${i}|sed s/_R1/_R2/) -x ../output/bowtie/refs/Apoc_ref.btindex -b 
done

echo "Reads aligned!" $(date)
```

Let's try running it again and seeing what happens. Submitted batch job 283318. Nope same error as above. I may need to add a bam file name. Added `-b ${i}`. Submitted batch job 283321. Didn't work, this is the error that I'm getting: `Warning: Output file '{i}' was specified without -S.  This will not work in future Bowtie 2 versions.  Please use -S instead.` `-S` means to write the output files as sam files, but I don't want them as sam files. I'm going to add the `-S` anyway and see what happens. I can always convert the sam to bam files. Submitted batch job 283322

Now I'm getting this error: 

```
Extra parameter(s) specified: "/data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/trimmed.AST-1412_R1_001.fastq.gz"
Note that if <mates> files are specified using -1/-2, a <singles> file cannot
also be specified.  Please run bowtie separately for mates and singles.
Error: Encountered internal Bowtie 2 exception (#1)
Command: /opt/software/Bowtie2/2.4.4-GCC-11.2.0/bin/bowtie2-align-s --wrapper basic-0 -x ../output/bowtie/refs/Apoc_ref.btindex -S /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/trimmed.AST-1412_R1_001.fastq.gz -1 /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/trimmed.AST-1412_R1_001.fastq.gz -2 /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/trimmed.AST-1412_R2_001.fastq.gz -b /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/trimmed.AST-1412_R1_001.fastq.gz 
(ERR): bowtie2-align exited with value 1
```

It looks like it thinks I'm trying to run a singles file but its not working...Going back! Changing the syntax of the array a little bit and removing the bam and sam commands

```
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=15
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/mRNA/scripts              
#SBATCH --error="bowtie_align_error" #if your job fails, the error report will be put in this file
#SBATCH --output="bowtie_align_output" #once your job is completed, any final job report comments will be put in this file

module load Bowtie2/2.4.4-GCC-11.2.0 

echo "Aligning reads" $(date)

cd /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/

array=($(ls /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/*_R1_001.fastq.gz)) # call the clean sequences - make an array to align

for i in ${array[@]}; do
    bowtie2 -1 ${i} -2 $(echo ${i}|sed s/_R1/_R2/) -x /data/putnamlab/jillashey/Astrangia2021/mRNA/output/bowtie/refs/Apoc_ref.btindex
done

echo "Reads aligned!" $(date)
```

Submitted batch job 283326. Now I got this error: 

```
Error, fewer reads in file specified with -1 than in file specified with -2
terminate called after throwing an instance of 'int'
(ERR): bowtie2-align died with signal 6 (ABRT) 
Unable to read file magic number
```

Unsure why bowtie is so angry at me. boo. Still need to troubleshoot as of 10/4/23. 

20240103

I may have found the source of the problem...when I trim my data, for some reason, the R1 files now have no data in them. Need to go back to trimming step. 


### Assemble and quantify transcripts

#### Hisat2

Going to be using stringtie to assemble and quantify transcripts 

```
cd /data/putnamlab/jillashey/Astrangia2021/mRNA
mkdir stringtie
cd scripts 
```

In the scripts folder: `nano assemble.sh`

```
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=15
#SBATCH --export=NONE
#SBATCH --mem=128GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/mRNA/output/hisat2/align              
#SBATCH --error="stringtie-hisat-assemble_error" #if your job fails, the error report will be put in this file
#SBATCH --output="stringtie-hisat-assemble_error" #once your job is completed, any final job report comments will be put in this file

module load StringTie/2.2.1-GCC-11.2.0

echo "Starting stringtie assembly" $(date)

# Transcript assembly: StringTie

array1=($(ls *.bam)) #Make an array of sequences to assemble

for i in ${array1[@]}; do
    sample_name=`echo $i| awk -F [_] '{print $1"_"$2"_"$3}'`
    stringtie -p 8 -e -B -G /data/putnamlab/jillashey/Astrangia_Genome/apoculata_v2.0.gff3 -A ${sample_name}.gene_abund.tab -o ${sample_name}.gtf ${i}
    echo "StringTie assembly for seq file ${i}" $(date)
done

echo "Stringtie assembly complete" $(date)
```

Submitted batch job 275232. Job was pending for about two days, but then ran in about an hour. 

Move GTF and gene abundance files into stringtie folder, as they are currently in the hisat2 align folder

```
cd /data/putnamlab/jillashey/Astrangia2021/mRNA/output/hisat2/align
mv *gtf ../../../stringtie/
mv *gene_abund.tab ../../../stringtie/
mv *.ctab ../../../stringtie/
mv stringtie-hisat-assemble_error ../../../scripts/
```

Make the gene count matrix! This step uses a script called prepDE.py that can be downloaded/copied from the [StringTie github repo](https://github.com/gpertea/stringtie/blob/master/prepDE.py). Copy the script into the scripts folder. 

```
cd /data/putnamlab/jillashey/Astrangia2021/mRNA/scripts
nano prepDE.py # copy script from github into here and save
```

In the scripts folder: `nano prepDE.sh`

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=200GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab              
#SBATCH --error="prepDE_error" #if your job fails, the error report will be put in this file
#SBATCH --output="prepDE_output" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/mRNA/stringtie

echo "Starting assembly and assembly QC" $(date)

# Load packages
module load Python/2.7.15-foss-2018b #Python
module load StringTie/2.2.1-GCC-11.2.0
module load GffCompare/0.12.6-GCC-11.2.0 #Transcript assembly QC: GFFCompare

# Make gtf_list.txt file
ls *.gtf > gtf_list.txt

stringtie --merge -e -p 8 -G /data/putnamlab/jillashey/Astrangia_Genome/apoculata_v2.0.gff3 -o Apoc_merged.gtf gtf_list.txt #Merge GTFs 
echo "Stringtie merge complete" $(date)

gffcompare -r /data/putnamlab/jillashey/Astrangia_Genome/apoculata_v2.0.gff3 -G -o merged Cvir_merged.gtf #Compute the accuracy 
echo "GFFcompare complete, Starting gene count matrix assembly..." $(date)

#Note: the merged part is actually redundant and unnecessary unless we perform the original stringtie step without the -e function and perform
#re-estimation with -e after stringtie --merge, but will redo the pipeline later and confirm that I get equal results.

#make gtf list text file
for filename in *bam.gtf; do echo $filename $PWD/$filename; done > listGTF.txt

python /data/putnamlab/jillashey/Astrangia2021/mRNA/scripts/prepDE.py -g Apoc_mRNA_count_matrix.csv -i listGTF.txt #Compile the gene count matrix

echo "Gene count matrix compiled." $(date)
```

The above code was written by Zoe in this [post](https://github.com/zdellaert/ZD_Putnam_Lab_Notebook/blob/master/_posts/2023-02-27-Point-Judith-RNAseq.md). Submitted batch job 275703





### 20240103 

I'm going to go in date order and then rewrite and make it pretty. Today, I figured out that the trimming step resulted in the R1 files having 0 reads. Let's try to fix that. When I look at the trimmed R1 files, it says: 

```
@HD     VN:1.0  SO:unsorted
@SQ     SN:chromosome_1 LN:21106950
@SQ     SN:chromosome_2 LN:21532011
@SQ     SN:chromosome_3 LN:22725890
@SQ     SN:chromosome_4 LN:27282911
@SQ     SN:chromosome_5 LN:29602567
@SQ     SN:chromosome_6 LN:30212294
@SQ     SN:chromosome_7 LN:30683508
@SQ     SN:chromosome_8 LN:29863146
@SQ     SN:chromosome_9 LN:31044737
@SQ     SN:chromosome_10        LN:33861578
@SQ     SN:chromosome_11        LN:33534404
@SQ     SN:chromosome_12        LN:42634786
@SQ     SN:chromosome_13        LN:42928715
@SQ     SN:chromosome_14        LN:58409227
@PG     ID:bowtie2      PN:bowtie2      VN:2.4.4        CL:"/opt/software/Bowtie2/2.4.4-GCC-11.2.0/bin/bowtie2-align-s --wrapper basic-0 -x ../output/bowtie/refs/Apoc_ref.btindex -1 /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/trimmed.AST-1065_R1_001.fastq.gz -2 /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/trimmed.AST-1065_R2_001.fastq.gz -b /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/trimmed.AST-1065_R1_001.fastq.gz"
```

Did I accidently replace my data with bowtie info? I'm going to check to make sure the raw reads look okay still and then rerun fastp. I'm going to try trimming just one sample first in interactive mode in the trim folder. 

```
module load fastp/0.19.7-foss-2018b

fastp --in1 ../raw/AST-1065_R1_001.fastq.gz \
--in2 ../raw/AST-1065_R2_001.fastq.gz \
--out1 test.AST-1065_R1_001.fastq.gz \
--out2 test.AST-1065_R2_001.fastq.gz \
--detect_adapter_for_pe \
--qualified_quality_phred 30 \
--unqualified_percent_limit 10 \
--length_required 100 \
--cut_right cut_right_window_size 5 cut_right_mean_quality 20
```

Didn't finish running, but it seems like this is the way to go. Rerun `fastp_QC.sh`. Submitted batch job 291968

Now I'm going to attempt to run bowtie2 again. I already made the bowtie2 genome index above so lets align the data. I need to be careful not to overwrite my trimmed files with bowtie info. First, I'm going to try with just one sample: 

```
interactive 

module load GCCcore/11.2.0 #I needed to add this to resolve conflicts between loaded GCCcore/9.3.0 and GCCcore/11.3.0
module load Bowtie2/2.4.4-GCC-11.2.0 

bowtie2 -1 trimmed.AST-1065_R1_001.fastq.gz -2 trimmed.AST-1065_R2_001.fastq.gz -x /data/putnamlab/jillashey/Astrangia2021/mRNA/output/bowtie/refs/Apoc_ref.btindex -S test
```

Batch job 291993. This appears to have been successful. Let's try to convert the test sam file into a bam file. 

```
module load SAMtools/1.9-foss-2018b
samtools sort -@ 8 -o test.bam test.sam
```

Gives me this error: `Illegal instruction (core dumped)`. 

### 20240104

Trimming finished last night and plots look good. Now I'm going to run bowtie2 to align all samples. Here's the bowtie2 code for all of the samples. In scripts: `nano bowtie2_align.sh`

```
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=15
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/mRNA/scripts              
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.out

module load Bowtie2/2.4.4-GCC-11.2.0 

echo "Aligning reads" $(date)

cd /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/

array=($(ls *R1_001.fastq.gz)) # call the clean sequences - make an array to align

for i in ${array[@]}; do
    bowtie2 -1 ${i} -2 $(echo ${i}|sed s/_R1/_R2/) -x /data/putnamlab/jillashey/Astrangia2021/mRNA/output/bowtie/refs/Apoc_ref.btindex -S -b --align-paired-reads
done

echo "Reads aligned!" $(date)
```

Submitted batch job 292070. failed immediately. Error was: `Warning: Output file '/data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/test.AST-1065_R1_001.fastq.gz.bam' was specified without -S.  This will not work in future Bowtie 2 versions.  Please use -S instead.`. Going to add `-S` even though I want a bam file. Submitted batch job 292080. Now this error: 

```
Extra parameter(s) specified: "/data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/test.AST-1065_R1_001.fastq.gz.bam"
Note that if <mates> files are specified using -1/-2, a <singles> file cannot
also be specified.  Please run bowtie separately for mates and singles.
Error: Encountered internal Bowtie 2 exception (#1)
```

I think it's angry that I specified a single output file for bam so I'm going to remove the file name. Submitted batch job 292081

Now it says this: 

```
Error: 0 mate files/sequences were specified with -1, but 1
mate files/sequences were specified with -2.  The same number of mate files/
sequences must be specified with -1 and -2.
Error: Encountered internal Bowtie 2 exception (#1)
Command: /opt/software/Bowtie2/2.4.4-GCC-11.2.0/bin/bowtie2-align-s --wrapper basic-0 -x /data/putnamlab/jillashey/Astrangia2021/mRNA/output/bowtie/refs/Apoc_ref.btindex -S -1 -2 /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/test.AST-1065_R2_001.fastq.gz -b --align-paired-reads /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/test.AST-1065_R1_001.fastq.gz 
(ERR): bowtie2-align exited with value 1
```

but there is a file in the -1 position? Maybe the array doesn't like having the full path name in the array itself? Going to remove that and rerun. Submitted batch job 292082. Same error as above... This is what it thinks my command is: `Command: /opt/software/Bowtie2/2.4.4-GCC-11.2.0/bin/bowtie2-align-s --wrapper basic-0 -x /data/putnamlab/jillashey/Astrangia2021/mRNA/output/bowtie/refs/Apoc_ref.btindex -S -1 -2 trimmed.AST-1065_R2_001.fastq.gz -b --align-paired-reads trimmed.AST-1065_R1_001.fastq.gz`. For some reason, it is moving the R1 to the end of the script. I'm going to space out using `\` linebreaks the script to see if that fixes it

```
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=15
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/mRNA/scripts              
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load Bowtie2/2.4.4-GCC-11.2.0 

echo "Aligning reads" $(date)

cd /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/

array=($(ls *R1_001.fastq.gz)) # call the clean sequences - make an array to align

for i in ${array[@]}; do
    bowtie2 -x /data/putnamlab/jillashey/Astrangia2021/mRNA/output/bowtie/refs/Apoc_ref.btindex \
	 -1 ${i} \
    -2 $(echo ${i}|sed s/_R1/_R2/) \
    -S \
    -b \
    --align-paired-reads
done

echo "Reads aligned!" $(date)
```


Submitted batch job 292083. SAME ERROR. Going to comment everything except for the -x, -1, and -2 out. Submitted batch job 292084. Appears to be running but the output seems like its going into the slurm output file...Canceled the job. I'm going to add the -b back in a specify a file name (`-b align.${i}`). This should output bam files with the specific file names. Submitted batch job 292085. Still seems to be outputting to the slurm output file...I'm going to stop bowtie and comment out the b and include the -S back in. Submitted batch job 292098

The aligned files are getting put in sam files but labeled as .gz. Will need to convert to sam files, then to bam files once alignment is finished. The output files are also labeled R1 but contain info for both R1 and R2 of that sample. 

### 20240105

Convert .gz files to proper .sam file name and convert to bam file. Test on one sample first in a test script. In the scripts folder: `nano test_sam_to_bam.sh`

```
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=15
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/mRNA/scripts              
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

cd /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim

#mv align.trimmed.AST-1065_R1_001.fastq.gz align.trimmed.AST-1065.sam

samtools sort -o align.trimmed.AST-1065.bam align.trimmed.AST-1065.sam
```

Submitted batch job 292161. Worked! Now I need to make do it for all the samples 

### 20240107 

First, I'm going to move the aligned reads to the bowtie output alignment folder: `mv /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/align* /data/putnamlab/jillashey/Astrangia2021/mRNA/output/bowtie/align`

In the scripts folder: `nano sam_to_bam.sh`

```
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=15
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/mRNA/scripts              
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

echo "Converting sam to bam" $(date)

cd /data/putnamlab/jillashey/Astrangia2021/mRNA/output/bowtie/align

array=($(ls *R1_001.fastq.gz)) # call the clean sequences - make an array to align

for i in ${array[@]}; do
		mv ${i} ${i}.sam
		samtools sort -o ${i}.bam ${i}.sam
		rm ${i}.sam
done

echo "Files converted to bam, sam files removed" $(date)
```

Submitted batch job 292223. Job has been pending for a while...cancelled the job. Adding `-@ 8` to the script to indicate to use 8 threads while running. Submitted batch job 292231

### 20240107 

Cancelled job 292231, as it was still pending. Editing the sbatch info for the script so that the time is 24 hrs, n tasks per node is 1 and mem is 100 GB. Submitted batch job 292244. Got an error saying it didn't recognize the `-@` argument. Removing that and resubmitting. Submitted batch job 292245. Now getting the error: `ls: cannot access *R1_001.fastq.gz: No such file or directory`. The aligned files are gone from the bowtie align folder...fml. Somehow, they must have gotten deleted from the sam to bam script? That's super annoying. 

I guess I need to rerun bowtie. Submitted batch job 292246

### 20240114 

Going to run a test sam to bam so that I can confirm that the files are not being deleted. 

In the scripts folder, modify the `test_sam_to_bam.sh` script: 

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
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

echo "Converting sam to bam for test sample" $(date)

cd /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim

mv align.test.AST-1065_R1_001.fastq.gz align.test.AST-1065.sam

samtools sort -o align.test.AST-1065.bam align.test.AST-1065.sam

echo "Sam to bam conversion complete for test sample" $(date)
```

Submitted batch job 292494. Took ~45 mins. Check % reads aligned

```
interactive 
module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools
samtools flagstat align.test.AST-1065.bam

57029408 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
21403407 + 0 mapped (37.53% : N/A)
57029408 + 0 paired in sequencing
28514704 + 0 read1
28514704 + 0 read2
15442422 + 0 properly paired (27.08% : N/A)
17601996 + 0 with itself and mate mapped
3801411 + 0 singletons (6.67% : N/A)
135006 + 0 with mate mapped to a different chr
64974 + 0 with mate mapped to a different chr (mapQ>=5)
```

37.53% mapped using bowtie2, whereas alignment with hisat2 for this sample was 47.17%. Even though the alignment with hisat2 was higher, I am more inclined to go with the bowtie2 alignment, as I want to be consistent in mapping algorithms for the mRNA and smRNA analysis. Other papers that have examined miRNAs in corals have also used bowtie for the mRNA analysis (i.e., Gajigan & Conaco 2017; Baumgarten et al. 2013, etc). 

Now let's do the sam to bam for all the samples. First, I'm going to move the aligned reads to the bowtie output alignment folder: `mv /data/putnamlab/jillashey/Astrangia2021/mRNA/data/trim/align* /data/putnamlab/jillashey/Astrangia2021/mRNA/output/bowtie/align`

In the scripts folder: `nano sam_to_bam.sh`

```
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=15
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/mRNA/scripts              
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

echo "Converting sam to bam" $(date)

cd /data/putnamlab/jillashey/Astrangia2021/mRNA/output/bowtie/align

array=($(ls *R1_001.fastq.gz)) # call the clean sequences - make an array to align

for i in ${array[@]}; do
		mv ${i} ${i}.sam
		samtools sort -o ${i}.bam ${i}.sam
		#rm ${i}.sam
done

echo "Files converted to bam" $(date)
```

Submitted batch job 292496. Took about 9.5 hours. Calculate % mapped from the bam files. In scripts, `nano map_percent.sh`

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
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

echo "Calculating the percentage of reads aligned" $(date)

cd /data/putnamlab/jillashey/Astrangia2021/mRNA/output/bowtie/align

for i in *.bam; do
    echo "${i}" >> mapped_reads_counts_Apoc
    samtools flagstat ${i} | grep "mapped (" >> mapped_reads_counts_Apoc
done

echo "Calculation complete" $(date)
```

Submitted batch job 292513



### 20240115

Assemble using stringtie. In scripts, `nano assemble.sh`

```
#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=128GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/mRNA/scripts              
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load StringTie/2.2.1-GCC-11.2.0

echo "Assembling transcripts using stringtie" $(date)

cd /data/putnamlab/jillashey/Astrangia2021/mRNA/output/bowtie/align

array1=($(ls *.bam)) #Make an array of sequences to assemble

for i in ${array1[@]}; do
    stringtie -p 8 -e -B -G /data/putnamlab/jillashey/Astrangia_Genome/apoculata_v2.0.gff3 -A ${i}.gene_abund.tab -o ${i}.gtf ${i}
done

echo "Assembly for each sample complete " $(date)
```









For lncRNA 

- CPC: https://github.com/biocoder/CPC2/blob/master/bin/CPC2.py 
- lncRNA discovery overview: https://github.com/urol-e5/deep-dive/blob/main/D-Apul/code/05.32-lncRNA-discovery-overview.Rmd 














NEED TO ASK HP: WAS THE SEQUENCING DONE RF, FF OR FR?