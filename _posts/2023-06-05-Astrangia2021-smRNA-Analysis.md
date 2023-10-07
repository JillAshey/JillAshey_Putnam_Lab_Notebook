---
layout: post
title: Astrangia 2021 small RNA analysis 
date: '2023-06-05'
categories: Analysis
tags: [Bioinformatics, Astrangia, small RNA]
projects: Astrangia 2021
---

## Astrangia 2021 small RNA analysis

These data came from my Astrangia 2021 experiment, during which adult Astrangia colonies were exposed to ambient and high temperatures for ~9 months. 

Files were downloaded to this location: `/data/putnamlab/KITT/hputnam/20230605_Astrangia_smallRNA`

### Set-up
#### Make new folders for this project in my directory 

```
cd /data/putnamlab/jillashey
mkdir Astrangia2021
cd Astrangia2021
mkdir smRNA
cd smRNA
mkdir data scripts fastqc
cd data
mkdir raw
cd raw
```

#### Copy the new files into raw data folder 

```
cp /data/putnamlab/KITT/hputnam/20230605_Astrangia_smallRNA/* .
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
zgrep -c "@GWNJ" *fastq.gz
AST-1065_R1_001.fastq.gz:18782226
AST-1065_R2_001.fastq.gz:18782226
AST-1105_R1_001.fastq.gz:18535712
AST-1105_R2_001.fastq.gz:18535712
AST-1147_R1_001.fastq.gz:43815757
AST-1147_R2_001.fastq.gz:43815757
AST-1412_R1_001.fastq.gz:17729353
AST-1412_R2_001.fastq.gz:17729353
AST-1560_R1_001.fastq.gz:19958419
AST-1560_R2_001.fastq.gz:19958419
AST-1567_R1_001.fastq.gz:18414936
AST-1567_R2_001.fastq.gz:18414936
AST-1617_R1_001.fastq.gz:17164109
AST-1617_R2_001.fastq.gz:17164109
AST-1722_R1_001.fastq.gz:17993435
AST-1722_R2_001.fastq.gz:17993435
AST-2000_R1_001.fastq.gz:18885883
AST-2000_R2_001.fastq.gz:18885883
AST-2007_R1_001.fastq.gz:17958643
AST-2007_R2_001.fastq.gz:17958643
AST-2302_R1_001.fastq.gz:17901570
AST-2302_R2_001.fastq.gz:17901570
AST-2360_R1_001.fastq.gz:17996561
AST-2360_R2_001.fastq.gz:17996561
AST-2398_R1_001.fastq.gz:18231685
AST-2398_R2_001.fastq.gz:18231685
AST-2404_R1_001.fastq.gz:17661430
AST-2404_R2_001.fastq.gz:17661430
AST-2412_R1_001.fastq.gz:18215455
AST-2412_R2_001.fastq.gz:18215455
AST-2512_R1_001.fastq.gz:17643371
AST-2512_R2_001.fastq.gz:17643371
AST-2523_R1_001.fastq.gz:17901421
AST-2523_R2_001.fastq.gz:17901421
AST-2563_R1_001.fastq.gz:18067665
AST-2563_R2_001.fastq.gz:18067665
AST-2729_R1_001.fastq.gz:18840062
AST-2729_R2_001.fastq.gz:18840062
AST-2755_R1_001.fastq.gz:18122482
AST-2755_R2_001.fastq.gz:18122482
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
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/smRNA/scripts              
#SBATCH --error="fastqc_raw_error" #if your job fails, the error report will be put in this file
#SBATCH --output="fastqc_raw_output" #once your job is completed, any final job report comments will be put in this file

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

for file in /data/putnamlab/jillashey/Astrangia2021/smRNA/data/raw/*fastq.gz
do 
fastqc $file --outdir /data/putnamlab/jillashey/Astrangia2021/smRNA/fastqc/raw
done

multiqc --interactive fastqc_results
```

Submitted batch job 261246

### Trim data 

I have not done trimming specific to small RNAs, but this [paper](https://www.sciencedirect.com/science/article/pii/S266591312030131X#fig1) gave a nice [workflow](https://www.dropbox.com/s/2t2gqvkqsr58mjv/PotlaP_miRNA_pipeline.zip?dl=0&file_subpath=%2Fpratibha-miRNApipeline%2Ffinal-pipeline.sh) for miRNA analysis. They suggested using cutadapt. I'm going to follow their code, which applies a minimum length of 18 and a max length of 30. It doesn't do any trimming of adapters, but we will see how the reads look after they go through this cutting. 

In scripts folder: `nano cutadapt.sh`

```
#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=15
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/smRNA/scripts              
#SBATCH --error="cutadapt_error" #if your job fails, the error report will be put in this file
#SBATCH --output="cutadapt_output" #once your job is completed, any final job report comments will be put in this file

module load cutadapt/3.5-GCCcore-11.2.0 

# Make array of sequences to cut 
cd /data/putnamlab/jillashey/Astrangia2021/smRNA/data/raw
array1=($(ls *R1_001.fastq.gz))

echo "Trimming reads so that min length is 18 bp and max length is 30 bp" $(date)

# cutadapt loop
for i in ${array1[@]}; do
	cutadapt --minimum-length=18 --maximum-length=30 -o /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/trimmed.${i} -p /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/trimmed.$(echo ${i}|sed s/_R1/_R2/) ${i} $(echo ${i}|sed s/_R1/_R2/)
done

echo "Trimming done!" $(date)
```

Submitted batch job 269180. Okay cutadapt not working. It is just cutting all of the reads because they are all long and the resulting trimmed file is just empty. Cancelling this job.

20230630
I should also ask Sam White/Javi about their trimming of miRNA data...

Questions for Javi/Sam
- details on seq for E5
- what trimming software did they use? 
- did they do something like setting the min to 18 and max to 30

From Hao et al. 2021: "Raw reads obtained from the sequencing machine were filtered to get clean tags according to the following rules: removing low quality reads containing more than one low quality (Q-value≤ 20) base or containing unknown nucleotides(N) to get the high-quality reads. Then, high-quality reads were filtered by removing reads without 3′ adapters, containing 5′ adapters, containing 3′ and 5′ adapters but no small RNA fragment between them, containing polyA in small RNA fragment and shorter than 18 nt to get clean tags. The clean tags were aligned with small RNAs in the GenBank database". This sounds like something i should try, but I'm not sure how to ID the 3' and 5' adapters in my sequences. 

trimming - try cutadapt, trimmomatic or trimgalore 

Looking at the adapter content MultiQC plot, it looks like the reads were processed using the illumina universal adapter and the illumina small rna 5' adapter. The R1 reads have the universal adapter and the R2 reads have the small rna 5' adapter. Not sure why that is. I looked the adapters sequences on Illumina and found this [post](https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/FastQC_Adapter_Kmer_files_fDG.htm) that says the Illumina Universal Adapter—AGATCGGAAGAG and Illumina Small RNA 5' Adapter—GATCGTCGGACT. I am unsure when this page wzs written, but I'm going to test them out. 

`nano cudadapt.sh`


```
#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=15
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/smRNA/scripts              
#SBATCH --error="cutadapt_error" #if your job fails, the error report will be put in this file
#SBATCH --output="cutadapt_output" #once your job is completed, any final job report comments will be put in this file

module load cutadapt/3.5-GCCcore-11.2.0 

# Make array of sequences to cut 
cd /data/putnamlab/jillashey/Astrangia2021/smRNA/data/raw
array1=($(ls *R1_001.fastq.gz))

echo "Trimming reads so that min length is 18 bp and max length is 30 bp" $(date)

echo "Starting to trim using the Illumina universal adapter" $(date)

for i in ${array1[@]}; do
	cutadapt -a AGATCGGAAGAG --minimum-length=18 --maximum-length=30 -o /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/trimmed.${i} -p /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/trimmed.$(echo ${i}|sed s/_R1/_R2/) ${i} $(echo ${i}|sed s/_R1/_R2/)
done

echo "Starting to trim using the Illumina small rna 5' adapter" $(date)

for i in ${array1[@]}; do
	cutadapt -a GATCGTCGGACT --minimum-length=18 --maximum-length=30 -o /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/trimmed.again.${i} -p /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/trimmed.again.$(echo ${i}|sed s/_R1/_R2/) ${i} $(echo ${i}|sed s/_R1/_R2/)
done

echo "Trimming done!" $(date)
```

Submitted batch job 275252. The trimming parameters seem to be too stringent. It is saying that the reads are either too long or too short. Not sure what this means. I'm going to try Sam's [code](https://robertslab.github.io/sams-notebook/2023/06/20/Trimming-and-QC-E5-Coral-sRNA-seq-Data-fro-A.pulchra-P.evermanni-and-P.meandrina-Using-FastQC-flexbar-and-MultiQC-on-Mox.html) where he trimmed some smRNAs using a software called flexbar. Need to ask Kevin Bryan to install flexbar. 

Flexbar installed! 

Let's try [Flexbar trimming code](https://robertslab.github.io/sams-notebook/posts/2023/2023-06-20-Trimming-and-QC---E5-Coral-sRNA-seq-Data-fro-A.pulchra-P.evermanni-and-P.meandrina-Using--FastQC-flexbar-and-MultiQC-on-Mox/) that Sam White wrote for the e5 small RNA analysis

First, I'm only going to try to trim one sample (2 reads) to see if flexbar works. 

```
cd /data/putnamlab/jillashey/Astrangia2021/smRNA/data/raw
```

First, make the NEB adapters fasta file.

```
nano NEB-adapters.fasta

>first
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>second
GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT
```

In scripts folder: `nano test_flexbar.sh`

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH --error="flexbar_raw_error" #if your job fails, the error report will be put in this file
#SBATCH --output="flexbar_raw_output" #once your job is completed, any final job report comments will be put in this file

module load Flexbar/3.5.0-foss-2018b  

cd /data/putnamlab/jillashey/Astrangia2021/smRNA/data/raw

R1_fastq=/data/putnamlab/jillashey/Astrangia2021/smRNA/data/raw/AST-1065_R1_001.fastq.gz 
R2_fastq=/data/putnamlab/jillashey/Astrangia2021/smRNA/data/raw/AST-1065_R2_001.fastq.gz 

flexbar \
-r ${R1_fastq} \
-p ${R2_fastq}  \
-a NEB-adapters.fasta \
-ap ON \
-qf i1.8 \
-qt 25 \
--post-trim-length 35 \
--target TEST_AST-1065 \
--zip-output GZ
```

Submitted batch job 284050. Finished in about 25 mins. 

Now I'm going to run fastqc on the test sample and see how it looks: 

```
module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

fastqc TEST_AST-1065_1.fastq.gz TEST_AST-1065_2.fastq.gz
multiqc *fastqc*```

ADD PLOTS 

The plots look good! I'm going to move forward w/ flexbar trimming. 

In scripts folder: `nano flexbar.sh`

```
#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=200GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH --error="flexbar_error" #if your job fails, the error report will be put in this file
#SBATCH --output="flexbar_output" #once your job is completed, any final job report comments will be put in this file

module load Flexbar/3.5.0-foss-2018b  
module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

cd /data/putnamlab/jillashey/Astrangia2021/smRNA/data/raw

echo "Trimming reads using flexbar" $(date)

array1=($(ls *R1_001.fastq.gz))

for i in ${array1[@]}; do
flexbar \
-r ${i} \
-p $(echo ${i}|sed s/_R1/_R2/) \
-a NEB-adapters.fasta \
-ap ON \
-qf i1.8 \
-qt 25 \
--post-trim-length 35 \
--target $(echo ${i}|sed s/_R1/_R2/) \
--zip-output GZ
done 

mv *trim ../trim

echo "Trimming complete, now running fastqc and multiqc" $(date)

cd /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim

for file in *fastq.gz
do 
fastqc $file 
done

multiqc *fastqc*

echo "QC complete!" $(date)
```

Submitted batch job 284064


