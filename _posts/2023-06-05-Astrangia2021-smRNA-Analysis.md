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
multiqc *fastqc*
```

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
```

Submitted batch job 284064

Move trimmed reads to trim folder

```
cd /data/putnamlab/jillashey/Astrangia2021/smRNA/data/raw
mv trim.AST-* ../trim/
```

### Trim QC

#### Run fastqc to quality check trim reads 

In scripts folder: `nano fastqc_trim.sh`

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
#SBATCH --error="fastqc_trim_error" #if your job fails, the error report will be put in this file
#SBATCH --output="fastqc_trim_output" #once your job is completed, any final job report comments will be put in this file

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

for file in /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/*fastq.gz
do 
fastqc $file --outdir /data/putnamlab/jillashey/Astrangia2021/smRNA/fastqc/trim
done

multiqc --interactive fastqc_results/trim
```

Submitted batch job 284426

#### Count number of reads per file 


```
zgrep -c "@GWNJ" *fastq.gz

trim.AST-1065_R2_001.fastq.gz_1.fastq.gz:17829111
trim.AST-1065_R2_001.fastq.gz_2.fastq.gz:17829111
trim.AST-1105_R2_001.fastq.gz_1.fastq.gz:17238126
trim.AST-1105_R2_001.fastq.gz_2.fastq.gz:17238126
trim.AST-1147_R2_001.fastq.gz_1.fastq.gz:40415224
trim.AST-1147_R2_001.fastq.gz_2.fastq.gz:40415224
trim.AST-1412_R2_001.fastq.gz_1.fastq.gz:16279555
trim.AST-1412_R2_001.fastq.gz_2.fastq.gz:16279555
trim.AST-1560_R2_001.fastq.gz_1.fastq.gz:17827024
trim.AST-1560_R2_001.fastq.gz_2.fastq.gz:17827024
trim.AST-1567_R2_001.fastq.gz_1.fastq.gz:16611397
trim.AST-1567_R2_001.fastq.gz_2.fastq.gz:16611397
trim.AST-1617_R2_001.fastq.gz_1.fastq.gz:16077717
trim.AST-1617_R2_001.fastq.gz_2.fastq.gz:16077717
trim.AST-1722_R2_001.fastq.gz_1.fastq.gz:16430221
trim.AST-1722_R2_001.fastq.gz_2.fastq.gz:16430221
trim.AST-2000_R2_001.fastq.gz_1.fastq.gz:17428854
trim.AST-2000_R2_001.fastq.gz_2.fastq.gz:17428854
trim.AST-2007_R2_001.fastq.gz_1.fastq.gz:16559551
trim.AST-2007_R2_001.fastq.gz_2.fastq.gz:16559551
trim.AST-2302_R2_001.fastq.gz_1.fastq.gz:16665370
trim.AST-2302_R2_001.fastq.gz_2.fastq.gz:16665370
trim.AST-2360_R2_001.fastq.gz_1.fastq.gz:16648356
trim.AST-2360_R2_001.fastq.gz_2.fastq.gz:16648356
trim.AST-2398_R2_001.fastq.gz_1.fastq.gz:16788208
trim.AST-2398_R2_001.fastq.gz_2.fastq.gz:16788208
trim.AST-2404_R2_001.fastq.gz_1.fastq.gz:16712903
trim.AST-2404_R2_001.fastq.gz_2.fastq.gz:16712903
trim.AST-2412_R2_001.fastq.gz_1.fastq.gz:17488508
trim.AST-2412_R2_001.fastq.gz_2.fastq.gz:17488508
trim.AST-2512_R2_001.fastq.gz_1.fastq.gz:16265716
trim.AST-2512_R2_001.fastq.gz_2.fastq.gz:16265716
trim.AST-2523_R2_001.fastq.gz_1.fastq.gz:16995265
trim.AST-2523_R2_001.fastq.gz_2.fastq.gz:16995265
trim.AST-2563_R2_001.fastq.gz_1.fastq.gz:17023002
trim.AST-2563_R2_001.fastq.gz_2.fastq.gz:17023002
trim.AST-2729_R2_001.fastq.gz_1.fastq.gz:17121869
trim.AST-2729_R2_001.fastq.gz_2.fastq.gz:17121869
trim.AST-2755_R2_001.fastq.gz_1.fastq.gz:16269823
trim.AST-2755_R2_001.fastq.gz_2.fastq.gz:16269823
```



20231130

Locations of cnidarian miRNA data 
- [Stylophora pistillata](https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR1130519&display=metadata)
- [Hydra](https://mirbase.org/browse/results/?organism=hma)
- [Nematostella](https://mirbase.org/browse/results/?organism=nve)
- [Acropora muricata, Montipora capricornis, Montipora foliosa, Pocillopora verrucosa](http://118.89.77.43:7081/miR/browse.php)

20240103
Should I retrim with fastp to keep it consistent with mRNA? Going to try it out and compare results from flexbar. In  trimmed smRNA folder, make new directory to put fastp trimmed reads

```
cd /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim
mkdir fastp
```

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
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/smRNA/scripts              
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.out

# Load modules needed 
module load fastp/0.19.7-foss-2018b
module load FastQC/0.11.8-Java-1.8
module load MultiQC/1.9-intel-2020a-Python-3.8.2

# Make an array of sequences to trim in raw data directory 

cd /data/putnamlab/jillashey/Astrangia2021/smRNA/data/raw
array1=($(ls *R1_001.fastq.gz))

echo "Read trimming of adapters started." $(date)

# fastp and fastqc loop 
for i in ${array1[@]}; do
    fastp --in1 ${i} \
        --in2 $(echo ${i}|sed s/_R1/_R2/)\
        --out1 /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/fastp/trimmed.${i} \
        --out2 /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/fastp/trimmed.$(echo ${i}|sed s/_R1/_R2/) \
        --detect_adapter_for_pe \
        --qualified_quality_phred 30 \
        --unqualified_percent_limit 10 \
        #--length_required 100 \
        --cut_right cut_right_window_size 5 cut_right_mean_quality 20
	 fastqc /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/fastp/trimmed.${i}
    fastqc /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/fastp/trimmed.$(echo ${i}|sed s/_R1/_R2/)
done

echo "Read trimming of adapters complete." $(date)

# Quality Assessment of Trimmed Reads
cd /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/fastp #go to output directory

# Compile MultiQC report from FastQC files 
multiqc --interactive ./  

echo "Cleaned MultiQC report generated." $(date)
```

Submitted batch job 291969. Took about 5 hours. The QC plots don't look amazing and the length is still at 100 bp. I'm going to rerun but adding the argument `--length_limit 30` for fastp. This means that reads longer than 30 bp will be discarded. Submitted batch job 291997. Downloaded the QC report and it looks okay. Low amount of reads but high duplication. Let's see how many counts are in each file: 

```
zgrep -c "@GWNJ" *fastq.gz

trimmed.AST-1065_R1_001.fastq.gz:7834561
trimmed.AST-1065_R2_001.fastq.gz:7834561
trimmed.AST-1105_R1_001.fastq.gz:8754651
trimmed.AST-1105_R2_001.fastq.gz:8754651
trimmed.AST-1147_R1_001.fastq.gz:22326820
trimmed.AST-1147_R2_001.fastq.gz:22326820
trimmed.AST-1412_R1_001.fastq.gz:8158050
trimmed.AST-1412_R2_001.fastq.gz:8158050
trimmed.AST-1560_R1_001.fastq.gz:8733402
trimmed.AST-1560_R2_001.fastq.gz:8733402
trimmed.AST-1567_R1_001.fastq.gz:9830273
trimmed.AST-1567_R2_001.fastq.gz:9830273
trimmed.AST-1617_R1_001.fastq.gz:8146294
trimmed.AST-1617_R2_001.fastq.gz:8146294
trimmed.AST-1722_R1_001.fastq.gz:9014021
trimmed.AST-1722_R2_001.fastq.gz:9014021
trimmed.AST-2000_R1_001.fastq.gz:10252309
trimmed.AST-2000_R2_001.fastq.gz:10252309
trimmed.AST-2007_R1_001.fastq.gz:9622779
trimmed.AST-2007_R2_001.fastq.gz:9622779
trimmed.AST-2302_R1_001.fastq.gz:8921101
trimmed.AST-2302_R2_001.fastq.gz:8921101
trimmed.AST-2360_R1_001.fastq.gz:8635502
trimmed.AST-2360_R2_001.fastq.gz:8635502
trimmed.AST-2398_R1_001.fastq.gz:9565008
trimmed.AST-2398_R2_001.fastq.gz:9565008
trimmed.AST-2404_R1_001.fastq.gz:8798584
trimmed.AST-2404_R2_001.fastq.gz:8798584
trimmed.AST-2412_R1_001.fastq.gz:7032530
trimmed.AST-2412_R2_001.fastq.gz:7032530
trimmed.AST-2512_R1_001.fastq.gz:7463762
trimmed.AST-2512_R2_001.fastq.gz:7463762
trimmed.AST-2523_R1_001.fastq.gz:8272150
trimmed.AST-2523_R2_001.fastq.gz:8272150
trimmed.AST-2563_R1_001.fastq.gz:8537237
trimmed.AST-2563_R2_001.fastq.gz:8537237
trimmed.AST-2729_R1_001.fastq.gz:7929477
trimmed.AST-2729_R2_001.fastq.gz:7929477
trimmed.AST-2755_R1_001.fastq.gz:9775688
trimmed.AST-2755_R2_001.fastq.gz:9775688
```

Flexbar keeps more reads but it seems like it is combining the two reads into one for each sample (which I don't want to do yet). Let's edit the flexbar code and see if I can fix that. For now, I just commented out the line `--target trim.$(echo ${i}|sed s/_R1/_R2/) \`, which names the files. Also changing max read length to 30. Submitted batch job 292029

20240104
Flexbar trimming finished last night but it just rewrote the files over one another so only the last sample has the files. Looking at the flexbar [documentation](https://github.com/seqan/flexbar/wiki/Options#output-selection), it seems like `--target` is the prefix for output file names or paths, where as `--output-reads` and `--output-reads2` is used for the output file for reads 1 and 2 instead of the target prefix usage. So I am just dumb and should've specified the output reads argument instead of target. 




Let's look at the flexbar script again: 

```
#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=200GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load Flexbar/3.5.0-foss-2018b  

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
-qt 30 \
-t /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar
--post-trim-length 30 \
--output-reads trim.${i} \
--output-reads2 trim.$(echo ${i}|sed s/_R1/_R2/) \
--zip-output GZ
done 
```

Made those changes for output reads to the script. Also changed script so that `-qt` was 30 (ie phred score of 30) and post trim length was 30 bp. Submitted batch job 292075. Kept failing, for some reason all my files were empty. I'm now re-copying them from the KITT directory (thank god for backups) and then will rerun flexbar. Added trim. prefix to the beginning of the output file names so that the original files don't get overwritten. Submitted batch job 292079. 

Checked back and it is still overwriting the files with the flexbar.log and flexbar fastq files...Also in the slurm error report, it tells me that the post trim length argument was not found. Looking at the code again, I didn't add a `\` after the -t line. Going to add it and rerun. Also going to try to take out the `--` and just do the `-`. Submitted batch job 292092

Job finished but the files are empty. Why does Flexbar hate me? The error file said it couldn't open the files that flexbar created but why did those need to be opened anyway? 

I'm going to go back to the original script that I ran because it seems to have worked, even though it named both files R2. 

In the scripts folder: `nano flexbar_og.sh`

```
#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=200GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load Flexbar/3.5.0-foss-2018b  

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
-qt 30 \
--post-trim-length 30 \
--target ${i} \
--zip-output GZ
done 
```

Submitted batch job 292095

While that is running, I'm going to figure out how to set up a conda environment using the conda unity [documentation](https://docs.unity.uri.edu/documentation/software/conda/). For whatever reason, Kevin Bryan recommended I use a conda environment. 

I'm going to make it in the putnam lab folder. 

First I need to load miniconda module: `module load Miniconda3/4.9.2`

Now I need to create a conda environment. 

```
conda create --prefix /data/putnamlab/miranda
```

If I need to update:

```
==> WARNING: A newer version of conda exists. <==
  current version: 4.9.2
  latest version: 23.11.0

Please update conda by running

    $ conda update -n base -c defaults conda
```

Now that the conda environment is created, I need to activate it.

```
conda activate /data/putnamlab/miranda
```

It told me my shell had not been properly configured. To do that: 

```
conda init

no change     /opt/software/Miniconda3/4.9.2/condabin/conda
no change     /opt/software/Miniconda3/4.9.2/bin/conda
no change     /opt/software/Miniconda3/4.9.2/bin/conda-env
no change     /opt/software/Miniconda3/4.9.2/bin/activate
no change     /opt/software/Miniconda3/4.9.2/bin/deactivate
no change     /opt/software/Miniconda3/4.9.2/etc/profile.d/conda.sh
no change     /opt/software/Miniconda3/4.9.2/etc/fish/conf.d/conda.fish
no change     /opt/software/Miniconda3/4.9.2/shell/condabin/Conda.psm1
no change     /opt/software/Miniconda3/4.9.2/shell/condabin/conda-hook.ps1
no change     /opt/software/Miniconda3/4.9.2/lib/python3.8/site-packages/xontrib/conda.xsh
no change     /opt/software/Miniconda3/4.9.2/etc/profile.d/conda.csh
modified      /home/jillashey/.bashrc

==> For changes to take effect, close and re-open your current shell. <==
```

I closed the terminal window and logged back in on a new window. Now let's activate!

```
conda activate /data/putnamlab/miranda
```

Now the conda environment is activated! My shell thing now looks like this: `(/data/putnamlab/miranda) [jillashey@ssh3 putnamlab]$`. To deactivate, do `conda deactivate`. 

Create a conda environment for mirdeep2 using the same steps above. 

I'm going to first work in the mirdeep2 environment. Activate the environment: `conda activate /data/putnamlab/mirdeep2`

Install mirdeep2 within the conda env: `conda install bioconda::mirdeep2`. This will take a few minutes to install and load the required packages. 

I'm going to try to run [mirDeep2](https://github.com/rajewsky-lab/mirdeep2) using code from the mirdeep2 github [tutorial](https://github.com/rajewsky-lab/mirdeep2/blob/master/TUTORIAL.md) and Sam White's [code](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/code/11-Peve-sRNAseq-miRdeep2.md#7-run-mirdeep2) from the E5 deep dive project. 

Before running any mirdeep2 modules, I need to upload some databases to the HPC and configure some of my files. Let's first configure my files. There are 2 files (R1 and R2) per sample, so I need to concatenate and collapse the reads. 

```
cat trimmed.AST-1065_R1_001.fastq.gz trimmed.AST-1065_R2_001.fastq.gz > cat.trimmed.AST-1065.fastq

module load FASTX-Toolkit/0.0.14-GCC-9.3.0 

# gunzip cat.trimmed.AST-1065.fastq.gz # files must be unzipped for collapsing; unzip if needed

fastx_collapser -v -i cat.trimmed.AST-1065.fastq -o collapse.cat.trimmed.AST-1065.fastq

head collapse.cat.trimmed.AST-1065.fastq 
>1-116635
TGGTCTATGGTGTAACTGGCAACACGTCTGT
>2-115039
ACAGACGTGTTGCCAGTTACACCATAGACCA
>3-104350
TGGTCTATGGTGTAACTGGCAACACGTCTGTT
>4-103158
AACAGACGTGTTGCCAGTTACACCATAGACCA
>5-71882
TGAAAATCTTTTCTCTGAAGTGGAA
```

As per the mirdeep2 documentation, The readID must end with _xNumber and is not allowed to contain whitespaces. So it has to have the format name_uniqueNumber_xnumber. 

```
sed '/^>/ s/-/_x/g' collapse.cat.trimmed.AST-1065.fastq \
| sed '/^>/ s/>/>seq_/' \
> collapse.cat.trimmed.AST-1065.fastq 

>seq_1_x116635
TGGTCTATGGTGTAACTGGCAACACGTCTGT
>seq_2_x115039
ACAGACGTGTTGCCAGTTACACCATAGACCA
>seq_3_x104350
TGGTCTATGGTGTAACTGGCAACACGTCTGTT
>seq_4_x103158
AACAGACGTGTTGCCAGTTACACCATAGACCA
>seq_5_x71882
TGAAAATCTTTTCTCTGAAGTGGAA
```

Next I need to reformat the genome fasta description lines. miRDeep2 can’t process genome FastAs with spaces in the description lines. I don't think the Apoc genome has any spaces but I'm going to double check. 

```
cd /data/putnamlab/jillashey/Astrangia_Genome/

grep "^>" apoculata.assembly.scaffolds_chromosome_level.fasta 
>chromosome_1
>chromosome_2
>chromosome_3
>chromosome_4
>chromosome_5
>chromosome_6
>chromosome_7
>chromosome_8
>chromosome_9
>chromosome_10
>chromosome_11
>chromosome_12
>chromosome_13
>chromosome_14
```

Nice, there are no spaces so I don't need to reformat. If I did, I would sub the spaces with underscores. 

Index the genome with bowtie (NOT bowtie2). In the scripts folder: `nano bowtie_build.sh`

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=15
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/smRNA/scripts              
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.out

module load GCCcore/11.3.0 #I needed to add this to resolve conflicts between loaded GCCcore/9.3.0 and GCCcore/11.3.0
module load Bowtie/1.3.1-GCC-11.3.0

# Index the reference genome for A. poculata 
bowtie-build /data/putnamlab/jillashey/Astrangia_Genome/apoculata.assembly.scaffolds_chromosome_level.fasta Apoc_ref.btindex

echo "Referece genome indexed!" $(date)
```

The indexed genome lives in the scripts folder for now. 

Load the miRbase mature miRNA fasta database onto the server. I downloaded it onto my computer (its not very large) and will now copy it to the server. I downloaded it on 1/3/24. It will live in the refs folder in the smRNA directory. 

```
cd /data/putnamlab/jillashey/Astrangia2021/smRNA
mkdir refs 
cd refs 
ls refs 
20240103_mature.fa

head 20240103_mature.fa 
>cel-let-7-5p MIMAT0000001 Caenorhabditis elegans let-7-5p
UGAGGUAGUAGGUUGUAUAGUU
>cel-let-7-3p MIMAT0015091 Caenorhabditis elegans let-7-3p
CUAUGCAAUUUUCUACCUUACC
>cel-lin-4-5p MIMAT0000002 Caenorhabditis elegans lin-4-5p
UCCCUGAGACCUCAAGUGUGA
>cel-lin-4-3p MIMAT0015092 Caenorhabditis elegans lin-4-3p
ACACCUGGGCUCUCCGGGUACC
>cel-miR-1-5p MIMAT0020301 Caenorhabditis elegans miR-1-5p
CAUACUUCCUUACAUGCCCAUA 
```

Check how many mature miRNA sequences there are in the file

```
zgrep -c ">" 20240103_mature.fa 
48885
```

I'm going to reformat the fasta header names here so there are no spaces 

```
sed '/^>/ s/ /_/g' 20240103_mature.fa \
| sed '/^>/ s/,//g' \
> 20240103_mature.fa

head 20240103_mature.fa 
>cel-let-7-5p_MIMAT0000001_Caenorhabditis_elegans_let-7-5p
UGAGGUAGUAGGUUGUAUAGUU
>cel-let-7-3p_MIMAT0015091_Caenorhabditis_elegans_let-7-3p
CUAUGCAAUUUUCUACCUUACC
>cel-lin-4-5p_MIMAT0000002_Caenorhabditis_elegans_lin-4-5p
UCCCUGAGACCUCAAGUGUGA
>cel-lin-4-3p_MIMAT0015092_Caenorhabditis_elegans_lin-4-3p
ACACCUGGGCUCUCCGGGUACC
>cel-miR-1-5p_MIMAT0020301_Caenorhabditis_elegans_miR-1-5p
CAUACUUCCUUACAUGCCCAUA
```

Do I need to change the U to T in `20240103_mature.fa`? Let's try mapping first and seeing. 

```
mapper.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/fastp/collapse.cat.trimmed.AST-1065.fastq -e -p /data/putnamlab/jillashey/Astrangia2021/smRNA/scripts/Apoc_ref.btindex -s reads_collapsed.fa -t reads_collapsed_vs_genome.arf -v  
```

Didn't take very long! Only a few seconds. It output this: 

```
discarding short reads
mapping reads to genome index
trimming unmapped nts in the 3' ends
Log file for this run is in mapper_logs and called mapper.log_10292
Mapping statistics

#desc	total	mapped	unmapped	%mapped	%unmapped
total: 4289826	379347	3910479	8.843	91.157
seq: 4289826	379347	3910479	8.843	91.157
```

Not a very high mapping rate but I wonder if that's normal. I may need to change the U to T in the miRbase file. It also produced a `bowtie.log` file:

```
less bowtie.log 

# reads processed: 1768
# reads with at least one reported alignment: 215 (12.16%)
# reads that failed to align: 1496 (84.62%)
# reads with alignments suppressed due to -m: 57 (3.22%)
Reported 613 alignments to 1 output stream(s)
```

The `reads_collapsed_vs_genome.arf` provides info about the sequences that did align: 

```
head reads_collapsed_vs_genome.arf 
seq_19_x21706	32	1	32	aacttttgacggtggatctcttggctcacgca	chromosome_2	32	42321	42352	aacttttgacggtggatctcttggctcacgca	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_19_x21706	32	1	32	aacttttgacggtggatctcttggctcacgca	chromosome_2	32	53070	53101	aacttttgacggtggatctcttggctcacgca	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_19_x21706	32	1	32	aacttttgacggtggatctcttggctcacgca	chromosome_2	32	20734	20765	aacttttgacggtggatctcttggctcacgca	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_19_x21706	32	1	32	aacttttgacggtggatctcttggctcacgca	chromosome_2	32	31478	31509	aacttttgacggtggatctcttggctcacgca	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_20_x21513	32	1	32	tgcgtgagccaagagatccaccgtcaaaagtt	chromosome_2	32	31478	31509	tgcgtgagccaagagatccaccgtcaaaagtt	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_20_x21513	32	1	32	tgcgtgagccaagagatccaccgtcaaaagtt	chromosome_2	32	20734	20765	tgcgtgagccaagagatccaccgtcaaaagtt	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_20_x21513	32	1	32	tgcgtgagccaagagatccaccgtcaaaagtt	chromosome_2	32	42321	42352	tgcgtgagccaagagatccaccgtcaaaagtt	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_20_x21513	32	1	32	tgcgtgagccaagagatccaccgtcaaaagtt	chromosome_2	32	53070	53101	tgcgtgagccaagagatccaccgtcaaaagtt	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_33_x17087	43	1	43	ttgctacgatcttctgagattaagcctttgttctaagatttgt	chromosome_2	43	879093	879135	ttgctacgatcttctgagattaagcctttgttctaagatttgt	+mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_33_x17087	43	1	43	ttgctacgatcttctgagattaagcctttgttctaagatttgt	chromosome_2	43	38360	38402	ttgctacgatcttctgagattaagcctttgttctaagatttgt	-mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
```

Not sure what all of this info means, but I will look into it. I also need to look into why/how the reads get collapsed because the collapse set left me with only 1796 sequences (compared to the 26561348 sequences I had in the cat file) and I want to make sure that's normal. 

Now lets run mirdeep2!!!!!!!!

```
miRDeep2.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/fastp/collapse.cat.trimmed.AST-1065.fastq /data/putnamlab/jillashey/Astrangia_Genome/apoculata.assembly.scaffolds_chromosome_level.fasta /data/putnamlab/jillashey/Astrangia2021/smRNA/reads_collapsed_vs_genome.arf /data/putnamlab/jillashey/Astrangia2021/smRNA/refs/20240103_mature.fa none none -t N.vectensis -P -v -g -1 2>report.log
```

I need to specify none 2x because I do not have the files for known miRNAs or known precursor miRNAs in this species. I got an error:

```
#Starting miRDeep2
/data/putnamlab/mirdeep2/bin/miRDeep2.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/fastp/collapse.cat.trimmed.AST-1065.fastq /data/putnamlab/jillashey/Astrangia_Genome/apoculata.assembly.scaffolds_chromosome_level.fasta /data/putnamlab/jillashey/Astrangia2021/smRNA/reads_collapsed_vs_genome.arf /data/putnamlab/jillashey/Astrangia2021/smRNA/refs/20240103_mature.fa none none -t N.vectensis -P -v -g -1

miRDeep2 started at 16:15:57


mkdir mirdeep_runs/run_04_01_2024_t_16_15_57

#testing input files
started: 16:16:06
sanity_check_mature_ref.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/refs/20240103_mature.fa


ended: 16:16:06
total:0h:0m:0s

sanity_check_reads_ready_file.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/fastp/collapse.cat.trimmed.AST-1065.fastq

started: 16:16:06
ESC[1;31mError: ESC[0mproblem with /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/fastp/collapse.cat.trimmed.AST-1065.fastq
Error in line 860: Either the sequence

GATGGAATTGTAGCAT


contains less than 17 characters or contains characters others than [acgtunACGTUN]

Please make sure that your file only comprises sequences that have at least 17 characters

containing letters [acgtunACGTUN]
```

My collapsed read file has some sequences that are <17 bp which mirdeep2 doesn't like. I need to remove sequences with <17 nts (or do it during the trimming step). Used chatgpt for the code below :)  

```
#!/bin/bash

# Define the input and output files
input_file="collapse.cat.trimmed.AST-1065.fastq"
output_file="17_collapse.cat.trimmed.AST-1065.fastq"

# Initialize the output file
> "$output_file"

# Use awk to process the sequences
awk '{
    if (substr($0, 1, 1) == ">") {
        header = $0
        getline
        sequence = $0
        if (length(sequence) >= 17) {
            print header >> "'$output_file'"
            print sequence >> "'$output_file'"
        }
    }
}' "$input_file"

```

Rerun mirdeep2 with the new file

```
miRDeep2.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/fastp/17_collapse.cat.trimmed.AST-1065.fastq /data/putnamlab/jillashey/Astrangia_Genome/apoculata.assembly.scaffolds_chromosome_level.fasta /data/putnamlab/jillashey/Astrangia2021/smRNA/reads_collapsed_vs_genome.arf /data/putnamlab/jillashey/Astrangia2021/smRNA/refs/20240103_mature.fa none none -t N.vectensis -P -v -g -1 2>report.log
```

Took about 2 mins and there were no alignments. I'm going to change the U to T in the `20240103_mature.fa` file, then rerun the mapping and mirdeep2 steps. 

```
#!/bin/bash

# Define the input and output files
input_file="20240103_mature.fa"  # Replace with your actual input file name
output_file="20240103_mature_T.fa"

# Initialize the output file
> "$output_file"

# Use awk to process the file
awk '{
    if (substr($0, 1, 1) == ">") {
        print $0 >> "'$output_file'"  # Print the identifier as is
    } else {
        gsub(/U/, "T", $0)  # Replace U with T in sequences
        print $0 >> "'$output_file'"
    }
}' "$input_file"

```

I guess I don't need to redo the mapping step because the mapping did not use the 20240103_mature.fa file. Therefore, I will proceed to the mirdeep2 step.

```
miRDeep2.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/fastp/17_collapse.cat.trimmed.AST-1065.fastq /data/putnamlab/jillashey/Astrangia_Genome/apoculata.assembly.scaffolds_chromosome_level.fasta /data/putnamlab/jillashey/Astrangia2021/smRNA/reads_collapsed_vs_genome.arf /data/putnamlab/jillashey/Astrangia2021/smRNA/refs/20240103_mature_T.fa none none -t N.vectensis -P -v -g -1 2>report.log
```

Still no alignment :'( I need to talk to Sam and ask him the following: 

- How many reads were left after he concatnated and collapsed his fastq files? 
- What was his alignment after the mapping step? 

Also what if I didn't collapse the reads? What if I just removed the heading, + sign, and quality scores and formatted it like mirdeep2 wants it? 

I just briefly reran the collapse step and now there are hundreds of thousands of sequences...may need to rerun mapping idk 

### 20240105

Flexbar finished running overnight, took about 9 hours. Now I need to QC it. First going to move it from the raw data folder to the trim data folder. 

In `/data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim`, I moved the old flexbar trimmed seqs to the foler `flexbar_old`. I moved the newly trimmed seqs from the raw data folder into `/data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar`. 

Now let's run fastqc on the newly trimmed samples using the `fastqc_trim.sh` from above, but changing the directories 

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
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.out

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

echo "QC for trimmed reads using flexbar with max length of 30 bp" $(date)

for file in /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar/*fastq.gz
do 
fastqc $file --outdir /data/putnamlab/jillashey/Astrangia2021/smRNA/fastqc/trim/flexbar
done

echo "QC complete, run multiqc" $(date)

multiqc --interactive fastqc_results/trim/flexbar
```

Submitted batch job 292159

If the data looks good, I will do another test run of mirdeep2. It finished in ~30 mins but it did not complete the multiQC. When I went into the folder to run the multiqc step, it appears to only have run it on the R2 files. I'm going to move the R1 and R2 fastqc info into separate folders and run the QC on them separately. There is probably a better way to do this idk. 

```
mkdir R1 R2
mv *_1_fastqc* R1
mv *_2_fastqc* R2
```

Now go into each folder and run multiqc separately. QC looks good for both reads. 

Next I need to cat, collapse and prep reads for mirdeep2. 

### 20240107

I'm going to write a script that will cat and collapse the reads. I'll do it on a test sample first. In the scripts folder: `nano test_cat_collapse.sh`. 

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
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load FASTX-Toolkit/0.0.14-GCC-9.3.0 

cd /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar

echo "Concatenating R1 and R2 for test sample" $(date)

cat AST-1065_R1_001.fastq.gz_1.fastq AST-1065_R1_001.fastq.gz_2.fastq > cat.AST-1065.fastq

echo "Collapsing redundant sequences with fastx collapse" $(date)

fastx_collapser -v -i cat.AST-1065.fastq -o collapse.cat.AST-1065.fastq

echo "Prep sequence IDs for mirdeep2 analysis" $(date)

sed '/^>/ s/-/_x/g' collapse.cat.AST-1065.fastq \
| sed '/^>/ s/>/>seq_/' \
> collapse.cat.AST-1065.fastq

echo "Done!" $(date)
```

Submitted batch job 292222. I did the genome and database prep already so I don't need to redo that. Took about 20 mins, but the collapsed file is empty...going to have the output file for the `sed` portion be `sed.collapse.cat.AST-1065.fastq` to see if the `sed` portion is what is happening to the files. Submitted batch job 292224. That iteration worked! 

Check how many sequences are in the collapsed file. 

```
zgrep -c ">" sed.collapse.cat.AST-1065.fastq 
11979585

head sed.collapse.cat.AST-1065.fastq 
>seq_1_x357414
TGGTCTATGGTGTAACTGGCAACACGTCTG
>seq_2_x138955
ACAGACGTGTTGCCAGTTACACCATAGACC
>seq_3_x125294
AACAGACGTGTTGCCAGTTACACCATAGAC
>seq_4_x98253
TGAAAATCTTTTCTCTGAAGTGGAA
>seq_5_x87633
TTCCACTTCAGAGAAAAGATTTTCA
```

Nice, almost 12 million sequences were retained. I was looking at the fastx collapse documentation and it said that the first number in the sequence id corresponded to a sequence and the second number corresponded to how many times that sequence appeared prior to the file being collapsed. So for example, `>seq_1_x357414` was the most represented sequence, as indicated by the 1, and it appeared 357414 times in the pre-collapsed file. 

Need to make sure that all of my sequences are >17 bp, as mirdeep2 does not run if sequences are present with <16 bp. 

My collapsed read file has some sequences that are <17 bp which mirdeep2 doesn't like. I need to remove sequences with <17 nts (or do it during the trimming step). Used chatgpt for the code below :)  

```
#!/bin/bash

# Define the input and output files
input_file="sed.collapse.cat.AST-1065.fastq"
output_file="17_sed.collapse.cat.AST-1065.fastq"

# Initialize the output file
> "$output_file"

# Use awk to process the sequences
awk '{
    if (substr($0, 1, 1) == ">") {
        header = $0
        getline
        sequence = $0
        if (length(sequence) >= 17) {
            print header >> "'$output_file'"
            print sequence >> "'$output_file'"
        }
    }
}' "$input_file"

zgrep -c ">" 17_sed.collapse.cat.AST-1065.fastq 
11979585
```

Retained all of the sequences. NOW lets attempt an mirdeep2 run. 

```
conda activate /data/putnamlab/mirdeep2
```

Map reads to genome first

```
mapper.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar/17_sed.collapse.cat.AST-1065.fastq -c -p /data/putnamlab/jillashey/Astrangia2021/smRNA/scripts/Apoc_ref.btindex -s 20240107_reads_collapsed.fa -t 20240107_reads_collapsed_vs_genome.arf -v  

discarding short reads
mapping reads to genome index
trimming unmapped nts in the 3' ends
Log file for this run is in mapper_logs and called mapper.log_56091
Mapping statistics

#desc	total	mapped	unmapped	%mapped	%unmapped
total: 35658222	3275049	32383173	9.185	90.815
seq: 35658222	3275049	32383173	9.185	90.815
```

Still got a pretty high % of unmapped reads...lets look at some of the files produced. 

```
head 20240107_reads_collapsed.fa 
>seq_1_x357414
TGGTCTATGGTGTAACTGGCAACACGTCTG
>seq_2_x138955
ACAGACGTGTTGCCAGTTACACCATAGACC
>seq_3_x125294
AACAGACGTGTTGCCAGTTACACCATAGAC
>seq_4_x98253
TGAAAATCTTTTCTCTGAAGTGGAA
>seq_5_x87633
TTCCACTTCAGAGAAAAGATTTTCA

zgrep -c ">" 20240107_reads_collapsed.fa 
11979585


head 20240107_reads_collapsed_vs_genome.arf 
seq_7_x81424	30	1	30	aacttttgacggtggatctcttggctcacg	chromosome_2	30	31480	31509	aacttttgacggtggatctcttggctcacg	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_7_x81424	30	1	30	aacttttgacggtggatctcttggctcacg	chromosome_2	30	42323	42352	aacttttgacggtggatctcttggctcacg	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_7_x81424	30	1	30	aacttttgacggtggatctcttggctcacg	chromosome_2	30	53072	53101	aacttttgacggtggatctcttggctcacg	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_7_x81424	30	1	30	aacttttgacggtggatctcttggctcacg	chromosome_2	30	20736	20765	aacttttgacggtggatctcttggctcacg	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_17_x38316	30	1	30	acaaatcttagaacaaaggcttaatctcag	chromosome_2	30	27510	27539	acaaatcttagaacaaaggcttaatctcag	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_17_x38316	30	1	30	acaaatcttagaacaaaggcttaatctcag	chromosome_2	30	49105	49134	acaaatcttagaacaaaggcttaatctcag	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_17_x38316	30	1	30	acaaatcttagaacaaaggcttaatctcag	chromosome_2	30	38360	38389	acaaatcttagaacaaaggcttaatctcag	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_17_x38316	30	1	30	acaaatcttagaacaaaggcttaatctcag	chromosome_2	30	879106	879135	acaaatcttagaacaaaggcttaatctcag	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_25_x32759	30	1	30	ttgctacgatcttctgagattaagcctttg	chromosome_2	30	879093	879122	ttgctacgatcttctgagattaagcctttg	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_25_x32759	30	1	30	ttgctacgatcttctgagattaagcctttg	chromosome_2	30	38373	38402	ttgctacgatcttctgagattaagcctttg	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

wc -l 20240107_reads_collapsed_vs_genome.arf 
1091161 20240107_reads_collapsed_vs_genome.arf
```

I'm not sure what the 20240107_reads_collapsed_vs_genome.arf file means. Let's see how many unique sequences are in that file 

```
cut -f1 20240107_reads_collapsed_vs_genome.arf | sort | uniq | wc -l 
523747
```

So ~500,000 unique sequences were mapped to the genome? That is how I am interpreting this. Let's try to run mirdeep2. 

```
miRDeep2.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar/17_sed.collapse.cat.AST-1065.fastq /data/putnamlab/jillashey/Astrangia_Genome/apoculata.assembly.scaffolds_chromosome_level.fasta /data/putnamlab/jillashey/Astrangia2021/smRNA/20240107_reads_collapsed_vs_genome.arf /data/putnamlab/jillashey/Astrangia2021/smRNA/refs/20240103_mature_T.fa none none -t N.vectensis -P -v -g -1 2>report.log
```

I'm going to let it run for ~10 mins then cut it off, as it probably takes a while (Sam White said his script took several days). I'm not sure if I can activate a conda env in a job script, emailed Kevin Bryan to ask. 

Cut the script off, this is as far as it got: 

```
#####################################
#                                   #
# miRDeep2.0.1.3                   #
#                                   #
# last change: 08/11/2019           #
#                                   #
#####################################

miRDeep2 started at 17:59:09


#Starting miRDeep2
#testing input files
#parsing genome mappings
#excising precursors
#preparing signature
#folding precursors
```

Yay, things are happening! Hopefully I can run this by the e5 meeting on Friday. 

### 20240107

Kevin Bryan confirmed that I can run a conda environment in a job script, I just need to add `-i` to the `#!/bin/bash` because of the way conda changes environments. Going to try this now on the test smRNA sample. In the scripts folder: `nano test_mirdeep2.sh`. 

```
#!/bin/bash -i
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/smRNA/scripts              
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load Miniconda3/4.9.2
conda activate /data/putnamlab/mirdeep2

echo "Starting mirdeep2 on test sample" $(date)

miRDeep2.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar/17_sed.collapse.cat.AST-1065.fastq /data/putnamlab/jillashey/Astrangia_Genome/apoculata.assembly.scaffolds_chromosome_level.fasta /data/putnamlab/jillashey/Astrangia2021/smRNA/20240107_reads_collapsed_vs_genome.arf /data/putnamlab/jillashey/Astrangia2021/smRNA/refs/20240103_mature_T.fa none none -t N.vectensis -P -v -g -1 2>report.log

echo "mirdeep2 concluded for test sample" $(date)

conda deactivate
```

Submitted batch job 292242. Job has been pending for a few mins and says its waiting for resources. 

### 20230109

After pending for about a day, mirdeep2 finally finished running on the test sample. It took about 3.5 hours to run. It created these folders/files in the scripts folder: 

```
drwxr-xr-x. 3 jillashey 4.0K Jan  9 07:49 mirdeep_runs
drwxr-xr-x. 2 jillashey 4.0K Jan  9 07:52 dir_prepare_signature1704804676
-rw-r--r--. 1 jillashey  377 Jan  9 11:23 error_09_01_2024_t_07_49_09.log
-rw-r--r--. 1 jillashey  63K Jan  9 11:24 result_09_01_2024_t_07_49_09.csv
-rw-r--r--. 1 jillashey 704K Jan  9 11:24 result_09_01_2024_t_07_49_09.html
drwxr-xr-x. 2 jillashey 4.0K Jan  9 11:28 pdfs_09_01_2024_t_07_49_09
-rw-r--r--. 1 jillashey  28K Jan  9 11:28 result_09_01_2024_t_07_49_09.bed
drwxr-xr-x. 2 jillashey 4.0K Jan  9 11:28 mirna_results_09_01_2024_t_07_49_09
-rw-r--r--. 1 jillashey  20K Jan  9 11:28 report.log
```

Let's look at them each. 

```
cd mirdeep_runs
cd run_09_01_2024_t_07_49_09 # folder 
ls
identified_precursors.fa  output.mrd  rfam_vs_precursor.bwt  run_09_01_2024_t_07_49_09_parameters  survey.csv
```

The identified precursors fasta file includes the precursor sequences (ie the part of the sequence that forms the pre-miRNA)

```
head identified_precursors.fa 
>chromosome_10_48090
ugauggagauggagaacgagaguggacuggacaguuuggcacugaagguucccuuuauaagcaguguuuuucuuucgacuacc
>M:chromosome_10_48090
TGTTTTTCTTTCGACTACC
>L:chromosome_10_48090
AGTGGACTGGACAGTTTGGCACTGAAGGTTCCCTTTATAAGCAG
>S:chromosome_10_48090
TGATGGAGATGGAGAACGAG
>chromosome_14_108653
cgcgcgcuauaguuacaguagcuauagcgcgcacuauaauuauagcagcuauagcgcacgcuauaguuagaaacuguagcgcgaguu

zgrep -c ">" identified_precursors.fa 
18695
```

Almost 19000 precursor sequences. 

In the output.mrd file, it has info on the different miRNAs identified I believe 

```
>chromosome_7_30929
score total                        2.4
score for star read(s)            -1.3
score for read counts    0
score for mfe                      2.1
score for randfold                 1.6
total read count                 13651
mature read count                13499
loop read count          0
star read count                    152
exp                                                                           fffffffffffffffffffMMMMMMMMMMMMMMMMMMMMlllllllllllllllSSSSSSSSSSSSSSSSSSSSffffffffffffffffffffffffff
ffffffffffff
obs                                                                           fffffffffffffffffffMMMMMMMMMMMMMMMMMMMMlllllllllllllSSSSSSSSSSSSSSSSSSSSSfffffffffffffffffffffffffff
ffffffffffff
pri_seq                                                                       cgcacugcaguugacgugaacccguagauccgaacuugugggauuuuucuccacaaguucggcuccaugguccacgugugcugugcucacaaacguugcu
acagcgugguca
pri_struct                                                                    .((((.(((....(((((...((((.((.((((((((((((((....)).)))))))))))).)).))))..))))).))).))))......((((((..
.)))))).....  #MM
seq_7183221_x1                                                                .................Uaacccguagauccgaacuugug............................................................
............  1
seq_2180915_x2                                                                ..................aacccguagauccgaacu................................................................
............  0
seq_9719163_x1                                                                ..................aacccguagauccgaUcuu...............................................................
............  1
seq_2267835_x2                                                                ..................Cacccguagauccgaacuu...............................................................
............  1
ola-miR-100_MIMAT0022614_Oryzias_latipes_miR-100                              ..................aacccguagauccgaacuu...............................................................
............  0
seq_126867_x21                                                                ..................aacccguagauccgaacuug..............................................................
............  0
seq_254962_x11                                                                ..................Cacccguagauccgaacuug..............................................................
............  1
sbo-miR-100_MIMAT0049501_Saimiri_boliviensis_miR-100                          ..................aacccguagauccgaacuugu.............................................................
............  0
dma-miR-100_MIMAT0049252_Daubentonia_madagascariensis_miR-100                 ..................aacccguagauccgaacuugu.............................................................
............  0
seq_199802_x14                                                                ..................aacccguagauccgaacuugC.............................................................
............  1
pmi-miR-100-5p_MIMAT0032156_Patiria_miniata_miR-100-5p                        ..................aacccguagauccgaacuugu.............................................................
............  0
seq_2153292_x2                                                                ..................aacccguagauccgaGcuugu.............................................................
............  1
seq_7747127_x1    

zgrep -c ">" output.mrd 
5231
```

The rfam vs precursor file includes information about where on the chromosomes the rRNAs and tRNAs are? 

```
head rfam_vs_precursor.bwt 
M:chromosome_13_66113	+	AM086652.1/1-576_RF00177;SSU_rRNA_5;	468	TGTTTCGGGATTGCAATG	IIIIIIIIIIIIIIIIII	3	
M:chromosome_13_66113	+	AF508778.1/21-597_RF00177;SSU_rRNA_5;	469	TGTTTCGGGATTGCAATG	IIIIIIIIIIIIIIIIII	3	
M:chromosome_13_66113	+	AJ310485.1/21-596_RF00177;SSU_rRNA_5;	468	TGTTTCGGGATTGCAATG	IIIIIIIIIIIIIIIIII	3	
M:chromosome_13_66113	+	DQ057346.1/21-597_RF00177;SSU_rRNA_5;	469	TGTTTCGGGATTGCAATG	IIIIIIIIIIIIIIIIII	3	
S:chromosome_13_66113	+	AACY021626480.1/155-83_RF00005;tRNA;	26	TTTGTTTCGTAAGCAAA	IIIIIIIIIIIIIIIII	2	7:T>C
S:chromosome_13_66113	+	AACY023301721.1/825-896_RF00005;tRNA;	25	TTTGTTTCGTAAGCAAA	IIIIIIIIIIIIIIIII	2	12:A>G
S:chromosome_13_66113	+	AACY022901721.1/116-188_RF00005;tRNA;	26	TTTGTTTCGTAAGCAAA	IIIIIIIIIIIIIIIII	2	12:A>G
M:chromosome_14_72583	+	CP000030.1/153914-153986_RF00005;tRNA;	3	CGGTTAGCTCAGTTGGTAGA	IIIIIIIIIIIIIIIIIIII	13	
M:chromosome_14_72583	+	AACY020037993.1/1235-1163_RF00005;tRNA;	3	CGGTTAGCTCAGTTGGTAGA	IIIIIIIIIIIIIIIIIIII	13	
M:chromosome_14_72583	+	AACY020166163.1/13-85_RF00005;tRNA;	3	CGGTTAGCTCAGTTGGTAGA	IIIIIIIIIIIIIIIIIIII	13

wc -l rfam_vs_precursor.bwt 
2241 rfam_vs_precursor.bwt
```

The run parameters file has the code specifics 

```
Start: 09_01_2024_t_07_49_09
Script  /data/putnamlab/mirdeep2/bin/miRDeep2.pl
args /data/putnamlab/mirdeep2/bin/miRDeep2.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar/17_sed.collapse.cat.AST-1065.fastq /data/putnamlab/jillashey/Astrangia_Genome/apoculata.assembly.scaffolds_chromosome_level.fasta /data/putnamlab/jillashey/Astrangia2021/smRNA/20240107_reads_collapsed_vs_genome.arf /data/putnamlab/jillashey/Astrangia2021/smRNA/refs/20240103_mature_T.fa none none -t N.vectensis -P -v -g -1

dir_with_tmp_files      dir_miRDeep2_09_01_2024_t_07_49_09
dir     /glfs/brick01/gv0/putnamlab/jillashey/Astrangia2021/smRNA/scripts
file_reads      /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar/17_sed.collapse.cat.AST-1065.fastq
file_genome     /data/putnamlab/jillashey/Astrangia_Genome/apoculata.assembly.scaffolds_chromosome_level.fasta
file_reads_vs_genome    /data/putnamlab/jillashey/Astrangia2021/smRNA/20240107_reads_collapsed_vs_genome.arf
file_mature_ref_this_species    /data/putnamlab/jillashey/Astrangia2021/smRNA/refs/20240103_mature_T.fa
file_mature_ref_other_species   none
option{t} =     N.vectensis
option{v} =     used
miRDeep runtime: 

started: 7:49:09
ended: 11:28:44
total:3h:39m:35s
```

The survey file inclues info about the mirdeep2 scores. This is the same info that is at the top of the hmtl and csv files. 

```
miRDeep2 score  novel miRNAs reported by miRDeep2       novel miRNAs, estimated false positives novel miRNAs, estimated true positives  known miRNAs in species known miRNAs in data    known miRNAs detected by miRDeep2       estimated signal-to-noise       excision gearing
10      49      3 +/- 2 46 +/- 2 (93 +/- 3%)    48885   83      1 (1%)  15.6    1
9       51      3 +/- 2 48 +/- 2 (93 +/- 3%)    48885   83      1 (1%)  15.1    1
8       56      4 +/- 2 52 +/- 2 (93 +/- 3%)    48885   83      1 (1%)  15.4    1
7       64      4 +/- 2 60 +/- 2 (94 +/- 3%)    48885   83      1 (1%)  16.3    1
6       68      4 +/- 2 64 +/- 2 (94 +/- 3%)    48885   83      1 (1%)  16      1
5       70      5 +/- 2 65 +/- 2 (93 +/- 3%)    48885   83      1 (1%)  15      1
4       72      5 +/- 2 67 +/- 2 (93 +/- 3%)    48885   83      1 (1%)  14      1
3       101     6 +/- 2 95 +/- 2 (94 +/- 2%)    48885   83      1 (1%)  16.2    1
2       183     9 +/- 3 174 +/- 3 (95 +/- 2%)   48885   83      72 (87%)        20.2    1
1       260     22 +/- 4        238 +/- 4 (91 +/- 2%)   48885   83      72 (87%)        11.6    1
0       315     53 +/- 6        262 +/- 6 (83 +/- 2%)   48885   83      72 (87%)        5.9     1
-1      365     97 +/- 8        268 +/- 8 (73 +/- 2%)   48885   83      72 (87%)        3.8     1
-2      462     150 +/- 10      312 +/- 10 (67 +/- 2%)  48885   83      72 (87%)        3.1     1
-3      707     228 +/- 14      479 +/- 14 (68 +/- 2%)  48885   83      72 (87%)        3.1     1
-4      1003    402 +/- 18      601 +/- 18 (60 +/- 2%)  48885   83      73 (88%)        2.5     1
-5      1167    750 +/- 23      417 +/- 23 (36 +/- 2%)  48885   83      73 (88%)        1.6     1
-6      1309    1227 +/- 35     82 +/- 35 (6 +/- 3%)    48885   83      73 (88%)        1.1     1
-7      1405    1733 +/- 42     0 +/- 0 (0 +/- 0%)      48885   83      73 (88%)        0.8     1
-8      1653    2208 +/- 45     0 +/- 0 (0 +/- 0%)      48885   83      73 (88%)        0.7     1
-9      2013    2615 +/- 49     0 +/- 0 (0 +/- 0%)      48885   83      73 (88%)        0.8     1
-10     2423    2951 +/- 53     0 +/- 0 (0 +/- 0%)      48885   83      73 (88%)        0.8     1
```

Going into the `dir_prepare_signature1704804676` from the scripts folder

```
cd dir_prepare_signature1704804676

ls
mature_vs_precursors.arf  precursors.ebwt.2.ebwt  precursors.ebwt.rev.1.ebwt  reads_vs_precursors.arf  signature_unsorted.arf.tmp
mature_vs_precursors.bwt  precursors.ebwt.3.ebwt  precursors.ebwt.rev.2.ebwt  reads_vs_precursors.bwt  signature_unsorted.arf.tmp2
precursors.ebwt.1.ebwt    precursors.ebwt.4.ebwt  precursors.fa               signature_unsorted.arf
```

Looked at the mature vs. precursors file

```
head mature_vs_precursors.arf 
hsa-miR-100-5p_MIMAT0000098_Homo_sapiens_miR-100-5p	22	1	22	aacccgtagatccgaacttgtg	chromosome_7_30929	22	19	40	aacccgtagatccgaacttgtg	+mmmmmmmmmmmmmmmmmmmmmm
hsa-miR-100-5p_MIMAT0000098_Homo_sapiens_miR-100-5p	22	1	22	aacccgtagatccgaacttgtg	chromosome_7_30930	22	69	90	aacccgtagatccgaacttgtg	+mmmmmmmmmmmmmmmmmmmmmm
mmu-miR-100-5p_MIMAT0000655_Mus_musculus_miR-100-5p	22	1	22	aacccgtagatccgaacttgtg	chromosome_7_30930	22	69	90	aacccgtagatccgaacttgtg	+mmmmmmmmmmmmmmmmmmmmmm
mmu-miR-100-5p_MIMAT0000655_Mus_musculus_miR-100-5p	22	1	22	aacccgtagatccgaacttgtg	chromosome_7_30929	22	19	40	aacccgtagatccgaacttgtg	+mmmmmmmmmmmmmmmmmmmmmm
rno-miR-100-5p_MIMAT0000822_Rattus_norvegicus_miR-100-5p	22	1	22	aacccgtagatccgaacttgtg	chromosome_7_30930	22	69	90	aacccgtagatccgaacttgtg	+	0	mmmmmmmmmmmmmmmmmmmmmm
rno-miR-100-5p_MIMAT0000822_Rattus_norvegicus_miR-100-5p	22	1	22	aacccgtagatccgaacttgtg	chromosome_7_30929	22	19	40	aacccgtagatccgaacttgtg	+	0	mmmmmmmmmmmmmmmmmmmmmm
gga-miR-100-5p_MIMAT0001178_Gallus_gallus_miR-100-5p	22	1	22	aacccgtagatccgaacttgtg	chromosome_7_30930	22	69	90	aacccgtagatccgaacttgtg	+mmmmmmmmmmmmmmmmmmmmmm
gga-miR-100-5p_MIMAT0001178_Gallus_gallus_miR-100-5p	22	1	22	aacccgtagatccgaacttgtg	chromosome_7_30929	22	19	40	aacccgtagatccgaacttgtg	+mmmmmmmmmmmmmmmmmmmmmm
aga-miR-100_MIMAT0001498_Anopheles_gambiae_miR-100	22	1	22	aacccgtagatccgaacttgtg	chromosome_7_30930	22	69	90	aacccgtagatccgaacttgtg	+mmmmmmmmmmmmmmmmmmmmmm
aga-miR-100_MIMAT0001498_Anopheles_gambiae_miR-100	22	1	22	aacccgtagatccgaacttgtg	chromosome_7_30929	22	19	40	aacccgtagatccgaacttgtg	+mmmmmmmmmmmmmmmmmmmmmm

wc -l mature_vs_precursors.arf 
191 mature_vs_precursors.arf
```

I'm not sure what this means...Need to look into this more. Is it providing info about the mirbase sequences in comparison to my own? It looks like most of them are related to miR-100, which makes sense as this is the only miRNA that is in bilaterians and cnidarians. 

Back in the scripts directory, look at the error file:

```
RNAfold: invalid option -- n
total number of rounds controls=100
1^M2^M3^M4^M5^M6^M7^M8^M9^M10^M11^M12^M13^M14^M15^M16^M17^M18^M19^M20^M21^M22^M23^M24^M25^M26^M27^M28^M29^M30^M31^M32^M33^M34^M35^M36^M37^M38^M39^M40^M41^M42^M43^M44^M45^M46^M47^M48^M49^M50^M51^M52^M53^M54^M55^M56^M57^M58^M59^M60^M61^M62^M63^M64^M65^M66^M67^M68^M69^M70^M71^M72^M73^M74^M75^M76^M77^M78^M79^M80^M81^M82^M83^M84^M85^M86^M87^M88^M89^M90^M91
^M92^M93^M94^M95^M96^M97^M98^M99^M100^Mcontrols performed
```

Not sure what this means either...some issue with the RNAfold option? But I got a randfold pvalue. 

The pdf folder includes a pdf file for each miRNA (?) identified and provides info about the scores and gives a nice graph about where the mature and star sequences are

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/astrangia2021_bioinf/smRNA/AST-1065_chromosome_11_50332.png)

In my first run of mirdeep2, 318 unique pdfs were produced. In the mirna_results_09_01_2024_t_07_49_09 results folder, there are .bed and .fa files. The bed files contain this info: 

```
less known_mature_09_01_2024_t_07_49_09_score-50_to_na.bed 

browser position chromosome_6:21902755-21902777
browser hide all
track name="notTrackname.known_miRNAs" description="known miRNAs detected by miRDeep2 for notTrackname" visibility=2
itemRgb="On";
chromosome_6    21902755        21902777        chromosome_6_22810      19.4    +       21902755        21902777        255,0,0
chromosome_7    21599230        21599250        chromosome_7_30929      2.4     -       21599230        21599250        0,0,255
chromosome_11   24411396        24411418        chromosome_11_54894     -3.5    -       24411396        24411418        0,0,255
```

and the fa files contain this info: 

```
less known_mature_09_01_2024_t_07_49_09_score-50_to_na.fa

>chromosome_6_22810
AAGAACACCCAAAATAGCTGAA
>chromosome_7_30929
ACCCGTAGATCCGAACTTGT
>chromosome_11_54894
GCGGGTGTGTGTGTGTGTGTGT
```

So I suppose its giving information about the known mature sequences and their location on the Astrangia genome. The other files in this folder contain bed and fa files for the known precursors, star and mature sequences, as well as the novel precursors, star and mature sequences. 

In the main scripts folder, there is also a file `result_09_01_2024_t_07_49_09` in .bed, .csv and .html format. It has the same info in all of them, but it is summarizing the parameters used, the survey info, the novel miRNAs and the known miRNAs. In the survey info, it says that 83 known miRNAs were identified in my data but only 3 are listed at the bottom? 

Anything that is listed as mirdeep2 results from Jan 9 are associated with AST-1065. 

I wonder if these results would change if I trimmed to 25 bp. I also wonder what would happen if I removed Nvectensis as my related species. Sam white put S.purpurtus as the related species, but Nvectensis is more closely related. I'm going to remove the related species and rerun mirdeep2. Submitted batch job 292316. Started running immediately, nice. 

I'm also going to rerun flexbar with trimming of 25 bp. In the trim data folder, make a folder for flexbar 25 bp

```
cd /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim
mkdir flexbar_25bp
```

In the scripts folder: `nano flexbar_25bp.sh`

```
#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=200GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load Flexbar/3.5.0-foss-2018b  

cd /data/putnamlab/jillashey/Astrangia2021/smRNA/data/raw

echo "Trimming reads to 25 bp using flexbar" $(date)

array1=($(ls *R1_001.fastq.gz))

for i in ${array1[@]}; do
flexbar \
-r ${i} \
-p $(echo ${i}|sed s/_R1/_R2/) \
-a NEB-adapters.fasta \
-ap ON \
-qf i1.8 \
-qt 30 \
-t /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar_25bp \
-k 25 \
-R trim.${i} \
-P trim.$(echo ${i}|sed s/_R1/_R2/) \
-z GZ
done

echo "Trimming complete" $(date)
```

Submitted batch job 292317. Immediately got an error saying: `ERROR: Could not open file trim.AST-2360_R1_001.fastq.gz.gz`. I commented out the `-z` command. Submitted batch job 292344

### 20230110

Hooray! Flexbar 25 bp and mirdeep2 test finished running. Let's look at the mirdeep2 results. As a reminder, I reran the mirdeep2 code but removed the specification of Nematostella as a related species. When removing this specification, I got marginally more (1-2 more) novel miRNAs predicted, but got the same number of known miRNAs identified. The mirdeep2 documentation states "it will in practice always improve miRDeep2 performance if miRNAs from some related species is input, even if it is not closely related." So it is likely best to keep Nematostella as a related species in the code. 

My next step is to run mirdeep2 on the newly trimmed (25 bp) reads. Again, I'm going to run it on a test sample. First, move the newly trimmed reads from the raw data folder to the flexbar 25 bp folder. 

```
cd /data/putnamlab/jillashey/Astrangia2021/smRNA/data/raw
mv trim* ../trim/flexbar_25bp/
```

Next run fastqc. In the fastqc folder, make a folder for the new flexbar 25bp QC data. 

```
cd /data/putnamlab/jillashey/Astrangia2021/smRNA/fastqc/trim
mkdir flexbar_25bp
cd flexbar_25bp
mkdir R1 R2
```

In the scripts folder, edit the `fastqc_trim.sh` script: 

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
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

echo "QC for trimmed reads using flexbar with max length of 25 bp" $(date)

for file in /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar_25bp/*fastq.gz
do
fastqc $file --outdir /data/putnamlab/jillashey/Astrangia2021/smRNA/fastqc/trim/flexbar_25bp
done

echo "FastQC complete" $(date)

#multiqc --interactive fastqc_results/trim/flexbar
```


### 20240112

Had to rerun the trimming bc I accidently left the 30 as the max length. Now im running the QC. Submitted batch job 292416

Once this has finished running, navigate to the flexbar 25bp fastqc folder and move the reads into R1 and R2 folders. 

```
cd /data/putnamlab/jillashey/Astrangia2021/smRNA/fastqc/trim/flexbar_25bp
mv *_R1* R1
mv *_R2* R2
```

Go into each folder and run MultiQC on each. I have to do this because for some reason, multiQC was just running QC stats on the R2 reads when the R1 and R2 reads were in the same folder. 

I think it would be a good idea to make a custom database that includes primarily cnidarian miRNAs, as there are not many in mirbase itself. [Baumgarten et al. (2017)](https://onlinelibrary.wiley.com/doi/10.1111/mec.14452) did something similar (under 2.3 miRNA annotation in their paper). 

- "miRNAs were then annotated using the miRDeep2 package with default settings (Friedlander et al., 2012). To identify putatively conserved miRNAs based on previous de novo annotations of other cnidarian genomes, we created a reference library of mature miRNA sequences from N. vectensis (Grimson et al., 2008; Moran et al., 2014), Hydra magnipapillata (Krishna et al., 2013) and Stylophora pistillata (Liew et al., 2014)." 

Therefore, I am going to make a custom database of the cnidarian miRNAs that I am aware of. 

- Stylophora pistillata - Table S11 from [Liew et al. 2014](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0091101)
- Aiptasia - Table S4 from [Baumgarten et al. 2017](https://onlinelibrary.wiley.com/doi/10.1111/mec.14452)
- Acropora digitifera - Table S5 from [Gajigan & Conaco 2017](https://onlinelibrary.wiley.com/doi/10.1111/mec.14130)
- the sea anemones Edwardsiella carnea, Scolanthus callimorphus, Metridium senile and Anemonia viridis, and the stony coral Acropora millepora - Tables S1 and S2 from [Praher et al. 2021](https://royalsocietypublishing.org/doi/full/10.1098/rspb.2020.3169#d1e1330)
- Hydra - Table S1 (SuppFile2) from [Krishna et al. 2012](https://academic.oup.com/nar/article/41/1/599/1146895#20334407)
- Nematostella - Table S1 from [Moran et al. 2014](https://genome.cshlp.org/content/24/4/651/suppl/DC1)
- Anemonia viridis - Table S2 from [Urbarova et al. 2018](https://academic.oup.com/gbe/article/10/2/410/4827693?login=true#supplementary-data)
- Nematostella + Hydra are also both on miRBase  

In google sheets, I gathered all the cnidarian miRNA sequences that I could find and made it into a csv file. This is the format: 

```
miRNA	Mature_miRNA_sequence	Species	Citation	Notes
spi-mir-temp-1	acccguagauccgaacuugugg	Stylophora pistillata	Liew et al. 2014	Matches miR-100 family.
spi-mir-temp-2	uaucgaauccgucaaaaagaga	Stylophora pistillata	Liew et al. 2014	NA
spi-mir-temp-3	ucagggauuguggugaguuaguu	Stylophora pistillata	Liew et al. 2014	NA
spi-mir-temp-4	aaagaaguacaagugguaggg	Stylophora pistillata	Liew et al. 2014	Exact match of nve-miR-2023.
spi-mir-temp-5	gagguccggaugguuga	Stylophora pistillata	Liew et al. 2014	NA
```

I downloaded the csv to my computer and manipulated it so that the first, third, fourth and fifth columns are the headers on one line and denoted with a ">". Then the sequence, in the second column, was put under the header. 

```
awk -F',' 'NR>1 {print ">"$1" "$2" "$3" "$4"\n"$5}' cnidarian_miRNAs.csv > cnidarian_miRNAs.fasta
```

Now go back to the miRbase fasta file and subset so that I make a file with Hydra and Nematostella sequences only 

```
awk '/>.*hma|>.*nve/ {print; getline; print}' 20240103_mature_T.fa > subset.fasta
```

I then copied the subset fasta info into the cnidarian_miRNAs.fasta. A complete cnidarian miRNA fasta! Reformat the fasta header names so there are no spaces. 

```
sed '/^>/ s/ /_/g' cnidarian_miRNAs.fasta \
| sed '/^>/ s/,//g' \
> cnidarian_miRNAs.fasta
```

Reformat sequences so that everything is uppercase

```
awk '/^>/ {print; getline; print toupper($0); next} {print}' cnidarian_miRNAs.fasta > cnidarian_miRNAs.fasta
```

I then copied the file to andromeda to `/data/putnamlab/jillashey/Astrangia2021/smRNA/refs`. Change the U to T in the fasta file

```
#!/bin/bash

# Define the input and output files
input_file="cnidarian_miRNA.fa"  # Replace with your actual input file name
output_file="cnidarian_miRNA_T.fa"

# Initialize the output file
> "$output_file"

# Use awk to process the file
awk '{
    if (substr($0, 1, 1) == ">") {
        print $0 >> "'$output_file'"  # Print the identifier as is
    } else {
        gsub(/U/, "T", $0)  # Replace U with T in sequences
        print $0 >> "'$output_file'"
    }
}' "$input_file"
```

Now I am going to run mirdeep2 with the cnidarian miRNA file! I'm going to modify the `nano test_mirdeep2.sh` so that the fasta file is the `cnidarian_miRNA_T.fa` file. Submitted batch job 292356. Didnt work. Removed the `-t` argument. Submitted batch job 292358...still not running. Getting this error: 

```
bash: cannot set terminal process group (-1): Function not implemented
bash: no job control in this shell
```

Need to troubleshoot this. 

When I look at the `report.log` file, it says: 

```
miRDeep2 started at 12:39:17


mkdir mirdeep_runs/run_12_01_2024_t_12_39_17

#testing input files
started: 12:39:23
sanity_check_mature_ref.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/refs/cnidarian_miRNA_T.fa

ESC[1;31mError: ESC[0mproblem with /data/putnamlab/jillashey/Astrangia2021/smRNA/refs/cnidarian_miRNA_T.fa
Error in line 64: The sequence

 SPIS


contains characters others than [acgtunACGTUN]

Please check your file for the following issues:

I.  Sequences are allowed only to comprise characters [ACGTNacgtn].
II. Identifiers are not allowed to have withespaces.
```

I looked at the `cnidarian_miRNA_T.fa` file and found that for a few sequences, it has Spis instead of the actual sequence: 

```
>spi-mir-temp-42 Stylophora pistillata Liew et al. 2014 NA
ugugcaagaauuugagucgcugg
>apa-mir-100 Exaiptasia pallida Baumgarten et al. 2017 "miR-100; Nve
 Spis
>apa-mir-2022a Exaiptasia pallida Baumgarten et al. 2017 "miR-2022; Nve
 Spis
>apa-mir-2023 Exaiptasia pallida Baumgarten et al. 2017 "miR-2023; Nve
 Spis
>apa-mir-2025 Exaiptasia pallida Baumgarten et al. 2017 "miR-2025; Nve
 Adi"
>apa-mir-2026 Exaiptasia pallida Baumgarten et al. 2017 miR-2026; Nve
aauuucaaauauccacugauug
>apa-mir-2030 Exaiptasia pallida Baumgarten et al. 2017 "miR-2030; Nve
 Spis
>apa-mir-2036 Exaiptasia pallida Baumgarten et al. 2017 "miR-2036; Nve
 Spis
>apa-mir-2037 Exaiptasia pallida Baumgarten et al. 2017 "miR-2037; Nve
 Spis"
>apa-mir-2050 Exaiptasia pallida Baumgarten et al. 2017 "miR-2050; Nve
 Spis
```

I'm thinking maybe it doesn't like the commas? Going back to csv and adding semi-colans instead of commas. 

I downloaded the csv to my computer and manipulated it so that the first, third, fourth and fifth columns are the headers on one line and denoted with a ">". Then the sequence, in the second column, was put under the header. 

```
awk -F',' 'NR>1 {print ">"$1" "$2" "$3" "$4"\n"$5}' cnidarian_miRNAs.csv > cnidarian_miRNAs.fasta
```

That seems to have fixed the problem. Reformat the fasta header names so there are no spaces. 

```
sed '/^>/ s/ /_/g' cnidarian_miRNAs.fasta \
| sed '/^>/ s/,//g' \
> cnidarian_miRNAs.fasta
```

Reformat sequences so that everything is uppercase

```
awk '/^>/ {print; getline; print toupper($0); next} {print}' cnidarian_miRNAs.fasta > cnidarian_miRNAs.fasta
```

I then copied the file to andromeda to `/data/putnamlab/jillashey/Astrangia2021/smRNA/refs`. Change the U to T in the fasta file

```
#!/bin/bash

# Define the input and output files
input_file="cnidarian_miRNAs.fa"  # Replace with your actual input file name
output_file="cnidarian_miRNAs_T.fa"

# Initialize the output file
> "$output_file"

# Use awk to process the file
awk '{
    if (substr($0, 1, 1) == ">") {
        print $0 >> "'$output_file'"  # Print the identifier as is
    } else {
        gsub(/U/, "T", $0)  # Replace U with T in sequences
        print $0 >> "'$output_file'"
    }
}' "$input_file"
```

Now concatenate the cnidarian miRNAs with the mature miRNA fasta from miRBase. 

```
cat 20240103_mature_T.fa cnidarian_miRNAs_T.fa > mature_mirbase_cnidarian_T.fa
```

Now edit the `test_mirdeep2.sh` so that `mature_mirbase_cnidarian_T.fa` is the input fasta file. Submitted batch job 292458. Failed. The report log file is telling me that I need to put none none before the `-t` flag. Changed that and resubmitted job. Submitted batch job 292459. Hooray appears to be running!

Sam White talked to Azenta (who did the sequencing for my project) and they recommended using Trimmomatic for trimming the small RNA reads. They also recommended tossing out read 2 and only using read 1 for analysis. I concatenate and collapse the reads anyway, so that shouldn't matter too much. Sam will get back to me about trimming info/code from Azenta soon. 

### 20240116

While I wait for Sam to get back to me about the trimming info, I'm going to run a mirdeep2 test on another sample (AST-2000). Unzip the files first. 

```
cd /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar
gunzip AST-2000_R1_001.fastq.*
```

Next, go to scripts folder and modify `test_cat_collapse.sh` so that the sample is AST-2000. Submitted batch job 292595. Took about 20 mins. 

```
head sed.collapse.cat.AST-2000.fastq
>seq_1_x719979
GCACTGGTGGTTCAGTGGTAGAATTCTC
>seq_2_x647300
GAGAATTCTACCACTGAACCACCAGTGC
>seq_3_x228711
GCACTGTGGTTCAGTGGTAGAATTCTC
>seq_4_x206224
GAGAATTCTACCACTGAACCACAGTGC
>seq_5_x161452
GCACTGGTGGTTCAGTGGTAGAATTCT

zgrep -c ">" sed.collapse.cat.AST-2000.fastq
10773578
```

Remove any sequences that are <17 nts. 

```
#!/bin/bash

# Define the input and output files
input_file="sed.collapse.cat.AST-2000.fastq"
output_file="17_sed.collapse.cat.AST-2000.fastq"

# Initialize the output file
> "$output_file"

# Use awk to process the sequences
awk '{
    if (substr($0, 1, 1) == ">") {
        header = $0
        getline
        sequence = $0
        if (length(sequence) >= 17) {
            print header >> "'$output_file'"
            print sequence >> "'$output_file'"
        }
    }
}' "$input_file"

zgrep -c ">" 17_sed.collapse.cat.AST-2000.fastq 
10773578
```

Looks like no sequences were removed. Now let's attempt the mirdeep2 run with AST-2000. Map reads to genome first.

```
conda activate /data/putnamlab/mirdeep2

mapper.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar/17_sed.collapse.cat.AST-2000.fastq -c -p /data/putnamlab/jillashey/Astrangia2021/smRNA/scripts/Apoc_ref.btindex -s 20240116_reads_collapsed.fa -t 20240116_reads_collapsed_vs_genome.arf -v  

discarding short reads
mapping reads to genome index
trimming unmapped nts in the 3' ends
Log file for this run is in mapper_logs and called mapper.log_5771
Mapping statistics

#desc	total	mapped	unmapped	%mapped	%unmapped
total: 34857708	5888651	28969057	16.893	83.107
seq: 34857708	5888651	28969057	16.893	83.107
```

Higher mapping % than AST-1065. Let's look at the files produced. 

```
head 20240116_reads_collapsed.fa
>seq_1_x719979
GCACTGGTGGTTCAGTGGTAGAATTCTC
>seq_2_x647300
GAGAATTCTACCACTGAACCACCAGTGC
>seq_3_x228711
GCACTGTGGTTCAGTGGTAGAATTCTC
>seq_4_x206224
GAGAATTCTACCACTGAACCACAGTGC
>seq_5_x161452
GCACTGGTGGTTCAGTGGTAGAATTCT

zgrep -c ">" 20240116_reads_collapsed.fa 
10773578

head 20240116_reads_collapsed_vs_genome.arf 
seq_7_x119692	30	1	30	aacttttgacggtggatctcttggctcacg	chromosome_2	30	31480	31509	aacttttgacggtggatctcttggctcacg	-	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_7_x119692	30	1	30	aacttttgacggtggatctcttggctcacg	chromosome_2	30	42323	42352	aacttttgacggtggatctcttggctcacg	-	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_7_x119692	30	1	30	aacttttgacggtggatctcttggctcacg	chromosome_2	30	53072	53101	aacttttgacggtggatctcttggctcacg	-	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_7_x119692	30	1	30	aacttttgacggtggatctcttggctcacg	chromosome_2	30	20736	20765	aacttttgacggtggatctcttggctcacg	-	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_15_x41699	30	1	30	tgcgtgagccaagagatccaccgtcaaaag	chromosome_2	30	31478	31507	tgcgtgagccaagagatccaccgtcaaaag	+	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_15_x41699	30	1	30	tgcgtgagccaagagatccaccgtcaaaag	chromosome_2	30	20734	20763	tgcgtgagccaagagatccaccgtcaaaag	+	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_15_x41699	30	1	30	tgcgtgagccaagagatccaccgtcaaaag	chromosome_2	30	42321	42350	tgcgtgagccaagagatccaccgtcaaaag	+	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_15_x41699	30	1	30	tgcgtgagccaagagatccaccgtcaaaag	chromosome_2	30	53070	53099	tgcgtgagccaagagatccaccgtcaaaag	+	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_27_x28319	30	1	30	tccgacactcagacagacatgctcctggga	chromosome_2	30	20610	20639	tccgacactcagacagacatgctcctggga	+	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_27_x28319	30	1	30	tccgacactcagacagacatgctcctggga	chromosome_2	30	42197	42226	tccgacactcagacagacatgctcctggga	+	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

cut -f1 20240116_reads_collapsed_vs_genome.arf | sort | uniq | wc -l 
1016532
```

Edit the `test_mirdeep2.sh` script to contain info for AST-2000

```
#!/bin/bash -i
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/lncRNA/scripts              
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

#module load Miniconda3/4.9.2
conda activate /data/putnamlab/mirdeep2

echo "Starting mirdeep2 on test sample AST-2000" $(date)

miRDeep2.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar/17_sed.collapse.cat.AST-2000.fastq /data/putnamlab/jillashey/Astrangia_Genome/apoculata.assembly.scaffolds_chromosome_level.fasta /data/putnamlab/jillashey/Astrangia2021/smRNA/20240116_reads_collapsed_vs_genome.arf /data/putnamlab/jillashey/Astrangia2021/smRNA/refs/mature_mirbase_cnidarian_T.fa none none -t N.vectensis -P -v -g -1 2>report.log

echo "mirdeep2 concluded for test sample AST-2000" $(date)

conda deactivate
```

Submitted batch job 292597. This took about 7 hours to run. 

I then compared the AST-1065 and AST-2000 samples in R to see if there were any overlapping sequences (code [here](https://github.com/JillAshey/Astrangia_repo/blob/main/scripts/compare_miRNAs.R)). I found 21 unique overlapping sequences between the two samples! Exciting. 

### 20240119

Hollie recommended that I compare the 25 bp vs 30 bp to see what the outcome was (ie did one trim length yield more miRNAs than the other). First, I need to prep the AST-2000 25bp sample for mirdeep2. 

Unzip the files. 

```
cd /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar_25bp
gunzip trim.AST-2000*
```

Next, go to scripts folder and modify `test_cat_collapse.sh` so that the sample is AST-2000 from the flexbar 25 bp folder. Submitted batch job 293019. Took about 15 mins

```
head sed.collapse.cat.AST-2000_25bp.fastq 
>seq_1_x1075819
GCACTGGTGGTTCAGTGGTAGAATT
>seq_2_x673045
GAGAATTCTACCACTGAACCACCAG
>seq_3_x339694
GCACTGTGGTTCAGTGGTAGAATTC
>seq_4_x212925
GAGAATTCTACCACTGAACCACAGT
>seq_5_x150098
AGAATTCTACCACTGAACCACCAGT

zgrep -c ">" sed.collapse.cat.AST-2000_25bp.fastq 
9235930
```

Remove any seqs that are <17 nts. 

```
#!/bin/bash

# Define the input and output files
input_file="sed.collapse.cat.AST-2000_25bp.fastq"
output_file="17_sed.collapse.cat.AST-2000_25bp.fastq"

# Initialize the output file
> "$output_file"

# Use awk to process the sequences
awk '{
    if (substr($0, 1, 1) == ">") {
        header = $0
        getline
        sequence = $0
        if (length(sequence) >= 17) {
            print header >> "'$output_file'"
            print sequence >> "'$output_file'"
        }
    }
}' "$input_file"

zgrep -c ">" 17_sed.collapse.cat.AST-2000_25bp.fastq 
9235930
```

No reads removed. Now map the reads to the genome using mirdeep2

```
conda activate /data/putnamlab/mirdeep2

mapper.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar_25bp/17_sed.collapse.cat.AST-2000_25bp.fastq -c -p /data/putnamlab/jillashey/Astrangia2021/smRNA/scripts/Apoc_ref.btindex -s 20240119_reads_collapsed.fa -t 20240119_reads_collapsed_vs_genome.arf -v  

discarding short reads
mapping reads to genome index
trimming unmapped nts in the 3' ends
Log file for this run is in mapper_logs and called mapper.log_51907
Mapping statistics

#desc	total	mapped	unmapped	%mapped	%unmapped
total: 34857708	6018696	28839012	17.266	82.734
seq: 34857708	6018696	28839012	17.266	82.734
```

Slightly higher mapping than the 30bp AST-2000. Let's look at the files 

```
head 20240119_reads_collapsed.fa
>seq_1_x1075819
GCACTGGTGGTTCAGTGGTAGAATT
>seq_2_x673045
GAGAATTCTACCACTGAACCACCAG
>seq_3_x339694
GCACTGTGGTTCAGTGGTAGAATTC
>seq_4_x212925
GAGAATTCTACCACTGAACCACAGT
>seq_5_x150098
AGAATTCTACCACTGAACCACCAGT

zgrep -c ">" 20240119_reads_collapsed.fa
9235930

head 20240119_reads_collapsed_vs_genome.arf
seq_6_x145559	25	1	25	aacttttgacggtggatctcttggc	chromosome_2	25	20741	20765	aacttttgacggtggatctcttggc	mmmmmmmmmmmmmmmmmmmmmmmmm
seq_6_x145559	25	1	25	aacttttgacggtggatctcttggc	chromosome_2	25	31485	31509	aacttttgacggtggatctcttggc	mmmmmmmmmmmmmmmmmmmmmmmmm
seq_6_x145559	25	1	25	aacttttgacggtggatctcttggc	chromosome_2	25	42328	42352	aacttttgacggtggatctcttggc	mmmmmmmmmmmmmmmmmmmmmmmmm
seq_6_x145559	25	1	25	aacttttgacggtggatctcttggc	chromosome_2	25	53077	53101	aacttttgacggtggatctcttggc	mmmmmmmmmmmmmmmmmmmmmmmmm
seq_18_x42856	25	1	25	tgcgtgagccaagagatccaccgtc	chromosome_2	25	42321	42345	tgcgtgagccaagagatccaccgtc	mmmmmmmmmmmmmmmmmmmmmmmmm
seq_18_x42856	25	1	25	tgcgtgagccaagagatccaccgtc	chromosome_2	25	53070	53094	tgcgtgagccaagagatccaccgtc	mmmmmmmmmmmmmmmmmmmmmmmmm
seq_18_x42856	25	1	25	tgcgtgagccaagagatccaccgtc	chromosome_2	25	31478	31502	tgcgtgagccaagagatccaccgtc	mmmmmmmmmmmmmmmmmmmmmmmmm
seq_18_x42856	25	1	25	tgcgtgagccaagagatccaccgtc	chromosome_2	25	20734	20758	tgcgtgagccaagagatccaccgtc	mmmmmmmmmmmmmmmmmmmmmmmmm
seq_27_x32954	25	1	25	tccgacactcagacagacatgctcc	chromosome_2	25	42197	42221	tccgacactcagacagacatgctcc	mmmmmmmmmmmmmmmmmmmmmmmmm
seq_27_x32954	25	1	25	tccgacactcagacagacatgctcc	chromosome_2	25	52946	52970	tccgacactcagacagacatgctcc	mmmmmmmmmmmmmmmmmmmmmmmmm

cut -f1 20240119_reads_collapsed_vs_genome.arf | sort | uniq | wc -l 
792395
```

Less reads collapsed vs genome for the 25 bp than the 30 bp AST-2000 sample. 

Edit the `test_mirdeep2.sh` script to contain info for AST-2000 25bp

```
#!/bin/bash -i
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/smRNA/scripts              
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

#module load Miniconda3/4.9.2
conda activate /data/putnamlab/mirdeep2

echo "Starting mirdeep2 on test sample AST-2000 trimmed to 25bp" $(date)

miRDeep2.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar_25bp/17_sed.collapse.cat.AST-2000_25bp.fastq /data/putnamlab/jillashey/Astrangia_Genome/apoculata.assembly.scaffolds_chromosome_level.fasta /data/putnamlab/jillashey/Astrangia2021/smRNA/20240119_reads_collapsed_vs_genome.arf /data/putnamlab/jillashey/Astrangia2021/smRNA/refs/mature_mirbase_cnidarian_T.fa none none -t N.vectensis -P -v -g -1 2>report.log

echo "mirdeep2 concluded for test sample AST-2000 trimmed to 25bp" $(date)

conda deactivate
```

Submitted batch job 293022

### 20240121

Took about 16 hours to run. Now doing to download the csv and compare the AST-2000 25 vs 30 bp results. 

AST-2000 25 bp
![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/astrangia2021_bioinf/smRNA/AST-2000_test_25bp.png)

AST-2000 30 bp
![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/astrangia2021_bioinf/smRNA/AST-2000_test_30bp.png)

The 25bp trimming had more miRNAs identified but more false positives. There were 111 unique sews when comparing 25bp v 30bp AST-2000. I think I am going to move forward with the 30bp trimming. 

Yay!!!!! Okay well now I can do the rest of the mirdeep2 runs for all the samples. Idk how I feel about putting all the samples in a loop because I don't want them to overwrite one another. I might do a separate script for each sample...If I do that, I'll need to make separate scripts for the mirdeep2 itself. 

In `/data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar`, gunzip all fastq files. 

I can write a single script to concatenate (with `cat` command) and collapse reads (with `fastx_collapser` from the [fastx toolkit](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html). In scripts folder: `nano cat_collapse.sh`. 

```
#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/smRNA/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load FASTX-Toolkit/0.0.14-GCC-9.3.0 

echo "Concatenate and collapse smRNA reads" $(date)

cd /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar

samples=$(ls *_R1_001.fastq.gz_1.fastq | sed 's/\(.*\)_R1_001.fastq.gz_1.fastq/\1/')

for sample in $samples
do
    # Concatenate paired-end reads
    cat "${sample}_R1_001.fastq.gz_1.fastq" "${sample}_R1_001.fastq.gz_2.fastq" > "cat.${sample}.fastq"
    echo "${sample} reads are concatenated into one file per sample" $(date)

    # Collapse concatenated reads
    fastx_collapser -v -i "cat.${sample}.fastq" -o "collapse.cat.${sample}.fastq"
    echo "${sample} reads collapsed" $(date)
done
```

Submitted batch job 293086

### 20240122

All files are now concatenated and collapsed. I need to now prep the sequence IDs for mirdeep2 with `sed`. I should've put this in the script. I'm also going to remove any reads <17nts. 

AST-1065

```
sed '/^>/ s/-/_x/g' collapse.cat.AST-1065.fastq \
| sed '/^>/ s/>/>seq_/' \
> sed.collapse.cat.AST-1065.fastq

#!/bin/bash

# Define the input and output files
input_file="sed.collapse.cat.AST-1065.fastq"
output_file="17_sed.collapse.cat.AST-1065.fastq"

# Initialize the output file
> "$output_file"

# Use awk to process the sequences
awk '{
    if (substr($0, 1, 1) == ">") {
        header = $0
        getline
        sequence = $0
        if (length(sequence) >= 17) {
            print header >> "'$output_file'"
            print sequence >> "'$output_file'"
        }
    }
}' "$input_file"

zgrep -c ">" 17_sed.collapse.cat.AST-1065.fastq 
11979585
```

Not doing AST-1105 bc of poor QC and mapping 

AST-1147

```
sed '/^>/ s/-/_x/g' collapse.cat.AST-1147.fastq \
| sed '/^>/ s/>/>seq_/' \
> sed.collapse.cat.AST-1147.fastq

#!/bin/bash

# Define the input and output files
input_file="sed.collapse.cat.AST-1147.fastq"
output_file="17_sed.collapse.cat.AST-1147.fastq"

# Initialize the output file
> "$output_file"

# Use awk to process the sequences
awk '{
    if (substr($0, 1, 1) == ">") {
        header = $0
        getline
        sequence = $0
        if (length(sequence) >= 17) {
            print header >> "'$output_file'"
            print sequence >> "'$output_file'"
        }
    }
}' "$input_file"

zgrep -c ">" 17_sed.collapse.cat.AST-1147.fastq 
18025765
```

AST-1412

```
sed '/^>/ s/-/_x/g' collapse.cat.AST-1412.fastq \
| sed '/^>/ s/>/>seq_/' \
> sed.collapse.cat.AST-1412.fastq

#!/bin/bash

# Define the input and output files
input_file="sed.collapse.cat.AST-1412.fastq"
output_file="17_sed.collapse.cat.AST-1412.fastq"

# Initialize the output file
> "$output_file"

# Use awk to process the sequences
awk '{
    if (substr($0, 1, 1) == ">") {
        header = $0
        getline
        sequence = $0
        if (length(sequence) >= 17) {
            print header >> "'$output_file'"
            print sequence >> "'$output_file'"
        }
    }
}' "$input_file"

zgrep -c ">" 17_sed.collapse.cat.AST-1412.fastq 
11251691
```

AST-1560

```
sed '/^>/ s/-/_x/g' collapse.cat.AST-1560.fastq \
| sed '/^>/ s/>/>seq_/' \
> sed.collapse.cat.AST-1560.fastq

#!/bin/bash

# Define the input and output files
input_file="sed.collapse.cat.AST-1560.fastq"
output_file="17_sed.collapse.cat.AST-1560.fastq"

# Initialize the output file
> "$output_file"

# Use awk to process the sequences
awk '{
    if (substr($0, 1, 1) == ">") {
        header = $0
        getline
        sequence = $0
        if (length(sequence) >= 17) {
            print header >> "'$output_file'"
            print sequence >> "'$output_file'"
        }
    }
}' "$input_file"

zgrep -c ">" 17_sed.collapse.cat.AST-1560.fastq 
11702120
```

AST-1567

```
sed '/^>/ s/-/_x/g' collapse.cat.AST-1567.fastq \
| sed '/^>/ s/>/>seq_/' \
> sed.collapse.cat.AST-1567.fastq

#!/bin/bash

# Define the input and output files
input_file="sed.collapse.cat.AST-1567.fastq"
output_file="17_sed.collapse.cat.AST-1567.fastq"

# Initialize the output file
> "$output_file"

# Use awk to process the sequences
awk '{
    if (substr($0, 1, 1) == ">") {
        header = $0
        getline
        sequence = $0
        if (length(sequence) >= 17) {
            print header >> "'$output_file'"
            print sequence >> "'$output_file'"
        }
    }
}' "$input_file"

zgrep -c ">" 17_sed.collapse.cat.AST-1567.fastq 
9871465
```

AST-1617

```
sed '/^>/ s/-/_x/g' collapse.cat.AST-1617.fastq \
| sed '/^>/ s/>/>seq_/' \
> sed.collapse.cat.AST-1617.fastq

#!/bin/bash

# Define the input and output files
input_file="sed.collapse.cat.AST-1617.fastq"
output_file="17_sed.collapse.cat.AST-1617.fastq"

# Initialize the output file
> "$output_file"

# Use awk to process the sequences
awk '{
    if (substr($0, 1, 1) == ">") {
        header = $0
        getline
        sequence = $0
        if (length(sequence) >= 17) {
            print header >> "'$output_file'"
            print sequence >> "'$output_file'"
        }
    }
}' "$input_file"

zgrep -c ">" 17_sed.collapse.cat.AST-1617.fastq 
8081542
```

AST-1722

```
sed '/^>/ s/-/_x/g' collapse.cat.AST-1722.fastq \
| sed '/^>/ s/>/>seq_/' \
> sed.collapse.cat.AST-1722.fastq

#!/bin/bash

# Define the input and output files
input_file="sed.collapse.cat.AST-1722.fastq"
output_file="17_sed.collapse.cat.AST-1722.fastq"

# Initialize the output file
> "$output_file"

# Use awk to process the sequences
awk '{
    if (substr($0, 1, 1) == ">") {
        header = $0
        getline
        sequence = $0
        if (length(sequence) >= 17) {
            print header >> "'$output_file'"
            print sequence >> "'$output_file'"
        }
    }
}' "$input_file"

zgrep -c ">" 17_sed.collapse.cat.AST-1722.fastq 
7574236
```

AST-2000

```
sed '/^>/ s/-/_x/g' collapse.cat.AST-2000.fastq \
| sed '/^>/ s/>/>seq_/' \
> sed.collapse.cat.AST-2000.fastq

#!/bin/bash

# Define the input and output files
input_file="sed.collapse.cat.AST-2000.fastq"
output_file="17_sed.collapse.cat.AST-2000.fastq"

# Initialize the output file
> "$output_file"

# Use awk to process the sequences
awk '{
    if (substr($0, 1, 1) == ">") {
        header = $0
        getline
        sequence = $0
        if (length(sequence) >= 17) {
            print header >> "'$output_file'"
            print sequence >> "'$output_file'"
        }
    }
}' "$input_file"

zgrep -c ">" 17_sed.collapse.cat.AST-2000.fastq 
10773578
```

AST-2007

```
sed '/^>/ s/-/_x/g' collapse.cat.AST-2007.fastq \
| sed '/^>/ s/>/>seq_/' \
> sed.collapse.cat.AST-2007.fastq

#!/bin/bash

# Define the input and output files
input_file="sed.collapse.cat.AST-2007.fastq"
output_file="17_sed.collapse.cat.AST-2007.fastq"

# Initialize the output file
> "$output_file"

# Use awk to process the sequences
awk '{
    if (substr($0, 1, 1) == ">") {
        header = $0
        getline
        sequence = $0
        if (length(sequence) >= 17) {
            print header >> "'$output_file'"
            print sequence >> "'$output_file'"
        }
    }
}' "$input_file"

zgrep -c ">" 17_sed.collapse.cat.AST-2007.fastq 
8653745
```

AST-2302

```
sed '/^>/ s/-/_x/g' collapse.cat.AST-2302.fastq \
| sed '/^>/ s/>/>seq_/' \
> sed.collapse.cat.AST-2302.fastq

#!/bin/bash

# Define the input and output files
input_file="sed.collapse.cat.AST-2302.fastq"
output_file="17_sed.collapse.cat.AST-2302.fastq"

# Initialize the output file
> "$output_file"

# Use awk to process the sequences
awk '{
    if (substr($0, 1, 1) == ">") {
        header = $0
        getline
        sequence = $0
        if (length(sequence) >= 17) {
            print header >> "'$output_file'"
            print sequence >> "'$output_file'"
        }
    }
}' "$input_file"

zgrep -c ">" 17_sed.collapse.cat.AST-2302.fastq 
11391780
```

AST-2360

```
sed '/^>/ s/-/_x/g' collapse.cat.AST-2360.fastq \
| sed '/^>/ s/>/>seq_/' \
> sed.collapse.cat.AST-2360.fastq

#!/bin/bash

# Define the input and output files
input_file="sed.collapse.cat.AST-2360.fastq"
output_file="17_sed.collapse.cat.AST-2360.fastq"

# Initialize the output file
> "$output_file"

# Use awk to process the sequences
awk '{
    if (substr($0, 1, 1) == ">") {
        header = $0
        getline
        sequence = $0
        if (length(sequence) >= 17) {
            print header >> "'$output_file'"
            print sequence >> "'$output_file'"
        }
    }
}' "$input_file"

zgrep -c ">" 17_sed.collapse.cat.AST-2360.fastq 
9737775
```

AST-2398

```
sed '/^>/ s/-/_x/g' collapse.cat.AST-2398.fastq \
| sed '/^>/ s/>/>seq_/' \
> sed.collapse.cat.AST-2398.fastq

#!/bin/bash

# Define the input and output files
input_file="sed.collapse.cat.AST-2398.fastq"
output_file="17_sed.collapse.cat.AST-2398.fastq"

# Initialize the output file
> "$output_file"

# Use awk to process the sequences
awk '{
    if (substr($0, 1, 1) == ">") {
        header = $0
        getline
        sequence = $0
        if (length(sequence) >= 17) {
            print header >> "'$output_file'"
            print sequence >> "'$output_file'"
        }
    }
}' "$input_file"

zgrep -c ">" 17_sed.collapse.cat.AST-2398.fastq
9507430
```

AST-2404

```
sed '/^>/ s/-/_x/g' collapse.cat.AST-2404.fastq \
| sed '/^>/ s/>/>seq_/' \
> sed.collapse.cat.AST-2404.fastq

#!/bin/bash

# Define the input and output files
input_file="sed.collapse.cat.AST-2404.fastq"
output_file="17_sed.collapse.cat.AST-2404.fastq"

# Initialize the output file
> "$output_file"

# Use awk to process the sequences
awk '{
    if (substr($0, 1, 1) == ">") {
        header = $0
        getline
        sequence = $0
        if (length(sequence) >= 17) {
            print header >> "'$output_file'"
            print sequence >> "'$output_file'"
        }
    }
}' "$input_file"

zgrep -c ">" 17_sed.collapse.cat.AST-2404.fastq
11599063
```

AST-2412

```
sed '/^>/ s/-/_x/g' collapse.cat.AST-2412.fastq \
| sed '/^>/ s/>/>seq_/' \
> sed.collapse.cat.AST-2412.fastq

#!/bin/bash

# Define the input and output files
input_file="sed.collapse.cat.AST-2412.fastq"
output_file="17_sed.collapse.cat.AST-2412.fastq"

# Initialize the output file
> "$output_file"

# Use awk to process the sequences
awk '{
    if (substr($0, 1, 1) == ">") {
        header = $0
        getline
        sequence = $0
        if (length(sequence) >= 17) {
            print header >> "'$output_file'"
            print sequence >> "'$output_file'"
        }
    }
}' "$input_file"

zgrep -c ">" 17_sed.collapse.cat.AST-2412.fastq
10485205
```

AST-2512

```
sed '/^>/ s/-/_x/g' collapse.cat.AST-2512.fastq \
| sed '/^>/ s/>/>seq_/' \
> sed.collapse.cat.AST-2512.fastq

#!/bin/bash

# Define the input and output files
input_file="sed.collapse.cat.AST-2512.fastq"
output_file="17_sed.collapse.cat.AST-2512.fastq"

# Initialize the output file
> "$output_file"

# Use awk to process the sequences
awk '{
    if (substr($0, 1, 1) == ">") {
        header = $0
        getline
        sequence = $0
        if (length(sequence) >= 17) {
            print header >> "'$output_file'"
            print sequence >> "'$output_file'"
        }
    }
}' "$input_file"

zgrep -c ">" 17_sed.collapse.cat.AST-2512.fastq
9710109
```

AST-2523

```
sed '/^>/ s/-/_x/g' collapse.cat.AST-2523.fastq \
| sed '/^>/ s/>/>seq_/' \
> sed.collapse.cat.AST-2523.fastq

#!/bin/bash

# Define the input and output files
input_file="sed.collapse.cat.AST-2523.fastq"
output_file="17_sed.collapse.cat.AST-2523.fastq"

# Initialize the output file
> "$output_file"

# Use awk to process the sequences
awk '{
    if (substr($0, 1, 1) == ">") {
        header = $0
        getline
        sequence = $0
        if (length(sequence) >= 17) {
            print header >> "'$output_file'"
            print sequence >> "'$output_file'"
        }
    }
}' "$input_file"

zgrep -c ">" 17_sed.collapse.cat.AST-2523.fastq
11748402
```

AST-2563

```
sed '/^>/ s/-/_x/g' collapse.cat.AST-2563.fastq \
| sed '/^>/ s/>/>seq_/' \
> sed.collapse.cat.AST-2563.fastq

#!/bin/bash

# Define the input and output files
input_file="sed.collapse.cat.AST-2563.fastq"
output_file="17_sed.collapse.cat.AST-2563.fastq"

# Initialize the output file
> "$output_file"

# Use awk to process the sequences
awk '{
    if (substr($0, 1, 1) == ">") {
        header = $0
        getline
        sequence = $0
        if (length(sequence) >= 17) {
            print header >> "'$output_file'"
            print sequence >> "'$output_file'"
        }
    }
}' "$input_file"

zgrep -c ">" 17_sed.collapse.cat.AST-2563.fastq
11689234
```

AST-2729

```
sed '/^>/ s/-/_x/g' collapse.cat.AST-2729.fastq \
| sed '/^>/ s/>/>seq_/' \
> sed.collapse.cat.AST-2729.fastq

#!/bin/bash

# Define the input and output files
input_file="sed.collapse.cat.AST-2729.fastq"
output_file="17_sed.collapse.cat.AST-2729.fastq"

# Initialize the output file
> "$output_file"

# Use awk to process the sequences
awk '{
    if (substr($0, 1, 1) == ">") {
        header = $0
        getline
        sequence = $0
        if (length(sequence) >= 17) {
            print header >> "'$output_file'"
            print sequence >> "'$output_file'"
        }
    }
}' "$input_file"

zgrep -c ">" 17_sed.collapse.cat.AST-2729.fastq
10703305
```

AST-2755

```
sed '/^>/ s/-/_x/g' collapse.cat.AST-2755.fastq \
| sed '/^>/ s/>/>seq_/' \
> sed.collapse.cat.AST-2755.fastq

#!/bin/bash

# Define the input and output files
input_file="sed.collapse.cat.AST-2755.fastq"
output_file="17_sed.collapse.cat.AST-2755.fastq"

# Initialize the output file
> "$output_file"

# Use awk to process the sequences
awk '{
    if (substr($0, 1, 1) == ">") {
        header = $0
        getline
        sequence = $0
        if (length(sequence) >= 17) {
            print header >> "'$output_file'"
            print sequence >> "'$output_file'"
        }
    }
}' "$input_file"

zgrep -c ">" 17_sed.collapse.cat.AST-2755.fastq
10834414
```

Okay everything is prepped! I'm going to remove the extra files made (ie cat, collapse, sed iterations of the file). Going back to the main folder, I'm going to make a `mirdeep2` folder that I will run all of my mirdeep2 code in so that the output is all in one place. 

```
cd /data/putnamlab/jillashey/Astrangia2021/smRNA
mkdir mirdeep2
cd mirdeep2
mkdir AST-1065 AST-1147 AST-1412 AST-1560 AST-1567 AST-1617 AST-1722 AST-2000 AST-2007 AST-2302 AST-2360 AST-2398 AST-2404 AST-2412 AST-2512 AST-2523 AST-2563 AST-2729 AST-2755
```

Now let's run mirdeep2 on all of the samples! I am paranoid that if I run the script in a loop for all the samples, something is going to overwrite. Therefore, I am going to run each sample individually. 

#### AST-1065

```
cd AST-1065
```

Run mapping step 

```
conda activate /data/putnamlab/mirdeep2

mapper.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar/17_sed.collapse.cat.AST-1065.fastq -c -p /data/putnamlab/jillashey/Astrangia2021/smRNA/scripts/Apoc_ref.btindex -s AST-1065_reads_collapsed.fa -t AST-1065_reads_collapsed_vs_genome.arf -v  

discarding short reads
mapping reads to genome index
trimming unmapped nts in the 3' ends
Log file for this run is in mapper_logs and called mapper.log_51087
Mapping statistics

#desc	total	mapped	unmapped	%mapped	%unmapped
total: 35658222	3275049	32383173	9.185	90.815
seq: 35658222	3275049	32383173	9.185	90.815
```

Look at the mapping results

```
head AST-1065_reads_collapsed.fa 
>seq_1_x357414
TGGTCTATGGTGTAACTGGCAACACGTCTG
>seq_2_x138955
ACAGACGTGTTGCCAGTTACACCATAGACC
>seq_3_x125294
AACAGACGTGTTGCCAGTTACACCATAGAC
>seq_4_x98253
TGAAAATCTTTTCTCTGAAGTGGAA
>seq_5_x87633
TTCCACTTCAGAGAAAAGATTTTCA

zgrep -c ">" AST-1065_reads_collapsed.fa 
11979585

head AST-1065_reads_collapsed_vs_genome.arf 
seq_7_x81424	30	1	30	aacttttgacggtggatctcttggctcacg	chromosome_2	30	31480	31509	aacttttgacggtggatctcttggctcacg	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_7_x81424	30	1	30	aacttttgacggtggatctcttggctcacg	chromosome_2	30	42323	42352	aacttttgacggtggatctcttggctcacg	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_7_x81424	30	1	30	aacttttgacggtggatctcttggctcacg	chromosome_2	30	53072	53101	aacttttgacggtggatctcttggctcacg	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_7_x81424	30	1	30	aacttttgacggtggatctcttggctcacg	chromosome_2	30	20736	20765	aacttttgacggtggatctcttggctcacg	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_17_x38316	30	1	30	acaaatcttagaacaaaggcttaatctcag	chromosome_2	30	27510	27539	acaaatcttagaacaaaggcttaatctcag	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_17_x38316	30	1	30	acaaatcttagaacaaaggcttaatctcag	chromosome_2	30	49105	49134	acaaatcttagaacaaaggcttaatctcag	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_17_x38316	30	1	30	acaaatcttagaacaaaggcttaatctcag	chromosome_2	30	38360	38389	acaaatcttagaacaaaggcttaatctcag	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_17_x38316	30	1	30	acaaatcttagaacaaaggcttaatctcag	chromosome_2	30	879106	879135	acaaatcttagaacaaaggcttaatctcag	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_25_x32759	30	1	30	ttgctacgatcttctgagattaagcctttg	chromosome_2	30	879093	879122	ttgctacgatcttctgagattaagcctttg	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_25_x32759	30	1	30	ttgctacgatcttctgagattaagcctttg	chromosome_2	30	38373	38402	ttgctacgatcttctgagattaagcctttg	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

cut -f1 AST-1065_reads_collapsed_vs_genome.arf | sort | uniq | wc -l 
523747

conda deactivate
```

Run mirdeep2. `nano AST-1065_mirdeep2.sh`

```
#!/bin/bash -i
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/smRNA/mirdeep2/AST-1065              
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

#module load Miniconda3/4.9.2
conda activate /data/putnamlab/mirdeep2

echo "Starting mirdeep2 on sample AST-1065 trimmed to 30bp" $(date)

miRDeep2.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar/17_sed.collapse.cat.AST-1065.fastq /data/putnamlab/jillashey/Astrangia_Genome/apoculata.assembly.scaffolds_chromosome_level.fasta /data/putnamlab/jillashey/Astrangia2021/smRNA/mirdeep2/AST-1065/AST-1065_reads_collapsed_vs_genome.arf /data/putnamlab/jillashey/Astrangia2021/smRNA/refs/mature_mirbase_cnidarian_T.fa none none -t N.vectensis -P -v -g -1 2>report.log

echo "mirdeep2 concluded for sample AST-1065 trimmed to 30bp" $(date)

conda deactivate
```

Submitted batch job 293148

#### AST-1147

```
cd /data/putnamlab/jillashey/Astrangia2021/smRNA/mirdeep2/AST-1147
```

Run mapping step 

```
conda activate /data/putnamlab/mirdeep2

mapper.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar/17_sed.collapse.cat.AST-1147.fastq -c -p /data/putnamlab/jillashey/Astrangia2021/smRNA/scripts/Apoc_ref.btindex -s AST-1147_reads_collapsed.fa -t AST-1147_reads_collapsed_vs_genome.arf -v  

discarding short reads
mapping reads to genome index
trimming unmapped nts in the 3' ends
Log file for this run is in mapper_logs and called mapper.log_53573
Mapping statistics

#desc	total	mapped	unmapped	%mapped	%unmapped
total: 80830448	23638038	57192410	29.244	70.756
seq: 80830448	23638038	57192410	29.244	70.756
```

Look at the mapping results

```
head AST-1147_reads_collapsed.fa 
>seq_1_x3434945
GCACTGGTGGTTCAGTGGTAGAATTCTCGC
>seq_2_x1817679
GGCGAGAATTCTACCACTGAACCACCAGTG
>seq_3_x926642
AGCGAGAATTCTACCACTGAACCACCAGTG
>seq_4_x602263
GCACTGTGGTTCAGTGGTAGAATTCTCGCC
>seq_5_x475440
GGCGAGAATTCTACCACTGAACCACAGTGC

zgrep -c ">" AST-1147_reads_collapsed.fa 
18025765

head AST-1147_reads_collapsed_vs_genome.arf 
seq_6_x427578	30	1	30	aacttttgacggtggatctcttggctcacg	chromosome_2	30	31480	31509	aacttttgacggtggatctcttggctcacg	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_6_x427578	30	1	30	aacttttgacggtggatctcttggctcacg	chromosome_2	30	42323	42352	aacttttgacggtggatctcttggctcacg	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_6_x427578	30	1	30	aacttttgacggtggatctcttggctcacg	chromosome_2	30	53072	53101	aacttttgacggtggatctcttggctcacg	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_6_x427578	30	1	30	aacttttgacggtggatctcttggctcacg	chromosome_2	30	20736	20765	aacttttgacggtggatctcttggctcacg	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_14_x202590	30	1	30	tgcgtgagccaagagatccaccgtcaaaag	chromosome_2	30	31478	31507	tgcgtgagccaagagatccaccgtcaaaag	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_14_x202590	30	1	30	tgcgtgagccaagagatccaccgtcaaaag	chromosome_2	30	20734	20763	tgcgtgagccaagagatccaccgtcaaaag	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_14_x202590	30	1	30	tgcgtgagccaagagatccaccgtcaaaag	chromosome_2	30	42321	42350	tgcgtgagccaagagatccaccgtcaaaag	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_14_x202590	30	1	30	tgcgtgagccaagagatccaccgtcaaaag	chromosome_2	30	53070	53099	tgcgtgagccaagagatccaccgtcaaaag	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_17_x127492	30	1	30	acaaatcttagaacaaaggcttaatctcag	chromosome_2	30	879106	879135	acaaatcttagaacaaaggcttaatctcag	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_17_x127492	30	1	30	acaaatcttagaacaaaggcttaatctcag	chromosome_2	30	38360	38389	acaaatcttagaacaaaggcttaatctcag	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

cut -f1 AST-1147_reads_collapsed_vs_genome.arf | sort | uniq | wc -l 
3227088

conda deactivate
```

Run mirdeep2. `nano AST-1147_mirdeep2.sh`

```
#!/bin/bash -i
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/smRNA/mirdeep2/AST-1147            
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

#module load Miniconda3/4.9.2
conda activate /data/putnamlab/mirdeep2

echo "Starting mirdeep2 on sample AST-1147 trimmed to 30bp" $(date)

miRDeep2.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar/17_sed.collapse.cat.AST-1147.fastq /data/putnamlab/jillashey/Astrangia_Genome/apoculata.assembly.scaffolds_chromosome_level.fasta /data/putnamlab/jillashey/Astrangia2021/smRNA/mirdeep2/AST-1147/AST-1147_reads_collapsed_vs_genome.arf /data/putnamlab/jillashey/Astrangia2021/smRNA/refs/mature_mirbase_cnidarian_T.fa none none -t N.vectensis -P -v -g -1 2>report.log

echo "mirdeep2 concluded for sample AST-1147 trimmed to 30bp" $(date)

conda deactivate
```

Submitted batch job 293149

#### AST-1412

```
cd /data/putnamlab/jillashey/Astrangia2021/smRNA/mirdeep2/AST-1412
```

Run mapping step 

```
conda activate /data/putnamlab/mirdeep2

mapper.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar/17_sed.collapse.cat.AST-1412.fastq -c -p /data/putnamlab/jillashey/Astrangia2021/smRNA/scripts/Apoc_ref.btindex -s AST-1412_reads_collapsed.fa -t AST-1412_reads_collapsed_vs_genome.arf -v  

discarding short reads
mapping reads to genome index
trimming unmapped nts in the 3' ends
Log file for this run is in mapper_logs and called mapper.log_56487
Mapping statistics

#desc	total	mapped	unmapped	%mapped	%unmapped
total: 32559110	3100518	29458592	9.523	90.477
seq: 32559110	3100518	29458592	9.523	90.477
```

Look at the mapping results

```
head AST-1412_reads_collapsed.fa 
>seq_1_x137015
GGAAGAGCACACGTCTGAACTCCAGTCACT
>seq_2_x100202
AACTTTTGACGGTGGATCTCTTGGCTCACG
>seq_3_x94707
TCGGACTGTAGAACTCTGAACGTGTAGATC
>seq_4_x78147
GCACTGGTGGTTCAGTGGTAGAATTCTCGC
>seq_5_x73890
TACTGGATAACTAAGGGAAAGTTTGGCTAA

zgrep -c ">" AST-1412_reads_collapsed.fa 
11251691

head AST-1412_reads_collapsed_vs_genome.arf 
seq_2_x100202	30	1	30	aacttttgacggtggatctcttggctcacg	chromosome_2	30	31480	31509	aacttttgacggtggatctcttggctcacg	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_2_x100202	30	1	30	aacttttgacggtggatctcttggctcacg	chromosome_2	30	42323	42352	aacttttgacggtggatctcttggctcacg	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_2_x100202	30	1	30	aacttttgacggtggatctcttggctcacg	chromosome_2	30	53072	53101	aacttttgacggtggatctcttggctcacg	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_2_x100202	30	1	30	aacttttgacggtggatctcttggctcacg	chromosome_2	30	20736	20765	aacttttgacggtggatctcttggctcacg	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_19_x39468	30	1	30	tgcgtgagccaagagatccaccgtcaaaag	chromosome_2	30	20734	20763	tgcgtgagccaagagatccaccgtcaaaag	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_19_x39468	30	1	30	tgcgtgagccaagagatccaccgtcaaaag	chromosome_2	30	42321	42350	tgcgtgagccaagagatccaccgtcaaaag	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_19_x39468	30	1	30	tgcgtgagccaagagatccaccgtcaaaag	chromosome_2	30	53070	53099	tgcgtgagccaagagatccaccgtcaaaag	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_19_x39468	30	1	30	tgcgtgagccaagagatccaccgtcaaaag	chromosome_2	30	31478	31507	tgcgtgagccaagagatccaccgtcaaaag	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_20_x31213	30	1	30	tccgacactcagacagacatgctcctggga	chromosome_2	30	52946	52975	tccgacactcagacagacatgctcctggga	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_20_x31213	30	1	30	tccgacactcagacagacatgctcctggga	chromosome_2	30	31354	31383	tccgacactcagacagacatgctcctggga	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

cut -f1 AST-1412_reads_collapsed_vs_genome.arf | sort | uniq | wc -l 
599077

conda deactivate
```

Run mirdeep2. `nano AST-1412_mirdeep2.sh`

```
#!/bin/bash -i
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/smRNA/mirdeep2/AST-1412           
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

#module load Miniconda3/4.9.2
conda activate /data/putnamlab/mirdeep2

echo "Starting mirdeep2 on sample AST-1412 trimmed to 30bp" $(date)

miRDeep2.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar/17_sed.collapse.cat.AST-1412.fastq /data/putnamlab/jillashey/Astrangia_Genome/apoculata.assembly.scaffolds_chromosome_level.fasta /data/putnamlab/jillashey/Astrangia2021/smRNA/mirdeep2/AST-1412/AST-1412_reads_collapsed_vs_genome.arf /data/putnamlab/jillashey/Astrangia2021/smRNA/refs/mature_mirbase_cnidarian_T.fa none none -t N.vectensis -P -v -g -1 2>report.log

echo "mirdeep2 concluded for sample AST-1412 trimmed to 30bp" $(date)

conda deactivate
```

Submitted batch job 293150

#### AST-1560

```
cd /data/putnamlab/jillashey/Astrangia2021/smRNA/mirdeep2/AST-1560
```

Run mapping step 

```
conda activate /data/putnamlab/mirdeep2

mapper.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar/17_sed.collapse.cat.AST-1560.fastq -c -p /data/putnamlab/jillashey/Astrangia2021/smRNA/scripts/Apoc_ref.btindex -s AST-1560_reads_collapsed.fa -t AST-1560_reads_collapsed_vs_genome.arf -v  

discarding short reads
mapping reads to genome index
trimming unmapped nts in the 3' ends
Log file for this run is in mapper_logs and called mapper.log_58565
Mapping statistics

#desc	total	mapped	unmapped	%mapped	%unmapped
total: 35654048	767007	34887041	2.151	97.849
seq: 35654048	767007	34887041	2.151	97.849
```

Look at the mapping results

```
head AST-1560_reads_collapsed.fa 
>seq_1_x404539
GACTTTGTAGCATAGGTAAGGTTAGTGCAT
>seq_2_x272137
TCAGATGCACTAACCTTACCTATGCTACAA
>seq_3_x78722
GGAAGAGCACACGTCTGAACTCCAGTCACG
>seq_4_x58838
TCGGACTGTAGAACTCTGAACGTGTAGATC
>seq_5_x46651
GCACTGATGGTTCAGTGGTAGAATTCTCGC

zgrep -c ">" AST-1560_reads_collapsed.fa 
11702120

head AST-1560_reads_collapsed_vs_genome.arf 
seq_22_x19193	30	1	30	caggtctgtgatgcccttagatgtccgggg	chromosome_2	30	32082	32111	caggtctgtgatgcccttagatgtacaggg	-	2	mmmmmmmmmmmmmmmmmmmmmmmmMmMmmm
seq_22_x19193	30	1	30	caggtctgtgatgcccttagatgtccgggg	chromosome_2	30	42925	42954	caggtctgtgatgcccttagatgtacaggg	-	2	mmmmmmmmmmmmmmmmmmmmmmmmMmMmmm
seq_22_x19193	30	1	30	caggtctgtgatgcccttagatgtccgggg	chromosome_2	30	21338	21367	caggtctgtgatgcccttagatgtacaggg	-	2	mmmmmmmmmmmmmmmmmmmmmmmmMmMmmm
seq_43_x14818	30	1	30	acaggtctgtgatgcccttagatgtccggg	chromosome_2	30	32083	32112	acaggtctgtgatgcccttagatgtacagg	-	2	mmmmmmmmmmmmmmmmmmmmmmmmmMmMmm
seq_43_x14818	30	1	30	acaggtctgtgatgcccttagatgtccggg	chromosome_2	30	42926	42955	acaggtctgtgatgcccttagatgtacagg	-	2	mmmmmmmmmmmmmmmmmmmmmmmmmMmMmm
seq_43_x14818	30	1	30	acaggtctgtgatgcccttagatgtccggg	chromosome_2	30	21339	21368	acaggtctgtgatgcccttagatgtacagg	-	2	mmmmmmmmmmmmmmmmmmmmmmmmmMmMmm
seq_72_x11510	30	1	30	aggtctgtgatgcccttagatgtccggggc	chromosome_2	30	32081	32110	aggtctgtgatgcccttagatgtacagggc	-	2	mmmmmmmmmmmmmmmmmmmmmmmMmMmmmm
seq_72_x11510	30	1	30	aggtctgtgatgcccttagatgtccggggc	chromosome_2	30	42924	42953	aggtctgtgatgcccttagatgtacagggc	-	2	mmmmmmmmmmmmmmmmmmmmmmmMmMmmmm
seq_72_x11510	30	1	30	aggtctgtgatgcccttagatgtccggggc	chromosome_2	30	21337	21366	aggtctgtgatgcccttagatgtacagggc	-	2	mmmmmmmmmmmmmmmmmmmmmmmMmMmmmm
seq_87_x10194	30	1	30	agccaggaatcctaaccgctagaccatctg	chromosome_12	30	41900395	41900424	agccaggaatcctaaccgctagaccatttg	-	1	mmmmmmmmmmmmmmmmmmmmmmmmmmmMmm

cut -f1 AST-1560_reads_collapsed_vs_genome.arf | sort | uniq | wc -l 
123911

conda deactivate
```

Run mirdeep2. `nano AST-1560_mirdeep2.sh`

```
#!/bin/bash -i
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/smRNA/mirdeep2/AST-1560          
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

#module load Miniconda3/4.9.2
conda activate /data/putnamlab/mirdeep2

echo "Starting mirdeep2 on sample AST-1560 trimmed to 30bp" $(date)

miRDeep2.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar/17_sed.collapse.cat.AST-1560.fastq /data/putnamlab/jillashey/Astrangia_Genome/apoculata.assembly.scaffolds_chromosome_level.fasta /data/putnamlab/jillashey/Astrangia2021/smRNA/mirdeep2/AST-1560/AST-1560_reads_collapsed_vs_genome.arf /data/putnamlab/jillashey/Astrangia2021/smRNA/refs/mature_mirbase_cnidarian_T.fa none none -t N.vectensis -P -v -g -1 2>report.log

echo "mirdeep2 concluded for sample AST-1560 trimmed to 30bp" $(date)

conda deactivate
```

Submitted batch job 293152

#### AST-1567

```
cd /data/putnamlab/jillashey/Astrangia2021/smRNA/mirdeep2/AST-1567
```

Run mapping step 

```
conda activate /data/putnamlab/mirdeep2

mapper.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar/17_sed.collapse.cat.AST-1567.fastq -c -p /data/putnamlab/jillashey/Astrangia2021/smRNA/scripts/Apoc_ref.btindex -s AST-1567_reads_collapsed.fa -t AST-1567_reads_collapsed_vs_genome.arf -v  

discarding short reads
mapping reads to genome index
trimming unmapped nts in the 3' ends
Log file for this run is in mapper_logs and called mapper.log_12237
Mapping statistics

#desc	total	mapped	unmapped	%mapped	%unmapped
total: 33222794	3149689	30073105	9.481	90.519
seq: 33222794	3149689	30073105	9.481	90.519
```

Look at the mapping results

```
head AST-1567_reads_collapsed.fa 
>seq_1_x330204
GCACTGGTGGTTCAGTGGTAGAATTCTCGC
>seq_2_x226911
GCACTGATGGTTCAGTGGTAGAATTCTCGC
>seq_3_x157338
GGCGAGAATTCTACCACTGAACCACCAGTG
>seq_4_x98181
AGCGAGAATTCTACCACTGAACCACCAGTG
>seq_5_x77909
GGAAGAGCACACGTCTGAACTCCAGTCACA

zgrep -c ">" AST-1567_reads_collapsed.fa 
9871465

head AST-1567_reads_collapsed_vs_genome.arf 
seq_9_x59542	30	1	30	aacttttgacggtggatctcttggctcacg	chromosome_2	30	53072	53101	aacttttgacggtggatctcttggctcacg	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_9_x59542	30	1	30	aacttttgacggtggatctcttggctcacg	chromosome_2	30	20736	20765	aacttttgacggtggatctcttggctcacg	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_9_x59542	30	1	30	aacttttgacggtggatctcttggctcacg	chromosome_2	30	31480	31509	aacttttgacggtggatctcttggctcacg	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_9_x59542	30	1	30	aacttttgacggtggatctcttggctcacg	chromosome_2	30	42323	42352	aacttttgacggtggatctcttggctcacg	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_24_x25786	30	1	30	tgcgtgagccaagagatccaccgtcaaaag	chromosome_2	30	20734	20763	tgcgtgagccaagagatccaccgtcaaaag	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_24_x25786	30	1	30	tgcgtgagccaagagatccaccgtcaaaag	chromosome_2	30	42321	42350	tgcgtgagccaagagatccaccgtcaaaag	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_24_x25786	30	1	30	tgcgtgagccaagagatccaccgtcaaaag	chromosome_2	30	53070	53099	tgcgtgagccaagagatccaccgtcaaaag	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_24_x25786	30	1	30	tgcgtgagccaagagatccaccgtcaaaag	chromosome_2	30	31478	31507	tgcgtgagccaagagatccaccgtcaaaag	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_37_x16700	29	1	29	aatggataaccctcaaccgtccggacctc	chromosome_14	29	6295369	6295397	aatggataaccctcaaccgtccggacctc	-	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmm
seq_54_x12863	30	1	30	tccgacactcagacagacatgctcctggga	chromosome_2	30	31354	31383	tccgacactcagacagacatgctcctggga	+	0	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

cut -f1 AST-1567_reads_collapsed_vs_genome.arf | sort | uniq | wc -l 
569922

conda deactivate
```

Run mirdeep2. `nano AST-1567_mirdeep2.sh`

```
#!/bin/bash -i
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/smRNA/mirdeep2/AST-1567        
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

#module load Miniconda3/4.9.2
conda activate /data/putnamlab/mirdeep2

echo "Starting mirdeep2 on sample AST-1567 trimmed to 30bp" $(date)

miRDeep2.pl /data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/flexbar/17_sed.collapse.cat.AST-1567.fastq /data/putnamlab/jillashey/Astrangia_Genome/apoculata.assembly.scaffolds_chromosome_level.fasta /data/putnamlab/jillashey/Astrangia2021/smRNA/mirdeep2/AST-1567/AST-1567_reads_collapsed_vs_genome.arf /data/putnamlab/jillashey/Astrangia2021/smRNA/refs/mature_mirbase_cnidarian_T.fa none none -t N.vectensis -P -v -g -1 2>report.log

echo "mirdeep2 concluded for sample AST-1567 trimmed to 30bp" $(date)

conda deactivate
```

Submitted batch job 293153

I'm going to call it quits for today. 3 samples are running through mirdeep2 now, and 2 samples are pending on the server. 







good resource for miranda 
- https://bioinformaticsworkbook.org/dataAnalysis/SmallRNA/Miranda_miRNA_Target_Prediction.html#gsc.tab=0 

Other potential miRNA target prediction programs: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7839038/#:~:text=TargetScan%20primarily%20predicts%20potential%20miRNA,to%20include%20only%20conserved%20sites. 