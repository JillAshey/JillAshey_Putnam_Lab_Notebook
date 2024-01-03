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

Submitted batch job 291969. Took about 5 hours. The QC plots don't look amazing and the length is still at 100 bp. I'm going to rerun but adding the argument `--length_limit 30` for fastp. This means that reads longer than 30 bp will be discarded. Submitted batch job 291997

















I'm going to try to run [mirDeep2](https://github.com/rajewsky-lab/mirdeep2) using code from the mirdeep2 github [tutorial](https://github.com/rajewsky-lab/mirdeep2/blob/master/TUTORIAL.md) and Sam White's [code](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/code/11-Peve-sRNAseq-miRdeep2.md#7-run-mirdeep2) from the E5 deep dive project. 

When I asked Kevin Bryan to install mirdeep2 and miranda on the server, he said: 

"Hi Jill, I’m out this week, but you should be able to install both of these into a conda environment. You can use Miniconda3/23.5.2-0 for that.

https://anaconda.org/bioconda/miranda
https://anaconda.org/bioconda/mirdeep2

Let me know if you run into any issues. You probably want to use --prefix /data/putnamlab/miranda or similar so it can be shared with the group."

I'll first try this and then if that doesn't work, I'll email him again to install. Use one of these commands to install: 

```
conda install -c bioconda mirdeep2
conda install -c "bioconda/label/cf201901" mirdeep2
```

Failed to install properly. When installing in an conda environment, it said I didn't have access to install this software. Emailed Kevin Bryan again to see if he can install them on the HPC. 









To run mirdeep2, I need the following inputs: 
- Collapsed reads (concatenated, unique reads)
- Genome fasta 
- [miRBase](https://mirbase.org/download/) mature miRNA fasta 

First, I need to concatenate (with `cat` command) and collapse my reads (with `fastx_collapser` from the [fastx toolkit](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html). I'm going to try with just one sample for now. In the `/data/putnamlab/jillashey/Astrangia2021/smRNA/data/trim/fastp`:

```
cat trimmed.AST-1065_R1_001.fastq.gz trimmed.AST-1065_R2_001.fastq.gz > cat.trimmed.AST-1065.fastq

module load FASTX-Toolkit/0.0.14-GCC-9.3.0 

# gunzip cat.trimmed.AST-1065.fastq.gz # files must be unzipped for collapsing; unzip if needed
fastx_collapser -v -i cat.trimmed.AST-1065.fastq -o collapse.cat.trimmed.AST-1065.fastq
```

There are 3 primary modules for mirdeep2: 

- miRDeep2 module
	- Input: reference genome, sequencing reads, positions of reads mapped against genome
	- Test format of input files so that any errors can be identified and corrected for
	- Potential miRNA precursors (ie pre-miRNAs) are excised from genome using mapped reads as guidelines
		- Takes the sequence that the mapped read covers and some additional sequence
	- Map the reads against the excised miRNA potential precursors using Bowtie
	- RNAfold tool is used to predict if the RNA secondary structures of each excised potential precursors resembles a typical miRNA hairpin structure
		- Structure must resemble a hairpin and the read must fall within the hairpin as would be expected from Dicer processing
		- Based on that info^^, the potential precursor is assigned a score that reflects how likely it is to be an actual miRNA
	- Output: info on every miRNA identified in the data (known or novel?)
- Mapper module
	- Input: reference genome, sequencing reads, adapter sequences (from library prep)
	- Reads are trimmed and adapters removed
	- Reads are mapped to genome with Bowtie
	- Output: file with processed reads, file with reads mapped against genome
- Quantifier module
	- Input: sequencing reads, known and precursor miRNAs from reference species
	- Reads and miRNA strands are mapped separately against precursors
	- Mappings of reads and miRNA strands are intersected —> reads that map to the same position as a given strand add to the read count of that miRNA strand
	- Output: file with read counts of all known miRNAs in data

Based on what I've read about mirdeep2, the mapper module should be run first and then mirdeep2 and quantifier modules. I may need to edit the file names before I run any mirdeep2 stuff because the documentation says: "The readID must end with _xNumber and is not allowed to contain whitespaces. has to have the format name_uniqueNumber_xnumber"

First, index the genome with bowtie (NOT bowtie2). In the scripts folder: `nano bowtie_build.sh`

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

Submitted batch job 291990. Getting some GCC conflict errors so added `module load GCCcore/11.3.0` to the script. Submitted batch job 291991


Run the mapper module 

```
mapper.pl reads.fa -c -j -p Apoc_ref.btindex -s COLLAPSED_READS.fa -t reads_collapsed_vs_genome.arf -v
```

Run the quantifier module 

```
quantifier.pl -m MIRBASE.MATURE -r COLLAPSED_READS.fa  
```

Run the mirdeep2 module 


