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

### 20241027 

Even though the QC looks horrible, I want to run short stack on this data to see what it looks like. In the scripts folder: `nano shortstack.sh`

```
#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

echo "Running short stack on mature trimmed miRNAs (R1) from Mcap DT project"

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim

# Load modules 
module load ShortStack/4.0.2-foss-2022a  
module load Kent_tools/442-GCC-11.3.0

# Run short stack
ShortStack \
--genomefile /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.assembly.fasta \
--readfile trim.fastp.31bp.10_small_RNA_S4_R1_001.fastq.gz \
trim.fastp.31bp.11_small_RNA_S5_R1_001.fastq.gz \
trim.fastp.31bp.13_S76_R1_001.fastq.gz \
trim.fastp.31bp.14_small_RNA_S6_R1_001.fastq.gz \
trim.fastp.31bp.23_S77_R1_001.fastq.gz \
trim.fastp.31bp.24_small_RNA_S7_R1_001.fastq.gz \
trim.fastp.31bp.26_small_RNA_S8_R1_001.fastq.gz \
trim.fastp.31bp.28_small_RNA_S9_R1_001.fastq.gz \
trim.fastp.31bp.35_S78_R1_001.fastq.gz \
trim.fastp.31bp.36_small_RNA_S10_R1_001.fastq.gz \
trim.fastp.31bp.37_small_RNA_S11_R1_001.fastq.gz \
trim.fastp.31bp.39_small_RNA_S12_R1_001.fastq.gz \
trim.fastp.31bp.47_small_RNA_S13_R1_001.fastq.gz \
trim.fastp.31bp.48_small_RNA_S14_R1_001.fastq.gz \
trim.fastp.31bp.51_small_RNA_S15_R1_001.fastq.gz \
trim.fastp.31bp.52_S79_R1_001.fastq.gz \
trim.fastp.31bp.60_S80_R1_001.fastq.gz \
trim.fastp.31bp.61_small_RNA_S16_R1_001.fastq.gz \
trim.fastp.31bp.62_small_RNA_S17_R1_001.fastq.gz \
trim.fastp.31bp.63_small_RNA_S18_R1_001.fastq.gz \
trim.fastp.31bp.6_small_RNA_S1_R1_001.fastq.gz \
trim.fastp.31bp.72_S81_R1_001.fastq.gz \
trim.fastp.31bp.73_small_RNA_S19_R1_001.fastq.gz \
trim.fastp.31bp.74_small_RNA_S20_R1_001.fastq.gz \
trim.fastp.31bp.75_small_RNA_S21_R1_001.fastq.gz \
trim.fastp.31bp.7_small_RNA_S2_R1_001.fastq.gz \
trim.fastp.31bp.85_S82_R1_001.fastq.gz \
trim.fastp.31bp.86_small_RNA_S22_R1_001.fastq.gz \
trim.fastp.31bp.87_small_RNA_S23_R1_001.fastq.gz \
trim.fastp.31bp.88_small_RNA_S24_R1_001.fastq.gz \
trim.fastp.31bp.8_small_RNA_S3_R1_001.fastq.gz \
trim.fastp.31bp.9_S75_R1_001.fastq.gz \
--known_miRNAs /data/putnamlab/jillashey/Astrangia2021/smRNA/refs/mature_mirbase_cnidarian_T.fa \
--outdir /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/shortstack \
--threads 10 \
--dn_mirna

echo "Short stack complete!"
```

Submitted batch job 345205. 

### 20241028

Shortstack finished running overnight but doesn't seem to have worked well. Getting this at the end of the error file:

```
Traceback (most recent call last):
  File "/opt/software/ShortStack/4.0.2-foss-2022a/ShortStack", line 3592, in <module>
    mir_qdata = mirna(args, merged_bam, fai, pmir_bedfile, read_count)
  File "/opt/software/ShortStack/4.0.2-foss-2022a/ShortStack", line 2259, in mirna
    dn_q_bedlines, dn_locus_bedlines = zip(*denovo_mloci1)
ValueError: not enough values to unpack (expected 2, got 0)
```

Also getting this in the error file for the samples: 

```
# reads processed: 18875806
# reads with at least one alignment: 999 (0.01%)
# reads that failed to align: 18874807 (99.99%)
Reported 2699 alignments
[bam_sort_core] merging from 0 files and 10 in-memory blocks...
# reads processed: 11804203
# reads with at least one alignment: 1216 (0.01%)
# reads that failed to align: 11802987 (99.99%)
Reported 2981 alignments
[bam_sort_core] merging from 0 files and 10 in-memory blocks...
# reads processed: 24558419
# reads with at least one alignment: 50 (0.00%)
# reads that failed to align: 24558369 (100.00%)
Reported 803 alignments
```

Vast majority of reads are not aligning which makes sense because I gave it very long reads (>150bp). This is probably why it failed. No results files or miRNA fasta files were generated either. Need to revisit trimming

### 20241105

Revisiting trimming with cutadapt. In the scripts folder, edit `cutadapt_trim.sh`

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
#module load cutadapt/3.5-GCCcore-11.2.0
module load FastQC/0.11.8-Java-1.8

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/

#echo "Count number of reads in each file" $(date)

#zgrep -c "@LH00" *.gz > smRNA_raw_read_count.txt

echo "Start read trimming with cutadapt, followed by QC" $(date)

# Make an array of sequences to trim
array1=($(ls *R1_001.fastq.gz))

# cutadapt loop
for i in ${array1[@]}; do
    cutadapt -a TGGAATTCTCGGGTGCCAAGG -u 3 -m 15 -q 20,20 --trim-n --max-n 0 -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim/trim_cutadapt_${i} $i
done

echo "Read trimming of adapters complete" $(date)
```

Submitted batch job 347763. Took almost 4 hours. I didn't put the QC code in here because I wanted to run it in interactive mode. 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim
interactive
module load FastQC/0.11.8-Java-1.8
fastqc *
```

Started running in interactive mode at 2:15pm

```
mv *qc ../../output/fastqc/trim
cd ../../output/fastqc/trim
module load MultiQC/1.9-intel-2020a-Python-3.8.2
multiqc *
```

Looks very similar to cutadapt attempt 1. The samples that were sequenced in the first batch were able to be successfully trimmed to 25 bp. Going to rerun cutadapt with some more stringent cutoffs and more adapter seqs that were identified in the individual fastqc outputs. In the scripts folder: `nano cutadapt_trim_stringent.sh`

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
#module load cutadapt/3.5-GCCcore-11.2.0
module load FastQC/0.11.8-Java-1.8

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/

#echo "Count number of reads in each file" $(date)

#zgrep -c "@LH00" *.gz > smRNA_raw_read_count.txt

echo "Start read trimming with cutadapt, followed by QC" $(date)

# Make an array of sequences to trim
array1=($(ls *R1_001.fastq.gz))

# cutadapt loop
for i in ${array1[@]}; do
cutadapt -a TGGAATTCTCGGGTGCCAAGG -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a GATCGGAAGAGCACACGTCTGAACTCCAGTCA -a GAAACGTTGGGTTGCGGTATTGGAAGAGCACACGTCTGAACTCCAGTCAC -a AGTGTCCTATCAGCTACCATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a CTCGGCATGGAAGCCGTGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a TAGCGTCGAACGGGCGCAATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a ACGTCTGAACTCCAGTCACTTACAGGAATCTGGGGGGGGGGGGGGGGGGG -a TATCGGTGAAACATCCTCATCGGAAGAGCACACGTCTGAACTCCAGTCAC -m 18 -M 30 -q 30,30 --trim-n --max-n 0 --discard-untrimmed -e 0.1 -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_${i} $i
done

## added adapters that were overrepresented sequences from the individual fastqc htmls 

echo "Read trimming of adapters complete" $(date)

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/
fastqc * -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/fastqc/trim_stringent
done

echo "Read trimming of adapters and QC complete. Starting multiqc" $(date)

module unload cutadapt/4.2-GCCcore-11.3.0
module load MultiQC/1.9-intel-2020a-Python-3.8.2

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/fastqc/trim_stringent
multiqc * 

echo "MultiQC complete" $(date)

echo "Count number of reads in each trimmed file" $(date)

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent
zgrep -c "@LH00" *.gz > smRNA_trim_stringent_read_count.txt
```

Submitted batch job 347786. Based on initial looks, a lot of the reads are getting tossed because they are too short. 

### 20241107

Took a while but ran. The lengths look much better, as does the GC content. There are some samples that have such a high duplication rate though...I know this is a result of PCR artifacts. I may run shortstack to see what the output looks like. I may need to do some adapter trimming as well with certain samples. 

next steps: 

- run short stack on stringent trimmed reads 
- run picard (picard/2.25.1-Java-11)
	- EstimateLibraryComplexity - Estimates the numbers of unique molecules in a sequencing library.  
	- ExtractIlluminaBarcodes - Tool determines the barcode for each read in an Illumina lane.  
	- MarkDuplicates - Identifies duplicate reads.  
USAGE: PicardCommandLine <program name> [-h]

Run short stack. 

In the scripts folder: `nano shortstack_trim_stringent.sh`

```
#!/bin/bash
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

echo "Running short stack on mature trimmed miRNAs from Mcap DT project"

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent

# Load modules 
module load ShortStack/4.0.2-foss-2022a  
module load Kent_tools/442-GCC-11.3.0

# Run short stack
ShortStack \
--genomefile /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.assembly.fasta \
--readfile trim_stringent_cutadapt_10_small_RNA_S4_R1_001.fastq.gz \
trim_stringent_cutadapt_11_small_RNA_S5_R1_001.fastq.gz \
trim_stringent_cutadapt_13_S76_R1_001.fastq.gz \
trim_stringent_cutadapt_14_small_RNA_S6_R1_001.fastq.gz \
trim_stringent_cutadapt_23_S77_R1_001.fastq.gz \
trim_stringent_cutadapt_24_small_RNA_S7_R1_001.fastq.gz \
trim_stringent_cutadapt_26_small_RNA_S8_R1_001.fastq.gz \
trim_stringent_cutadapt_28_small_RNA_S9_R1_001.fastq.gz \
trim_stringent_cutadapt_35_S78_R1_001.fastq.gz \
trim_stringent_cutadapt_36_small_RNA_S10_R1_001.fastq.gz \
trim_stringent_cutadapt_37_small_RNA_S11_R1_001.fastq.gz \
trim_stringent_cutadapt_39_small_RNA_S12_R1_001.fastq.gz \
trim_stringent_cutadapt_47_small_RNA_S13_R1_001.fastq.gz \
trim_stringent_cutadapt_48_small_RNA_S14_R1_001.fastq.gz \
trim_stringent_cutadapt_51_small_RNA_S15_R1_001.fastq.gz \
trim_stringent_cutadapt_52_S79_R1_001.fastq.gz \
trim_stringent_cutadapt_60_S80_R1_001.fastq.gz \
trim_stringent_cutadapt_61_small_RNA_S16_R1_001.fastq.gz \
trim_stringent_cutadapt_62_small_RNA_S17_R1_001.fastq.gz \
trim_stringent_cutadapt_63_small_RNA_S18_R1_001.fastq.gz \
#trim_stringent_cutadapt_6_small_RNA_S1_R1_001.fastq.gz \
#trim_stringent_cutadapt_72_S81_R1_001.fastq.gz \
#trim_stringent_cutadapt_73_small_RNA_S19_R1_001.fastq.gz \
#trim_stringent_cutadapt_74_small_RNA_S20_R1_001.fastq.gz \
#trim_stringent_cutadapt_75_small_RNA_S21_R1_001.fastq.gz \
#trim_stringent_cutadapt_7_small_RNA_S2_R1_001.fastq.gz \
#trim_stringent_cutadapt_85_S82_R1_001.fastq.gz \
#trim_stringent_cutadapt_86_small_RNA_S22_R1_001.fastq.gz\
#trim_stringent_cutadapt_87_small_RNA_S23_R1_001.fastq.gz \
#trim_stringent_cutadapt_88_small_RNA_S24_R1_001.fastq.gz \
#trim_stringent_cutadapt_8_small_RNA_S3_R1_001.fastq.gz \
#trim_stringent_cutadapt_9_S75_R1_001.fastq.gz \
--known_miRNAs /data/putnamlab/jillashey/Astrangia2021/smRNA/refs/mature_mirbase_cnidarian_T.fa \
--outdir /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/shortstack_trim_stringent \
--threads 10 \
--dn_mirna

echo "Short stack complete!"
```

Submitted batch job 348221 - only doing a subset for now







Zoe cutadapt trimming: https://github.com/zdellaert/LaserCoral/blob/1b776313512822d368e957591890b2228c590bd5/code/RNA-seq-bioinf.md 






code to extract fasta info from shortstack: https://github.com/urol-e5/deep-dive/blob/main/F-Pmea/code/13.2.1.1-Pmea-sRNAseq-ShortStack-FastA-extraction.md 

