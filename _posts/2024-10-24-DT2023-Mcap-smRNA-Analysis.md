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

### 20241108

Took a while but ran. The lengths look much better, as does the GC content. There are some samples that have such a high duplication rate though...I know this is a result of PCR artifacts. I may run shortstack to see what the output looks like. I may need to do some adapter trimming as well with certain samples. 

next steps: 

- run short stack on stringent trimmed reads 
- run picard (picard/2.25.1-Java-11)
	- EstimateLibraryComplexity - Estimates the numbers of unique molecules in a sequencing library.  
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
# Run short stack
ShortStack \
--genomefile /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.assembly.fasta \
--readfile trim_stringent_cutadapt_10_small_RNA_S4_R1_001.fastq.gz \
trim_stringent_cutadapt_11_small_RNA_S5_R1_001.fastq.gz \
trim_stringent_cutadapt_13_S76_R1_001.fastq.gz \
trim_stringent_cutadapt_63_small_RNA_S18_R1_001.fastq.gz \
--known_miRNAs /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/mature_mirbase_cnidarian_T.fa \
--outdir /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/shortstack_trim_stringent \
--threads 10 \
--dn_mirna

echo "Short stack complete!"
```

Submitted batch job 348221 - only doing a subset for now. Ran but ended with this error: `/var/spool/slurmd/job348221/slurm_script: line 56: --known_miRNAs: command not found`, even though that is a command...Maybe copy `mature_mirbase_cnidarian_T.fa` to this folder? 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data
cp /data/putnamlab/jillashey/Astrangia2021/smRNA/refs/mature_mirbase_cnidarian_T.fa .
```

Edit the `shortstack_trim_stringent.sh` script so that the known miRNAs path is the new path (`/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/mature_mirbase_cnidarian_T.fa`). Submitted batch job 348256. Got the same error. Trying on a subset and making sure that the backslashes are correct. Submitted batch job 348269

Let's try running picard in interactive mode

```
interactive 
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/shortstack_trim_stringent/ShortStack_1731009886
module load picard/2.25.1-Java-11
To execute picard run: java -jar $EBROOTPICARD/picard.jar
```

I just want to test picard on one of the bam files (`trim_stringent_cutadapt_51_small_RNA_S15_R1_001.bam`). 

```
# Run EstimateLibraryComplexity
java -jar $EBROOTPICARD/picard.jar EstimateLibraryComplexity \
    I=trim_stringent_cutadapt_24_small_RNA_S7_R1_001.bam \
    O=trim_stringent_cutadapt_24_small_RNA_complexity_metrics.txt
    
INFO	2024-11-08 09:24:40	EstimateLibraryComplexity	

********** NOTE: Picard's command line syntax is changing.
**********
********** For more information, please see:
********** https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition)
**********
********** The command line looks like this in the new syntax:
**********
**********    EstimateLibraryComplexity -I trim_stringent_cutadapt_51_small_RNA_S15_R1_001.bam -O trim_stringent_cutadapt_51_small_RNA_complexity_metrics.txt
**********


09:24:42.296 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/software/picard/2.25.1-Java-11/picard.jar!/com/intel/gkl/native/libgkl_compression.so
[Fri Nov 08 09:24:42 EST 2024] EstimateLibraryComplexity INPUT=[trim_stringent_cutadapt_51_small_RNA_S15_R1_001.bam] OUTPUT=trim_stringent_cutadapt_51_small_RNA_complexity_metrics.txt    MIN_IDENTICAL_BASES=5 MAX_DIFF_RATE=0.03 MIN_MEAN_QUALITY=20 MAX_GROUP_RATIO=500 MAX_READ_LENGTH=0 MIN_GROUP_COUNT=2 READ_NAME_REGEX=<optimized capture of last three ':' separated fields as numeric values> OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 MAX_OPTICAL_DUPLICATE_SET_SIZE=300000 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=33053659 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false
[Fri Nov 08 09:24:42 EST 2024] Executing as jillashey@n063.cluster.com on Linux 3.10.0-1160.119.1.el7.x86_64 amd64; OpenJDK 64-Bit Server VM 11.0.2+9; Deflater: Intel; Inflater: Intel; Provider GCS is not available; Picard version: 2.25.1
INFO	2024-11-08 09:24:42	EstimateLibraryComplexity	Will store 33053659 read pairs in memory before sorting.
INFO	2024-11-08 09:24:45	EstimateLibraryComplexity	Finished reading - read 0 records - moving on to scanning for duplicates.
[Fri Nov 08 09:24:45 EST 2024] picard.sam.markduplicates.EstimateLibraryComplexity done. Elapsed time: 0.05 minutes.
Runtime.totalMemory()=2035417088
```

didn't really seem to run successfully...I got an output file but no information in the output file. Let's try running MarkDuplicates

```
java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I trim_stringent_cutadapt_51_small_RNA_S15_R1_001.bam -O marked_dups_trim_stringent_cutadapt_51_small_RNA_S15_R1_001.bam -M marked_dup_51_small_RNA_metrics.txt

09:32:36.015 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/opt/software/picard/2.25.1-Java-11/picard.jar!/com/intel/gkl/native/libgkl_compression.so
[Fri Nov 08 09:32:36 EST 2024] MarkDuplicates INPUT=[trim_stringent_cutadapt_51_small_RNA_S15_R1_001.bam] OUTPUT=marked_dups_trim_stringent_cutadapt_51_small_RNA_S15_R1_001.bam METRICS_FILE=marked_dup_51_small_RNA_metrics.txt    MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 SORTING_COLLECTION_SIZE_RATIO=0.25 TAG_DUPLICATE_SET_MEMBERS=false REMOVE_SEQUENCING_DUPLICATES=false TAGGING_POLICY=DontTag CLEAR_DT=true DUPLEX_UMI=false ADD_PG_TAG_TO_READS=true REMOVE_DUPLICATES=false ASSUME_SORTED=false DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_NAME=MarkDuplicates READ_NAME_REGEX=<optimized capture of last three ':' separated fields as numeric values> OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 MAX_OPTICAL_DUPLICATE_SET_SIZE=300000 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false
[Fri Nov 08 09:32:36 EST 2024] Executing as jillashey@n063.cluster.com on Linux 3.10.0-1160.119.1.el7.x86_64 amd64; OpenJDK 64-Bit Server VM 11.0.2+9; Deflater: Intel; Inflater: Intel; Provider GCS is not available; Picard version: 2.25.1
INFO	2024-11-08 09:32:36	MarkDuplicates	Start of doWork freeMemory: 2018640208; totalMemory: 2035417088; maxMemory: 31136546816
INFO	2024-11-08 09:32:36	MarkDuplicates	Reading input file and constructing read end information.
INFO	2024-11-08 09:32:36	MarkDuplicates	Will retain up to 112813575 data points before spilling to disk.
INFO	2024-11-08 09:32:37	MarkDuplicates	Read 102803 records. 0 pairs never matched.
INFO	2024-11-08 09:32:38	MarkDuplicates	After buildSortedReadEndLists freeMemory: 1285863816; totalMemory: 2214113280; maxMemory: 31136546816
INFO	2024-11-08 09:32:38	MarkDuplicates	Will retain up to 973017088 duplicate indices before spilling to disk.
Killed
```

Failed bleh. Also no miRNA loci found in short stack...may need to adjust trimming parameters. 

### 20241112

Need to retrim but first going to run short stack with the samples from the initial sequencing run. These have lower duplication rates and I do not think that they have as much PCR duplication. In the scripts folder: `nano shortstack_trim_stringent_sub.sh`

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
# Run short stack
ShortStack \
--genomefile /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.assembly.fasta \
--readfile trim_stringent_cutadapt_13_S76_R1_001.fastq.gz \
trim_stringent_cutadapt_23_S77_R1_001.fastq.gz \
trim_stringent_cutadapt_35_S78_R1_001.fastq.gz \
trim_stringent_cutadapt_52_S79_R1_001.fastq.gz \
trim_stringent_cutadapt_60_S80_R1_001.fastq.gz \
trim_stringent_cutadapt_72_S81_R1_001.fastq.gz \
trim_stringent_cutadapt_85_S82_R1_001.fastq.gz \
trim_stringent_cutadapt_9_S75_R1_001.fastq.gz \
--known_miRNAs /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/mature_mirbase_cnidarian_T.fa \
--threads 10 \
--dn_mirna

echo "Short stack complete!"
```

Submitted batch job 348560. Ran in a few hours. Output said it found 3 miRNA loci but not sure if this includes known miRNAs...

```
awk '$3 == "mature_miRNA"' Results.gff3 > filtered_mature_miRNA.gff
head filtered_mature_miRNA.gff
Montipora_capitata_HIv3___Scaffold_2    ShortStack      mature_miRNA    46381100        46381121        476     +       .       ID=Cluster_2999.mature;Parent=Cluster_2999
Montipora_capitata_HIv3___Scaffold_8    ShortStack      mature_miRNA    22810894        22810915        204     +       .       ID=Cluster_12334.mature;Parent=Cluster_12334
Montipora_capitata_HIv3___Scaffold_14   ShortStack      mature_miRNA    2986538 2986560 1       +       .       ID=Cluster_27569.mature;Parent=Cluster_27569

awk -F'\t' '$20 == "Y"' Results.txt > Results_miRNA.txt
head Results_miRNA.txt
Montipora_capitata_HIv3___Scaffold_2:46381078-46381169  Cluster_2999    Montipora_capitata_HIv3___Scaffold_2    46381078        46381169        92      1012    59      1.0     +       ACAAUGUUUCGGCUUGUUCCCG  476     108     8       169     705     17      5       22      Y       spi-mir-temp-14_Stylophora_pistillata_Liew_et_al._2014_NA
Montipora_capitata_HIv3___Scaffold_8:22810843-22810935  Cluster_12334   Montipora_capitata_HIv3___Scaffold_8    22810843        22810935        93      4483    95      1.0     +       AUGCAGCGGAAAUCGAACCUGGG 3053    154     15      43      256     3224    791     23      Y       NA
Montipora_capitata_HIv3___Scaffold_14:2986487-2986580   Cluster_27569   Montipora_capitata_HIv3___Scaffold_14   2986487 2986580 94      49      7       1.0     +       ACGGUGAAAGUCGUCUCAAUAUUCA       35      2       35      0       1       1       10      N       Y       spi-mir-temp-30_Stylophora_pistillata_Liew_et_al._2014_Close_match_of_nve-miR-2036.;eca-nve-F-miR-2036_Edwardsiella_carnea_Praher_et_al._2021_Transcriptome-level;Adi-Mir-2036_3p_Acropora_digitifera_Gajigan_&_Conaco_2017_nve-miR-2036-3p;_nve-miR-2036-3p;_spi-miR-temp-30;ami-nve-F-miR-2036-3p_Acropora_millepora_Praher_et_al._2021_NA;spi-nve-F-miR-2036_Stylophora_pistillata_Praher_et_al._2021_NA
```

Clearly didn't pick up all miRNAs. Maybe trying running mirdeep2?? I did a lot of trouble shooting with AST and mirdeep2 (see post [here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/76d256ef61ada1c96b99fd09d07bfca3df781bd2/_posts/2023-06-05-Astrangia2021-smRNA-Analysis.md)). First, I need to map the samples to the genome with the `mapper.pl` script. Genome must first be indexed by bowtie-build (NOT bowtie2). Because I have multiple samples, I need to make a config file that contains fq file locations and a unique 3 letter code (see mirdeep2 [documentation](https://www.mdc-berlin.de/content/mirdeep2-documentation?mdcbl%5B0%5D=%2Fn-rajewsky%23t-data%2Csoftware%26resources&mdcbv=71nDTh7VzOJOW6SFGuFySs4mus4wnovu-t2LZzV2dL8&mdcot=6&mdcou=20738&mdctl=0)). I will also be unzipping the files so remove the .gz from file name. 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent

nano config.txt

/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_10_small_RNA_S4_R1_001.fastq s10
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_11_small_RNA_S5_R1_001.fastq s11
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_13_S76_R1_001.fastq s13
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_14_small_RNA_S6_R1_001.fastq s14
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_23_S77_R1_001.fastq s23
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_24_small_RNA_S7_R1_001.fastq s24
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_26_small_RNA_S8_R1_001.fastq s26
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_28_small_RNA_S9_R1_001.fastq s28
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_35_S78_R1_001.fastq s35
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_36_small_RNA_S10_R1_001.fastq s36
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_37_small_RNA_S11_R1_001.fastq s37
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_39_small_RNA_S12_R1_001.fastq s39
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_47_small_RNA_S13_R1_001.fastq s47
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_48_small_RNA_S14_R1_001.fastq s48
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_51_small_RNA_S15_R1_001.fastq s51
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_52_S79_R1_001.fastq s52
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_60_S80_R1_001.fastq s60
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_61_small_RNA_S16_R1_001.fastq s61
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_62_small_RNA_S17_R1_001.fastq s62
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_63_small_RNA_S18_R1_001.fastq s63
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_6_small_RNA_S1_R1_001.fastq s06
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_72_S81_R1_001.fastq s72
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_73_small_RNA_S19_R1_001.fastq s73
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_74_small_RNA_S20_R1_001.fastq s74
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_75_small_RNA_S21_R1_001.fastq s75
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_7_small_RNA_S2_R1_001.fastq s07
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_85_S82_R1_001.fastq s85
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_86_small_RNA_S22_R1_001.fastq s86
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_87_small_RNA_S23_R1_001.fastq s87
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_88_small_RNA_S24_R1_001.fastq s88
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_8_small_RNA_S3_R1_001.fastq s08
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_9_S75_R1_001.fastq s09
```

In the `config.txt` file, I have the path and file name, as well as the unique 3 letter code (s followed by sample number).

In the scripts folder: `nano mapper_mirdeep2.sh`

```
#!/bin/bash -i
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

module load GCCcore/11.3.0 #I needed to add this to resolve conflicts between loaded GCCcore/9.3.0 and GCCcore/11.3.0
module load Bowtie/1.3.1-GCC-11.3.0

echo "Index Mcap genome" $(date)

# Index the reference genome for Mcap 
bowtie-build /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.assembly.fasta /data/putnamlab/jillashey/genome/Mcap/V3/Mcap_ref.btindex

echo "Referece genome indexed!" $(date)
echo "Unload unneeded packages and run mapper script for trimmed stringent reads" $(date)

module unload module load GCCcore/11.3.0 
module unload Bowtie/1.3.1-GCC-11.3.0

conda activate /data/putnamlab/mirdeep2

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent

gunzip *.fastq.gz

mapper.pl config.txt -e -d -h -j -l 18 -m -p /data/putnamlab/jillashey/genome/Mcap/V3/Mcap_ref.btindex -s mapped_reads.fa -t mapped_reads_vs_genome.arf

echo "Mapping complete for trimmed stringent reads" $(date)

conda deactivate 
```

Submitted batch job 348612. The flags mean the following: 

- `-e`- use fastq files as input
- `-d` - input is config file 
- `-h` - convert to fasta format 
- `-i` - convert RNA to DNA alphabet 
- `-j` - remove entries with non-canonical letters
- `-l` - discard reads shorter than 18 nt
- `-m` - collapse reads 
- `-p` - map to genome--genome has to be indexed by bowtie build 
- `-s` - mapped reads fasta output
- `-t` - read mappings v genome output

Also consider adding -q which allows mapping with one mismatch in seed but mapping takes longer. Ran in about an hour. Look at output

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent

zgrep -c ">" mapped_reads.fa
21804779

wc -l mapped_reads_vs_genome.arf
8098667 mapped_reads_vs_genome.arf

less bowtie.log 
# reads processed: 2215013
# reads with at least one reported alignment: 416126 (18.79%)
# reads that failed to align: 1500959 (67.76%)
# reads with alignments suppressed due to -m: 297928 (13.45%)
Reported 841100 alignments to 1 output stream(s)
```

I can now run mirdeep2 but this is a problem for tomorrow!

### 20241113

Short stack and mirdeep2 output has been in this folder: `/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent`. Going to move to output folder. 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output
mkdir mirdeep2_trim_stringent

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent
mv ShortStack_1731* ../../output/shortstack_trim_stringent/
mv mapp* ../../output/mirdeep2_trim_stringent/
mv bowtie.log ../../output/mirdeep2_trim_stringent/
mv config.txt ../../output/mirdeep2_trim_stringent/
```

The output files from the `mapper.pl` portion of mirdeep2 will be used as input for the `mirdeep2.pl` portion of the pipeline. In the scripts folder: `nano mirna_predict_mirdeep2.sh`

```
#!/bin/bash -i
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

#module load Miniconda3/4.9.2
conda activate /data/putnamlab/mirdeep2

echo "Starting mirdeep2 with reads that were trimmed stringently via cutadapt" $(date)

miRDeep2.pl /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_trim_stringent/mapped_reads.fa /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.assembly.fasta /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_trim_stringent/mapped_reads_vs_genome.arf /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/mature_mirbase_cnidarian_T.fa none none -P -v -g -1 2>report.log

echo "mirdeep2 complete" $(date)

conda deactivate
```

Submitted batch job 348663. In the past, mirdeep2 has taken a few days to process my AST samples so I am guessing it will take a similar amount of time. The flags mean the following: 

- `mapped_reads.fa` - fasta file with sequencing reads 
- `Montipora_capitata_HIv3.assembly.fasta` - genome fasta 
- `mapped_reads_vs_genome.arf` - mapped reads to the genome in miRDeep2 arf format
- `mature_mirbase_cnidarian_T.fa` - known miRNAs 
- `none` - no fasta file for precursor sequences. I accidently added two `none` flags in there, hopefully that wont mess things up too much. 
- `-P` - report novel miRNAs 
- `-v` - output verbose mode 
- `-g -1` - report all potential miRNAs regardless of low score
- `2>report.log` - redirects errors to log file 

### 20241118

Took about 2 days to run and most of output was put in scripts folder. Downloaded output to my computer. 

### 20241121

I need to revisit trimming. I downloaded `slurm-347786.out` (most recent trimming iteration output) to my computer to look at sample adapter information. 

```
# M10
Total reads processed:              19,773,964
Reads with adapters:                12,283,878 (62.1%)

== Read fate breakdown ==
Reads that were too short:           3,481,028 (17.6%)
Reads that were too long:            7,532,588 (38.1%)
Reads with too many N:                  33,125 (0.2%)
Reads discarded as untrimmed:           25,736 (0.1%)
Reads written (passing filters):     8,701,487 (44.0%)

Total basepairs processed: 2,985,868,564 bp
Quality-trimmed:              49,791,977 bp (1.7%)
Total written (filtered):    175,293,612 bp (5.9%)

# M11
Total reads processed:              12,310,413
Reads with adapters:                10,586,911 (86.0%)

== Read fate breakdown ==
Reads that were too short:           6,691,286 (54.4%)
Reads that were too long:            1,817,067 (14.8%)
Reads with too many N:                  14,264 (0.1%)
Reads discarded as untrimmed:           13,200 (0.1%)
Reads written (passing filters):     3,774,596 (30.7%)

Total basepairs processed: 1,858,872,363 bp
Quality-trimmed:              32,849,546 bp (1.8%)
Total written (filtered):     79,147,991 bp (4.3%)

# M13
Total reads processed:              25,732,779
Reads with adapters:                25,613,472 (99.5%)

== Read fate breakdown ==
Reads that were too short:           7,832,390 (30.4%)
Reads that were too long:            4,570,545 (17.8%)
Reads with too many N:                   6,929 (0.0%)
Reads discarded as untrimmed:           48,788 (0.2%)
Reads written (passing filters):    13,274,127 (51.6%)

Total basepairs processed: 3,885,649,629 bp
Quality-trimmed:             133,678,654 bp (3.4%)
Total written (filtered):    328,216,706 bp (8.4%)

# M14
Total reads processed:              35,754,439
Reads with adapters:                29,105,761 (81.4%)

== Read fate breakdown ==
Reads that were too short:          26,385,621 (73.8%)
Reads that were too long:            6,927,638 (19.4%)
Reads with too many N:                   8,413 (0.0%)
Reads discarded as untrimmed:           15,538 (0.0%)
Reads written (passing filters):     2,417,229 (6.8%)

Total basepairs processed: 5,398,920,289 bp
Quality-trimmed:              90,439,586 bp (1.7%)
Total written (filtered):     53,521,904 bp (1.0%)

# M23
Total reads processed:              28,608,340
Reads with adapters:                28,488,320 (99.6%)

== Read fate breakdown ==
Reads that were too short:           8,300,653 (29.0%)
Reads that were too long:            5,340,020 (18.7%)
Reads with too many N:                   8,149 (0.0%)
Reads discarded as untrimmed:           50,412 (0.2%)
Reads written (passing filters):    14,909,106 (52.1%)

Total basepairs processed: 4,319,859,340 bp
Quality-trimmed:             121,050,542 bp (2.8%)
Total written (filtered):    381,785,695 bp (8.8%)

# M24 
Total reads processed:              16,185,558
Reads with adapters:                15,051,209 (93.0%)

== Read fate breakdown ==
Reads that were too short:           1,238,159 (7.6%)
Reads that were too long:            1,298,021 (8.0%)
Reads with too many N:                  53,492 (0.3%)
Reads discarded as untrimmed:           17,232 (0.1%)
Reads written (passing filters):    13,578,654 (83.9%)

Total basepairs processed: 2,444,019,258 bp
Quality-trimmed:              51,942,286 bp (2.1%)
Total written (filtered):    258,626,182 bp (10.6%)

# M26
Total reads processed:              45,463,982
Reads with adapters:                33,200,737 (73.0%)

== Read fate breakdown ==
Reads that were too short:          26,954,608 (59.3%)
Reads that were too long:           12,621,446 (27.8%)
Reads with too many N:                  21,951 (0.0%)
Reads discarded as untrimmed:           29,590 (0.1%)
Reads written (passing filters):     5,836,387 (12.8%)

Total basepairs processed: 6,865,061,282 bp
Quality-trimmed:             114,030,951 bp (1.7%)
Total written (filtered):    124,668,281 bp (1.8%)

# M28 
Total reads processed:              29,330,724
Reads with adapters:                22,529,784 (76.8%)

== Read fate breakdown ==
Reads that were too short:          20,934,300 (71.4%)
Reads that were too long:            6,864,087 (23.4%)
Reads with too many N:                   5,150 (0.0%)
Reads discarded as untrimmed:           36,179 (0.1%)
Reads written (passing filters):     1,491,008 (5.1%)

Total basepairs processed: 4,428,939,324 bp
Quality-trimmed:              57,180,891 bp (1.3%)
Total written (filtered):     33,212,070 bp (0.7%)

# M35
Total reads processed:              27,415,125
Reads with adapters:                27,297,487 (99.6%)

== Read fate breakdown ==
Reads that were too short:           8,746,075 (31.9%)
Reads that were too long:            4,714,853 (17.2%)
Reads with too many N:                   7,301 (0.0%)
Reads discarded as untrimmed:           46,476 (0.2%)
Reads written (passing filters):    13,900,420 (50.7%)

Total basepairs processed: 4,139,683,875 bp
Quality-trimmed:             135,071,303 bp (3.3%)
Total written (filtered):    350,393,120 bp (8.5%)

# M36
Total reads processed:              21,669,140
Reads with adapters:                15,229,442 (70.3%)

== Read fate breakdown ==
Reads that were too short:           3,530,235 (16.3%)
Reads that were too long:            6,397,943 (29.5%)
Reads with too many N:                  45,695 (0.2%)
Reads discarded as untrimmed:           37,563 (0.2%)
Reads written (passing filters):    11,657,704 (53.8%)

Total basepairs processed: 3,272,040,140 bp
Quality-trimmed:              51,150,713 bp (1.6%)
Total written (filtered):    228,309,042 bp (7.0%)

# M37
Total reads processed:              20,247,832
Reads with adapters:                19,802,633 (97.8%)

== Read fate breakdown ==
Reads that were too short:           7,675,729 (37.9%)
Reads that were too long:              434,204 (2.1%)
Reads with too many N:                  45,858 (0.2%)
Reads discarded as untrimmed:           20,730 (0.1%)
Reads written (passing filters):    12,071,311 (59.6%)

Total basepairs processed: 3,057,422,632 bp
Quality-trimmed:              62,466,338 bp (2.0%)
Total written (filtered):    246,713,865 bp (8.1%)

# M39
Total reads processed:              20,298,000
Reads with adapters:                14,625,467 (72.1%)

== Read fate breakdown ==
Reads that were too short:          13,581,604 (66.9%)
Reads that were too long:            5,768,966 (28.4%)
Reads with too many N:                   3,465 (0.0%)
Reads discarded as untrimmed:           11,546 (0.1%)
Reads written (passing filters):       932,419 (4.6%)

Total basepairs processed: 3,064,998,000 bp
Quality-trimmed:              41,859,848 bp (1.4%)
Total written (filtered):     19,320,488 bp (0.6%)

# M47
Total reads processed:              61,683,862
Reads with adapters:                56,145,283 (91.0%)

== Read fate breakdown ==
Reads that were too short:          13,015,507 (21.1%)
Reads that were too long:            5,696,266 (9.2%)
Reads with too many N:                 166,609 (0.3%)
Reads discarded as untrimmed:          153,126 (0.2%)
Reads written (passing filters):    42,652,354 (69.1%)

Total basepairs processed: 9,314,263,162 bp
Quality-trimmed:             181,161,945 bp (1.9%)
Total written (filtered):    846,330,853 bp (9.1%)

# M48
Total reads processed:              26,207,354
Reads with adapters:                 5,677,260 (21.7%)

== Read fate breakdown ==
Reads that were too short:           3,896,981 (14.9%)
Reads that were too long:           20,434,620 (78.0%)
Reads with too many N:                   6,817 (0.0%)
Reads discarded as untrimmed:           90,086 (0.3%)
Reads written (passing filters):     1,778,850 (6.8%)

Total basepairs processed: 3,957,310,454 bp
Quality-trimmed:              46,499,881 bp (1.2%)
Total written (filtered):     35,591,311 bp (0.9%)

# M51
Total reads processed:              20,348,304
Reads with adapters:                 7,710,258 (37.9%)

== Read fate breakdown ==
Reads that were too short:           7,042,815 (34.6%)
Reads that were too long:           12,673,695 (62.3%)
Reads with too many N:                   2,273 (0.0%)
Reads discarded as untrimmed:           14,418 (0.1%)
Reads written (passing filters):       615,103 (3.0%)

Total basepairs processed: 3,072,593,904 bp
Quality-trimmed:              48,509,512 bp (1.6%)
Total written (filtered):     12,632,989 bp (0.4%)

# M52
Total reads processed:              29,279,123
Reads with adapters:                29,137,374 (99.5%)

== Read fate breakdown ==
Reads that were too short:           7,298,239 (24.9%)
Reads that were too long:            6,020,859 (20.6%)
Reads with too many N:                   8,657 (0.0%)
Reads discarded as untrimmed:           61,778 (0.2%)
Reads written (passing filters):    15,889,590 (54.3%)

Total basepairs processed: 4,421,147,573 bp
Quality-trimmed:             151,159,079 bp (3.4%)
Total written (filtered):    408,145,425 bp (9.2%)

# M60
Total reads processed:              24,332,450
Reads with adapters:                24,230,741 (99.6%)

== Read fate breakdown ==
Reads that were too short:           9,734,570 (40.0%)
Reads that were too long:            3,070,229 (12.6%)
Reads with too many N:                   5,938 (0.0%)
Reads discarded as untrimmed:           38,908 (0.2%)
Reads written (passing filters):    11,482,805 (47.2%)

Total basepairs processed: 3,674,199,950 bp
Quality-trimmed:             119,459,545 bp (3.3%)
Total written (filtered):    282,122,042 bp (7.7%)

# M61
Total reads processed:              24,731,279
Reads with adapters:                 6,861,730 (27.7%)

== Read fate breakdown ==
Reads that were too short:           6,034,816 (24.4%)
Reads that were too long:           17,805,944 (72.0%)
Reads with too many N:                   2,980 (0.0%)
Reads discarded as untrimmed:           76,252 (0.3%)
Reads written (passing filters):       811,287 (3.3%)

Total basepairs processed: 3,734,423,129 bp
Quality-trimmed:              46,085,583 bp (1.2%)
Total written (filtered):     16,211,139 bp (0.4%)

# M62
Total reads processed:              15,145,513
Reads with adapters:                 7,633,247 (50.4%)

== Read fate breakdown ==
Reads that were too short:           4,771,918 (31.5%)
Reads that were too long:            7,508,877 (49.6%)
Reads with too many N:                  11,032 (0.1%)
Reads discarded as untrimmed:           45,692 (0.3%)
Reads written (passing filters):     2,807,994 (18.5%)

Total basepairs processed: 2,286,972,463 bp
Quality-trimmed:              37,090,659 bp (1.6%)
Total written (filtered):     55,944,873 bp (2.4%)

# M63
Total reads processed:              21,903,406
Reads with adapters:                10,379,795 (47.4%)

== Read fate breakdown ==
Reads that were too short:           9,069,045 (41.4%)
Reads that were too long:           11,612,371 (53.0%)
Reads with too many N:                   4,485 (0.0%)
Reads discarded as untrimmed:           28,940 (0.1%)
Reads written (passing filters):     1,188,565 (5.4%)

Total basepairs processed: 3,307,414,306 bp
Quality-trimmed:              46,947,183 bp (1.4%)
Total written (filtered):     25,709,935 bp (0.8%)

# M6 
Total reads processed:              24,211,196
Reads with adapters:                15,778,017 (65.2%)

== Read fate breakdown ==
Reads that were too short:          15,305,643 (63.2%)
Reads that were too long:            8,428,574 (34.8%)
Reads with too many N:                   1,536 (0.0%)
Reads discarded as untrimmed:           43,373 (0.2%)
Reads written (passing filters):       432,070 (1.8%)

Total basepairs processed: 3,655,890,596 bp
Quality-trimmed:              52,333,701 bp (1.4%)
Total written (filtered):      9,537,955 bp (0.3%)

# M72 
Total reads processed:              25,266,689
Reads with adapters:                25,144,539 (99.5%)

== Read fate breakdown ==
Reads that were too short:           7,803,487 (30.9%)
Reads that were too long:            3,726,653 (14.7%)
Reads with too many N:                   7,388 (0.0%)
Reads discarded as untrimmed:           51,139 (0.2%)
Reads written (passing filters):    13,678,022 (54.1%)

Total basepairs processed: 3,815,270,039 bp
Quality-trimmed:             137,376,345 bp (3.6%)
Total written (filtered):    328,694,194 bp (8.6%)

# M73 
Total reads processed:               5,419,554
Reads with adapters:                 5,395,120 (99.5%)

== Read fate breakdown ==
Reads that were too short:           3,459,671 (63.8%)
Reads that were too long:              320,201 (5.9%)
Reads with too many N:                   6,668 (0.1%)
Reads discarded as untrimmed:            4,293 (0.1%)
Reads written (passing filters):     1,628,721 (30.1%)

Total basepairs processed:   818,352,654 bp
Quality-trimmed:              17,216,446 bp (2.1%)
Total written (filtered):     32,665,596 bp (4.0%)

# M74
Total reads processed:              23,599,752
Reads with adapters:                 6,346,294 (26.9%)

== Read fate breakdown ==
Reads that were too short:           5,447,500 (23.1%)
Reads that were too long:           17,208,033 (72.9%)
Reads with too many N:                   3,295 (0.0%)
Reads discarded as untrimmed:           55,067 (0.2%)
Reads written (passing filters):       885,857 (3.8%)

Total basepairs processed: 3,563,562,552 bp
Quality-trimmed:              44,030,663 bp (1.2%)
Total written (filtered):     18,225,152 bp (0.5%)

# M75
Total reads processed:              39,936,022
Reads with adapters:                31,308,251 (78.4%)

== Read fate breakdown ==
Reads that were too short:          26,334,128 (65.9%)
Reads that were too long:            9,283,704 (23.2%)
Reads with too many N:                  15,318 (0.0%)
Reads discarded as untrimmed:           53,526 (0.1%)
Reads written (passing filters):     4,249,346 (10.6%)

Total basepairs processed: 6,030,339,322 bp
Quality-trimmed:             104,351,331 bp (1.7%)
Total written (filtered):     94,582,130 bp (1.6%)

# M7
Total reads processed:              30,275,553
Reads with adapters:                24,226,251 (80.0%)

== Read fate breakdown ==
Reads that were too short:          11,063,587 (36.5%)
Reads that were too long:            6,232,765 (20.6%)
Reads with too many N:                  48,605 (0.2%)
Reads discarded as untrimmed:           51,715 (0.2%)
Reads written (passing filters):    12,878,881 (42.5%)

Total basepairs processed: 4,571,608,503 bp
Quality-trimmed:              88,504,444 bp (1.9%)
Total written (filtered):    255,714,395 bp (5.6%)

# M85
Total reads processed:              21,835,402
Reads with adapters:                21,711,666 (99.4%)

== Read fate breakdown ==
Reads that were too short:           7,627,502 (34.9%)
Reads that were too long:            2,616,361 (12.0%)
Reads with too many N:                   5,995 (0.0%)
Reads discarded as untrimmed:           48,112 (0.2%)
Reads written (passing filters):    11,537,432 (52.8%)

Total basepairs processed: 3,297,145,702 bp
Quality-trimmed:             123,423,025 bp (3.7%)
Total written (filtered):    268,088,833 bp (8.1%)

# M86
Total reads processed:              17,514,990
Reads with adapters:                 8,429,216 (48.1%)

== Read fate breakdown ==
Reads that were too short:           5,771,993 (33.0%)
Reads that were too long:            9,103,226 (52.0%)
Reads with too many N:                  10,039 (0.1%)
Reads discarded as untrimmed:           38,292 (0.2%)
Reads written (passing filters):     2,591,440 (14.8%)

Total basepairs processed: 2,644,763,490 bp
Quality-trimmed:              46,364,967 bp (1.8%)
Total written (filtered):     52,301,229 bp (2.0%)

# M87
Total reads processed:              35,917,496
Reads with adapters:                30,802,070 (85.8%)

== Read fate breakdown ==
Reads that were too short:          27,183,568 (75.7%)
Reads that were too long:            5,601,464 (15.6%)
Reads with too many N:                  11,049 (0.0%)
Reads discarded as untrimmed:           64,118 (0.2%)
Reads written (passing filters):     3,057,297 (8.5%)

Total basepairs processed: 5,423,541,896 bp
Quality-trimmed:              82,622,227 bp (1.5%)
Total written (filtered):     68,971,401 bp (1.3%)

# M88
Total reads processed:              28,792,588
Reads with adapters:                11,393,165 (39.6%)

== Read fate breakdown ==
Reads that were too short:           9,948,192 (34.6%)
Reads that were too long:           17,523,034 (60.9%)
Reads with too many N:                   4,898 (0.0%)
Reads discarded as untrimmed:           30,002 (0.1%)
Reads written (passing filters):     1,286,462 (4.5%)

Total basepairs processed: 4,347,680,788 bp
Quality-trimmed:              58,097,687 bp (1.3%)
Total written (filtered):     27,541,768 bp (0.6%)

# M8
Total reads processed:              37,732,448
Reads with adapters:                21,548,435 (57.1%)

== Read fate breakdown ==
Reads that were too short:          19,322,853 (51.2%)
Reads that were too long:           16,300,593 (43.2%)
Reads with too many N:                   7,780 (0.0%)
Reads discarded as untrimmed:           49,950 (0.1%)
Reads written (passing filters):     2,051,272 (5.4%)

Total basepairs processed: 5,697,599,648 bp
Quality-trimmed:              82,384,975 bp (1.4%)
Total written (filtered):     45,156,438 bp (0.8%)

# M9
Total reads processed:              25,911,154
Reads with adapters:                25,800,281 (99.6%)

== Read fate breakdown ==
Reads that were too short:           8,632,667 (33.3%)
Reads that were too long:            4,601,815 (17.8%)
Reads with too many N:                   6,587 (0.0%)
Reads discarded as untrimmed:           45,543 (0.2%)
Reads written (passing filters):    12,624,542 (48.7%)

Total basepairs processed: 3,912,584,254 bp
Quality-trimmed:             153,459,563 bp (3.9%)
Total written (filtered):    305,187,572 bp (7.8%)
```

A lot of samples had a high percentage of reads that were too short for the 18 bp cutoff. There were also some samples that did not have many reads with adapters, which means I need to go back to the fastqc files to identify adapters/sequences that were overrepresented in these samples and add them to the cut adapt code. By doing this, I will capture more reads with adapters for trimming. 

Next trimming iteration: decrease min length to 15 bp, look at fastqc for more adapter sequences 

### 20241203

I am going to look at the outputs from the [individual sample fastQC](https://github.com/JillAshey/DevelopmentalTimeseries/tree/main/data/Molecular/smRNA/trim) to see what adapters I am missing. I will focus primarily on the samples that did not have high adapter content in my trimmed stringent iteration. I am hoping there is no limit to number of adapters. 

Most of the samples I did in the initial sequencing batch have no adaptere hits for overrepresented sequences. Using the individual fastQC files, I pulled all the adapters that had >1% of possible contamination and got over 400 unique adapter sequences lol. Instead, I am going to go through the individual fastQCs again while looking at the cutadapt output to see which samples did not have many reads with adapters. 

Adapters to include: 

```
TGGAATTCTCGGGTGCCAAGG 
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA 
GATCGGAAGAGCACACGTCTGAACTCCAGTCA 
GAAACGTTGGGTTGCGGTATTGGAAGAGCACACGTCTGAACTCCAGTCAC 
AGTGTCCTATCAGCTACCATCGGAAGAGCACACGTCTGAACTCCAGTCAC 
CTCGGCATGGAAGCCGTGATCGGAAGAGCACACGTCTGAACTCCAGTCAC 
TAGCGTCGAACGGGCGCAATCGGAAGAGCACACGTCTGAACTCCAGTCAC 
ACGTCTGAACTCCAGTCACTTACAGGAATCTGGGGGGGGGGGGGGGGGGG 
TATCGGTGAAACATCCTCATCGGAAGAGCACACGTCTGAACTCCAGTCAC 
ACCGTTGTTCGAGCAGGAATCGGAAGAGCACACGTCTGAACTCCAGTCAC # M10
TCGGAAGAGCACACGTCTGAACTCCAGTCACTACCGAGGATCTGGGGGGG # M11
TCGGAAGAGCACACGTCTGAACTCCAGTCACTCGTAGTGATCTGGGGGGG # M14
TCGGAAGAGCACACGTCTGAACTCCAGTCACCGTTAGAAATCTGGGGGGG # M26
TCGGAAGAGCACACGTCTGAACTCCAGTCACCTACGACAATCTGGGGGGG # M28
CTCGGCATGGAAGCCGTGATCGGAAGAGCACACGTCTGAACTCCAGTCAC # M36
AAGAGCACACGTCTGAACTCCAGTCACTAAGTGGTATCTGGGGGGGGGGG # M39
ACGTCTGAACTCCAGTCACCGGACAACATCTGGGGGGGGGGGGGGGGGGG # M48
TCGGAAGAGCACACGTCTGAACTCCAGTCACATATGGATATCTGGGGGGG # M51
TCGGAAGAGCACACGTCTGAACTCCAGTCACGCGCAAGCATCTAGGGGGG # M61
GAAGAGCACACGTCTGAACTCCAGTCACCCGTGAAGATCTGGGGGGGGGG # M62
TCGGAAGAGCACACGTCTGAACTCCAGTCACAAGATACTATCTGGGGGGG # M63
TCGGAAGAGCACACGTCTGAACTCCAGTCACATGAGGCCATCTGGGGGGG # M6
ACGTCTGAACTCCAGTCACGGAGCGTCATCTGGGGGGGGGGGGGGGGGGG # M74
TCGGAAGAGCACACGTCTGAACTCCAGTCACATGGCATGATCTGGGGGGG # M75
CTCGGCATGGAAGCCGTGATCGGAAGAGCACACGTCTGAACTCCAGTCAC # M7
CTCGGCATGGAAGCCGTGATCGGAAGAGCACACGTCTGAACTCCAGTCAC # M86
TCGGAAGAGCACACGTCTGAACTCCAGTCACGCAATGCAATCTGGGGGGG # M87
TCGGAAGAGCACACGTCTGAACTCCAGTCACGTTCCAATATCTAGGGGGG # M88
TCGGAAGAGCACACGTCTGAACTCCAGTCACGATTCTGCATCTGGGGGGG # M8
```

Make new folders for output

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data
mkdir trim_stringent_take2
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/fastqc
mkdir trim_stringent_take2
```

My last cutadapt run (job 347786) took about a day to run. I am going to set the time higher this time since I am adding so many more adapters. In the scripts folder: `nano cutadapt_trim_stringent_take2.sh`

```
#!/bin/bash
#SBATCH -t 100:00:00
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
cutadapt -b TGGAATTCTCGGGTGCCAAGG -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -b GATCGGAAGAGCACACGTCTGAACTCCAGTCA -b GAAACGTTGGGTTGCGGTATTGGAAGAGCACACGTCTGAACTCCAGTCAC -b AGTGTCCTATCAGCTACCATCGGAAGAGCACACGTCTGAACTCCAGTCAC -b CTCGGCATGGAAGCCGTGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -b TAGCGTCGAACGGGCGCAATCGGAAGAGCACACGTCTGAACTCCAGTCAC -b ACGTCTGAACTCCAGTCACTTACAGGAATCTGGGGGGGGGGGGGGGGGGG -b TATCGGTGAAACATCCTCATCGGAAGAGCACACGTCTGAACTCCAGTCAC -b ACCGTTGTTCGAGCAGGAATCGGAAGAGCACACGTCTGAACTCCAGTCAC -b TCGGAAGAGCACACGTCTGAACTCCAGTCACTACCGAGGATCTGGGGGGG -b TCGGAAGAGCACACGTCTGAACTCCAGTCACTCGTAGTGATCTGGGGGGG -b TCGGAAGAGCACACGTCTGAACTCCAGTCACCGTTAGAAATCTGGGGGGG -b TCGGAAGAGCACACGTCTGAACTCCAGTCACCTACGACAATCTGGGGGGG -b CTCGGCATGGAAGCCGTGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -b AAGAGCACACGTCTGAACTCCAGTCACTAAGTGGTATCTGGGGGGGGGGG -b ACGTCTGAACTCCAGTCACCGGACAACATCTGGGGGGGGGGGGGGGGGGG -b TCGGAAGAGCACACGTCTGAACTCCAGTCACATATGGATATCTGGGGGGG -b TCGGAAGAGCACACGTCTGAACTCCAGTCACGCGCAAGCATCTAGGGGGG -b GAAGAGCACACGTCTGAACTCCAGTCACCCGTGAAGATCTGGGGGGGGGG -b TCGGAAGAGCACACGTCTGAACTCCAGTCACAAGATACTATCTGGGGGGG -b TCGGAAGAGCACACGTCTGAACTCCAGTCACATGAGGCCATCTGGGGGGG -b ACGTCTGAACTCCAGTCACGGAGCGTCATCTGGGGGGGGGGGGGGGGGGG -b TCGGAAGAGCACACGTCTGAACTCCAGTCACATGGCATGATCTGGGGGGG -b CTCGGCATGGAAGCCGTGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -b CTCGGCATGGAAGCCGTGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -b TCGGAAGAGCACACGTCTGAACTCCAGTCACGCAATGCAATCTGGGGGGG -b TCGGAAGAGCACACGTCTGAACTCCAGTCACGTTCCAATATCTAGGGGGG -b TCGGAAGAGCACACGTCTGAACTCCAGTCACGATTCTGCATCTGGGGGGG -m 15 -M 35 -q 30,30 --trim-n --max-n 0 --discard-untrimmed -e 0.1 -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent_take2/trim_stringent_take2_cutadapt_${i} $i
done

## added even more adapters that were overrepresented sequences from the individual fastqc htmls 

echo "Read trimming of adapters complete" $(date)

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent_take2/
fastqc * -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/fastqc/trim_stringent_take2
done

echo "Read trimming of adapters and QC complete. Starting multiqc" $(date)

module unload cutadapt/4.2-GCCcore-11.3.0
module load MultiQC/1.9-intel-2020a-Python-3.8.2

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/fastqc/trim_stringent_take2
multiqc * 

echo "MultiQC complete" $(date)

echo "Count number of reads in each trimmed file" $(date)

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent_take2
zgrep -c "@LH00" *.gz > smRNA_trim_stringent_take2_read_count.txt
```
 
Submitted batch job 352143. In this iteration, I used cutadapt with -b: For this type of adapter, the sequence is specified with -b ADAPTER (or use the longer spelling --anywhere ADAPTER). The adapter may appear in the beginning (even degraded), within the read, or at the end of the read (even partially). The decision which part of the read to remove is made as follows: If there is at least one base before the found adapter, then the adapter is considered to be a 3’ adapter and the adapter itself and everything following it is removed. Otherwise, the adapter is considered to be a 5’ adapter and it is removed from the read, but the sequence after it remains.

discussion w. zoe
- do they all have adapters? 
- blast overrepresented sequences 
- reverse complement of sequence 
- what would data look like if i dont throw out anything >30bp
- sequences of index primers 

Bleh trimming is hard and confusing. 

### 20241206

While trimming takes place, I am going to run the quantifier module to obtain counts for the miRNAs. From the mirdeep2 output, I identified putative novel and known miRNAs (see [code](https://github.com/JillAshey/DevelopmentalTimeseries/blob/main/scripts/miRNA_discovery.Rmd)). For the novel miRNAs, I filtered the list so that mirdeep2 score >10, no rfam info, at least 10 reads in mature read count, and significant randfold p-value; 52 putative novel miRNAs were identified. For the known miRNAs, I filtered the list so that mirdeep2 score >0, no rfam info, at least 10 reads in mature read count, and significant ranfold p-value; 7 putative known miRNAs were identified. 59 miRNAs total! 

```
# novel 
Montipora_capitata_HIv3___Scaffold_2_98529
Montipora_capitata_HIv3___Scaffold_8_488597
Montipora_capitata_HIv3___Scaffold_4_264216
Montipora_capitata_HIv3___Scaffold_9_587940
Montipora_capitata_HIv3___Scaffold_9_625013
Montipora_capitata_HIv3___Scaffold_2_71488
Montipora_capitata_HIv3___Scaffold_9_613239
Montipora_capitata_HIv3___Scaffold_11_759461
Montipora_capitata_HIv3___Scaffold_12_859123
Montipora_capitata_HIv3___Scaffold_10_631026
Montipora_capitata_HIv3___Scaffold_8_585867
Montipora_capitata_HIv3___Scaffold_2_124999
Montipora_capitata_HIv3___Scaffold_8_578139
Montipora_capitata_HIv3___Scaffold_11_826703
Montipora_capitata_HIv3___Scaffold_10_729473
Montipora_capitata_HIv3___Scaffold_10_667068
Montipora_capitata_HIv3___Scaffold_151_1075760
Montipora_capitata_HIv3___Scaffold_14_1017829
Montipora_capitata_HIv3___Scaffold_5_309429
Montipora_capitata_HIv3___Scaffold_8_508612
Montipora_capitata_HIv3___Scaffold_8_553085
Montipora_capitata_HIv3___Scaffold_7_456843
Montipora_capitata_HIv3___Scaffold_4_291272
Montipora_capitata_HIv3___Scaffold_6_400996
Montipora_capitata_HIv3___Scaffold_13_954675
Montipora_capitata_HIv3___Scaffold_8_542727
Montipora_capitata_HIv3___Scaffold_14_999529
Montipora_capitata_HIv3___Scaffold_14_994729
Montipora_capitata_HIv3___Scaffold_9_617967
Montipora_capitata_HIv3___Scaffold_5_322659
Montipora_capitata_HIv3___Scaffold_6_419669
Montipora_capitata_HIv3___Scaffold_10_645863
Montipora_capitata_HIv3___Scaffold_2_87665
Montipora_capitata_HIv3___Scaffold_1_35393
Montipora_capitata_HIv3___Scaffold_4_232369
Montipora_capitata_HIv3___Scaffold_5_375083
Montipora_capitata_HIv3___Scaffold_2_64557
Montipora_capitata_HIv3___Scaffold_3_149822
Montipora_capitata_HIv3___Scaffold_9_606053
Montipora_capitata_HIv3___Scaffold_11_801837
Montipora_capitata_HIv3___Scaffold_8_508613
Montipora_capitata_HIv3___Scaffold_8_535535
Montipora_capitata_HIv3___Scaffold_11_796076
Montipora_capitata_HIv3___Scaffold_11_787579
Montipora_capitata_HIv3___Scaffold_11_816162
Montipora_capitata_HIv3___Scaffold_13_949160
Montipora_capitata_HIv3___Scaffold_11_841599
Montipora_capitata_HIv3___Scaffold_5_310446
Montipora_capitata_HIv3___Scaffold_6_439094
Montipora_capitata_HIv3___Scaffold_5_365006
Montipora_capitata_HIv3___Scaffold_2_67225
Montipora_capitata_HIv3___Scaffold_5_307262

# known 
Montipora_capitata_HIv3___Scaffold_14_1018693
Montipora_capitata_HIv3___Scaffold_2_96328
Montipora_capitata_HIv3___Scaffold_8_565039
Montipora_capitata_HIv3___Scaffold_8_586427
Montipora_capitata_HIv3___Scaffold_14_992284
Montipora_capitata_HIv3___Scaffold_12_910953
Montipora_capitata_HIv3___Scaffold_1_4183
```

I now want to do 2 things - quantify the miRNAs and use miranda to assess miRNA binding to 3'UTR. Let's quantify first. 

All output files from mirdeep2 were put in the scripts folder. Move them to the proper output folder. 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/scripts
mv dir_prepare_signature1731506059/ ../output/mirdeep2_trim_stringent/
mv result_13_11_2024_t_08_48_44.* ../output/mirdeep2_trim_stringent/
mv mirna_results_13_11_2024_t_08_48_44/ ../output/mirdeep2_trim_stringent/
mv report.log ../output/mirdeep2_trim_stringent/
mv pdfs_13_11_2024_t_08_48_44/ ../output/mirdeep2_trim_stringent/
mv mirdeep_runs/ ../output/mirdeep2_trim_stringent/
```

In this folder (`/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_trim_stringent/mirna_results_13_11_2024_t_08_48_44`), there are bed and fasta files for the known and novel mature, star and precursor sequences. I need to filter them by the miRNAs that I identified. Make a text file with the novel and known names of putative miRNAs.

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_trim_stringent/mirna_results_13_11_2024_t_08_48_44
nano putative_miRNA_list.txt 
Montipora_capitata_HIv3___Scaffold_2_98529
Montipora_capitata_HIv3___Scaffold_8_488597
Montipora_capitata_HIv3___Scaffold_4_264216
Montipora_capitata_HIv3___Scaffold_9_587940
Montipora_capitata_HIv3___Scaffold_9_625013
Montipora_capitata_HIv3___Scaffold_2_71488
Montipora_capitata_HIv3___Scaffold_9_613239
Montipora_capitata_HIv3___Scaffold_11_759461
Montipora_capitata_HIv3___Scaffold_12_859123
Montipora_capitata_HIv3___Scaffold_10_631026
Montipora_capitata_HIv3___Scaffold_8_585867
Montipora_capitata_HIv3___Scaffold_2_124999
Montipora_capitata_HIv3___Scaffold_8_578139
Montipora_capitata_HIv3___Scaffold_11_826703
Montipora_capitata_HIv3___Scaffold_10_729473
Montipora_capitata_HIv3___Scaffold_10_667068
Montipora_capitata_HIv3___Scaffold_151_1075760
Montipora_capitata_HIv3___Scaffold_14_1017829
Montipora_capitata_HIv3___Scaffold_5_309429
Montipora_capitata_HIv3___Scaffold_8_508612
Montipora_capitata_HIv3___Scaffold_8_553085
Montipora_capitata_HIv3___Scaffold_7_456843
Montipora_capitata_HIv3___Scaffold_4_291272
Montipora_capitata_HIv3___Scaffold_6_400996
Montipora_capitata_HIv3___Scaffold_13_954675
Montipora_capitata_HIv3___Scaffold_8_542727
Montipora_capitata_HIv3___Scaffold_14_999529
Montipora_capitata_HIv3___Scaffold_14_994729
Montipora_capitata_HIv3___Scaffold_9_617967
Montipora_capitata_HIv3___Scaffold_5_322659
Montipora_capitata_HIv3___Scaffold_6_419669
Montipora_capitata_HIv3___Scaffold_10_645863
Montipora_capitata_HIv3___Scaffold_2_87665
Montipora_capitata_HIv3___Scaffold_1_35393
Montipora_capitata_HIv3___Scaffold_4_232369
Montipora_capitata_HIv3___Scaffold_5_375083
Montipora_capitata_HIv3___Scaffold_2_64557
Montipora_capitata_HIv3___Scaffold_3_149822
Montipora_capitata_HIv3___Scaffold_9_606053
Montipora_capitata_HIv3___Scaffold_11_801837
Montipora_capitata_HIv3___Scaffold_8_508613
Montipora_capitata_HIv3___Scaffold_8_535535
Montipora_capitata_HIv3___Scaffold_11_796076
Montipora_capitata_HIv3___Scaffold_11_787579
Montipora_capitata_HIv3___Scaffold_11_816162
Montipora_capitata_HIv3___Scaffold_13_949160
Montipora_capitata_HIv3___Scaffold_11_841599
Montipora_capitata_HIv3___Scaffold_5_310446
Montipora_capitata_HIv3___Scaffold_6_439094
Montipora_capitata_HIv3___Scaffold_5_365006
Montipora_capitata_HIv3___Scaffold_2_67225
Montipora_capitata_HIv3___Scaffold_5_307262
Montipora_capitata_HIv3___Scaffold_14_1018693
Montipora_capitata_HIv3___Scaffold_2_96328
Montipora_capitata_HIv3___Scaffold_8_565039
Montipora_capitata_HIv3___Scaffold_8_586427
Montipora_capitata_HIv3___Scaffold_14_992284
Montipora_capitata_HIv3___Scaffold_12_910953
Montipora_capitata_HIv3___Scaffold_1_4183
```

Cat the known and novel miRNA fasta together and filter fasta so that I keep only the miRNAs in `putative_miRNA_list.txt`. 

```
cat novel_mature_13_11_2024_t_08_48_44_score-50_to_na.fa known_mature_13_11_2024_t_08_48_44_score-50_to_na.fa > putative_miRNAs.fa

grep -F -f putative_miRNA_list.txt putative_miRNAs.fa -A 1 --no-group-separator > putative_miRNAs_filt.fa
```

Do the same with the precursor sequences. 

```
cat novel_pres_13_11_2024_t_08_48_44_score-50_to_na.fa known_pres_13_11_2024_t_08_48_44_score-50_to_na.fa > putative_precursors.fa

grep -F -f putative_miRNA_list.txt putative_precursors.fa -A 1 --no-group-separator > putative_precursors_filt.fa
```

Similar to the mapper module, I can use a `config.txt` file to specify all the samples. I also need to use the `mapped_reads.fa` as the collapsed reads. 

In the scripts folder: `nano quantifier_mirdeep2.sh`

```
#!/bin/bash -i
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

echo "Quantifying smRNA counts" $(date)

conda activate /data/putnamlab/mirdeep2

quantifier.pl -r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_trim_stringent/mapped_reads.fa -p /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_trim_stringent/mirna_results_13_11_2024_t_08_48_44/putative_precursors_filt.fa -m /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_trim_stringent/mirna_results_13_11_2024_t_08_48_44/putative_miRNAs_filt.fa -c /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_trim_stringent/config.txt 

echo "Quantifying complete!" $(date)

conda deactivate
```

Submitted batch job 352797. Success! But...

```
# reads processed: 21804779
# reads with at least one reported alignment: 5795 (0.03%)
# reads that failed to align: 21798984 (99.97%)
Reported 6985 alignments to 1 output stream(s)
```

Only ~5700 reads aligned...I am going to mess around with the cutoffs. Going to add the following: 

- `-k` - also considers precursor-mature mappings that have different IDs
- `-g 3` - number of allowed mismatches when mapping reads to precursors, default is 1
- `-e 4` - number of nt upstream of mature sequence to consider, default is 2 
- `-f 7` - number of nt downstream of mature sequence to consider, default is 5 
- `-w` - considers whole precursor as mature sequence 

Submitted batch job 352801

Identify 3'UTRs. First, identify counts of each feature from gff file 

```
cd /data/putnamlab/jillashey/genome/Mcap/V3

module load BEDTools/2.30.0-GCC-11.3.0

grep -v '^#' Montipora_capitata_HIv3.genes.gff3 | cut -s -f 3 | sort | uniq -c | sort -rn > all_features_mcap.txt
256031 exon
 256031 CDS
  54384 transcript
```

Generate individual gff for transcripts 

```
grep $'\ttranscript\t' Montipora_capitata_HIv3.genes.gff3 > Mcap_transcript.gtf
```

Extract scaffold lengths 

```
cat is Montipora_capitata_HIv3.assembly.fasta | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Mcap_Chromosome_lengths.txt

wc -l Mcap_Chromosome_lengths.txt 
1699 Mcap_Chromosome_lengths.txt
```

Extract scaffold names 

```
awk -F" " '{print $1}' Mcap_Chromosome_lengths.txt > Mcap_Chromosome_names.txt
```

Sort gtf 

```
bedtools sort -faidx Mcap_Chromosome_names.txt -i Mcap_transcript.gtf > Mcap_transcript_sorted.gtf
```

Create 1kb 3'UTR in gff 

```
bedtools flank -i Mcap_transcript_sorted.gtf -g Mcap_Chromosome_lengths.txt -l 0 -r 1000 -s | awk '{gsub("transcript","3prime_UTR",$3); print $0 }' | awk '{if($5-$4 > 3)print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | tr ' ' '\t' > Mcap_3UTR_1kb.gtf
```

Subtract portions of 3'UTRs that overlap with nearby genes 

```
bedtools subtract -a Mcap_3UTR_1kb.gtf -b Mcap_transcript_sorted.gtf > Mcap_3UTR_1kb_corrected.gtf
```

Extract 3'UTR sequences from genome 

```
awk '{print $1 "\t" $4-1 "\t" $5 "\t" $9 "\t" "." "\t" $7}' Mcap_3UTR_1kb_corrected.gtf | sed 's/"//g' > Mcap_3UTR_1kb_corrected.bed

bedtools getfasta -fi Montipora_capitata_HIv3.assembly.fasta -bed Mcap_3UTR_1kb_corrected.bed -fo Mcap_3UTR_1kb.fasta -name

grep -c ">" Mcap_3UTR_1kb.fasta
55148
```

Run miranda for the Mcap data. In the scripts folder: `nano miranda_strict_all_1kb_mcap.sh`

```
#!/bin/bash -i
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

echo "Mcap starting miranda run with all genes and miRNAs with energy cutoff <-20 and strict binding invoked"$(date)

module load Miniconda3/4.9.2
conda activate /data/putnamlab/conda/miranda 

miranda /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_trim_stringent/mirna_results_13_11_2024_t_08_48_44/putative_miRNAs_filt.fa /data/putnamlab/jillashey/genome/Mcap/V3/Mcap_3UTR_1kb.fasta -en -20 -strict -out /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/miranda_trim_stringent/miranda_strict_all_1kb_mcap.tab

conda deactivate

echo "miranda run finished!"$(date)
echo "counting number of interactions possible" $(date)

zgrep -c "Performing Scan" /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/miranda_trim_stringent/miranda_strict_all_1kb_mcap.tab

echo "Parsing output" $(date)
grep -A 1 "Scores for this hit:" /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/miranda_trim_stringent/miranda_strict_all_1kb_mcap.tab | sort | grep '>' > /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/miranda_trim_stringent/miranda_strict_all_1kb_parsed_mcap.txt

echo "counting number of putative interactions predicted" $(date)
wc -l /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/miranda_trim_stringent/miranda_strict_all_1kb_parsed_mcap.txt

echo "Mcap DT miranda script complete" $(date)
```

Submitted batch job 352829. And the quantifier script finished running!

```
# reads processed: 21804779
# reads with at least one reported alignment: 20407 (0.09%)
# reads that failed to align: 21784372 (99.91%)
Reported 24370 alignments to 1 output stream(s)
```

Still very low alignment rate...and when looking at the csv, only the samples that were sequenced in the first batch got read counts. The others have 0s for basically all miRNAs. Need to look back at trimming. 

Miranda script ran! 

```
counting number of interactions possible Fri Dec 6 16:07:31 EST 2024
3253732
Parsing output Fri Dec 6 16:07:43 EST 2024
counting number of putative interactions predicted Fri Dec 6 16:07:44 EST 2024
25457 /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/miranda_trim_stringent/miranda_strict_all_1kb_parsed_mcap.txt
```

Copy onto local computer. Rerunning the QC for raw data (Submitted batch job 352838). 

### 20241230

Back at it! Holidays and travel have prevented me from working more on this. I reran the QC for the raw data and put it in my repo (link [here](https://github.com/JillAshey/DevelopmentalTimeseries/tree/main/data/Molecular/smRNA/raw_qc)). There is a lot of duplication. Even in the individual fastqc files, there is high duplication for each sample. My sequence duplication plotes look exactly like example 4 (see below) from this [page](https://proteo.me.uk/2013/09/a-new-way-to-look-at-duplication-in-fastqc-v0-11/): 

![](https://proteo.me.uk/wp-content/uploads/2013/09/high_duplication.png)

On the page, the author writes about this example: "This was about the worst example I could find.  This was an RNA-Seq library which had some nasty contamination and PCR duplication.  you can see that in the raw data the majority of reads come from sequences which are present tens of thousands of times in the library, but that there is also strong representation from reads which are present hundreds or thousands of times, suggesting that the duplication is widespreads. In this case deduplicating may help the experiment, but it causes the loss of 93% of the original data so it’s a pretty drastic step to take." Need to discuss with Hollie. 

I am going to rerun fastp on R1. Make new folder for them. 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data
mkdir trim_fastp
```

I already have a script (`fastp_trim_R1.sh`) that will run this. I am not going to set a length limit for now. I am also changing `--trim_poly_x` to `--trim_poly_g`, which is what was used in the [e5 deep dive code](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/code/08.1-Apul-sRNAseq-trimming-R1-only.Rmd). Submitted batch job 354306

### 20241231

Talked to Hollie briefly about my concerns. My data looks relatively similar to the deep dive [smRNA QC](https://gannet.fish.washington.edu/Atumefaciens/20230620-E5_coral-fastqc-flexbar-multiqc-sRNAseq/A_pulchra/raw_fastqc/multiqc_report.html). We also know we have a small diversity library (ie there are not that many unique miRNAs) and we sequenced to 28 million reads, which is a lot of sequencing to throw at a low diversity library. The Zymo miRNA library prep kit recommended to only sequence 75bp single end, but we did 150bp paired end. Following Sam's [notebook](https://robertslab.github.io/sams-notebook/posts/2023/2023-06-20-Trimming-and-QC---E5-Coral-sRNA-seq-Data-fro-A.pulchra-P.evermanni-and-P.meandrina-Using--FastQC-flexbar-and-MultiQC-on-Mox/), I am going to use flexbar to trim on R1 only. I am first only going to trim to 75 bp and then see how that looks. 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data
mkdir flexbar_75bp
```

In the scripts folder: `nano flexbar_75bp.sh`

```
#!/bin/bash 
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

echo "Flexbar trimming to 75bp" $(date)

module load Flexbar/3.5.0-foss-2018b  
module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw

# Make an array of sequences to trim 
array1=($(ls *R1_001.fastq.gz))

for i in ${array1[@]}; do
flexbar \
-r ${i} \
-a /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
-qf i1.8 \
-qt 25 \
--post-trim-length 75 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.${i} \
--zip-output GZ
done 

echo "Flexbar trimming complete, run QC" $(date)

for file in /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/*fastq.gz
do 
fastqc $file --outdir /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/fastqc/flexbar_75bp
done

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/fastqc/flexbar_75bp
multiqc *

echo "QC complete" $(date)
```

Submitted batch job 354418

### 20250101

Happy new year! Flexbar ran in about 9 hours. The QC data is [here](https://github.com/JillAshey/DevelopmentalTimeseries/tree/main/data/Molecular/smRNA/flexbar_75bp_qc). The length of samples are from 20 - 70bp, quite a range. The adapter content looks good though. Still a lot of overrepresented sequences and sequence duplication but we are going to move forward. Now I need to run the mapper portion of mirdeep2. According to the mirdeep2 [documentation](https://www.mdc-berlin.de/content/mirdeep2-documentation?mdcbl%5B0%5D=%2Fn-rajewsky%23t-data%2Csoftware%26resources&mdcbv=71nDTh7VzOJOW6SFGuFySs4mus4wnovu-t2LZzV2dL8&mdcot=6&mdcou=20738&mdctl=0) because I have multiple samples, I need to make a config file that contains the fq file locations and a unique 3 letter code. I will also be unzipping the files so remove the .gz from file name. 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp

nano config.txt 

/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.10_small_RNA_S4_R1_001.fastq.gz.fastq s10
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.11_small_RNA_S5_R1_001.fastq.gz.fastq s11
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.13_S76_R1_001.fastq.gz.fastq s13
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.14_small_RNA_S6_R1_001.fastq.gz.fastq s14
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.23_S77_R1_001.fastq.gz.fastq s23
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.24_small_RNA_S7_R1_001.fastq.gz.fastq s24
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.26_small_RNA_S8_R1_001.fastq.gz.fastq s26
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.28_small_RNA_S9_R1_001.fastq.gz.fastq s28
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.35_S78_R1_001.fastq.gz.fastq s35
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.36_small_RNA_S10_R1_001.fastq.gz.fastq s36
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.37_small_RNA_S11_R1_001.fastq.gz.fastq s37
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.39_small_RNA_S12_R1_001.fastq.gz.fastq s39
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.47_small_RNA_S13_R1_001.fastq.gz.fastq s47
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.48_small_RNA_S14_R1_001.fastq.gz.fastq s48
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.51_small_RNA_S15_R1_001.fastq.gz.fastq s51
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.52_S79_R1_001.fastq.gz.fastq s52
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.60_S80_R1_001.fastq.gz.fastq s60
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.61_small_RNA_S16_R1_001.fastq.gz.fastq s61
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.62_small_RNA_S17_R1_001.fastq.gz.fastq s62
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.63_small_RNA_S18_R1_001.fastq.gz.fastq s63
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.6_small_RNA_S1_R1_001.fastq.gz.fastq s06
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.72_S81_R1_001.fastq.gz.fastq s72
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.73_small_RNA_S19_R1_001.fastq.gz.fastq s73
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.74_small_RNA_S20_R1_001.fastq.gz.fastq s74
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.75_small_RNA_S21_R1_001.fastq.gz.fastq s75
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.7_small_RNA_S2_R1_001.fastq.gz.fastq s07
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.85_S82_R1_001.fastq.gz.fastq s85
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.86_small_RNA_S22_R1_001.fastq.gz.fastq s86
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.87_small_RNA_S23_R1_001.fastq.gz.fastq s87
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.88_small_RNA_S24_R1_001.fastq.gz.fastq s88
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.8_small_RNA_S3_R1_001.fastq.gz.fastq s08
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/trim.flexbar.75bp.9_S75_R1_001.fastq.gz.fastq s09
```

In the config.txt file, I have the path and file name, as well as the unique 3 letter code (s followed by sample number). In the scripts folder, I already have a mapper script (`mapper_mirdeep2.sh`) that I am going to edit with the updated paths. 

```
#!/bin/bash -i
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

#module load GCCcore/11.3.0 #I needed to add this to resolve conflicts between loaded GCCcore/9.3.0 and GCCcore/11.3.0
#module load Bowtie/1.3.1-GCC-11.3.0

#echo "Index Mcap genome" $(date)

# Index the reference genome for Mcap 
#bowtie-build /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.assembly.fasta /data/putnamlab/jillashey/genome/Mcap/V3/Mcap_ref.btindex

#module unload module load GCCcore/11.3.0 
#module unload Bowtie/1.3.1-GCC-11.3.0

#echo "Map trimmed reads with mirdeep2 mapper script" $(date)

conda activate /data/putnamlab/mirdeep2

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp

gunzip *.fastq.gz

mapper.pl config.txt -e -d -h -j -l 18 -m -q -p /data/putnamlab/jillashey/genome/Mcap/V3/Mcap_ref.btindex -s mapped_reads.fa -t mapped_reads_vs_genome.arf

echo "Mapping complete for trimmed stringent reads" $(date)

conda deactivate
```

Submitted batch job 354426. The flags mean the following: 

- `-e`- use fastq files as input
- `-d` - input is config file 
- `-h` - convert to fasta format 
- `-i` - convert RNA to DNA alphabet 
- `-j` - remove entries with non-canonical letters
- `-l` - discard reads shorter than 18 nt
- `-m` - collapse reads 
- `-p` - map to genome--genome has to be indexed by bowtie build 
- `-s` - mapped reads fasta output
- `-t` - read mappings v genome output
- `-q` - map with one mismatch in the seed 

### 20250102

Here is the mapping results: 

```
#desc   total   mapped  unmapped        %mapped %unmapped
total: 628271581        61907926        566363655       9.854   90.146
s06: 16391063   52942   16338121        0.323   99.677
s07: 26105703   4952449 21153254        18.971  81.029
s08: 20914298   102971  20811327        0.492   99.508
s09: 17402129   3044859 14357270        17.497  82.503
s10: 17319666   1252959 16066707        7.234   92.766
s11: 8702825    573453  8129372 6.589   93.411
s13: 18013957   4252214 13761743        23.605  76.395
s14: 22983933   205680  22778253        0.895   99.105
s23: 20418458   5323812 15094646        26.074  73.926
s24: 15238105   569210  14668895        3.735   96.265
s26: 29779739   1057567 28722172        3.551   96.449
s28: 20579195   104219  20474976        0.506   99.494
s35: 18792516   4519750 14272766        24.051  75.949
s36: 20592535   3164178 17428357        15.366  84.634
s37: 18149246   1286810 16862436        7.090   92.910
s39: 15148458   118789  15029669        0.784   99.216
s47: 57405680   13777537        43628143        24.000  76.000
s48: 23820154   361072  23459082        1.516   98.484
s51: 15857877   55736   15802141        0.351   99.649
s52: 22109116   5553213 16555903        25.117  74.883
s60: 14726932   3047585 11679347        20.694  79.306
s61: 21442064   142987  21299077        0.667   99.333
s62: 12666509   226583  12439926        1.789   98.211
s63: 14685987   125076  14560911        0.852   99.148
s72: 17626646   2760425 14866221        15.661  84.339
s73: 5016164    2083222 2932942 41.530  58.470
s74: 20074738   130657  19944081        0.651   99.349
s75: 23958496   411982  23546514        1.720   98.280
s85: 14366598   1732149 12634449        12.057  87.943
s86: 13903555   500450  13403105        3.599   96.401
s87: 21859884   216463  21643421        0.990   99.010
s88: 22219355   200927  22018428        0.904   99.096
```

The samples that were trimmed to be 30-20bp mapped better than the samples that were trimmed to 69-75bp. Should I trim now so that the rest of the samples are also ~30bp? Lets rerun the trimming but trim to 35bp max. 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data
mkdir flexbar_35bp
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/fastqc
mkdir flexbar_35bp
```

In the scripts folder: `nano flexbar_35bp.sh`

```
#!/bin/bash 
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

echo "Flexbar trimming to 35bp" $(date)

module load Flexbar/3.5.0-foss-2018b  
module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw

# Make an array of sequences to trim 
array1=($(ls *R1_001.fastq.gz))

for i in ${array1[@]}; do
flexbar \
-r ${i} \
-a /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
-qf i1.8 \
-qt 25 \
--post-trim-length 35 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.${i} \
--zip-output GZ
done 

echo "Flexbar trimming complete, run QC" $(date)

for file in /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/*fastq.gz
do 
fastqc $file --outdir /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/fastqc/flexbar_35bp
done

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/fastqc/flexbar_35bp
multiqc *

echo "QC complete" $(date)
```

Submitted batch job 354480. While this runs, going to run mirdeep2 miRNA predictions on the reads trimmed to 75bp. First, make new output folders and move mapper mirdeep2 output from data folder. 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output
mkdir mirdeep2_flexbar_75bp mirdeep2_flexbar_35bp
cd ../data/flexbar_75bp/
mv mapp* ../../output/mirdeep2_flexbar_75bp/
mv bowtie.log ../../output/mirdeep2_flexbar_75bp/
```

In the scripts folder: `nano mirdeep2_flexbar_75bp.sh`

```
#!/bin/bash -i
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

#module load Miniconda3/4.9.2
conda activate /data/putnamlab/mirdeep2

echo "Starting mirdeep2 with reads trimmed with flexbar to max 75bp" $(date)

miRDeep2.pl /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_flexbar_75bp/mapped_reads.fa /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.assembly.fasta /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_flexbar_75bp/mapped_reads_vs_genome.arf /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/mature_mirbase_cnidarian_T.fa none none -P -v -g -1 2>report.log

echo "mirdeep2 complete" $(date)

conda deactivate
```

Submitted batch job 354482. This will take about 2 days. 

### 20250103

Flexbar 35bp trimming finished running (QC data [here](https://github.com/JillAshey/DevelopmentalTimeseries/tree/main/data/Molecular/smRNA/flexbar_35bp_qc)). Higher amounts of duplication than the 75bp trimming iteration. Going to run mapping on these reads. Because I have multiple samples, I need to make a config file that contains the fq file locations and a unique 3 letter code. I will also be unzipping the files so remove the .gz from file name. 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp

nano config.txt 

/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.10_small_RNA_S4_R1_001.fastq.gz.fastq s10
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.11_small_RNA_S5_R1_001.fastq.gz.fastq s11
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.13_S76_R1_001.fastq.gz.fastq s13
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.14_small_RNA_S6_R1_001.fastq.gz.fastq s14
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.23_S77_R1_001.fastq.gz.fastq s23
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.24_small_RNA_S7_R1_001.fastq.gz.fastq s24
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.26_small_RNA_S8_R1_001.fastq.gz.fastq s26
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.28_small_RNA_S9_R1_001.fastq.gz.fastq s28
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.35_S78_R1_001.fastq.gz.fastq s35
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.36_small_RNA_S10_R1_001.fastq.gz.fastq s36
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.37_small_RNA_S11_R1_001.fastq.gz.fastq s37
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.39_small_RNA_S12_R1_001.fastq.gz.fastq s39
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.47_small_RNA_S13_R1_001.fastq.gz.fastq s47
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.48_small_RNA_S14_R1_001.fastq.gz.fastq s48
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.51_small_RNA_S15_R1_001.fastq.gz.fastq s51
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.52_S79_R1_001.fastq.gz.fastq s52
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.60_S80_R1_001.fastq.gz.fastq s60
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.61_small_RNA_S16_R1_001.fastq.gz.fastq s61
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.62_small_RNA_S17_R1_001.fastq.gz.fastq s62
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.63_small_RNA_S18_R1_001.fastq.gz.fastq s63
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.6_small_RNA_S1_R1_001.fastq.gz.fastq s06
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.72_S81_R1_001.fastq.gz.fastq s72
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.73_small_RNA_S19_R1_001.fastq.gz.fastq s73
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.74_small_RNA_S20_R1_001.fastq.gz.fastq s74
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.75_small_RNA_S21_R1_001.fastq.gz.fastq s75
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.7_small_RNA_S2_R1_001.fastq.gz.fastq s07
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.85_S82_R1_001.fastq.gz.fastq s85
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.86_small_RNA_S22_R1_001.fastq.gz.fastq s86
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.87_small_RNA_S23_R1_001.fastq.gz.fastq s87
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.88_small_RNA_S24_R1_001.fastq.gz.fastq s88
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.8_small_RNA_S3_R1_001.fastq.gz.fastq s08
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/trim.flexbar.35bp.9_S75_R1_001.fastq.gz.fastq s09
```

In the config.txt file, I have the path and file name, as well as the unique 3 letter code (s followed by sample number). In the scripts folder: `nano mapper_mirdeep2_flexbar_35bp.sh`

```
#!/bin/bash -i
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

#module load GCCcore/11.3.0 #I needed to add this to resolve conflicts between loaded GCCcore/9.3.0 and GCCcore/11.3.0
#module load Bowtie/1.3.1-GCC-11.3.0

#echo "Index Mcap genome" $(date)

# Index the reference genome for Mcap 
#bowtie-build /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.assembly.fasta /data/putnamlab/jillashey/genome/Mcap/V3/Mcap_ref.btindex

#module unload module load GCCcore/11.3.0 
#module unload Bowtie/1.3.1-GCC-11.3.0

echo "Map trimmed reads (flexbar 35bp) with mirdeep2 mapper script" $(date)

conda activate /data/putnamlab/mirdeep2

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp

gunzip *.fastq.gz

mapper.pl config.txt -e -d -h -j -l 18 -m -q -p /data/putnamlab/jillashey/genome/Mcap/V3/Mcap_ref.btindex -s mapped_reads.fa -t mapped_reads_vs_genome.arf

echo "Mapping complete for trimmed reads (flexbar 35bp)" $(date)

conda deactivate
```

Submitted batch job 354532

### 20250105

mirdeep2 miRNA prediction finished running early this morning. 





https://dnatech.ucdavis.edu/faqs/should-i-remove-pcr-duplicates-from-my-rna-seq-data 

Next steps for trimming ????

- Look at each raw QC file and trim each one using adapters and index primers 
- blast overrepresented seqs -- maybe they are something ?
- reverse complement of adapters and primers 
- what would data look like if i dont throw out anything >30bp 

Zoe cutadapt trimming: https://github.com/zdellaert/LaserCoral/blob/1b776313512822d368e957591890b2228c590bd5/code/RNA-seq-bioinf.md 


