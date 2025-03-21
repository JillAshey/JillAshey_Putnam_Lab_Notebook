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

Submitted batch job 354532. Move mapper files into proper folder once this is done running 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp
mv mapp* ../../output/mirdeep2_flexbar_35bp/
mv bowtie.log ../../output/mirdeep2_flexbar_35bp/
```

### 20250105

mirdeep2 miRNA prediction finished running early this morning. 

### 20250122

Ended up focusing on other things for sicb presentation but coming back to this analysis now. As a reminder, I ran two different trimming iterations with flexbar: one that trimmed data to 75bp max (`/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp`) and one that trimmed data to 35bp max (`/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp`). Here is what the mapping percentages looked like for each iteration: 

```
## 75bp trim: /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_flexbar_75bp
# reads processed: 2850328
# reads with at least one reported alignment: 1319485 (46.29%)
# reads that failed to align: 726277 (25.48%)
# reads with alignments suppressed due to -m: 804566 (28.23%)
Reported 2569169 alignments to 1 output stream(s)

## 35bp trim: /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_flexbar_35bp
# reads processed: 2781873
# reads with at least one reported alignment: 1318985 (47.41%)
# reads that failed to align: 685674 (24.65%)
# reads with alignments suppressed due to -m: 777214 (27.94%)
Reported 2567887 alignments to 1 output stream(s)
```

Kept and aligned slightly more reads using the 35bp trimming. 

I also ran the miRNA prediction on the 75bp trimmed reads. Move output into correct folder

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/scripts
mv report.log ../output/mirdeep2_flexbar_75bp/
mv *02_01_2025_t_14_40_55* ../output/mirdeep2_flexbar_75bp/
```

I downloaded the pdfs and html outputs to my computer and a quick look at them shows me they look as expected for miRNAs! Hooray! 

For the 35bp trimmed reads, I have run the mapper module and now need to run the mirdeep2 module. In the scripts folder: `nano mirdeep2_flexbar_35bp.sh`

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

echo "Starting mirdeep2 with reads trimmed with flexbar to max 35bp" $(date)

miRDeep2.pl /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_flexbar_35bp/mapped_reads.fa /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.assembly.fasta /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_flexbar_35bp/mapped_reads_vs_genome.arf /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/mature_mirbase_cnidarian_T.fa none none -P -v -g -1 2>report.log

echo "mirdeep2 on 35bp max reads complete" $(date)

conda deactivate
```

Submitted batch job 356032. While this is running, identify the putative miRNAs from the 75bp max mirdeep2 run using the [miRNA discovery code](https://github.com/JillAshey/DevelopmentalTimeseries/blob/main/scripts/miRNA_discovery.Rmd). 

In this folder (`/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_flexbar_75bp/mirna_results_02_01_2025_t_14_40_55`), there are bed and fasta files for the known and novel mature, star and precursor sequences. I need to filter them by the miRNAs that I identified. Make a text file with the novel and known names of putative miRNAs.

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_flexbar_75bp/mirna_results_02_01_2025_t_14_40_55
nano flexbar_75bp_putative_miRNA_list.txt
Montipora_capitata_HIv3___Scaffold_14_866635
Montipora_capitata_HIv3___Scaffold_8_480075
Montipora_capitata_HIv3___Scaffold_8_498263
Montipora_capitata_HIv3___Scaffold_12_774709
Montipora_capitata_HIv3___Scaffold_1_3519
Montipora_capitata_HIv3___Scaffold_8_414648
Montipora_capitata_HIv3___Scaffold_2_83947
Montipora_capitata_HIv3___Scaffold_8_415215
Montipora_capitata_HIv3___Scaffold_9_531157
Montipora_capitata_HIv3___Scaffold_2_60896
Montipora_capitata_HIv3___Scaffold_9_521005
Montipora_capitata_HIv3___Scaffold_8_497799
Montipora_capitata_HIv3___Scaffold_2_119703
Montipora_capitata_HIv3___Scaffold_2_106363
Montipora_capitata_HIv3___Scaffold_8_491119
Montipora_capitata_HIv3___Scaffold_11_702925
Montipora_capitata_HIv3___Scaffold_10_620169
Montipora_capitata_HIv3___Scaffold_12_773601
Montipora_capitata_HIv3___Scaffold_1_48204
Montipora_capitata_HIv3___Scaffold_14_865867
Montipora_capitata_HIv3___Scaffold_151_915764
Montipora_capitata_HIv3___Scaffold_5_262937
Montipora_capitata_HIv3___Scaffold_8_432272
Montipora_capitata_HIv3___Scaffold_7_388195
Montipora_capitata_HIv3___Scaffold_4_247620
Montipora_capitata_HIv3___Scaffold_6_340912
Montipora_capitata_HIv3___Scaffold_8_461265
Montipora_capitata_HIv3___Scaffold_14_850199
Montipora_capitata_HIv3___Scaffold_9_525121
Montipora_capitata_HIv3___Scaffold_11_696080
Montipora_capitata_HIv3___Scaffold_14_846123
Montipora_capitata_HIv3___Scaffold_10_585455
Montipora_capitata_HIv3___Scaffold_4_197697
Montipora_capitata_HIv3___Scaffold_5_318975
Montipora_capitata_HIv3___Scaffold_3_127504
Montipora_capitata_HIv3___Scaffold_9_515007
Montipora_capitata_HIv3___Scaffold_11_681829
Montipora_capitata_HIv3___Scaffold_8_432273
Montipora_capitata_HIv3___Scaffold_8_455105
Montipora_capitata_HIv3___Scaffold_13_807296
Montipora_capitata_HIv3___Scaffold_5_263768
Montipora_capitata_HIv3___Scaffold_11_634165
Montipora_capitata_HIv3___Scaffold_9_530478
Montipora_capitata_HIv3___Scaffold_6_373318
Montipora_capitata_HIv3___Scaffold_2_57253
```

Cat the known and novel miRNA fasta together and filter fasta so that I keep only the miRNAs in `flexbar_75bp_putative_miRNA_list.txt`.

```
cat known_mature_02_01_2025_t_14_40_55_score-50_to_na.fa novel_mature_02_01_2025_t_14_40_55_score-50_to_na.fa > flexbar_75bp_putative_miRNAs.fa

grep -F -f flexbar_75bp_putative_miRNA_list.txt flexbar_75bp_putative_miRNAs.fa -A 1 --no-group-separator > flexbar_75bp_putative_miRNAs_filt.fa
```

Do the same with the precursor sequences

```
cat known_pres_02_01_2025_t_14_40_55_score-50_to_na.fa novel_pres_02_01_2025_t_14_40_55_score-50_to_na.fa > flexbar_75bp_putative_pres.fa

grep -F -f flexbar_75bp_putative_miRNA_list.txt flexbar_75bp_putative_pres.fa -A 1 --no-group-separator > flexbar_75bp_putative_pres_filt.fa
```

Similar to the mapper module, I can use a `config.txt` file to specify all the samples. I also need to use the `mapped_reads.fa` as the collapsed reads.

In the scripts folder: `nano quant_mirdeep2_flexbar_75bp.sh`

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

echo "Quantifying smRNA counts for flexbar 75bp trimmed reads" $(date)

conda activate /data/putnamlab/mirdeep2

quantifier.pl -r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_flexbar_75bp/mapped_reads.fa -p /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_flexbar_75bp/mirna_results_02_01_2025_t_14_40_55/flexbar_75bp_putative_pres_filt.fa -m /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_flexbar_75bp/mirna_results_02_01_2025_t_14_40_55/flexbar_75bp_putative_miRNAs_filt.fa -c /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_75bp/config.txt -k -g 3 -e 4 -f 7 -w

echo "Quantifying complete for flexbar 75bp trimmed reads!" $(date)

conda deactivate
```

Submitted batch job 356039

### 20250210

I am back and ready to trim. Let's look at the output of the quant script I ran with the 75bp trim.

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/scripts
less slurm-356039.error

#desc   total   mapped  unmapped        %mapped %unmapped
total: 628271581        259598  628011983       0.041   99.959
s06: 16391063   72      16390991        0.000   100.000
s07: 26105703   0       26105703        0.000   100.000
s08: 20914298   247     20914051        0.001   99.999
s09: 17402129   16948   17385181        0.097   99.903
s10: 17319666   20      17319646        0.000   100.000
s11: 8702825    0       8702825 0.000   100.000
s13: 18013957   27107   17986850        0.150   99.850
s14: 22983933   416     22983517        0.002   99.998
s23: 20418458   27534   20390924        0.135   99.865
s24: 15238105   1055    15237050        0.007   99.993
s26: 29779739   2       29779737        0.000   100.000
s28: 20579195   13      20579182        0.000   100.000
s35: 18792516   36247   18756269        0.193   99.807
s36: 20592535   0       20592535        0.000   100.000
s37: 18149246   0       18149246        0.000   100.000
s39: 15148458   0       15148458        0.000   100.000
s47: 57405680   5       57405675        0.000   100.000
s48: 23820154   0       23820154        0.000   100.000
s51: 15857877   4       15857873        0.000   100.000
s52: 22109116   35915   22073201        0.162   99.838
s60: 14726932   36548   14690384        0.248   99.752
s61: 21442064   3       21442061        0.000   100.000
s62: 12666509   0       12666509        0.000   100.000
s63: 14685987   2139    14683848        0.015   99.985
s72: 17626646   38663   17587983        0.219   99.781
s73: 5016164    0       5016164 0.000   100.000
s74: 20074738   1       20074737        0.000   100.000
s75: 23958496   849     23957647        0.004   99.996
s85: 14366598   35275   14331323        0.246   99.754
s86: 13903555   8       13903547        0.000   100.000
s87: 21859884   500     21859384        0.002   99.998
s88: 22219355   27      22219328        0.000   100.000
```

This is showing the number of reads mapped to the miRNAs identified. The samples that I sent for test sequencing look good but the others not so much. I also ran mirdeep2 for the 35 bp trimmed reads. Move all of this output to the correct output folder 

```
# Data from mirdeep2 quant with 75bp trim
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/scripts/expression_analyses
mv expression_analyses_1737578046 ../../output/mirdeep2_flexbar_75bp
cd ../
mv *1737578046* ../output/mirdeep2_flexbar_75bp/

# Data from mirdeep2 with 35bp trim
mv *22_01_2025_t_14_50_27* ../output/mirdeep2_flexbar_35bp/
mv mirdeep_runs ../output/mirdeep2_flexbar_35bp/
mv dir_prepare_signature1737576068 ../output/mirdeep2_flexbar_35bp/
```

It doesn't look like the 75bp trimming did much to improve mapping rates. My next step is to identify the miRNAs from the mirdeep2 output with the 35bp reads. Copy `result_22_01_2025_t_14_50_27.csv` from `/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_flexbar_35bp` onto local computer to ID miRNAs. In `/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_flexbar_35bp/mirna_results_22_01_2025_t_14_50_27`, there are bed and fasta files for the known and novel mature, star and precursor sequences. I need to filter them by the miRNAs that I identified. Make a text file with the novel and known names of putative miRNAs.

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_flexbar_35bp/mirna_results_22_01_2025_t_14_50_27
nano flexbar_35bp_putative_miRNA_list.txt
Montipora_capitata_HIv3___Scaffold_1157_961045
Montipora_capitata_HIv3___Scaffold_8_414648
Montipora_capitata_HIv3___Scaffold_2_83947
Montipora_capitata_HIv3___Scaffold_8_415215
Montipora_capitata_HIv3___Scaffold_9_531157
Montipora_capitata_HIv3___Scaffold_2_60896
Montipora_capitata_HIv3___Scaffold_9_521005
Montipora_capitata_HIv3___Scaffold_8_497799
Montipora_capitata_HIv3___Scaffold_2_119703
Montipora_capitata_HIv3___Scaffold_2_106363
Montipora_capitata_HIv3___Scaffold_8_491119
Montipora_capitata_HIv3___Scaffold_11_702925
Montipora_capitata_HIv3___Scaffold_10_620169
Montipora_capitata_HIv3___Scaffold_12_773601
Montipora_capitata_HIv3___Scaffold_151_915764
Montipora_capitata_HIv3___Scaffold_5_262937
Montipora_capitata_HIv3___Scaffold_8_432272
Montipora_capitata_HIv3___Scaffold_7_388195
Montipora_capitata_HIv3___Scaffold_4_247620
Montipora_capitata_HIv3___Scaffold_6_340912
Montipora_capitata_HIv3___Scaffold_8_461265
Montipora_capitata_HIv3___Scaffold_14_850199
Montipora_capitata_HIv3___Scaffold_9_525121
Montipora_capitata_HIv3___Scaffold_14_846123
Montipora_capitata_HIv3___Scaffold_4_197697
Montipora_capitata_HIv3___Scaffold_5_318975
Montipora_capitata_HIv3___Scaffold_3_127504
Montipora_capitata_HIv3___Scaffold_9_515007
Montipora_capitata_HIv3___Scaffold_11_681829
Montipora_capitata_HIv3___Scaffold_8_432273
Montipora_capitata_HIv3___Scaffold_8_455105
Montipora_capitata_HIv3___Scaffold_13_807296
Montipora_capitata_HIv3___Scaffold_5_263768
Montipora_capitata_HIv3___Scaffold_6_373318
Montipora_capitata_HIv3___Scaffold_9_530478
Montipora_capitata_HIv3___Scaffold_2_57253
Montipora_capitata_HIv3___Scaffold_5_310384
Montipora_capitata_HIv3___Scaffold_14_866635
Montipora_capitata_HIv3___Scaffold_8_480075
Montipora_capitata_HIv3___Scaffold_8_498263
Montipora_capitata_HIv3___Scaffold_12_774709
Montipora_capitata_HIv3___Scaffold_1_3519
```

Cat known and novel miRNA fasta files together and filter fasta so that I keep only miRNAs in 

```
cat known_mature_22_01_2025_t_14_50_27_score-50_to_na.fa novel_mature_22_01_2025_t_14_50_27_score-50_to_na.fa > flexbar_35bp_putative_miRNAs.fa

grep -F -f flexbar_35bp_putative_miRNA_list.txt flexbar_35bp_putative_miRNAs.fa -A 1 --no-group-separator > flexbar_35bp_putative_miRNAs_filt.fa
```

Do the same with the precursor sequences

```
cat known_pres_22_01_2025_t_14_50_27_score-50_to_na.fa novel_pres_22_01_2025_t_14_50_27_score-50_to_na.fa > flexbar_35bp_putative_pres.fa
grep -F -f flexbar_35bp_putative_miRNA_list.txt flexbar_35bp_putative_pres.fa -A 1 --no-group-separator > flexbar_35bp_putative_pres_filt.fa
```

Similar to the mapper module, I can use a `config.txt` file to specify all the samples. I also need to use the mapped_reads.fa as the collapsed reads.

In the scripts folder: `nano quant_mirdeep2_flexbar_35bp.sh`

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

echo "Quantifying smRNA counts for flexbar 35bp trimmed reads" $(date)

conda activate /data/putnamlab/mirdeep2

quantifier.pl -r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_flexbar_35bp/mapped_reads.fa -p /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_flexbar_35bp/mirna_results_22_01_2025_t_14_50_27/flexbar_35bp_putative_pres_filt.fa -m /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_flexbar_35bp/mirna_results_22_01_2025_t_14_50_27/flexbar_35bp_putative_miRNAs_filt.fa -c /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/config.txt -k -g 3 -e 4 -f 7 -w

echo "Quantifying complete for flexbar 35bp trimmed reads!" $(date)

conda deactivate
```

Submitted batch job 360856. Looking at the mirdeep2 output from the identification module, the counts in the results files are still quite low so not sure this trimming was as effective as I had hoped. The quantification is complete, let's look at how it did: 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/scripts
less slurm-360856.error

#desc   total   mapped  unmapped        %mapped %unmapped
total: 628271581        235065  628036516       0.037   99.963
s06: 16391063   36      16391027        0.000   100.000
s07: 26105703   0       26105703        0.000   100.000
s08: 20914298   247     20914051        0.001   99.999
s09: 17402129   15012   17387117        0.086   99.914
s10: 17319666   13      17319653        0.000   100.000
s11: 8702825    0       8702825 0.000   100.000
s13: 18013957   17840   17996117        0.099   99.901
s14: 22983933   418     22983515        0.002   99.998
s23: 20418458   22185   20396273        0.109   99.891
s24: 15238105   4312    15233793        0.028   99.972
s26: 29779739   0       29779739        0.000   100.000
s28: 20579195   12      20579183        0.000   100.000
s35: 18792516   26114   18766402        0.139   99.861
s36: 20592535   0       20592535        0.000   100.000
s37: 18149246   0       18149246        0.000   100.000
s39: 15148458   0       15148458        0.000   100.000
s47: 57405680   0       57405680        0.000   100.000
s48: 23820154   0       23820154        0.000   100.000
s51: 15857877   1       15857876        0.000   100.000
s52: 22109116   29771   22079345        0.135   99.865
s60: 14726932   31455   14695477        0.214   99.786
s61: 21442064   0       21442064        0.000   100.000
s62: 12666509   0       12666509        0.000   100.000
s63: 14685987   2139    14683848        0.015   99.985
s72: 17626646   45491   17581155        0.258   99.742
s73: 5016164    0       5016164 0.000   100.000
s74: 20074738   0       20074738        0.000   100.000
s75: 23958496   813     23957683        0.003   99.997
s85: 14366598   38691   14327907        0.269   99.731
s86: 13903555   6       13903549        0.000   100.000
s87: 21859884   482     21859402        0.002   99.998
s88: 22219355   27      22219328        0.000   100.000
```

Bleh...basically the same as before. Back to the drawing board for trimming. I may need to trim each sample a very specific way. 

### 20250212

Digging in deep. I am thinking that I have a lot of PCR duplication but I am not sure if it is PCR or biological. UMIs are unique molecular identifiers that can be used to distinguish PCR dups from biologically real dups. I know that I put adapters and Illumina sequences but not sure if I or the sequencing facility added UMIs...I do not think that I did in my library prep protocol


Let me look at the raw sequences for one sample (M6). I think I will only be using R1. 

```
# Raw M6 fastq
@LH00260:104:22GLL5LT4:2:1101:7522:1056 1:N:0:ATGAGGCCAT+CAATTAACGT
CNCACGTCTGAACTCCAGTCACATGAAGCCATCTAGGGGGGTGGGTGGGGGTTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
+
I#IIIII9II-IIIIIIIIIIIIII9-IIIII-----IIII-9I999-9----999999--999999I999I9II9-IIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@LH00260:104:22GLL5LT4:2:1101:11148:1056 1:N:0:ATGAGTCCAT+CNATTAACGT
ANCGGAGCGCACACGTCTGAACTCCAGTCACATAAGGCCATCTGGGGGGGGGGGGTGGGGGGGGGGGGGGGGGGGGGGGGAGGGGGGGGGGGGGGGGGGGGGTGGGGGGGGGTTGTAGGGGGGGGGGGGGGGGGGGTGGGGGGGGGGTAGG
+
9#9I9I-99IIIII9IIII999IIIIIIIII99-9I9II9999II-9III9II-9III-I-----I---II-I9IIII-I-99I9-9II9I--I-9-9I99I-9I99-9999-----9--99-9-9--9999-9-9--99---I-99-9--
@LH00260:104:22GLL5LT4:2:1101:14635:1070 1:N:0:ATGAGGCCAT+CAATTAACGT
ANATACACGTCTGAACTCCAGTCACATGAGGCCATCTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
+
I#IIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIII9III--IIIIII9I-II9-9II9I-II-IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

A lot of Gs present in the latter half of the read. Let's look at the M6sample trimmed to 35bp. While this is the sample sample, these are not the same reads as the ones above. 

```
@LH00260:104:22GLL5LT4:2:1101:4018:1098 1:N:0:ATGAGGCCAT+CAATTAACGT
CACACGTCTGAACTCCAGTCACATGAGGCCATCTG
+
IIIIII9II9IIIIIIIIIIIIII9IIIII9IIII
@LH00260:104:22GLL5LT4:2:1101:6001:1112 1:N:0:ATGAGGCCGT+CAATTAACGT
ACACGTCTGAACTCCAGTCACATGAGGCCATCTGG
+
III9IIIIIIIIII-IIIIIIIIIIII99IIII9I
@LH00260:104:22GLL5LT4:2:1101:6422:1112 1:N:0:ATGAGGCCAT+CAATTAACGT
GAGCACACGTCTGACCTCCAGTCACATGAGGCCAT
+
II-I9IIIIIIIIIIIIIIIII9III9IIIIII9I
```

Interesting - the I symbol corresponds to a phred score of 40 and a 0.0001 probability of an incorrect base call. 9 corresponds to a phred score of 24 with a 0.004 probablity of an incorrect base call. When I trim with flexbar, I am setting the minimum quality threshold to 25, althought the default is 20. 

I found the same read in both raw and trimmed data as an example:

```
# Raw read 
@LH00260:104:22GLL5LT4:2:1101:6001:1112 1:N:0:ATGAGGCCGT+CAATTAACGT
ACACGTCTGAACTCCAGTCACATGAGGCCATCTGGGGGGGGTTTTTTTGGGGGGGCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
+
III9IIIIIIIIII-IIIIIIIIIIII99IIII9I9-9II-----9-9----99-----99-99--99-9999999I9---99--99I-II--9-I-III99II9IIIIIIIII-I9IIIIIIII-IIIIIII9IIIIII9III-99IIII

# Flexbar 35bp trim read 
@LH00260:104:22GLL5LT4:2:1101:6001:1112 1:N:0:ATGAGGCCGT+CAATTAACGT
ACACGTCTGAACTCCAGTCACATGAGGCCATCTGG
+
III9IIIIIIIIII-IIIIIIIIIIII99IIII9I
```

Okay now what? I also emailed Zymo to see if they have any recommendations. I want to convert one of the fastq samples into a fasta file and then blast it against the Mcap genome. In the scripts folder: `nano fasta_blast_M6.sh`

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

module load GCC/9.3.0 
module load seqtk/1.3-GCC-9.3.0
module load BLAST+/2.9.0-iimpi-2019b

echo "Convert trimmed M6 fastq file to fasta file" $(date)

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/

seqtk seq -a trim.flexbar.35bp.6_small_RNA_S1_R1_001.fastq.gz.fastq > trim.flexbar.35bp.6_small_RNA_S1_R1_001.fasta

echo "Fastq to fastq complete for M6" $(date)
echo "Make blast db from Mcap genome" $(date)

makeblastdb -in /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.assembly.fasta -out /data/putnamlab/jillashey/genome/Mcap/V3/Mcap_genome_V3 -dbtype nucl

echo "Blast M6 fasta against Mcap genome" $(date)

blastn -query trim.flexbar.35bp.6_small_RNA_S1_R1_001.fasta -db /data/putnamlab/jillashey/genome/Mcap/V3/Mcap_genome_V3 -out M6_genome_blast_results.txt -outfmt 0

blastn -query trim.flexbar.35bp.6_small_RNA_S1_R1_001.fasta -db /data/putnamlab/jillashey/genome/Mcap/V3/Mcap_genome_V3 -out M6_genome_blast_results_tab.txt -outfmt 6 -max_target_seqs 3

echo "Blast complete" $(date)
```

Submitted batch job 361447. Idk how this will help me but I want to see that the sequences actually map to the genome. 

I could also collapse the reads...The blast script above finished running in about an hour. Look at the data

```
head M6_genome_blast_results_tab.txt
LH00260:104:22GLL5LT4:2:1101:34494:1154	Montipora_capitata_HIv3___Scaffold_1699	100.000	34	0	0	1	34	695	662	4.97e-10	63.9
LH00260:104:22GLL5LT4:2:1101:34494:1154	Montipora_capitata_HIv3___Scaffold_1288	100.000	34	0	0	1	34	4524	4557	4.97e-10	63.9
LH00260:104:22GLL5LT4:2:1101:34494:1154	Montipora_capitata_HIv3___Scaffold_1157	100.000	34	0	0	1	34	4349	4316	4.97e-10	63.9
LH00260:104:22GLL5LT4:2:1101:34494:1154	Montipora_capitata_HIv3___Scaffold_1157	100.000	34	0	0	1	34	19191	19158	4.97e-10	63.9
LH00260:104:22GLL5LT4:2:1101:27186:1364	Montipora_capitata_HIv3___Scaffold_1699	100.000	30	0	0	1	30	3489	3460	5.29e-08	56.5
LH00260:104:22GLL5LT4:2:1101:27186:1364	Montipora_capitata_HIv3___Scaffold_1288	100.000	30	0	0	1	30	1731	1760	5.29e-08	56.5
LH00260:104:22GLL5LT4:2:1101:27186:1364	Montipora_capitata_HIv3___Scaffold_1288	100.000	30	0	0	1	30	11607	11636	5.29e-08	56.5
LH00260:104:22GLL5LT4:2:1101:27186:1364	Montipora_capitata_HIv3___Scaffold_1157	100.000	30	0	0	1	30	7145	7116	5.29e-08	56.5
LH00260:104:22GLL5LT4:2:1101:44965:1546	Montipora_capitata_HIv3___Scaffold_839	100.000	32	0	0	2	33	7302	7333	5.84e-09	60.2
LH00260:104:22GLL5LT4:2:1101:44965:1546	Montipora_capitata_HIv3___Scaffold_478	100.000	32	0	0	2	33	25150	25181	5.84e-09	60.2

# Number of unique terms in first column - ie the sequencing reads 
cut -f1 M6_genome_blast_results_tab.txt | sort -u | wc -l
45855

# Number of unique terms in first column - ie the genome
cut -f2 M6_genome_blast_results_tab.txt | sort -u | wc -l
1136

# Number of lines in file 
wc -l M6_genome_blast_results_tab.txt
219193
```

- Total Number of Reads: 16391063
- Unique Reads that Mapped (from BLAST results): 45855
- Total Alignments (from BLAST results): 219193

So this information indicates that only 0.28% of the total reads in this sample uniquely mapped to the genome (with the caveat that this was blast not an actual aligner). This suggests a very high level of duplication. Going to collapse the reads to see how many unique reads we actually have. 

```
module load FASTX-Toolkit/0.0.14-GCC-9.3.0 
fastx_collapser -v -i trim.flexbar.35bp.6_small_RNA_S1_R1_001.fastq.gz.fastq -o collapse.trim.flexbar.35bp.6_small_RNA_S1_R1_001.fastq
```

Here are the collapsed read information:

```
Input: 16391063 sequences (representing 16391063 reads)
Output: 179728 sequences (representing 16391063 reads)

head collapse.trim.flexbar.35bp.6_small_RNA_S1_R1_001.fastq 
>1-2633861
CACACGTCTGAACTCCAGTCACATGAGGCCATCTG
>2-2598169
GCACACGTCTGAACTCCAGTCACATGAGGCCATCT
>3-1314236
AGCACACGTCTGAACTCCAGTCACATGAGGCCATC
>4-1219064
GAGCACACGTCTGAACTCCAGTCACATGAGGCCAT
>5-710687
AGAGCACACGTCTGAACTCCAGTCACATGAGGCCA

grep -c ">" collapse.trim.flexbar.35bp.6_small_RNA_S1_R1_001.fastq 
179728

grep -c ">" trim.flexbar.35bp.6_small_RNA_S1_R1_001.fasta
16391063
```

Okay so our library is really not complex at all. Duplication is v high. I am going to use bowtie to map the collapsed reads to the genome and mirbase. In the scripts folder: `nano bowtie_M6.sh`

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

module load Bowtie/1.3.1-GCC-11.3.0
#other options if needed: Bowtie/1.2.3-GCC-8.3.0, Bowtie/1.3.1-GCC-11.2.0, Bowtie/1.2.2-foss-2018b 

echo "Use bowtie to align collapsed reads to genome for M6" $(date)

# genome already has bowtie index built, do not need to rebuild

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/

bowtie -a --verbose --best --strata -n 1 -x /data/putnamlab/jillashey/genome/Mcap/V3/Mcap_ref.btindex -f /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/collapse.trim.flexbar.35bp.6_small_RNA_S1_R1_001.fastq -S /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/aligned_genome_M6.sam

echo "Bowtie alignment complete for M6 to genome" $(date)
#echo "Build mirbase bowtie index" $(date)

bowtie-build /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/mature_mirbase_cnidarian_T.fa /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/mirbase_ref

echo "Build complte, align collapsed reads to mirbase for M6" $(date)

bowtie -a --verbose --best --strata -n 1 -x /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/mirbase_ref -f /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/collapse.trim.flexbar.35bp.6_small_RNA_S1_R1_001.fastq -S /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/aligned_mirbase_M6.sam

echo "Alignments complete!" $(date)
```

There are a lot of different alignment options with bowtie in terms of mismatches and what not, so I can play with settings if needed. Submitted batch job 361504. Did the genome in about an hour. Needed to fix the mirbase ref file path. Submitted batch job 361505

I am looking at the genome aligned file and I put some of the sequences that got hits on the genome into blast and they are turning up as rRNA...

### 20250213

Downloading `aligned_genome_M6.sam` and `aligned_mirbase_M6.sam` to my computer to look at them on IGB. Looking at the sam file aligned to the genome...doesn't look great. Here is an example of one of the chromosomes: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/refs/heads/master/images/IGB_mcap_miRNA_genome_example_chrom.png)

There is a lot of sequences aligning all over the place but its hard to know if any are biologically meaningful. I also looked at one of the sequences close up:

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/refs/heads/master/images/IGB_mcap_miRNA_genome_example_seq.png)

This read is just GAGAGAGA repeating. When looking at the sequences aligned to the mirbase db, I see several sequences aligning to nve-miR-9415 from Nematostella but no sequences are aligning to any other miRNAs, which is strange to me. Let's go back to the server and bowtie align the collapsed M6 file against the putative miRNAs identified from the Mcap data. 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_flexbar_35bp/mirna_results_22_01_2025_t_14_50_27

grep -c ">" flexbar_35bp_putative_miRNAs_filt.fa
42
```

The `flexbar_35bp_putative_miRNAs_filt.fa` file contains the putative miRNAs that were identified from the 35bp trimmed data. When I quantified the counts of these miRNAs for each sample, a lot of the re-amped samples got 0 counts for many of the putative miRNAs. Build a bowtie index for this fasta file 

```
interactive
module load Bowtie/1.3.1-GCC-11.3.0

bowtie-build flexbar_35bp_putative_miRNAs_filt.fa putative_miRNA_ref
```

Align the collapsed M6 file to the `putative_miRNA_ref` index.

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/

bowtie -a --verbose --best --strata -n 1 -x /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_flexbar_35bp/mirna_results_22_01_2025_t_14_50_27/putative_miRNA_ref -f collapse.trim.flexbar.35bp.6_small_RNA_S1_R1_001.fastq -S aligned_putative_mirna_M6.sam

# reads processed: 179728
# reads with at least one alignment: 8 (0.00%)
# reads that failed to align: 179720 (100.00%)
Reported 8 alignments
```

8 alignments reported LOL. When I look through the sam file itself, I can't even find the aligned reads. THIS SUCKS I SUCK!!!!!!

I am going to collapse one of the samples that went through the first sequencing iteration. There were 8 samples that we sent (n=1 per time point) initially to be sequenced to check if the library prep worked. I did not re-amp these samples and they turned out with much less duplication than the rest of the samples. First, collapse M9 fastq

```
module load FASTX-Toolkit/0.0.14-GCC-9.3.0 
fastx_collapser -v -i trim.flexbar.35bp.9_S75_R1_001.fastq.gz.fastq -o collapse.trim.flexbar.35bp.9_S75_R1_001.fastq
Input: 17402129 sequences (representing 17402129 reads)
Output: 2781873 sequences (representing 17402129 reads)

head collapse.trim.flexbar.35bp.9_S75_R1_001.fastq
>1-1339235
AAACCTTTGTTCTAAGATCGA
>2-969023
AATGTCTGTCTGAGGGTCGAA
>3-715497
AATCTGAGATTCAACCTTTGTTCTAAGATCGA
>4-607457
ACCCGGGAGCATGTCTGTCTGAGGGTCGAA
>5-451812
AGTCTGTCTGAGGGTCGAA

grep -c ">" collapse.trim.flexbar.35bp.9_S75_R1_001.fastq
2781873

grep -c "@LH" trim.flexbar.35bp.9_S75_R1_001.fastq.gz.fastq
17402129
```

Library still isn't super complex, but more so than the M6 sample. 


I am going to use bowtie to map the collapsed reads to the genome and mirbase. In the scripts folder: `nano bowtie_M9.sh`

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

module load Bowtie/1.3.1-GCC-11.3.0
#other options if needed: Bowtie/1.2.3-GCC-8.3.0, Bowtie/1.3.1-GCC-11.2.0, Bowtie/1.2.2-foss-2018b 

echo "Use bowtie to align collapsed reads to genome for M9" $(date)

# all indices already built, no need to redo 

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_35bp/

bowtie -a --verbose --best --strata -n 1 -x /data/putnamlab/jillashey/genome/Mcap/V3/Mcap_ref.btindex -f collapse.trim.flexbar.35bp.9_S75_R1_001.fastq -S aligned_genome_M9.sam

echo "Bowtie alignment complete for M9 to genome" $(date)
echo "Use bowtie to align collapsed reads to mirbase for M9" $(date)

bowtie -a --verbose --best --strata -n 1 -x /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/mirbase_ref -f collapse.trim.flexbar.35bp.9_S75_R1_001.fastq -S aligned_mirbase_M9.sam

echo "Bowtie alignment complete for M9 to mirbase" $(date)
echo "Use bowtie to align collapsed reads to putative miRNAs for M9" $(date)

bowtie -a --verbose --best --strata -n 1 -x /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_flexbar_35bp/mirna_results_22_01_2025_t_14_50_27/putative_miRNA_ref -f collapse.trim.flexbar.35bp.9_S75_R1_001.fastq -S aligned_putative_mirna_M9.sam

echo "Alignments complete!" $(date)
```

### 20250214

Talked w/ Hollie yesterday about the data. Our take-aways: Chapter 3 the problem child - samples that were sequenced in the second batch have high duplication and a lot of adapter/index sequences. We are expecting to have a lot of duplication anyway because these libraries are not very complex. We suspect that this is due to the re-amplification that I had to do with these samples, as the Taq premix didn't initially work. I aligned one sample yesterday and very few reads aligned to the genome. I am going to use the index sequences as inputs for the trimming step to see if I can rid the samples of their trash sequences. If that does not work, we will use n=1 of our initial miRNA sequencing samples and assume that the miRNA pattern is the same across all time points. In this case, we will only look at the genes that those miRNA are binding and look at expression through time.

It always comes down to more trimming. Taking another look at `/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent_take2/` this trimming iteration. I looked at the multiqc and it definitely didn't do what I wanted. still high levels of duplication. Okay...should I run each sample individually?? 

### 20250217

Talked with Zoe about the trimming. She went on a trimming [saga](https://github.com/zdellaert/LaserCoral/blob/main/code/notes/Trim_Saga.md) of her own for WGBS data. We talked about adjusting trimming parameters for flexbar. Specifically, I'm going to play around with `--adapter-trim-end`, `--adapter-trim-end`, and `--adapter-error-rate`. I'm going to try some different iterations on M6. 

Iterations that I will try today: 

- `--adapter-trim-end` to ANY - this will remove adapter sequence in any location of the read and the longer side of the read will remain after removal of overlap
- `--adapter-min-overlap` to 4 - this will allow more potential matches and removal 
- `--adapter-error-rate` to 0.2 
- `--htrim-right G` - trim any poly Gs on the right side of the read 

In the scripts folder: `nano flexbar_M6_20250217.sh`

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
#SBATCH --array=0-3

echo "Flexbar trimming iterations" $(date)

module load Flexbar/3.5.0-foss-2018b  
module load FastQC/0.11.9-Java-11

output_dir="/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250217"
mkdir -p $output_dir

# Define argument variations
arg_types=("trim_any" "trim_any_ao4" "trim_any_ae0.2" "trim_any_htrim_g")

cmds=(
"--adapter-trim-end ANY"
"--adapter-trim-end ANY --adapter-min-overlap 4"
"--adapter-trim-end ANY --adapter-error-rate 0.2"
"--adapter-trim-end ANY --htrim-right G"
)

# Get the trimming command for this array task
arg_type=${arg_types[$SLURM_ARRAY_TASK_ID]}
cmd=${cmds[$SLURM_ARRAY_TASK_ID]}

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapter-preset SmallRNA \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
 ${cmd} \
--zip-output GZ \
--threads 16 \
-t "${output_dir}/${arg_type}_M6"

# Run FastQC on the output
fastqc "${output_dir}/${arg_type}_M6.fastq.gz" -o $output_dir

echo "Finished processing ${arg_type}"
```

Submitted batch job 361795-98. No output though...Got this in the error files: `Adapter trim-end should be RIGHT for adapter presets.` Editing the script to use the illumina fasta sequences that I have on the server already. These should be the exact same as the presets. Submitted batch job 361799_*. Large output to one file: XXX. But not fastq or fastqc files produced. Got this error: `/var/spool/slurmd/job361800/slurm_script: line 39: -a: command not found`. Going to comment out the fastqc portion and change the `-a` to `--adapters`. Submitted batch job 361809_0-3. This is annoying. Going to comment out the cmds and arg lines for now and just run them one at a time. 

In the scripts folder: `nano flexbar_M6_20250217.sh`

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

echo "Flexbar trimming iterations" $(date)

module load Flexbar/3.5.0-foss-2018b  
module load FastQC/0.11.9-Java-11

output_dir="/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250217"

echo "Flexbar trimming iterations - trim ANY" $(date)

# Define argument variations
arg_type=("trim_any")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
#--adapter-preset SmallRNA \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--adapter-trim-end ANY \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
# ${cmd} \
--zip-output GZ \
--threads 16 \
-t "${output_dir}/${arg_type}_M6"

# Run FastQC on the output
fastqc "${output_dir}/${arg_type}_M6.fastq.gz" -o $output_dir

echo "Flexbar trimming iterations - trim ANY with min overlap of 4" $(date)

# Define argument variations
arg_type=("trim_any_ao4")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
#--adapter-preset SmallRNA \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--adapter-trim-end ANY \
--adapter-min-overlap 4 \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
# ${cmd} \
--zip-output GZ \
--threads 16 \
-t "${output_dir}/${arg_type}_M6"

# Run FastQC on the output
fastqc "${output_dir}/${arg_type}_M6.fastq.gz" -o $output_dir

echo "Flexbar trimming iterations - trim ANY with error rate of 0.2" $(date)

# Define argument variations
arg_type=("trim_any_ae0.2")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
#--adapter-preset SmallRNA \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--adapter-trim-end ANY \
--adapter-error-rate 0.2 \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
# ${cmd} \
--zip-output GZ \
--threads 16 \
-t "${output_dir}/${arg_type}_M6"

# Run FastQC on the output
fastqc "${output_dir}/${arg_type}_M6.fastq.gz" -o $output_dir

echo "Flexbar trimming iterations - trim ANY with polyG trimming right" $(date)

# Define argument variations
arg_type=("trim_any_htrim_g")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
#--adapter-preset SmallRNA \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--adapter-trim-end ANY \
--htrim-right G \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
# ${cmd} \
--zip-output GZ \
--threads 16 \
-t "${output_dir}/${arg_type}_M6"

# Run FastQC on the output
fastqc "${output_dir}/${arg_type}_M6.fastq.gz" -o $output_dir

echo "Flexbar trimming iterations complete" $(date)
```

Submitted batch job 361815. Okay reads are still getting written out to `flexbarOut.fastq` so I am stopping this job. Going to edit the code so that it is `--target` instead of `-t`. Actually first I'm going to run in interactive mode to see what is going on. Lol there was an s at the end of `arg_type` when I set the variable but not when I wrote it to the output. Remove the s, set `--target` and rerun. Submitted batch job 361821. STILL writing to flexbarOut.fastq!!!!! Letting it run for 5 mins and then checking. Okay cancelling. Trying to run ONE iteration. In the scripts folder: `nano flexbar_M6_trim_any.sh`

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

echo "Flexbar trimming iterations" $(date)

module load Flexbar/3.5.0-foss-2018b  
module load FastQC/0.11.9-Java-11

#output_dir="/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250217"

echo "Flexbar trimming iterations - trim ANY" $(date)

# Define argument variations
arg_type=("trim_any")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
#--adapter-preset SmallRNA \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--adapter-trim-end ANY \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
# ${cmd} \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250217/${arg_type}_M6

# Run FastQC on the output
fastqc "${output_dir}/${arg_type}_M6.fastq.gz" -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250217/
```

Submitted batch job 361822. STILL WRITING OUT TO THE FLEXBAR FILE. This is annoying. Got this error: 

```
/var/spool/slurmd/job361822/slurm_script: line 28: --adapters: command not found
/var/spool/slurmd/job361822/slurm_script: line 33: --zip-output: command not found
Skipping '/trim_any_M6.fastq.gz' which didn't exist, or couldn't be read
```

???? Why do you hate me flexbar? Going into interactive mode and running this: 

```
flexbar -r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz --adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta --adapter-trim-end ANY --qtrim-format i1.8 --qtrim-threshold 25 --zip-output GZ
```

Ran for a little then stopped it. Appears to be working correctly. Zoe recommended I check to make sure there are no weird spaces and to remove the commented out lines from my script (ie adapter preset and cmd). Submitted batch job 361824. Okay this seems to be working! Cancelling this job so I can add all of the iterations to this script. Submitted batch job 361825. Ran in 6 mins! But only the first iteration ran. Going to change the `--target` so that the full path is there for all iterations. Submitted batch job 361826. Download fastqc files from here (`/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250217`) onto local computer (see qc files [here](https://github.com/JillAshey/DevelopmentalTimeseries/tree/main/data/Molecular/smRNA/flexbar_iterations)). The one that looks the best is the `trim_any_htrim_g` iteration. Here, the vast majority of the adapter content is removed and the polyGs are gone. There are still some adapter overrepresented sequences but at a much lower percentage than before. My next step is to redo this trimming iteration with increased numbers of mismatches

### 20250218

More trimming today. First, going to try to set the `--adapter-min-overlap` to 3 and 2 to see what happens in conjunction with `--adapter-trim-end ANY` and `--htrim-right G`. In the scripts folder: `nano flexbar_M6_20250218_part1.sh`

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

echo "Flexbar trimming iterations" $(date)

module load Flexbar/3.5.0-foss-2018b  
module load FastQC/0.11.9-Java-11

mkdir /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218
output_dir = /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218

echo "Flexbar trimming iterations - trim ANY with polyG trimming right and min overlap of 3" $(date)

# Define argument variations
arg_type=("trim_any_htrim_g_ao3")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--adapter-trim-end ANY \
--htrim-right G \
--adapter-min-overlap 3 \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6

# Run FastQC on the output
fastqc "${output_dir}/${arg_type}_M6.fastq.gz" -o $output_dir

echo "Flexbar trimming iterations - trim ANY with polyG trimming right and min overlap of 2" $(date)

# Define argument variations
arg_type=("trim_any_htrim_g_ao2")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--adapter-trim-end ANY \
--htrim-right G \
--adapter-min-overlap 2 \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6

# Run FastQC on the output
fastqc "${output_dir}/${arg_type}_M6.fastq.gz" -o $output_dir

echo "Flexbar trimming iterations complete" $(date)
```

Submitted batch job 361892. Ran in about 6 mins. This is what the log files look like: 

```
# trim_any_htrim_g_ao3
Filtering statistics
====================
Processed reads                   24211196
  skipped due to uncalled bases     465454
  short prior to adapter removal         0
  finally skipped short reads     19978118
Discarded reads overall           20443572
Remaining reads                    3767624   (15%)

Processed bases   3655890596
Remaining bases    127805448   (3% of input)

# trim_any_htrim_g_ao2
Filtering statistics
====================
Processed reads                   24211196
  skipped due to uncalled bases     465454
  short prior to adapter removal         0
  finally skipped short reads     19978796
Discarded reads overall           20444250
Remaining reads                    3766946   (15%)

Processed bases   3655890596
Remaining bases    127020871   (3% of input)
```

Pretty similar results. Compared to the only htrim from yesterday: 

```
# trim_any_htrim_g
Filtering statistics
====================
Processed reads                   24211196
  skipped due to uncalled bases     465454
  short prior to adapter removal         0
  finally skipped short reads     19978118
Discarded reads overall           20443572
Remaining reads                    3767624   (15%)

Processed bases   3655890596
Remaining bases    127805448   (3% of input)
```

The `trim_any_htrim_g` and `trim_any_htrim_g_ao3` appear to be identical. Looking at the fastqc, I am still getting some adapter content in the overrepresented sequences. Let's try increasing the error rate. In the scripts folder: `nano flexbar_M6_20250218_part2.sh`

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

echo "Flexbar trimming iterations" $(date)

module load Flexbar/3.5.0-foss-2018b  
module load FastQC/0.11.9-Java-11

echo "Flexbar trimming iterations - trim ANY with polyG trimming right, min overlap of 3, error rate of 0.3" $(date)

# Define argument variations
arg_type=("trim_any_htrim_g_ao3_ae0.3")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--adapter-trim-end ANY \
--htrim-right G \
--adapter-min-overlap 3 \
--adapter-error-rate 0.3 \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6

# Run FastQC on the output
fastqc "/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6.fastq.gz" -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/

echo "Flexbar trimming iterations - trim ANY with polyG trimming right, min overlap of 3, error rate of 0.5" $(date)

# Define argument variations
arg_type=("trim_any_htrim_g_ao3_ae0.5")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--adapter-trim-end ANY \
--htrim-right G \
--adapter-min-overlap 3 \
--adapter-error-rate 0.5 \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6

# Run FastQC on the output
fastqc "/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6.fastq.gz" -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/

echo "Flexbar trimming iterations complete" $(date)
```

Submitted batch job 361895. Let's look at the log files: 

```
# trim_any_htrim_g_ao3_ae0.3
Filtering statistics
====================
Processed reads                   24211196
  skipped due to uncalled bases     465454
  short prior to adapter removal         0
  finally skipped short reads     22995599
Discarded reads overall           23461053
Remaining reads                     750143   (3%)

Processed bases   3655890596
Remaining bases     29142835   (0% of input)

# trim_any_htrim_g_ao3_ae0.5
Filtering statistics
====================
Processed reads                   24211196
  skipped due to uncalled bases     465454
  short prior to adapter removal         0
  finally skipped short reads     23018985
Discarded reads overall           23484439
Remaining reads                     726757   (3%)

Processed bases   3655890596
Remaining bases     27824550   (0% of input)
```

When looking at the fastqc, it looks like I got rid of the adapter content but there are some sequences that are overrepresented, such as `CATGAGGCCATCTGGGGGGGGGGGTTTTTTTTT`. The QC output says no hit but it looks like it might be an adapter? Instead of messing around with the error rate (though might come back to this), let's try to add the M6 index primer sequences to the adapter fasta file. The primers for M6 are `CAAGCAGAAGACGGCATACGAGATGGCCTCATGTGACTGGAGTTCAGACGTGT` for i7 and `AATGATACGGCGACCACCGAGATCTACACGTTAATTGACACTCTTTCCCTACACGAC` for i5. 

In the scripts folder: `nano flexbar_M6_20250218_part3.sh`

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

echo "Flexbar trimming iterations" $(date)

module load Flexbar/3.5.0-foss-2018b  
module load FastQC/0.11.9-Java-11

echo "Flexbar trimming iterations - trim ANY with polyG trimming right and min overlap of 3" $(date)

# Define argument variations
arg_type=("trim_any_htrim_g_ao3_index_primers")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--adapter-trim-end ANY \
--htrim-right G \
--adapter-min-overlap 3 \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6

# Run FastQC on the output
fastqc "/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6.fastq.gz" -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/

echo "Flexbar trimming iterations - trim ANY with polyG trimming right and min overlap of 2" $(date)

# Define argument variations
arg_type=("trim_any_htrim_g_ao2_index_primers")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--adapter-trim-end ANY \
--htrim-right G \
--adapter-min-overlap 2 \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6

# Run FastQC on the output
fastqc "/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6.fastq.gz" -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/

echo "Flexbar trimming iterations complete" $(date)
```

Submitted batch job 361899. Looking at the QC, this didn't do anything to improve from the only polyG trimming. Also looking at the QC again, the decrease of adapter overlap to 2 and 3 did not improve the QC relative to the iteration where only polyGs were trimmed.  I'm not sure what adapter matches, mismatches, or gaps mean or if it would help me...

Let's try to include `-ac, --adapter-revcomp STRING - Include reverse complements of adapters. One of ON and ONLY.` and max length of 35. 

In the scripts folder: `nano flexbar_M6_20250218_part4.sh`

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

echo "Flexbar trimming iterations" $(date)

module load Flexbar/3.5.0-foss-2018b  
module load FastQC/0.11.9-Java-11

echo "Flexbar trimming iterations - trim ANY with polyG trimming right and reverse complement of adapters" $(date)

# Define argument variations
arg_type=("trim_any_htrim_g_reverse")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--adapter-trim-end ANY \
--htrim-right G \
--adapter-revcomp ON \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6

# Run FastQC on the output
fastqc "/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6.fastq.gz" -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/

echo "Flexbar trimming iterations - trim ANY with polyG trimming right and max length of 35bp" $(date)

# Define argument variations
arg_type=("trim_any_htrim_g_max35")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--adapter-trim-end ANY \
--htrim-right G \
--qtrim-format i1.8 \
--post-trim-length 35 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6

# Run FastQC on the output
fastqc "/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6.fastq.gz" -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/

echo "Flexbar trimming iterations complete" $(date)
```

Submitted batch job 361914. The reverse did not change anything. The trim down to max 35bp also did not change the adapter content in the overrepresented sequences. Going to try removing the `--adapter-trim-end ANY` so that it trims from the right (default) and increasing the error rate. 

In the scripts folder: `nano flexbar_M6_20250218_part5.sh`

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

echo "Flexbar trimming iterations" $(date)

module load Flexbar/3.5.0-foss-2018b  
module load FastQC/0.11.9-Java-11

echo "Flexbar trimming iterations - polyG trimming right and default error rate" $(date)

# Define argument variations
arg_type=("trim_htrim_g")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--htrim-right G \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6

# Run FastQC on the output
fastqc "/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6.fastq.gz" -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/

echo "Flexbar trimming iterations - polyG trimming right and error rate of 0.3" $(date)

# Define argument variations
arg_type=("trim_htrim_g_ae0.3")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--htrim-right G \
--adapter-error-rate 0.3 \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6

# Run FastQC on the output
fastqc "/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6.fastq.gz" -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/

echo "Flexbar trimming iterations - polyG trimming right and error rate of 0.5" $(date)

# Define argument variations
arg_type=("trim_htrim_g_ae0.5")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--htrim-right G \
--adapter-error-rate 0.5 \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6

# Run FastQC on the output
fastqc "/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6.fastq.gz" -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/
```

Submitted batch job 361949. Still got some adapter in there. Going to try the above script again but replace the illumina fasta with the preset adapters that come in the software. They should be the same, I just want to check. 

`nano flexbar_M6_20250218_part6.sh`

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

echo "Flexbar trimming iterations" $(date)

module load Flexbar/3.5.0-foss-2018b  
module load FastQC/0.11.9-Java-11

echo "Flexbar trimming iterations - polyG trimming right, default error rate, preset adapters" $(date)

# Define argument variations
arg_type=("trim_htrim_g_preset_adapters")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapter-preset SmallRNA \
--htrim-right G \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6

# Run FastQC on the output
fastqc "/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6.fastq.gz" -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/

echo "Flexbar trimming iterations - polyG trimming right, error rate of 0.3, preset adapters" $(date)

# Define argument variations
arg_type=("trim_htrim_g_ae0.3_preset_adapters")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapter-preset SmallRNA \
--htrim-right G \
--adapter-error-rate 0.3 \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6

# Run FastQC on the output
fastqc "/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6.fastq.gz" -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/

echo "Flexbar trimming iterations - polyG trimming right, error rate of 0.5, preset adapters" $(date)

# Define argument variations
arg_type=("trim_htrim_g_ae0.5_preset_adapters")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapter-preset SmallRNA \
--htrim-right G \
--adapter-error-rate 0.5 \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6

# Run FastQC on the output
fastqc "/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6.fastq.gz" -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/
```

Submitted batch job 361966. okay so no dice there. I think I need to start playing around with the alignment score...but I don't understand what that means. Here are the options in flexbar: 

```
-am, --adapter-match INTEGER
      Alignment match score. Default: 1.
-ai, --adapter-mismatch INTEGER
      Alignment mismatch score. Default: -1.
-ag, --adapter-gap INTEGER
      Alignment gap score. Default: -6.
```

In the manual, it says: "The scoring scheme can be adjusted separately for detection of barcodes and adapters. This includes alignment match, mismatch and gap scores, see program options page. For example, it could make sense to specify a larger score for gaps, when the data's sequencing platform has high indel error rates: `flexbar -r reads.fasta -a adapters.fasta --adapter-gap -4`. If the gap score is set to a value of -4, the score of a gap corresponds to 4 mismatches and can be compensated by 4 matches." Not sure where the program options page is...This is a nice simple illustration of what these things mean: 

![](https://miro.medium.com/v2/resize:fit:954/1*2Wh0jTmRhXLcJQsHhirviw.png)

It appears that alignment scores (commonly using the [Smith & Waterman algorithm](https://www.sciencedirect.com/science/article/pii/0022283681900875)) are computed in flexbar. Although not flexbar, this [page](https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/1010/index.php?manual=Alignment_scoring_match_thresholds.html) was helpful in understanding how these things are computed. In our case, a match in the adapter sequence gives a score of 1, a mismatch gives a score of -1, and a gap gives a score of -6. I don't know if there is a minimum alignment score that needs to happen for each adapter? Here are my hypotheses about what will happen if I change the following: 

- Increase `--adapter-match` to 5 - matches more favorable, more conservative trimming 
- Increase `--adapter-mismatch` from -1 to -0.5 - mismatches less costly, potentially leading to more adapter detection and trimming but potential for false positives
- Increase `--adapter-gap` from -6 to -3 - gaps less costly, potentially leading to more adapter detection and trimming but potential for false positives

I am going to increase `--adapter-mismatch` from -1 to -0.5, increase `--adapter-gap` from -6 to -3, and combine both. In the scripts folder: `nano flexbar_M6_20250218_part7.sh`

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

echo "Flexbar trimming iterations" $(date)

module load Flexbar/3.5.0-foss-2018b  
module load FastQC/0.11.9-Java-11

echo "Flexbar trimming iterations - polyG trimming right, adapter mismatch at 0" $(date)

# Define argument variations
arg_type=("trim_htrim_g_mismatch0")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--htrim-right G \
--adapter-mismatch 0 \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6

# Run FastQC on the output
fastqc "/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6.fastq.gz" -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/

echo "Flexbar trimming iterations - polyG trimming right, adapter gap at -3" $(date)

# Define argument variations
arg_type=("trim_htrim_g_gap3")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--htrim-right G \
--adapter-gap -3 \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6

# Run FastQC on the output
fastqc "/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6.fastq.gz" -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/

echo "Flexbar trimming iterations - polyG trimming right, mismatch at 0, adapter gap at -3" $(date)

# Define argument variations
arg_type=("trim_htrim_g_mismatch0_gap3")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--htrim-right G \
--adapter-mismatch 0 \
--adapter-gap -3 \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6

# Run FastQC on the output
fastqc "/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6.fastq.gz" -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/

##### ADD ANY back in for adapters

echo "Flexbar trimming iterations - any, polyG trimming right, adapter mismatch at 0" $(date)

# Define argument variations
arg_type=("trim_any_htrim_g_mismatch0")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--adapter-trim-end ANY \
--htrim-right G \
--adapter-mismatch 0 \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6

# Run FastQC on the output
fastqc "/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6.fastq.gz" -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/

echo "Flexbar trimming iterations - any, polyG trimming right, adapter gap at -3" $(date)

# Define argument variations
arg_type=("trim_any_htrim_g_gap3")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--adapter-trim-end ANY \
--htrim-right G \
--adapter-gap -3 \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6

# Run FastQC on the output
fastqc "/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6.fastq.gz" -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/

echo "Flexbar trimming iterations - any, polyG trimming right, mismatch at 0, adapter gap at -3" $(date)

# Define argument variations
arg_type=("trim_any_htrim_g_mismatch0_gap3")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--adapter-trim-end ANY \
--htrim-right G \
--adapter-mismatch 0 \
--adapter-gap -3 \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6

# Run FastQC on the output
fastqc "/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6.fastq.gz" -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/

```

Submitted batch job 361968. Turns out that I need to give the mismatch option an integer, not a decimal. Going to change to 0 and rerun. Submitted batch job 361969. Using ANY is definitely the way to go. I got about 2 million reads when using `--adapter-trim-end ANY --htrim-right G --adapter-gap -3` and 1% as the highest overrepresented sequences, which were TruSeq adapters. When I added `--adapter-gap -3` to the function, the number of reads retained increased to 3.8 million but higher precentage of overrepresented TruSeq adapters (2.1%). What about if I set the `--adapter-gap` to -2, -1, and 0?

In the scripts folder: `nano flexbar_M6_20250218_part8.sh`

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

echo "Flexbar trimming iterations" $(date)

module load Flexbar/3.5.0-foss-2018b  
module load FastQC/0.11.9-Java-11

echo "Flexbar trimming iterations - any, polyG trimming right, adapter gap at -2" $(date)

# Define argument variations
arg_type=("trim_any_htrim_g_gap2")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--adapter-trim-end ANY \
--htrim-right G \
--adapter-gap -2 \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6

# Run FastQC on the output
fastqc "/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6.fastq.gz" -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/

echo "Flexbar trimming iterations - any, polyG trimming right, adapter gap at -1" $(date)

# Define argument variations
arg_type=("trim_any_htrim_g_gap1")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--adapter-trim-end ANY \
--htrim-right G \
--adapter-gap -1 \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6

# Run FastQC on the output
fastqc "/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6.fastq.gz" -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/

echo "Flexbar trimming iterations - any, polyG trimming right, adapter gap at 0" $(date)

# Define argument variations
arg_type=("trim_any_htrim_g_gap0")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--adapter-trim-end ANY \
--htrim-right G \
--adapter-gap 0 \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6

# Run FastQC on the output
fastqc "/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6.fastq.gz" -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/
```

Submitted batch job 361980. Output for -2 and -1 gap similar to -3 gap but less sequences are retained (1.6 v 2.1 million, respectively). When gap was 0, more sequences were kept (4 million) but more adapter content remained in overrepresented sequences. I also saw some overrepresented seqs that had polyT tails...should I add this to trimming?

Let's evaluate what I have done so far: 

- I get down to ~3.7 million sequences when using `--adapter-trim-end ANY --htrim-right G`, `--adapter-trim-end ANY --htrim-right G --adapter-min-overlap 2`, `--adapter-trim-end ANY --htrim-right G --adapter-min-overlap 3`, `--adapter-trim-end ANY --htrim-right G --adapter-min-overlap 2 INDEX Primers`, `--adapter-trim-end ANY --htrim-right G --adapter-min-overlap 3 INDEX Primers`, `--adapter-trim-end ANY --htrim-right G --adapter-revcomp ON`, `--adapter-trim-end ANY --htrim-right G --post-trim-length 35`, and `--adapter-trim-end ANY --htrim-right G --adapter-mismatch 0 --adapter-gap -3`. For all of these iterations, I am still getting some residual adapter content (2.1% for highest overrepresented seq, which is typically a TruSeq Adapter index 1)
- Trimming with `--adapter-trim-end ANY` tends to remove more adapters than without it 
- Using the preset adapters does not work as well as the illumina adapters that I have in a fasta 
- The only time I got no hits was `--adapter-trim-end ANY --htrim-right G --adapter-min-overlap 3 --adapter-error rate 0.3` and `--adapter-trim-end ANY --htrim-right G --adapter-min-overlap 3 --adapter-error rate 0.5`. But this also left <1 million seqs and sequences with polyTs
- The lowest that I have been able to get the overrepresented sequences percentage is with `--adapter-trim-end ANY --htrim-right G --adapter-gap -3`, where it was 1%

I feel pretty good about the 1% adapter content, I think this might be as good as it gets. First, lets trim to 75 bp. In the scripts folder: `nano flexbar_M6_20250218_part9.sh`

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

echo "Flexbar trimming iterations" $(date)

module load Flexbar/3.5.0-foss-2018b  
module load FastQC/0.11.9-Java-11

echo "Flexbar trimming iterations - any, polyG trimming right, adapter gap at -3, 75bp max length" $(date)

# Define argument variations
arg_type=("trim_any_htrim_g_gap3_max75")

flexbar \
-r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--adapter-trim-end ANY \
--htrim-right G \
--post-trim-length 75 \
--adapter-gap -3 \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6

# Run FastQC on the output
fastqc "/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/${arg_type}_M6.fastq.gz" -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/

```

Submitted batch job 361992. Let's also do 50 max bp. Submitted batch job 361993. Okay LETS MAP!!!!!!!! Using the mirdeep2 mapper module. Make an output folder

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output
mkdir mirdeep2_trim_any_htrim_g_gap3_max50
```

In the scripts folder: `nano mapper_mirdeep2_trim_any_htrim_g_gap3_max50.sh`

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

echo "Map trimmed reads (trim_any_htrim_g_gap3_max50) with mirdeep2 mapper script" $(date)

conda activate /data/putnamlab/mirdeep2

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218

gunzip trim_any_htrim_g_gap3_max50_M6.fastq.gz

mapper.pl /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/trim_any_htrim_g_gap3_max50_M6.fastq -e -h -j -l 18 -m -q -p /data/putnamlab/jillashey/genome/Mcap/V3/Mcap_ref.btindex -s /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_trim_any_htrim_g_gap3_max50/mapped_reads.fa -t /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_trim_any_htrim_g_gap3_max50/mapped_reads_vs_genome.arf

echo "Mapping complete for trimmed reads (trim_any_htrim_g_gap3_max50)" $(date)

conda deactivate
```

Submitted batch job 361994. Got this error: `prefix 1:N:0:ATGAGGCCAT+CAATTAACGT does not contain exactly three alphabet letters`. Why do you care? Removed `-d` argument, which says that input file is config file (it usually is a config file but only running one sample rn). Submitted batch job 361995. 

```
Mapping statistics
#desc   total   mapped  unmapped        %mapped %unmapped
total: 2137429  15302   2122127 0.716   99.284
seq: 2137429    15302   2122127 0.716   99.284
```

Still very low...let's try running the mirdeep2 module, which does the miRNA predictions

In the scripts folder: `nano mirdeep2_trim_any_htrim_g_gap3_max50.sh`

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

echo "Starting mirdeep2 with trimmed reads (trim_any_htrim_g_gap3_max50)" $(date)

miRDeep2.pl /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_trim_any_htrim_g_gap3_max50/mapped_reads.fa /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.assembly.fasta /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_trim_any_htrim_g_gap3_max50/mapped_reads_vs_genome.arf /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/mature_mirbase_cnidarian_T.fa none none -P -v -g -1 2>report.log

echo "mirdeep2 on trimmed reads (trim_any_htrim_g_gap3_max50)" $(date)

conda deactivate
```

Submitted batch job 361997. Cool cool nothing was predicted. 

```
Mapping mature,star and loop sequences against index
# reads processed: 21
# reads with at least one reported alignment: 4 (19.05%)
# reads that failed to align: 17 (80.95%)
Reported 32 alignments to 1 output stream(s)
```

I want to also map using bowtie (though I think mirdeep2 uses bowtie). In the scripts folder: `nano bowtie_trim_any_htrim_g_gap3_max50.sh`

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

module load Bowtie/1.3.1-GCC-11.3.0
#other options if needed: Bowtie/1.2.3-GCC-8.3.0, Bowtie/1.3.1-GCC-11.2.0, Bowtie/1.2.2-foss-2018b 

echo "Use bowtie to align collapsed reads to genome for trim_any_htrim_g_gap3_max50_M6 fastq" $(date)

# all indices already built, no need to redo 

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/

bowtie -a --verbose --seedmms 3 -x /data/putnamlab/jillashey/genome/Mcap/V3/Mcap_ref.btindex -q trim_any_htrim_g_gap3_max50_M6.fastq -S aligned_genome_trim_any_htrim_g_gap3_max50_M6.sam

echo "Bowtie alignment complete for trim_any_htrim_g_gap3_max50_M6 fastq to genome" $(date)
```

Submitted batch job 362001. maybe also try running shortstack?? In the scripts folder: `nano shortstack_trim_any_htrim_g_gap3_max50.sh`

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

echo "Running short stack using trim_any_htrim_g_gap3_max50"

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/

# Load modules 
module load ShortStack/4.0.2-foss-2022a  
module load Kent_tools/442-GCC-11.3.0

# Run short stack
# Run short stack
ShortStack \
--genomefile /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.assembly.fasta \
--readfile trim_any_htrim_g_gap3_max50_M6.fastq \
--known_miRNAs /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/mature_mirbase_cnidarian_T.fa \
--threads 10 \
--dn_mirna

echo "Short stack complete!"
```

Submitted batch job 362003. Failed, no miRNAs IDed. I feel like I have so many things to do...where to start? okay lets trim the good batch of reads with the updating trimming parameters. I'm going to move the raw reads of the good batch into its own folder. 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw
mkdir first_batch
mv 13_S76* first_batch/
mv 23_S77_R* first_batch/
mv 35_S78_R* first_batch/
mv 52_S79_R* first_batch/
mv 60_S80_R* first_batch/
mv 72_S81_R* first_batch/
mv 85_S82_R* first_batch/
mv 9_S75_R* first_batch/
cd ../
mkdir flexbar_first_batch
```

Going to trim the first batch using the following: `--adapter-trim-end ANY --htrim-right G --adapter-gap -3`. In the scripts folder: `nano flexbar_trim_any_htrim_g_gap3_max75_first_batch.sh`

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

echo "Flexbar trimming iterations" $(date)

module load Flexbar/3.5.0-foss-2018b  
module load FastQC/0.11.9-Java-11

echo "Flexbar trimming first batch - any, polyG trimming right, adapter gap at -3, 75bp max length" $(date)

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/first_batch

# Make an array of sequences to trim 
array1=($(ls *R1_001.fastq.gz))

for i in ${array1[@]}; do
flexbar \
-r ${i} \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--adapter-trim-end ANY \
--htrim-right G \
--post-trim-length 75 \
--adapter-gap -3 \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max75_${i}
done

echo "Flexbar trimming complete, run QC" $(date)

for file in /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/*fastq.gz
do 
fastqc $file 
done
```

Submitted batch job 362019. Also going to run one that does the same thing as above but trims to 35bp max. Submitted batch job 362025

Now I am going to collapse some of the raw reads: M6 and M9. I will then use the collapsed reads to make a blast db which I can blast miR-100 against. I want to ground truth that the universal miRNA is in these reads. Use the miR-100 sequences from other species that I have obtained.

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/
nano mir100.fasta

>Apul mir100 short stack from deep dive expression
UCCCGUAGAUCCGAACUUGU
>Peve mir100 short stack from deep dive expression
UCCCGUAGAUCCGAACUUGU
>Ptuh mir100 short stack from deep dive expression
ACCCGUAGAUCCGAACUUGU
>AST mir100 short stack from AST 2021 experiment
ACCCGUAGAUCCGAACUUGU
```

In the scripts folder: `nano collapse_raw_blast_mir100_M6_M9.sh`

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

echo "unzip reads" $(date)

gunzip /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq.gz 
gunzip /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/first_batch/9_S75_R1_001.fastq.gz

echo "unzipping complete, collapse raw reads" $(date)

module load FASTX-Toolkit/0.0.14-GCC-9.3.0 

fastx_collapser -v -i /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/6_small_RNA_S1_R1_001.fastq -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/collapse_6_small_RNA_S1_R1_001.fasta

fastx_collapser -v -i /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/first_batch/9_S75_R1_001.fastq -o /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/first_batch/collapse_9_S75_R1_001.fasta

echo "collapsing complete, make blast dbs" $(date)

module purge 
module load BLAST+/2.9.0-iimpi-2019b

makeblastdb -in data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/collapse_6_small_RNA_S1_R1_001.fasta -out data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/collapse_6_raw -dbtype nucl

makeblastdb -in /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/first_batch/collapse_9_S75_R1_001.fasta -out /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/first_batch/collapse_9_raw -dbtype nucl

echo "blast dbs complete, blast mir100 against dbs" $(date)

blastn -query /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/mir100.fasta -db data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/collapse_6_raw -out data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/mir100_M6_raw_blast_results_tab.txt -outfmt 6 -max_target_seqs 3

blastn -query /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/mir100.fasta -db /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/first_batch/collapse_9_raw -out data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/first_batch/mir100_M9_raw_blast_results_tab.txt -outfmt 6 -max_target_seqs 3

wc -l data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/mir100_M6_raw_blast_results_tab.txt
wc -l data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/first_batch/mir100_M9_raw_blast_results_tab.txt

echo "Blast complete" $(date)
```

Submitted batch job 362020. Accidently deleted the output files so gotta run again. Submitted batch job 362027

### 20250219 

Okay so the things collapsed but it does not look like the dbs were made or anything. Commenting out the collapse lines and running again. Submitted batch job 362040. Being weird and making new files and folders within existing files and folders. Going into interactive mode and running. `mir100_M6_raw_blast_results_tab.txt` got 0 hits. `mir100_M9_raw_blast_results_tab.txt` also got 0 hits...

Using the trimmed reads from the first batch, run the mapper.pl portion of mirdeep2 with both the 35 and 75bp trimmed reads. Make config files for both trim options. 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/

nano config_35bp.txt
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max35_9_S75_R1$.fastq s09
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max35_13_S76_R1_001.fastq.gz.fastq s13
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max35_23_S77_R1_001.fastq.gz.fastq s23
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max35_35_S78_R1_001.fastq.gz.fastq s35
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max35_52_S79_R1_001.fastq.gz.fastq s52
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max35_60_S80_R1_001.fastq.gz.fastq s60
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max35_72_S81_R1_001.fastq.gz.fastq s72
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max35_85_S82_R1_001.fastq.gz.fastq s85

nano config_75bp.txt
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max75_9_S75_R1_001.fastq.fastq s09
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max75_13_S76_R1_001.fastq.gz.fastq s13
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max75_23_S77_R1_001.fastq.gz.fastq s23
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max75_35_S78_R1_001.fastq.gz.fastq s35
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max75_52_S79_R1_001.fastq.gz.fastq s52
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max75_60_S80_R1_001.fastq.gz.fastq s60
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max75_72_S81_R1_001.fastq.gz.fastq s72
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max75_85_S82_R1_001.fastq.gz.fastq s85
```

In the scripts folder: `nano mapper_mirdeep2_first_batch.sh`

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

#echo "Referece genome indexed!" $(date)
#echo "Unload unneeded packages and run mapper script for trimmed stringent reads" $(date)

#module unload module load GCCcore/11.3.0 
#module unload Bowtie/1.3.1-GCC-11.3.0

conda activate /data/putnamlab/mirdeep2

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/

gunzip *.fastq.gz

echo "Start mapping for 35bp trimmed reads from first batch" $(date)
mapper.pl config_35bp.txt -e -d -h -j -l 18 -m -p /data/putnamlab/jillashey/genome/Mcap/V3/Mcap_ref.btindex -s mapped_reads_35bp.fa -t mapped_reads_vs_genome_35bp.arf

echo "Start mapping for 75bp trimmed reads from first batch" $(date)
mapper.pl config_75bp.txt -e -d -h -j -l 18 -m -p /data/putnamlab/jillashey/genome/Mcap/V3/Mcap_ref.btindex -s mapped_reads_75bp.fa -t mapped_reads_vs_genome_75bp.arf

echo "Mapping complete for trimmed reads" $(date)

conda deactivate 
```

Submitted batch job 362044. Cancelling bowtie run (`362001`). Also try mapping the first batch of reads using short stack. In the scripts folder: `nano shortstack_first_batch.sh`

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

echo "Running short stack on trimmed miRNAs from first seq batch"

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/

# Load modules 
module load ShortStack/4.0.2-foss-2022a  
module load Kent_tools/442-GCC-11.3.0

echo "Running short stack on 35bp trimmed miRNAs from first seq batch"

# Run short stack
ShortStack \
--genomefile /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.assembly.fasta \
--readfile trim_any_htrim_g_gap3_max35_9_S75_R1$.fastq \
trim_any_htrim_g_gap3_max35_13_S76_R1_001.fastq.gz.fastq \
trim_any_htrim_g_gap3_max35_23_S77_R1_001.fastq.gz.fastq \
trim_any_htrim_g_gap3_max35_35_S78_R1_001.fastq.gz.fastq \
trim_any_htrim_g_gap3_max75_52_S79_R1_001.fastq.gz.fastq \
trim_any_htrim_g_gap3_max35_60_S80_R1_001.fastq.gz.fastq \
trim_any_htrim_g_gap3_max35_72_S81_R1_001.fastq.gz.fastq \
trim_any_htrim_g_gap3_max35_85_S82_R1_001.fastq.gz.fastq \
--known_miRNAs /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/mature_mirbase_cnidarian_T.fa \
--threads 10 \
--dn_mirna

echo "Short stack complete for 35bp first batch, run shortstack on 75bp trimmed reads from first seq batch"

# Run short stack
ShortStack \
--genomefile /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.assembly.fasta \
--readfile trim_any_htrim_g_gap3_max75_9_S75_R1$.fastq \
trim_any_htrim_g_gap3_max75_13_S76_R1_001.fastq.gz.fastq \
trim_any_htrim_g_gap3_max75_23_S77_R1_001.fastq.gz.fastq \
trim_any_htrim_g_gap3_max75_35_S78_R1_001.fastq.gz.fastq \
trim_any_htrim_g_gap3_max75_52_S79_R1_001.fastq.gz.fastq \
trim_any_htrim_g_gap3_max75_60_S80_R1_001.fastq.gz.fastq \
trim_any_htrim_g_gap3_max75_72_S81_R1_001.fastq.gz.fastq \
trim_any_htrim_g_gap3_max75_85_S82_R1_001.fastq.gz.fastq \
--known_miRNAs /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/mature_mirbase_cnidarian_T.fa \
--threads 10 \
--dn_mirna

echo "Shortstack complete!" $(date)
```

Submitted batch job 362047

Let's try to use BWA to align the trimmed M6 sample, as BWA (as well as STAR) can map partially aligned reads. In the scripts folder: `nano bwa_trim_any_htrim_g_gap3_max50_M6.sh`

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

module load BWA/0.7.17-foss-2018b
module load SAMtools/1.9-foss-2018b

echo "Index Mcap genome for BWA" $(date)

cd /data/putnamlab/jillashey/genome/Mcap/V3/

bwa index -p Mcap_bwa Montipora_capitata_HIv3.assembly.fasta

echo "index complete, BWA alignment with trim_any_htrim_g_gap3_max50_M6.fastq" $(date)

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/

bwa mem -t 8 -a -k 15 -T 20 -L 5,5 /data/putnamlab/jillashey/genome/Mcap/V3/Mcap_bwa trim_any_htrim_g_gap3_max50_M6.fastq | samtools view -b > bwa_trim_any_htrim_g_gap3_max50_M6.bam

echo "BWA alignment complete!" $(date)
```

Submitted batch job 362046. Trying to look at the bam file on IGB but it says an index is required. I downloaded the index but still not working. Calculate some alignment stats. 

```
module load SAMtools/1.9-foss-2018b

# Only mapped reads
samtools view -c -F 260 bwa_trim_any_htrim_g_gap3_max50_M6.bam 
107053

samtools stats bwa_trim_any_htrim_g_gap3_max50_M6.bam
# This file was produced by samtools stats (1.9+htslib-1.9) and can be plotted using plot-bamstats
# This file contains statistics for all reads.
# The command line was:  stats bwa_trim_any_htrim_g_gap3_max50_M6.bam
# CHK, Checksum	[2]Read Names	[3]Sequences	[4]Qualities
# CHK, CRC32 of reads which passed filtering followed by addition (32bit overflow)
CHK	00b80d6f	e489230f	340b6f95
# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.
SN	raw total sequences:	2137429
SN	filtered sequences:	0
SN	sequences:	2137429
SN	is sorted:	0
SN	1st fragments:	2137429
SN	last fragments:	0
SN	reads mapped:	104385
SN	reads mapped and paired:	0	# paired-end technology bit set + both mates mapped
SN	reads unmapped:	2033044
SN	reads properly paired:	0	# proper-pair bit set
SN	reads paired:	0	# paired-end technology bit set
SN	reads duplicated:	0	# PCR or optical duplicate bit set
SN	reads MQ0:	70943	# mapped and MQ=0
SN	reads QC failed:	0
SN	non-primary alignments:	635928
SN	total length:	69516885	# ignores clipping
SN	total first fragment length:	69516885	# ignores clipping
SN	total last fragment length:	0	# ignores clipping
SN	bases mapped:	4199002	# ignores clipping
SN	bases mapped (cigar):	2641333	# more accurate
SN	bases trimmed:	0
SN	bases duplicated:	0
SN	mismatches:	39905	# from NM fields
SN	error rate:	1.510790e-02	# mismatches / bases mapped (cigar)
SN	average length:	32
SN	average first fragment length:	33
SN	average last fragment length:	0
SN	maximum length:	50
SN	maximum first fragment length:	0
SN	maximum last fragment length:	0
SN	average quality:	34.6
SN	insert size average:	0.0
SN	insert size standard deviation:	0.0
SN	inward oriented pairs:	0
SN	outward oriented pairs:	0
SN	pairs with other orientation:	0
SN	pairs on different chromosomes:	0
SN	percentage of properly paired reads (%):	0.0
```

Okay this is a lot of good info. Only 104,385 reads were aligned to the genome. There are a high number of non-primary alignments (635928), meaning many reads are aligning to multiple places in the genome. I can maybe use htseq count to quantify?

```
# sort 
samtools sort bwa_trim_any_htrim_g_gap3_max50_M6.bam -o bwa_trim_any_htrim_g_gap3_max50_M6.sorted.bam

# bam to sam
samtools view -h bwa_trim_any_htrim_g_gap3_max50_M6.sorted.bam > bwa_trim_any_htrim_g_gap3_max50_M6.sorted.sam

module load HTSeq/0.11.2-foss-2018b-Python-3.6.6 
htseq-count -f bam -r pos -i ID -t transcript bwa_trim_any_htrim_g_gap3_max50_M6.sorted.bam /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.genes.gff3 > counts_50bp_M6.txt
```

Didn't work, only the transcript info is in there. Converted bam to sam and looked at it in IGB. Looks similar to the previous one (image above). 

Run mirdeep2 with 35bp first batch samples. In the scripts folder: `nano mirdeep2_trim_any_htrim_g_gap3_35bp.sh`

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

echo "Starting mirdeep2 with trimmed reads (trim_any_htrim_g_gap3_max35)" $(date)

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch

miRDeep2.pl mapped_reads_35bp.fa /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.assembly.fasta mapped_reads_vs_genome_35bp.arf /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/mature_mirbase_cnidarian_T.fa none none -P -v -g -1 2>report.log

echo "mirdeep2 on trimmed reads (trim_any_htrim_g_gap3_max35)" $(date)

conda deactivate
```

Submitted batch job 362086. No miRNAs IDed...

Run flexbar so it trims to 30bp. `nano flexbar_trim_any_htrim_g_gap3_max30_first_batch.sh`

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

echo "Flexbar trimming iterations" $(date)

module load Flexbar/3.5.0-foss-2018b  
module load FastQC/0.11.9-Java-11

echo "Flexbar trimming first batch - any, polyG trimming right, adapter gap at -3, 30bp max length" $(date)

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/first_batch

# Make an array of sequences to trim 
array1=($(ls *R1_001.fastq.gz))

for i in ${array1[@]}; do
flexbar \
-r ${i} \
--adapters /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/illumina_adapters.fasta \
--adapter-trim-end ANY \
--htrim-right G \
--post-trim-length 30 \
--adapter-gap -3 \
--qtrim-format i1.8 \
--qtrim-threshold 25 \
--zip-output GZ \
--threads 16 \
--target /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max30_${i}
done

echo "Flexbar trimming complete, run QC" $(date)

for file in /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max30*fastq.gz
do 
fastqc $file 
done
```

Submitted batch job 362093

On the successful run previously, I used cutadapt. Maybe I should go back to these reads and map...make a config file from the reads I want to use. 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent
nano config_trim_stringent.txt

/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_9_S75_R1_001.fastq s09
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_13_S76_R1_001.fastq s13
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_23_S77_R1_001.fastq s23
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_35_S78_R1_001.fastq s35
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_52_S79_R1_001.fastq s52
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_60_S80_R1_001.fastq s60
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_72_S81_R1_001.fastq s72
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/trim_stringent_cutadapt_85_S82_R1_001.fastq s85
```

In the scripts folder: `nano mapper_mirdeep2_trim_stringent_first_batch.sh`

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

#echo "Referece genome indexed!" $(date)
#echo "Unload unneeded packages and run mapper script for trimmed stringent reads" $(date)

#module unload module load GCCcore/11.3.0 
#module unload Bowtie/1.3.1-GCC-11.3.0

conda activate /data/putnamlab/mirdeep2

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent

#gunzip *.fastq.gz

echo "Start mapping for 30bp cutadapt trimmed reads from first batch" $(date)
mapper.pl config_trim_stringent.txt -e -d -h -j -l 18 -m -p /data/putnamlab/jillashey/genome/Mcap/V3/Mcap_ref.btindex -s mapped_reads_trim_stringent_first_batch.fa -t mapped_reads_vs_genome_trim_stringent_first_batch.arf

echo "Mapping complete for trimmed reads" $(date)

conda deactivate 
```

Submitted batch job 362089. Look at the output

```
# Slurm output error file 
#desc   total   mapped  unmapped        %mapped %unmapped
total: 107296044        8932283 98363761        8.325   91.675
s09: 12624542   833865  11790677        6.605   93.395
s13: 13274127   1065833 12208294        8.029   91.971
s23: 14909106   1540921 13368185        10.335  89.665
s35: 13900420   1330386 12570034        9.571   90.429
s52: 15889590   1626703 14262887        10.238  89.762
s60: 11482805   973997  10508808        8.482   91.518
s72: 13678022   917316  12760706        6.706   93.294
s85: 11537432   643262  10894170        5.575   94.425
```

Mapping pretty good! This is what I expected. Move the output to new folder in output directory. 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent
mkdir ../../output/mirdeep2_trim_stringent_first_batch
mv mappe* ../../output/mirdeep2_trim_stringent_first_batch
mv bowtie.log ../../output/mirdeep2_trim_stringent_first_batch/
```

Run mirdeep2 with the trim stringent first batch samples. In the scripts folder: `nano mirdeep2_trim_stringent_first_batch.sh`

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

echo "Starting mirdeep2 with cutadapt trimmed reads - first batch" $(date)

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_trim_stringent_first_batch

miRDeep2.pl mapped_reads_trim_stringent_first_batch.fa /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.assembly.fasta mapped_reads_vs_genome_trim_stringent_first_batch.arf /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/mature_mirbase_cnidarian_T.fa none none -P -v -g -1 2>report.log

echo "mirdeep2 on cutadapt trimmed reads - first batch" $(date)

conda deactivate
```

Submitted batch job 362092. 

Flexbar trimming to 30bp finished running. Make a config file as input for mirdeep2

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch
nano config_trim_any_htrim_g_gap3_max30.txt

/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max30_9_S75_R1_001.fastq.gz.fastq s09
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max30_13_S76_R1_001.fastq.gz.fastq s13
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max30_23_S77_R1_001.fastq.gz.fastq s23
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max30_35_S78_R1_001.fastq.gz.fastq s35
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max30_52_S79_R1_001.fastq.gz.fastq s52
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max30_60_S80_R1_001.fastq.gz.fastq s60
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max30_72_S81_R1_001.fastq.gz.fastq s72
/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch/trim_any_htrim_g_gap3_max30_85_S82_R1_001.fastq.gz.fastq s85
```

In the scripts folder: `nano mapper_mirdeep2_trim_any_htrim_g_gap3_max30_first_batch.sh`

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

#echo "Referece genome indexed!" $(date)
#echo "Unload unneeded packages and run mapper script for trimmed stringent reads" $(date)

#module unload module load GCCcore/11.3.0 
#module unload Bowtie/1.3.1-GCC-11.3.0

conda activate /data/putnamlab/mirdeep2

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_first_batch

gunzip *.fastq.gz

echo "Start mapping for 30bp cutadapt trimmed reads from first batch" $(date)
mapper.pl config_trim_any_htrim_g_gap3_max30.txt -e -d -h -j -l 18 -m -p /data/putnamlab/jillashey/genome/Mcap/V3/Mcap_ref.btindex -s mapped_reads_trim_any_htrim_g_gap3_max30_first_batch.fa -t mapped_reads_vs_genome_trim_any_htrim_g_gap3_max30_first_batch.arf

echo "Mapping complete for trimmed reads" $(date)

conda deactivate 
```

Submitted batch job 362118. 

```
Mapping statistics

#desc   total   mapped  unmapped        %mapped %unmapped
total: 184302709        32565   184270144       0.018   99.982
s09: 24398852   1860    24396992        0.008   99.992
s13: 22050120   4304    22045816        0.020   99.980
s23: 20245901   6886    20239015        0.034   99.966
s35: 22019933   4771    22015162        0.022   99.978
s52: 29206532   5201    29201331        0.018   99.982
s60: 24265740   5454    24260286        0.022   99.978
s72: 23100469   1754    23098715        0.008   99.992
s85: 19015162   2335    19012827        0.012   99.988
```

Maybe flexbar just sucks?????

### 20250223 

mirdeep2 using the trimmed stringent first batch reads finished running! Output is here: `/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_trim_stringent_first_batch`. Copy the csv and html file onto my local computer. Identified 57 putative miRNAs using [this script](https://github.com/JillAshey/DevelopmentalTimeseries/blob/main/scripts/miRNA_discovery.Rmd) - 50 novel, 7 known including miRNA-100.

```
# novel 
Montipora_capitata_HIv3___Scaffold_2_97469
Montipora_capitata_HIv3___Scaffold_8_483169
Montipora_capitata_HIv3___Scaffold_9_581562
Montipora_capitata_HIv3___Scaffold_9_618173
Montipora_capitata_HIv3___Scaffold_9_606515
Montipora_capitata_HIv3___Scaffold_10_624116
Montipora_capitata_HIv3___Scaffold_14_1001044
Montipora_capitata_HIv3___Scaffold_14_1001054
Montipora_capitata_HIv3___Scaffold_2_123667
Montipora_capitata_HIv3___Scaffold_8_571913
Montipora_capitata_HIv3___Scaffold_11_817971
Montipora_capitata_HIv3___Scaffold_10_721713
Montipora_capitata_HIv3___Scaffold_10_659884
Montipora_capitata_HIv3___Scaffold_151_1064446
Montipora_capitata_HIv3___Scaffold_14_1007097
Montipora_capitata_HIv3___Scaffold_5_306063
Montipora_capitata_HIv3___Scaffold_8_503090
Montipora_capitata_HIv3___Scaffold_8_547143
Montipora_capitata_HIv3___Scaffold_7_451889
Montipora_capitata_HIv3___Scaffold_4_288176
Montipora_capitata_HIv3___Scaffold_11_742252
Montipora_capitata_HIv3___Scaffold_6_396712
Montipora_capitata_HIv3___Scaffold_8_536899
Montipora_capitata_HIv3___Scaffold_14_989061
Montipora_capitata_HIv3___Scaffold_14_984317
Montipora_capitata_HIv3___Scaffold_9_611209
Montipora_capitata_HIv3___Scaffold_10_638801
Montipora_capitata_HIv3___Scaffold_6_415143
Montipora_capitata_HIv3___Scaffold_2_86747
Montipora_capitata_HIv3___Scaffold_5_319191
Montipora_capitata_HIv3___Scaffold_1_35067
Montipora_capitata_HIv3___Scaffold_4_229985
Montipora_capitata_HIv3___Scaffold_5_371097
Montipora_capitata_HIv3___Scaffold_2_63883
Montipora_capitata_HIv3___Scaffold_3_148198
Montipora_capitata_HIv3___Scaffold_9_599461
Montipora_capitata_HIv3___Scaffold_11_793245
Montipora_capitata_HIv3___Scaffold_8_503091
Montipora_capitata_HIv3___Scaffold_8_529775
Montipora_capitata_HIv3___Scaffold_11_787538
Montipora_capitata_HIv3___Scaffold_11_779141
Montipora_capitata_HIv3___Scaffold_11_807462
Montipora_capitata_HIv3___Scaffold_13_939196
Montipora_capitata_HIv3___Scaffold_11_832609
Montipora_capitata_HIv3___Scaffold_9_617400
Montipora_capitata_HIv3___Scaffold_6_434354
Montipora_capitata_HIv3___Scaffold_5_361076
Montipora_capitata_HIv3___Scaffold_186_1068044
Montipora_capitata_HIv3___Scaffold_2_66515
Montipora_capitata_HIv3___Scaffold_5_303930

# known
Montipora_capitata_HIv3___Scaffold_14_1007949
Montipora_capitata_HIv3___Scaffold_2_95308
Montipora_capitata_HIv3___Scaffold_8_558953
Montipora_capitata_HIv3___Scaffold_8_580085
Montipora_capitata_HIv3___Scaffold_14_981912
Montipora_capitata_HIv3___Scaffold_12_901357
Montipora_capitata_HIv3___Scaffold_1_4143
```

I now want to do 2 things - quantify the miRNAs and use miranda to assess miRNA binding to 3'UTR. Let's quantify first. In this folder (`/data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_trim_stringent_first_batch/mirna_results_19_02_2025_t_14_06_12`), there are bed and fasta files for the known and novel mature, star and precursor sequences. I need to filter them by the miRNAs that I identified. Make a text file with the novel and known names of putative miRNAs.

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_trim_stringent_first_batch/mirna_results_19_02_2025_t_14_06_12

nano putative_miRNA_list.txt
Montipora_capitata_HIv3___Scaffold_2_97469
Montipora_capitata_HIv3___Scaffold_8_483169
Montipora_capitata_HIv3___Scaffold_9_581562
Montipora_capitata_HIv3___Scaffold_9_618173
Montipora_capitata_HIv3___Scaffold_9_606515
Montipora_capitata_HIv3___Scaffold_10_624116
Montipora_capitata_HIv3___Scaffold_14_1001044
Montipora_capitata_HIv3___Scaffold_14_1001054
Montipora_capitata_HIv3___Scaffold_2_123667
Montipora_capitata_HIv3___Scaffold_8_571913
Montipora_capitata_HIv3___Scaffold_11_817971
Montipora_capitata_HIv3___Scaffold_10_721713
Montipora_capitata_HIv3___Scaffold_10_659884
Montipora_capitata_HIv3___Scaffold_151_1064446
Montipora_capitata_HIv3___Scaffold_14_1007097
Montipora_capitata_HIv3___Scaffold_5_306063
Montipora_capitata_HIv3___Scaffold_8_503090
Montipora_capitata_HIv3___Scaffold_8_547143
Montipora_capitata_HIv3___Scaffold_7_451889
Montipora_capitata_HIv3___Scaffold_4_288176
Montipora_capitata_HIv3___Scaffold_11_742252
Montipora_capitata_HIv3___Scaffold_6_396712
Montipora_capitata_HIv3___Scaffold_8_536899
Montipora_capitata_HIv3___Scaffold_14_989061
Montipora_capitata_HIv3___Scaffold_14_984317
Montipora_capitata_HIv3___Scaffold_9_611209
Montipora_capitata_HIv3___Scaffold_10_638801
Montipora_capitata_HIv3___Scaffold_6_415143
Montipora_capitata_HIv3___Scaffold_2_86747
Montipora_capitata_HIv3___Scaffold_5_319191
Montipora_capitata_HIv3___Scaffold_1_35067
Montipora_capitata_HIv3___Scaffold_4_229985
Montipora_capitata_HIv3___Scaffold_5_371097
Montipora_capitata_HIv3___Scaffold_2_63883
Montipora_capitata_HIv3___Scaffold_3_148198
Montipora_capitata_HIv3___Scaffold_9_599461
Montipora_capitata_HIv3___Scaffold_11_793245
Montipora_capitata_HIv3___Scaffold_8_503091
Montipora_capitata_HIv3___Scaffold_8_529775
Montipora_capitata_HIv3___Scaffold_11_787538
Montipora_capitata_HIv3___Scaffold_11_779141
Montipora_capitata_HIv3___Scaffold_11_807462
Montipora_capitata_HIv3___Scaffold_13_939196
Montipora_capitata_HIv3___Scaffold_11_832609
Montipora_capitata_HIv3___Scaffold_9_617400
Montipora_capitata_HIv3___Scaffold_6_434354
Montipora_capitata_HIv3___Scaffold_5_361076
Montipora_capitata_HIv3___Scaffold_186_1068044
Montipora_capitata_HIv3___Scaffold_2_66515
Montipora_capitata_HIv3___Scaffold_5_303930
Montipora_capitata_HIv3___Scaffold_14_1007949
Montipora_capitata_HIv3___Scaffold_2_95308
Montipora_capitata_HIv3___Scaffold_8_558953
Montipora_capitata_HIv3___Scaffold_8_580085
Montipora_capitata_HIv3___Scaffold_14_981912
Montipora_capitata_HIv3___Scaffold_12_901357
Montipora_capitata_HIv3___Scaffold_1_4143
```

Cat the known and novel miRNA fasta together and filter fasta so that I keep only the miRNAs in `putative_miRNA_list.txt`.

```
cat novel_mature_19_02_2025_t_14_06_12_score-50_to_na.fa known_mature_19_02_2025_t_14_06_12_score-50_to_na.fa > putative_miRNAs.fa

grep -F -f putative_miRNA_list.txt putative_miRNAs.fa -A 1 --no-group-separator > putative_miRNAs_filt.fa

grep -c ">" putative_miRNAs_filt.fa 
57
```

Do the same with the precursor sequences. 

```
cat novel_pres_19_02_2025_t_14_06_12_score-50_to_na.fa known_pres_19_02_2025_t_14_06_12_score-50_to_na.fa > putative_precursors.fa

grep -F -f putative_miRNA_list.txt putative_precursors.fa -A 1 --no-group-separator > putative_precursors_filt.fa

grep -c ">" putative_precursors_filt.fa 
57
```

Similar to the mapper module, I can use a `config.txt` file to specify all the samples. I also need to use the `mapped_reads.fa` as the collapsed reads.

In the scripts folder: `nano quantifier_mirdeep2_trim_stringent_first_batch.sh`


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

echo "Quantifying smRNA counts for trim stringent first batch" $(date)

conda activate /data/putnamlab/mirdeep2

quantifier.pl -r /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_trim_stringent_first_batch/mapped_reads_trim_stringent_first_batch.fa -p /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_trim_stringent_first_batch/mirna_results_19_02_2025_t_14_06_12/putative_precursors_filt.fa -m /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_trim_stringent_first_batch/mirna_results_19_02_2025_t_14_06_12/putative_miRNAs_filt.fa -c /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/trim_stringent/config_trim_stringent.txt

echo "Quantifying complete for trim stringent first batch!" $(date)

conda deactivate
```

Submitted batch job 362306. Output: 

```
#desc   total   mapped  unmapped        %mapped %unmapped
total: 107296044        39310   107256734       0.037   99.963
s09: 12624542   1669    12622873        0.013   99.987
s13: 13274127   2281    13271846        0.017   99.983
s23: 14909106   3916    14905190        0.026   99.974
s35: 13900420   5233    13895187        0.038   99.962
s52: 15889590   6293    15883297        0.040   99.960
s60: 11482805   6215    11476590        0.054   99.946
s72: 13678022   7534    13670488        0.055   99.945
s85: 11537432   6169    11531263        0.053   99.947
```

Kinda low but not super unusual. Move the output, which is currently in the scripts folder to `mirdeep2_trim_stringent_first_batch`.

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/scripts
mv miRNAs_expressed_all_samples_1740330261
cd expression_analyses
mv expression_analyses_1740330261 ../../output/mirdeep2_trim_stringent_first_batch/
```

Copy output files to local computer as well. Run miranda on this miRNA data as well. Make miranda output folder. 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/
mkdir miranda_trim_stringent_first_batch
```

In the scripts folder: `nano miranda_strict_all_1kb_mcap_trim_stringent_first_batch.sh`

```
#!/bin/bash -i
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

echo "Mcap starting miranda run with all genes and miRNAs with energy cutoff <-20 and strict binding invoked"$(date)
echo "miRNAs generated from cutadapt trim stringent first batch"$(date)


module load Miniconda3/4.9.2
conda activate /data/putnamlab/conda/miranda 

miranda /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_trim_stringent_first_batch/mirna_results_19_02_2025_t_14_06_12/putative_miRNAs_filt.fa /data/putnamlab/jillashey/genome/Mcap/V3/Mcap_3UTR_1kb.fasta -en -20 -strict -out /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/miranda_trim_stringent_first_batch/miranda_strict_all_1kb_mcap_trim_stringent_first_batch.tab

conda deactivate

echo "miranda run finished!"$(date)
echo "counting number of interactions possible" $(date)

zgrep -c "Performing Scan" /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/miranda_trim_stringent_first_batch/miranda_strict_all_1kb_mcap_trim_stringent_first_batch.tab

echo "Parsing output" $(date)
grep -A 1 "Scores for this hit:" /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/miranda_trim_stringent_first_batch/miranda_strict_all_1kb_mcap_trim_stringent_first_batch.tab | sort | grep '>' > /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/miranda_trim_stringent_first_batch/miranda_strict_all_1kb_mcap_trim_stringent_first_batch_parsed.txt

echo "counting number of putative interactions predicted" $(date)
wc -l /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/miranda_trim_stringent_first_batch/miranda_strict_all_1kb_mcap_trim_stringent_first_batch_parsed.txt

echo "Mcap DT miranda script complete" $(date)
```

Submitted batch job 362307. Finished running in about 45 mins. 

```
counting number of interactions possible Sun Feb 23 13:16:24 EST 2025
3143436
Parsing output Sun Feb 23 13:16:40 EST 2025
counting number of putative interactions predicted Sun Feb 23 13:16:42 EST 2025
24162 /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/miranda_trim_stringent_first_batch/miranda_strict_all_1kb_mcap_trim_stringent_first_batch_parsed.txt
```

Lots of potential interactions! Copy miranda text file to local computer. 

### 20250304

Run miranda on coding sequences as well as 3UTR. 

In the scripts folder: `nano miranda_strict_all_cds_mcap_trim_stringent_first_batch.sh`

```
#!/bin/bash -i
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

echo "Mcap starting miranda run with all genes (CDS sequences) and miRNAs with energy cutoff <-20 and strict binding invoked"$(date)
echo "miRNAs generated from cutadapt trim stringent first batch"$(date)


module load Miniconda3/4.9.2
conda activate /data/putnamlab/conda/miranda 

miranda /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/mirdeep2_trim_stringent_first_batch/mirna_results_19_02_2025_t_14_06_12/putative_miRNAs_filt.fa /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.genes.cds.fna -en -20 -strict -out /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/miranda_trim_stringent_first_batch/miranda_strict_all_cds_mcap_trim_stringent_first_batch.tab

conda deactivate

echo "miranda run finished!"$(date)
echo "counting number of interactions possible" $(date)

zgrep -c "Performing Scan" /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/miranda_trim_stringent_first_batch/miranda_strict_all_cds_mcap_trim_stringent_first_batch.tab

echo "Parsing output" $(date)
grep -A 1 "Scores for this hit:" /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/miranda_trim_stringent_first_batch/miranda_strict_all_cds_mcap_trim_stringent_first_batch.tab | sort | grep '>' > /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/miranda_trim_stringent_first_batch/miranda_strict_all_cds_mcap_trim_stringent_first_batch_parsed.txt

echo "counting number of putative interactions predicted" $(date)
wc -l /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/output/miranda_trim_stringent_first_batch/miranda_strict_all_cds_mcap_trim_stringent_first_batch_parsed.txt

echo "Mcap DT miranda script complete" $(date)
```

Submitted batch job 362621













to do - lncRNA mRNA blast or something 
















mapped using bwa and was trying to quantify 

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

module load kallisto/0.48.0-gompi-2022a
cd /data/putnamlab/jillashey/genome/Mcap/V3/
kallisto index -i Mcap_V3.idx Montipora_capitata_HIv3.assembly.fasta

cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/flexbar_iterations_20250218/
kallisto quant -i index_name.idx -o output_directory -b 100 --single -l 200 -s 20 read.fastq
```








Things to try / think about 

- Aligner - how is it scoring the alignments? what about mismatches 
- Blast mir100 against raw/trimmed reads - running on raw reads 2/18/25
- Run another bad sample to see if I can trim and align
- Run all good samples to move forward with n=1 - flexbar trimming 2/18/25
- Use bwa, star or mapper that aligns part of the read 
- Use the truseq single index [adapters](https://support-docs.illumina.com/SHARE/AdapterSequences/Content/SHARE/AdapterSeq/TruSeq/SingleIndexes.html)

