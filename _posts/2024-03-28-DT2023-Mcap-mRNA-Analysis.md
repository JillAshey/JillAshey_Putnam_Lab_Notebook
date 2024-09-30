---
layout: post
title: Developmental 2023 Timeseries mRNA analysis (initial)
date: '2024-03-28'
categories: Analysis
tags: [Bioinformatics, mRNA, Montipora capitata]
projects: Developmental Timseries 2023
---

## Developmental 2023 Timeseries mRNA analysis (initial)

These data came from my developmental timeseries experiment in 2023 with *Montipora capitata* in Hawaii. In this experiment, *Montipora capitata* embryos and larvae were exposed to ambient and heat stress over 72 hours from embryo to swimming larvae. Samples were collected at 8 time points: 1, 4, 9, 14, 22, 28, 48 and 72 hpf. The github for this project is [here](https://github.com/JillAshey/Hawaii_Developmental_TimeSeries). 

To see which time points we want to sequence further, we ran an intial library prep and sequencing run on 8 samples (n=1 from each time point at ambient). The library prep was done by me using the [Zymo Ribofree](https://www.zymoresearch.com/pages/ngs-library-prep-kits#page-2) library prep kit (see my notebook [posts](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook) for more information on library prep). The sequencing was done through Oklahoma Medical Research Foundation NGS Core, which I found on Genohub. Here is the sequencing information: 

- Instrument: Illumina NovaSeq X Plus - 10B - PE 150 Cycle
- Read length: 2 x 150bp (Paired End)
- Pricing unit: per lane
- Number of samples: 8 libraries
- Guaranteed number of pass filter PE reads/sample: 20M (10M in each direction)
- Deliverables: FastQ files uploaded to Genohub project bucket

Files were downloaded to this location on Andromeda: `/data/putnamlab/KITT/hputnam/20240328_Mcap_RNASeq_Devo/`. Time to analyze! I'm going to write my notebook in chronological order by date and then will reorganize once the workflow is complete.

### 20240328

Sequences were received from sequencer today! Make a new directory in my own folder on Andromeda for this project: 

```
cd /data/putnamlab/jillashey
mkdir DT_Mcap_2023
cd DT_Mcap_2023
mkdir data scripts output
cd data 
mkdir raw trim 
cd raw
```

Copy the files into the raw data folder

```
cp /data/putnamlab/KITT/hputnam/20240328_Mcap_RNASeq_Devo/* .
```

Check md and initial raw QC was already done by Hollie. Raw QC is [here](https://github.com/JillAshey/Hawaii_Developmental_TimeSeries/blob/main/2023/data/Molecular/mRNA/20240328__Mcap_RNASeq_Devo_fastqc_multiqc_report.html). Overall, QC looks pretty good. The phred scores are high across the board. There is a decent amount of duplication but this is not unusual in early life stage samples. There is quite a bit of adapter content present so that definitely needs to come out. I'm going to use fastp to trim the reads. In the scripts folder: `nano fastp_qc.sh`

```
#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/DT_Mcap_2023/scripts            
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

# Load modules 
module load fastp/0.19.7-foss-2018b
module load FastQC/0.11.8-Java-1.8
module load MultiQC/1.9-intel-2020a-Python-3.8.2

cd /data/putnamlab/jillashey/DT_Mcap_2023/data/raw/

echo "Count number of reads in each file" $(date)

zgrep -c "@LH00" *.gz > raw_read_count.txt

echo "Start read trimming with fastp, followed by QC" $(date)

# Make an array of sequences to trim 
array1=($(ls *R1_001.fastq.gz))

# fastp and fastqc loop 
for i in ${array1[@]}; do
    fastp --in1 ${i} \
        --in2 $(echo ${i}|sed s/_R1/_R2/)\
        --out1 /data/putnamlab/jillashey/DT_Mcap_2023/data/trim/trim.${i} \
        --out2 /data/putnamlab/jillashey/DT_Mcap_2023/data/trim/trim.$(echo ${i}|sed s/_R1/_R2/) \
        --detect_adapter_for_pe \
        --qualified_quality_phred 30 \
        --unqualified_percent_limit 10 \
        --length_required 100 \
        --cut_right cut_right_window_size 5 cut_right_mean_quality 20 \
	 fastqc /data/putnamlab/jillashey/DT_Mcap_2023/output/fastqc/trim/trim.${i} \
    fastqc /data/putnamlab/jillashey/DT_Mcap_2023/output/fastqc/trim/trim.$(echo ${i}|sed s/_R1/_R2/)
done

echo "Read trimming of adapters and QC complete. Starting multiqc" $(date)

cd /data/putnamlab/jillashey/DT_Mcap_2023/output/fastqc/trim
multiqc *

echo "MultiQC complete" $(date)
```

Submitted batch job 310216

### 20240328

For whatever reason, the QC portion of the script didn't run, so I'm going to run a QC trim script. In the scripts folder: `nano fastqc_trim.sh`

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/DT_Mcap_2023/scripts            
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

echo "Count number of reads in each trimmed file" $(date)

zgrep -c "@LH00" /data/putnamlab/jillashey/DT_Mcap_2023/data/trim/*.gz > /data/putnamlab/jillashey/DT_Mcap_2023/data/trim/trim_read_count.txt

echo "Start fastqc" $(date)

for file in /data/putnamlab/jillashey/DT_Mcap_2023/data/trim/*fastq.gz
do 
fastqc $file --outdir /data/putnamlab/jillashey/DT_Mcap_2023/output/fastqc/trim
done

echo "Fastqc complete, start multiqc" $(date)

cd /data/putnamlab/jillashey/DT_Mcap_2023/output/fastqc/trim

multiqc *

echo "multiqc complete!" $(date)
```

Submitted batch job 310233. Here are the raw read counts: 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/data/raw
less raw_read_count.txt
13_S31_R1_001.fastq.gz:25264563
13_S31_R2_001.fastq.gz:25264563
23_S32_R1_001.fastq.gz:36786710
23_S32_R2_001.fastq.gz:36786710
35_S33_R1_001.fastq.gz:27430001
35_S33_R2_001.fastq.gz:27430001
52_S34_R1_001.fastq.gz:41778167
52_S34_R2_001.fastq.gz:41778167
60_S35_R1_001.fastq.gz:43121328
60_S35_R2_001.fastq.gz:43121328
72_S36_R1_001.fastq.gz:33397282
72_S36_R2_001.fastq.gz:33397282
85_S37_R1_001.fastq.gz:24484998
85_S37_R2_001.fastq.gz:24484998
9_S30_R1_001.fastq.gz:30559603
9_S30_R2_001.fastq.gz:30559603
```

Here are the trimmed read counts: 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/data/raw
less trim_read_count.txt
/data/putnamlab/jillashey/DT_Mcap_2023/data/trim/trim.13_S31_R1_001.fastq.gz:17766477
/data/putnamlab/jillashey/DT_Mcap_2023/data/trim/trim.13_S31_R2_001.fastq.gz:17766477
/data/putnamlab/jillashey/DT_Mcap_2023/data/trim/trim.23_S32_R1_001.fastq.gz:24346211
/data/putnamlab/jillashey/DT_Mcap_2023/data/trim/trim.23_S32_R2_001.fastq.gz:24346211
/data/putnamlab/jillashey/DT_Mcap_2023/data/trim/trim.35_S33_R1_001.fastq.gz:17755579
/data/putnamlab/jillashey/DT_Mcap_2023/data/trim/trim.35_S33_R2_001.fastq.gz:17755579
/data/putnamlab/jillashey/DT_Mcap_2023/data/trim/trim.52_S34_R1_001.fastq.gz:26921564
/data/putnamlab/jillashey/DT_Mcap_2023/data/trim/trim.52_S34_R2_001.fastq.gz:26921564
/data/putnamlab/jillashey/DT_Mcap_2023/data/trim/trim.60_S35_R1_001.fastq.gz:27316402
/data/putnamlab/jillashey/DT_Mcap_2023/data/trim/trim.60_S35_R2_001.fastq.gz:27316402
/data/putnamlab/jillashey/DT_Mcap_2023/data/trim/trim.72_S36_R1_001.fastq.gz:20445985
/data/putnamlab/jillashey/DT_Mcap_2023/data/trim/trim.72_S36_R2_001.fastq.gz:20445985
/data/putnamlab/jillashey/DT_Mcap_2023/data/trim/trim.85_S37_R1_001.fastq.gz:16126033
/data/putnamlab/jillashey/DT_Mcap_2023/data/trim/trim.85_S37_R2_001.fastq.gz:16126033
/data/putnamlab/jillashey/DT_Mcap_2023/data/trim/trim.9_S30_R1_001.fastq.gz:21574873
/data/putnamlab/jillashey/DT_Mcap_2023/data/trim/trim.9_S30_R2_001.fastq.gz:21574873
```

A decent amount of reads were removed from the samples. Between 60-70% of reads were retained in these samples. This is slightly unusual to me, as the raw QC looked decent (esp in terms of high phred scores) and in other datasets, I usually retain 80-90% of the reads. This is the first time I'm using genohub/OK sequencing center so maybe thats it? There is also a high duplication rate, which may be contributing. The QC for trimmed reads ran in about an hour. Trimmed QC is [here](https://github.com/JillAshey/Hawaii_Developmental_TimeSeries/blob/main/2023/data/Molecular/mRNA/trim_multiqc_report.html). The adapter content is all gone which is great. Quality still looks good. In samples 13 and 30, there is a small spike in the GC normal distribution but overall, still looks normally distributed. 

Align reads to Mcap genome using hisat2. I'm using Mcap genome V3, which I obtained from the Rutgers [website](http://cyanophora.rutgers.edu/montipora/) on 3/28/24. In the scripts folder: `nano align.sh`

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/DT_Mcap_2023/scripts            
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

# load modules needed
module load HISAT2/2.2.1-foss-2019b #Alignment to reference genome: HISAT2
module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

cd /data/putnamlab/jillashey/DT_Mcap_2023/data/trim/

echo "Building genome reference" $(date)

# index the reference genome for Mcap output index to working directory
hisat2-build -f /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.assembly.fasta Mcap_ref
echo "Referece genome indexed. Starting alingment" $(date)

# This script exports alignments as bam files
# sorts the bam file because Stringtie takes a sorted file for input (--dta)
# removes the sam file because it is no longer needed

array=($(ls *_R1_001.fastq.gz)) # call the clean sequences - make an array to align

for i in ${array[@]}; do
        sample_name=`echo $i| awk -F [.] '{print $2}'`
        hisat2 -p 8 --rna-strandness RF --dta -x Mcap_ref -1 ${i} -2 $(echo ${i}|sed s/_R1/_R2/) -S ${sample_name}.sam
        samtools sort -@ 8 -o ${sample_name}.bam ${sample_name}.sam
                echo "${i} bam-ified!"
        rm ${sample_name}.sam
done

echo "Alignment complete!" $(date)
echo "View mapping percentages " $(date)

for i in *.bam; do
     echo "${i}" >> mapped_reads_counts_Mcap.txt
     samtools flagstat ${i} | grep "mapped (" >> mapped_reads_counts_Mcap.txt
 done
``` 

Submitted batch job 310234. Ran in 2 hours. Here's the mapping information: 

```
13_S31_R1_001.bam
17303971 + 0 mapped (45.92% : N/A)
23_S32_R1_001.bam
39709957 + 0 mapped (68.10% : N/A)
35_S33_R1_001.bam
27912358 + 0 mapped (66.19% : N/A)
52_S34_R1_001.bam
45865671 + 0 mapped (74.17% : N/A)
60_S35_R1_001.bam
46178976 + 0 mapped (73.88% : N/A)
72_S36_R1_001.bam
36113018 + 0 mapped (76.53% : N/A)
85_S37_R1_001.bam
32246414 + 0 mapped (79.69% : N/A)
9_S30_R1_001.bam
26154642 + 0 mapped (56.50% : N/A)
```

These mapping percentages look pretty good. I'm not surprised that sample 13 is so low. This sample is from 4 hpf where I messed up the sampling protocol and added too much salt water (with larvae) to shield. While there was RNA in the samples, the RNA appeared relatively degraded on the gel (see gel [here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-02-10-MiniprepPlus-DNA%3ARNA-extractions-McapLarvae-DT.md)). Move the bam files + mapped percentages text file to the hisat2 output folder. 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/data/trim
mv *bam ../../output/hisat2/
mv mapped_reads_counts_Mcap.txt ../../output/hisat2/
```

Assemble reads using stringtie. In the scripts folder: `nano assemble.sh`

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/DT_Mcap_2023/scripts             
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load StringTie/2.2.1-GCC-11.2.0

echo "Assembling transcripts using stringtie" $(date)

cd /data/putnamlab/jillashey/DT_Mcap_2023/output/hisat2

array1=($(ls *.bam)) #Make an array of sequences to assemble

for i in ${array1[@]}; do
    stringtie -p 8 -e -B -G /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.genes.gff3 -A ${i}.gene_abund.tab -o ${i}.gtf ${i}
done

echo "Assembly for each sample complete " $(date)
```

Submitted batch job 310259. Move the gtf and tab files to the stringtie output folder.

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/output/hisat2
mv *gtf ../stringtie/
mv *tab ../stringtie/
cd ../stringtie 
```

Make a list of gtfs 

```
ls *.gtf > gtf_list.txt
```

Merge gtfs into a single gtf with stringtie

```
module purge
module load StringTie/2.1.4-GCC-9.3.0

stringtie --merge -e -p 8 -G /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.genes.gff3 -o Mcap_merged.gtf gtf_list.txt #Merge GTFs 

wc -l Mcap_merged.gtf 
310417 Mcap_merged.gtf
```

Assess accuracy of merged assembly

```
module purge
module load GffCompare/0.12.1-GCCcore-8.3.0

gffcompare -r /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.genes.gff3 -G -o merged Mcap_merged.gtf #Compute the accuracy 

  54384 reference transcripts loaded.
  2 duplicate reference transcripts discarded.
  54384 query transfrags loaded.
```

Check out the merged file 

```
less merged.stats

#= Summary for dataset: Mcap_merged.gtf 
#     Query mRNAs :   54384 in   54185 loci  (36023 multi-exon transcripts)
#            (141 multi-transcript loci, ~1.0 transcripts per locus)
# Reference mRNAs :   54382 in   54185 loci  (36023 multi-exon)
# Super-loci w/ reference transcripts:    54185
#-----------------| Sensitivity | Precision  |
        Base level:   100.0     |   100.0    |
        Exon level:   100.0     |   100.0    |
      Intron level:   100.0     |   100.0    |
Intron chain level:   100.0     |   100.0    |
  Transcript level:   100.0     |   100.0    |
       Locus level:   100.0     |   100.0    |

     Matching intron chains:   36023
       Matching transcripts:   54363
              Matching loci:   54183

          Missed exons:       0/256029  (  0.0%)
           Novel exons:       0/256028  (  0.0%)
        Missed introns:       0/201643  (  0.0%)
         Novel introns:       0/201643  (  0.0%)
           Missed loci:       0/54185   (  0.0%)
            Novel loci:       0/54185   (  0.0%)

 Total union super-loci across all input datasets: 54185 
54384 out of 54384 consensus transcripts written in merged.annotated.gtf (0 discarded as redundant)
```

Make gtf list text file for gene count matrix creation 

```
for filename in *bam.gtf; do echo $filename $PWD/$filename; done > listGTF.txt
```

Download the prepDE.py script from [here](https://github.com/gpertea/stringtie/blob/master/prepDE.py) and put it in the stringtie output folder. Load python and compile the gene count matrix

```
module purge
module load Python/2.7.18-GCCcore-9.3.0

python prepDE.py -g Mcap_gene_count_matrix.csv -i listGTF.txt

wc -l Mcap_gene_count_matrix.csv 
54385 Mcap_gene_count_matrix.csv
```

Copy the matrices onto my local computer. Woohoo! 






QC for Mcap DT 2023 data that came back from sequencer 9/19/24

```
less 10_S12_R1_001.fastq.gz 
@LH00260:104:22GLL5LT4:2:1101:13405:1070 1:N:0:ACTAAGATAT+CCGCGGTTGT
CNTTGTGATCAGTCCCGGGTGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTAAGATATCTGGGGGGGCGTGTTTTTTTTTTGAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
+
I#IIIIII99999II9IIIIIIIIIII--IIII9IIIIIIIIIIIIIIIIIIII-I-II9-I99--9-I9I9----999I999-9------9-999--999999999999999IIIII9IIIIIIIIIIIIIIIIIIIIII9IIIIIIIII

less 10_S12_R2_001.fastq.gz 
@LH00260:104:22GLL5LT4:2:1101:13405:1070 2:N:0:ACTAAGATAT+CCGCGGTTGT
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
+
IIIIIIIII-IIIIIIIII9III999I9IIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

less 10_S28_R1_001.fastq.gz 
@LH00260:104:22GLL5LT4:3:1101:3039:1028 1:N:0:ACTAAGATAT+CNGCGGTTGT
CNCCCCTCTTGCCAGATACTTGTTGTCAGATGTAACAGCTAGGCACAGGACATGTCCCATGTGGCCCTTATGTCCATCATTAGCCTTGTGGCAACCCGCTATCTTTCCTTGTTTGCAGCCAGTTTCAATATCCCACTTGATAATAGAACAG
+
I#9I9IIIIIIIIIIIII-9IIIIIIIIIIIIII9IIIIIIIIIIIIIIIIIII-IIIII9I9II-III9I-IIIIII9IIIII99I-I9II9I--II--II999IIIIIIIIII9-III9I9IIIIIIIII9-I9I99IIIII-9-II99

less 10_S28_R2_001.fastq.gz
@LH00260:104:22GLL5LT4:3:1101:3039:1028 2:N:0:ACTAAGATAT+CNGCGGTTGT
TGCTGTCAGAGTGCTCAAGGGTCATCAGCTGTCTGTCACCTGTCTCTGTGTGTCCCCTGCTGACAAGTTTGTTTTCTCTGGGTCTAAGGCCTGTTCCATTATCAAGTGGGATATTGAAACTGGCTGCAAACAAGGAAAGATATCGGGTTGC
+
I9IIIIIIIIIIIIIIIIIII-9IIIIIIIIIIIIIIII9IIIII9IIII9IIIII-II9III9I9III9IIII99IIIII9IIIIIII-II-IIII99IIII9IIIIIIIIIIIIII9III9IIIIII9IIIIIIIIIIIIIII9IIIII
```

which is small RNA sequence and which is RNA sequence? 

nano fastqc_all_raw.sh

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/DT_Mcap_2023/test            
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

echo "Start fastqc" $(date)

cd /data/putnamlab/jillashey/DT_Mcap_2023/test/

for file in *fastq.gz
do 
fastqc $file
done

echo "Fastqc complete, start multiqc" $(date)

multiqc *

echo "multiqc complete!" $(date)
```

Submitted batch job 338507

### 20240924

```
md5sum *RNA*.fastq.gz
702e86f018e45af059dcf669dd9e9b69  10_RNAseq_S28_R1_001.fastq.gz
683c7346bc1540d5a9b683fcafe8d839  10_RNAseq_S28_R2_001.fastq.gz
6a8e0308d34a7fe8e1dd950093f2d8cc  10_small_RNA_S4_R1_001.fastq.gz
5df1e17c449a0b2ca1709d23d9218a56  10_small_RNA_S4_R2_001.fastq.gz
4f46b70db8b4fd7f4c5f7c6d238ff672  11_RNAseq_S29_R1_001.fastq.gz
2c9db533ddd6fd9d1c0bf317ef64f512  11_RNAseq_S29_R2_001.fastq.gz
cbc48eb41ca388fb271a134d1c91cf06  11_small_RNA_S5_R1_001.fastq.gz
8e38afa5a4fa8e7aba3a66673dda99de  11_small_RNA_S5_R2_001.fastq.gz
9c23113bfd69cef0307438ca7279b117  14_RNAseq_S30_R1_001.fastq.gz
d27a1c25717ed55a2f19c279ed856509  14_RNAseq_S30_R2_001.fastq.gz
bb71fd3c990ce6f43be33a9e64520a7a  14_small_RNA_S6_R1_001.fastq.gz
99463142ad98ae02877b22ab4f227657  14_small_RNA_S6_R2_001.fastq.gz
bf857b42e61e2acfe5683537ccaebc8f  24_RNAseq_S31_R1_001.fastq.gz
2caea6f8376cb1def3f4d56c4a34794f  24_RNAseq_S31_R2_001.fastq.gz
f633b41622ecb8e3c046183dfd991ed0  24_small_RNA_S7_R1_001.fastq.gz
b38dee47e573c9a3eadcf79462bfb539  24_small_RNA_S7_R2_001.fastq.gz
466170e4bcfd6212622f32c554d52310  26_RNAseq_S32_R1_001.fastq.gz
550f784d7010a90a7979261acc4c247b  26_RNAseq_S32_R2_001.fastq.gz
2fcfbaa6de7007ada3af07b428dbd808  26_small_RNA_S8_R1_001.fastq.gz
9580845c7ca925b374459502f02fcc88  26_small_RNA_S8_R2_001.fastq.gz
8c0f6b8c1eb6544f22837be614d986da  28_RNAseq_S33_R1_001.fastq.gz
d3287f2a30504dfc950003d205cb648f  28_RNAseq_S33_R2_001.fastq.gz
3a63dab0d402c9c7b16db1370217ffa6  28_small_RNA_S9_R1_001.fastq.gz
0c56c73926283b43e6869996211be968  28_small_RNA_S9_R2_001.fastq.gz
5ecd3e61c2a8ed327b8115ae1e3c39eb  36_RNAseq_S34_R1_001.fastq.gz
37213fa3c4293381337465cad6d8ac9b  36_RNAseq_S34_R2_001.fastq.gz
84d42290a59f7f0cd1ce0214fd2031c6  36_small_RNA_S10_R1_001.fastq.gz
2275f89f69dec7f9f46f7969d91dd211  36_small_RNA_S10_R2_001.fastq.gz
300af022d27f8c796ddae7ad6cfd09b5  37_RNAseq_S35_R1_001.fastq.gz
7f0aa8a5fad7f1ef4e3d0b362debd9be  37_RNAseq_S35_R2_001.fastq.gz
83bdde7d4200b56eb81a80c371ff8530  37_small_RNA_S11_R1_001.fastq.gz
7c2bd1fc45277927f41eb771844e96a1  37_small_RNA_S11_R2_001.fastq.gz
8b68ef85d8a3427ea28b9b26602f6f08  39_RNAseq_S36_R1_001.fastq.gz
f619aa13b71e0fda37028361296c6890  39_RNAseq_S36_R2_001.fastq.gz
9712e83fccdee9c882721ce04a207162  39_small_RNA_S12_R1_001.fastq.gz
e5ead367c5104dfbbe7cecdc41589cb7  39_small_RNA_S12_R2_001.fastq.gz
6539a11de0504b9cd27ca0e99105eb31  47_RNAseq_S37_R1_001.fastq.gz
bcb55aadf544806e64973c8787270433  47_RNAseq_S37_R2_001.fastq.gz
e4f5279eef748705acc52b768d1ea06b  47_small_RNA_S13_R1_001.fastq.gz
35f57befe5167b7badfec76ed76164f1  47_small_RNA_S13_R2_001.fastq.gz
deeabbc714346a2be395ffee05534954  48_RNAseq_S38_R1_001.fastq.gz
3e57052be925bcf99a98d4bb4e2d0ef0  48_RNAseq_S38_R2_001.fastq.gz
f49407dd83da7ad68d6a353b1208f615  48_small_RNA_S14_R1_001.fastq.gz
995f6c2dbae731b5ae27e12f77966a4c  48_small_RNA_S14_R2_001.fastq.gz
81dd42bd91c0489520d5d6f43e93f7a7  51_RNAseq_S39_R1_001.fastq.gz
97005a5febfcd778642d296cfcff9804  51_RNAseq_S39_R2_001.fastq.gz
df7108ac712799380dc588e461d8faec  51_small_RNA_S15_R1_001.fastq.gz
15dd7cef395d93ad66412dba8576e022  51_small_RNA_S15_R2_001.fastq.gz
a17dfb969e4a6ac0a2b56aa037835ebd  61_RNAseq_S40_R1_001.fastq.gz
ef4dd53f6703e684bb048c6564c354cf  61_RNAseq_S40_R2_001.fastq.gz
3445966f7329f286a98a522835e398f8  61_small_RNA_S16_R1_001.fastq.gz
641c4aab05aead0e6493d54c906fc93e  61_small_RNA_S16_R2_001.fastq.gz
06111c15e28d3f0f27ecccc62906720d  62_RNAseq_S41_R1_001.fastq.gz
daced4a41b537bfad3788acc74465a15  62_RNAseq_S41_R2_001.fastq.gz
54468dc20565a6c032b6d73e6469fb46  62_small_RNA_S17_R1_001.fastq.gz
e28992a3905a6f178391d2dedc1f822f  62_small_RNA_S17_R2_001.fastq.gz
1c27d12f33601d51ea106d8ed48979e4  63_RNAseq_S42_R1_001.fastq.gz
16c999fa3820f42463791dcb03ec742b  63_RNAseq_S42_R2_001.fastq.gz
3a72bf0f147ba3d76251b62d2ea3c1dc  63_small_RNA_S18_R1_001.fastq.gz
6c243256ab289457fda664e86fc37d76  63_small_RNA_S18_R2_001.fastq.gz
26573c077d8d8da839ae711cedf42169  6_RNAseq_S25_R1_001.fastq.gz
3b7192772b515e3a2042f2dbf93c49f7  6_RNAseq_S25_R2_001.fastq.gz
652a9bbfd2eea5f6386be36ecb3fe7f7  6_small_RNA_S1_R1_001.fastq.gz
032ec5fd341c19df7d459a9509b64530  6_small_RNA_S1_R2_001.fastq.gz
b919bc773e19bc85c474d9b4d186a034  73_RNAseq_S43_R1_001.fastq.gz
799498b20adb92837e50f5c895caa54e  73_RNAseq_S43_R2_001.fastq.gz
bbf88e6ae4e152a0df2f7aad92108235  73_small_RNA_S19_R1_001.fastq.gz
da5e7f96e85eb981aaa76c5352c1a49b  73_small_RNA_S19_R2_001.fastq.gz
c87909a4b6cab7c6e8c77d89d32c9c1c  74_RNAseq_S44_R1_001.fastq.gz
8df0ef1d222ddde86c4597f03b0460e6  74_RNAseq_S44_R2_001.fastq.gz
b922733a75397fc7fd813404b373495c  74_small_RNA_S20_R1_001.fastq.gz
4b23d7649be0cfa2717f378939c655b5  74_small_RNA_S20_R2_001.fastq.gz
6fef936167af9a4721e22e5d811cc736  75_RNAseq_S45_R1_001.fastq.gz
8d4eaa76e3c47a3912aea942a7aa1c4b  75_RNAseq_S45_R2_001.fastq.gz
155b7e19e5d7a0fd8741256b3439909e  75_small_RNA_S21_R1_001.fastq.gz
2585139012a94991a2d928c3d943181c  75_small_RNA_S21_R2_001.fastq.gz
3af0d57415887b98c7ba8664f4be2175  7_RNAseq_S26_R1_001.fastq.gz
e82390f5537be0bec7eac9a6494dd042  7_RNAseq_S26_R2_001.fastq.gz
53605664bacdbc1e37b3120e2ac1e715  7_small_RNA_S2_R1_001.fastq.gz
7e7f794900855f697e760346fd88eb4b  7_small_RNA_S2_R2_001.fastq.gz
1090d5de0a6f71add02897147e640f3f  86_RNAseq_S46_R1_001.fastq.gz
87d08fd60ca278fe9c146ed525964da7  86_RNAseq_S46_R2_001.fastq.gz
46f3c89fa284d6a66d7a20814fe1f0fe  86_small_RNA_S22_R1_001.fastq.gz
c1618b96765199bc4db850b6805f64df  86_small_RNA_S22_R2_001.fastq.gz
fffe1bb3b1f92327fab218cf707ebdaf  87_RNAseq_S47_R1_001.fastq.gz
74a26b03b2527e0f3c51ffc555ef2455  87_RNAseq_S47_R2_001.fastq.gz
1aa383be2213f2f7982e8d88eefafca0  87_small_RNA_S23_R1_001.fastq.gz
b3e1c3a9314aec52f69dd88c6f699d60  87_small_RNA_S23_R2_001.fastq.gz
b28bdacd246571e26e248c858c07721d  88_RNAseq_S48_R1_001.fastq.gz
8795588fb6d7b74b162735241750060d  88_RNAseq_S48_R2_001.fastq.gz
9b4562d4eb6f1d662b8a04d1293a9053  88_small_RNA_S24_R1_001.fastq.gz
9bb1259ee25a5d410a8fae5b7f007c9a  88_small_RNA_S24_R2_001.fastq.gz
70f7d1fb97f559af56814d7f4b3c02c8  8_RNAseq_S27_R1_001.fastq.gz
57cf892971851db9f64d51f24acaab6b  8_RNAseq_S27_R2_001.fastq.gz
2b43ef5206c66fba3499301bbf00d963  8_small_RNA_S3_R1_001.fastq.gz
f458050f088dcfda08219a4773f21092  8_small_RNA_S3_R2_001.fastq.gz
```

### 20240926 

nano fastqc_RNAseq_raw.sh

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/scripts            
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

echo "Start fastqc" $(date)

cd /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/data/raw

for file in *fastq.gz
do 
fastqc $file
done

echo "Fastqc complete, start multiqc" $(date)

multiqc *

echo "multiqc complete!" $(date)
```

Submitted batch job 340171

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

Submitted batch job 340170
