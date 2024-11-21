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
mkdir mRNA
cd mRNA
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

### 20240919

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

Submitted batch job 340171. Raw QC for RNAseq data [here](https://github.com/JillAshey/DevelopmentalTimeseries/blob/main/data/Molecular/mRNA/multiqc_report_RNA_raw.html). 

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


### 20241023 

Now that I have all my data and have QCed the data, I am going to trim. For the RNAseq data, there are definitely some batch effects between the sequencing runs. May need to remove the samples that were sequenced earlier (March 2024) downstream. 

Trim using fastp. In the scripts folder: `nano XXX`

```
#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/scripts            
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

# Load modules 
module load fastp/0.19.7-foss-2018b
module load FastQC/0.11.8-Java-1.8
module load MultiQC/1.9-intel-2020a-Python-3.8.2

cd /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/data/raw/

echo "Count number of reads in each file" $(date)

zgrep -c "@LH00" *.gz > raw_read_count.txt

echo "Start read trimming with fastp, followed by QC" $(date)

# Make an array of sequences to trim 
array1=($(ls *R1_001.fastq.gz))

# fastp and fastqc loop 
for i in ${array1[@]}; do
    fastp --in1 ${i} \
        --in2 $(echo ${i}|sed s/_R1/_R2/)\
        --out1 /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/data/trim/trim.${i} \
        --out2 /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/data/trim/trim.$(echo ${i}|sed s/_R1/_R2/) \
        --detect_adapter_for_pe \
        --qualified_quality_phred 30 \
        --unqualified_percent_limit 10 \
        --length_required 100 \
        --cut_right cut_right_window_size 5 cut_right_mean_quality 20 \
	 fastqc /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output/fastqc/trim/trim.${i} \
    fastqc /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output/fastqc/trim/trim.$(echo ${i}|sed s/_R1/_R2/)
done

echo "Read trimming of adapters and QC complete. Starting multiqc" $(date)

cd /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output/fastqc/trim
multiqc *

echo "MultiQC complete" $(date)

echo "Count number of reads in each file" $(date)

zgrep -c "@LH00" *.gz > trim_read_count.txt
```

Submitted batch job 344847. QC looks good! Adapter content gone and phred score high across bases. Trimmed QC for RNAseq data [here](https://github.com/JillAshey/DevelopmentalTimeseries/blob/main/data/Molecular/mRNA/multiqc_report_RNA_trim.html). 

Align to Mcap V3 genome. In the scripts folder: `nano align.sh`

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/scripts            
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

# load modules needed
module load HISAT2/2.2.1-foss-2019b #Alignment to reference genome: HISAT2
module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

cd /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/data/trim/

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

Submitted batch job 344964. Mapping percentages: 

```
10_RNAseq_S28_R1_001.bam
4239251 + 0 mapped (42.36% : N/A)
11_RNAseq_S29_R1_001.bam
5991800 + 0 mapped (33.98% : N/A)
13_S31_R1_001.bam
17304441 + 0 mapped (45.92% : N/A)
14_RNAseq_S30_R1_001.bam
8085710 + 0 mapped (45.42% : N/A)
23_S32_R1_001.bam
39709953 + 0 mapped (68.10% : N/A)
24_RNAseq_S31_R1_001.bam
13563672 + 0 mapped (62.72% : N/A)
26_RNAseq_S32_R1_001.bam
14010844 + 0 mapped (61.96% : N/A)
28_RNAseq_S33_R1_001.bam
11577688 + 0 mapped (64.26% : N/A)
35_S33_R1_001.bam
27912360 + 0 mapped (66.19% : N/A)
36_RNAseq_S34_R1_001.bam
15565022 + 0 mapped (64.14% : N/A)
37_RNAseq_S35_R1_001.bam
16439270 + 0 mapped (66.31% : N/A)
39_RNAseq_S36_R1_001.bam
15933330 + 0 mapped (70.04% : N/A)
47_RNAseq_S37_R1_001.bam
20838605 + 0 mapped (81.81% : N/A)
48_RNAseq_S38_R1_001.bam
20274902 + 0 mapped (80.73% : N/A)
51_RNAseq_S39_R1_001.bam
19029047 + 0 mapped (72.68% : N/A)
52_S34_R1_001.bam
45865662 + 0 mapped (74.17% : N/A)
60_S35_R1_001.bam
46178975 + 0 mapped (73.88% : N/A)
61_RNAseq_S40_R1_001.bam
12591636 + 0 mapped (63.30% : N/A)
62_RNAseq_S41_R1_001.bam
14910190 + 0 mapped (67.58% : N/A)
63_RNAseq_S42_R1_001.bam
20914310 + 0 mapped (79.37% : N/A)
6_RNAseq_S25_R1_001.bam
11665744 + 0 mapped (62.94% : N/A)
72_S36_R1_001.bam
36113020 + 0 mapped (76.53% : N/A)
73_RNAseq_S43_R1_001.bam
13295206 + 0 mapped (73.43% : N/A)
74_RNAseq_S44_R1_001.bam
13728145 + 0 mapped (68.90% : N/A)
75_RNAseq_S45_R1_001.bam
8505675 + 0 mapped (68.58% : N/A)
7_RNAseq_S26_R1_001.bam
11945532 + 0 mapped (58.78% : N/A)
85_S37_R1_001.bam
32246414 + 0 mapped (79.69% : N/A)
86_RNAseq_S46_R1_001.bam
15919552 + 0 mapped (75.01% : N/A)
87_RNAseq_S47_R1_001.bam
11527670 + 0 mapped (68.50% : N/A)
88_RNAseq_S48_R1_001.bam
13406820 + 0 mapped (71.74% : N/A)
8_RNAseq_S27_R1_001.bam
10389582 + 0 mapped (58.35% : N/A)
9_S30_R1_001.bam
26154642 + 0 mapped (56.50% : N/A)
```

Most samples ranged from 56-80% mapping, which is good. A couple of samples (10, 11, 13, 14) has <50% mapping. This is not unexpected, as these samples were all 4 hpf, which is a period of maternal/zygotic turnover and duplication. 4 hpf is also the timepoint where I messed up the sampling protocol and added too much salt water (with larvae) to shield. While there was RNA in the samples, the RNA appeared relatively degraded on the gel (see 13 in gel picture [here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-02-10-MiniprepPlus-DNA%3ARNA-extractions-McapLarvae-DT.md)). 

Move the bam files + mapped percentages text file to the hisat2 output folder. 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/data/trim
mv *bam ../../output/hisat2/
mv mapped_reads_counts_Mcap.txt ../../output/hisat2/
```

Assemble reads using stringtie using the updated gff from AH (see code [here](https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/scripts/bioinformatics/fix_gff_format.Rmd)). In the scripts folder: `nano assemble.sh`

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

module load StringTie/2.2.1-GCC-11.2.0

echo "Assembling transcripts using stringtie" $(date)

cd /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output/hisat2

array1=($(ls *.bam)) #Make an array of sequences to assemble

for i in ${array1[@]}; do
    stringtie -p 8 -e -B -G /data/putnamlab/ashuffmyer/mcap-2023-rnaseq/references/Montipora_capitata_HIv3.genes_fixed.gff3 -A ${i}.gene_abund.tab -o ${i}.gtf ${i}
done

echo "Assembly for each sample complete " $(date)
```

Submitted batch job 344979. Move the gtf and tab files to the stringtie output folder.

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output/hisat2
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
interactive 
module load StringTie/2.1.4-GCC-9.3.0

stringtie --merge -e -p 8 -G /data/putnamlab/ashuffmyer/mcap-2023-rnaseq/references/Montipora_capitata_HIv3.genes_fixed.gff3 -o Mcap_merged.gtf gtf_list.txt #Merge GTFs 

wc -l Mcap_merged.gtf 
310417 Mcap_merged.gtf
```

Assess accuracy of merged assembly with gffcompare 

```
module purge
module load GffCompare/0.12.1-GCCcore-8.3.0

gffcompare -r /data/putnamlab/ashuffmyer/mcap-2023-rnaseq/references/Montipora_capitata_HIv3.genes_fixed.gff3 -G -o merged Mcap_merged.gtf #Compute the accuracy 
  54384 reference transcripts loaded.
  2 duplicate reference transcripts discarded.
  54384 query transfrags loaded.
```

Look at merge stats 

```
less merged.stats

# gffcompare v0.12.1 | Command line was:
#gffcompare -r /data/putnamlab/ashuffmyer/mcap-2023-rnaseq/references/Montipora_capitata_HIv3.genes_fixed.gff3 -G -o merged Mcap_merged.gtf
#

#= Summary for dataset: Mcap_merged.gtf 
#     Query mRNAs :   54384 in   54185 loci  (36023 multi-exon transcripts)
#            (141 multi-transcript loci, ~1.0 transcripts per locus)
# Reference mRNAs :   54382 in   54185 loci  (36023 multi-exon)
# Super-loci w/ reference transcripts:    54185
#-----------------| Sensitivity | Precision  |
        Base level:   100.0     |   100.0    |
        Exon level:    99.9     |    99.9    |
      Intron level:   100.0     |   100.0    |
Intron chain level:   100.0     |   100.0    |
  Transcript level:    99.9     |    99.9    |
       Locus level:   100.0     |   100.0    |

     Matching intron chains:   36023
       Matching transcripts:   54309
              Matching loci:   54173

          Missed exons:       0/256029  (  0.0%)
           Novel exons:       0/256026  (  0.0%)
        Missed introns:       0/201643  (  0.0%)
         Novel introns:       0/201643  (  0.0%)
           Missed loci:       0/54185   (  0.0%)
            Novel loci:       0/54185   (  0.0%)

 Total union super-loci across all input datasets: 54185 
54384 out of 54384 consensus transcripts written in merged.annotated.gtf (0 discarded as redundant)

head merged.tracking 
TCONS_00000001	XLOC_000001	Montipora_capitata_HIv3___RNAseq.g4584.t1|Montipora_capitata_HIv3___RNAseq.g4584.t1	=	q1:MSTRG.4|Montipora_capitata_HIv3___RNAseq.g4584.t1|13|0.000000|0.000000|1.125630|1854
TCONS_00000002	XLOC_000002	Montipora_capitata_HIv3___RNAseq.g4586.t1|Montipora_capitata_HIv3___RNAseq.g4586.t1	=	q1:MSTRG.30|Montipora_capitata_HIv3___RNAseq.g4586.t1|22|0.000000|0.000000|0.144260|2124
TCONS_00000003	XLOC_000003	Montipora_capitata_HIv3___TS.g26272.t1|Montipora_capitata_HIv3___TS.g26272.t1	=	q1:MSTRG.30|Montipora_capitata_HIv3___TS.g26272.t1|1|0.000000|0.000000|0.073721|589
TCONS_00000004	XLOC_000004	Montipora_capitata_HIv3___TS.g26273.t1|Montipora_capitata_HIv3___TS.g26273.t1	=	q1:MSTRG.5|Montipora_capitata_HIv3___TS.g26273.t1|1|0.000000|0.000000|0.042288|612
TCONS_00000005	XLOC_000005	Montipora_capitata_HIv3___RNAseq.g4588.t1|Montipora_capitata_HIv3___RNAseq.g4588.t1	=	q1:MSTRG.6|Montipora_capitata_HIv3___RNAseq.g4588.t1|1|0.000000|0.000000|0.016101|1101
TCONS_00000006	XLOC_000006	Montipora_capitata_HIv3___RNAseq.g4589.t1|Montipora_capitata_HIv3___RNAseq.g4589.t1	=	q1:MSTRG.7|Montipora_capitata_HIv3___RNAseq.g4589.t1|92|0.000000|0.000000|0.033519|10527
TCONS_00000007	XLOC_000007	Montipora_capitata_HIv3___TS.g26276.t1|Montipora_capitata_HIv3___TS.g26276.t1	=	q1:MSTRG.8|Montipora_capitata_HIv3___TS.g26276.t1|1|0.000000|0.000000|0.191457|261
TCONS_00000008	XLOC_000008	Montipora_capitata_HIv3___TS.g26277.t1|Montipora_capitata_HIv3___TS.g26277.t1	=	q1:MSTRG.9|Montipora_capitata_HIv3___TS.g26277.t1|1|0.000000|0.000000|0.085318|273
TCONS_00000009	XLOC_000009	Montipora_capitata_HIv3___RNAseq.g4592.t1|Montipora_capitata_HIv3___RNAseq.g4592.t1	=	q1:MSTRG.12|Montipora_capitata_HIv3___RNAseq.g4592.t1|2|0.000000|0.000000|1.483990|267
TCONS_00000010	XLOC_000010	Montipora_capitata_HIv3___TS.g26284.t1|Montipora_capitata_HIv3___TS.g26284.t1	=	q1:MSTRG.19|Montipora_capitata_HIv3___TS.g26284.t1|6|0.000000|0.000000|0.017805|945
```

Looking at the gffcompare [documentation](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml), the tracking file matches up transcripts between samples. Includes the MSTRG ID? Interesting. Is this file matching up the MSTRG to the gene ids? 

Make gtf list text file for gene count matrix creation 

```
for filename in *bam.gtf; do echo $filename $PWD/$filename; done > listGTF.txt
```

Download the prepDE.py script from [here](https://github.com/gpertea/stringtie/blob/master/prepDE.py) and put it in the scripts folder (I already had it in here). Load python and compile the gene count matrix

```
module purge
module load Python/2.7.18-GCCcore-9.3.0

python /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/scripts/prepDE.py -g Mcap_gene_count_matrix.csv -i listGTF.txt

wc -l Mcap_gene_count_matrix.csv 
54385 Mcap_gene_count_matrix.csv
```

Copy to local computer. Hooray! 

### 20241108 

Identify lncRNAs using my [Astrangia code](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-01-16-Astrangia2021-lncRNA-Analysis.md) and Zach's [oyster code](https://github.com/zbengt/oyster-lnc/blob/main/code/10-count-matrices-DESeq2-final.Rmd) as reference. 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output/stringtie

awk '$3 == "transcript" && $1 !~ /^#/' merged.annotated.gtf | awk '($5 - $4 > 199) || ($4 - $5 > 199)' > Mcap_lncRNA_candidates.gtf

head Mcap_lncRNA_candidates.gtf 
Montipora_capitata_HIv3___Scaffold_1	StringTie	transcript	82397	95409	.	+	.	transcript_id "Montipora_capitata_HIv3___RNAseq.g4584.t1"; gene_id "MSTRG.4"; gene_name "Montipora_capitata_HIv3___RNAseq.g4584.t1"; xloc "XLOC_000001"; ref_gene_id "Montipora_capitata_HIv3___RNAseq.g4584.t1"; cmp_ref "Montipora_capitata_HIv3___RNAseq.g4584.t1"; class_code "="; tss_id "TSS1";
Montipora_capitata_HIv3___Scaffold_1	StringTie	transcript	109801	163388	.	+	.	transcript_id "Montipora_capitata_HIv3___RNAseq.g4586.t1"; gene_id "MSTRG.30"; gene_name "Montipora_capitata_HIv3___RNAseq.g4586.t1"; xloc "XLOC_000002"; ref_gene_id "Montipora_capitata_HIv3___RNAseq.g4586.t1"; cmp_ref "Montipora_capitata_HIv3___RNAseq.g4586.t1"; class_code "="; tss_id "TSS2";
Montipora_capitata_HIv3___Scaffold_1	StringTie	transcript	162568	163156	.	+	.	transcript_id "Montipora_capitata_HIv3___TS.g26272.t1"; gene_id "MSTRG.30"; gene_name "Montipora_capitata_HIv3___TS.g26272.t1"; xloc "XLOC_000003"; ref_gene_id "Montipora_capitata_HIv3___TS.g26272.t1"; cmp_ref "Montipora_capitata_HIv3___TS.g26272.t1"; class_code "="; tss_id "TSS3";
Montipora_capitata_HIv3___Scaffold_1	StringTie	transcript	169950	170561	.	+	.	transcript_id "Montipora_capitata_HIv3___TS.g26273.t1"; gene_id "MSTRG.5"; gene_name "Montipora_capitata_HIv3___TS.g26273.t1"; xloc "XLOC_000004"; ref_gene_id "Montipora_capitata_HIv3___TS.g26273.t1"; cmp_ref "Montipora_capitata_HIv3___TS.g26273.t1"; class_code "="; tss_id "TSS4";
Montipora_capitata_HIv3___Scaffold_1	StringTie	transcript	170982	172082	.	+	.	transcript_id "Montipora_capitata_HIv3___RNAseq.g4588.t1"; gene_id "MSTRG.6"; gene_name "Montipora_capitata_HIv3___RNAseq.g4588.t1"; xloc "XLOC_000005"; ref_gene_id "Montipora_capitata_HIv3___RNAseq.g4588.t1"; cmp_ref "Montipora_capitata_HIv3___RNAseq.g4588.t1"; class_code "="; tss_id "TSS5";
Montipora_capitata_HIv3___Scaffold_1	StringTie	transcript	176400	276376	.	+	.	transcript_id "Montipora_capitata_HIv3___RNAseq.g4589.t1"; gene_id "MSTRG.7"; gene_name "Montipora_capitata_HIv3___RNAseq.g4589.t1"; xloc "XLOC_000006"; ref_gene_id "Montipora_capitata_HIv3___RNAseq.g4589.t1"; cmp_ref "Montipora_capitata_HIv3___RNAseq.g4589.t1"; class_code "="; tss_id "TSS6";
Montipora_capitata_HIv3___Scaffold_1	StringTie	transcript	204191	204451	.	+	.	transcript_id "Montipora_capitata_HIv3___TS.g26276.t1"; gene_id "MSTRG.8"; gene_name "Montipora_capitata_HIv3___TS.g26276.t1"; xloc "XLOC_000007"; ref_gene_id "Montipora_capitata_HIv3___TS.g26276.t1"; cmp_ref "Montipora_capitata_HIv3___TS.g26276.t1"; class_code "="; tss_id "TSS7";
Montipora_capitata_HIv3___Scaffold_1	StringTie	transcript	223366	223638	.	+	.	transcript_id "Montipora_capitata_HIv3___TS.g26277.t1"; gene_id "MSTRG.9"; gene_name "Montipora_capitata_HIv3___TS.g26277.t1"; xloc "XLOC_000008"; ref_gene_id "Montipora_capitata_HIv3___TS.g26277.t1"; cmp_ref "Montipora_capitata_HIv3___TS.g26277.t1"; class_code "="; tss_id "TSS8";
Montipora_capitata_HIv3___Scaffold_1	StringTie	transcript	330242	330991	.	+	.	transcript_id "Montipora_capitata_HIv3___RNAseq.g4592.t1"; gene_id "MSTRG.12"; gene_name "Montipora_capitata_HIv3___RNAseq.g4592.t1"; xloc "XLOC_000009"; ref_gene_id "Montipora_capitata_HIv3___RNAseq.g4592.t1"; cmp_ref "Montipora_capitata_HIv3___RNAseq.g4592.t1"; class_code "="; tss_id "TSS9";
Montipora_capitata_HIv3___Scaffold_1	StringTie	transcript	396628	400449	.	+	.	transcript_id "Montipora_capitata_HIv3___TS.g26284.t1"; gene_id "MSTRG.19"; gene_name "Montipora_capitata_HIv3___TS.g26284.t1"; xloc "XLOC_000010"; ref_gene_id "Montipora_capitata_HIv3___TS.g26284.t1"; cmp_ref "Montipora_capitata_HIv3___TS.g26284.t1"; class_code "="; tss_id "TSS10";

wc -l Mcap_lncRNA_candidates.gtf
54314 Mcap_lncRNA_candidates.gtf
```

Use bedtools to extract sequence data from the reference genome based on the coordinates from the GTF. The resulting seqs represent potential lncRNA candidates. 

```
interactive 
module load BEDTools/2.30.0-GCC-11.3.0 
bedtools getfasta -fi /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.assembly.fasta -bed Mcap_lncRNA_candidates.gtf -fo Mcap_lncRNA_candidates.fasta -fullHeader -split

zgrep -c ">" Mcap_lncRNA_candidates.fasta
```

Run CPC2. 

```
module load CPC2/1.0.1-foss-2022a
CPC2.py -i /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output/stringtie/Mcap_lncRNA_candidates.fasta

awk '$8 != "coding"' cpc2output.txt.txt > Mcap_noncoding_transcripts_info.txt
head Mcap_noncoding_transcripts_info.txt
#ID	transcript_length	peptide_length	Fickett_score	pI	ORF_integrity	coding_probability	label
transcript::Montipora_capitata_HIv3___Scaffold_1:204190-204451	261	87	0.35153	10.518226432800294	1	0.10615	noncoding
transcript::Montipora_capitata_HIv3___Scaffold_1:223365-223638	273	91	0.45305	9.851301002502442	1	0.437989	noncoding
transcript::Montipora_capitata_HIv3___Scaffold_1:330241-330991	750	53	0.44283	4.5459486007690435	1	0.100817	noncoding
transcript::Montipora_capitata_HIv3___Scaffold_1:499117-506138	7021	102	0.29595000000000005	9.135637474060061	1	0.0984293	noncoding
transcript::Montipora_capitata_HIv3___Scaffold_1:581275-583890	2615	84	0.39626000000000006	9.359922981262208	1	0.154755	noncoding
transcript::Montipora_capitata_HIv3___Scaffold_1:600468-626017	25549	131	0.2921	10.076875877380367	1	0.19239	noncoding
transcript::Montipora_capitata_HIv3___Scaffold_1:696118-715336	19218	107	0.29367000000000004	10.136251258850098	1	0.108674	noncoding
transcript::Montipora_capitata_HIv3___Scaffold_1:1058733-1069888	11155	118	0.29367000000000004	5.867625617980956	1	0.3621	noncoding
transcript::Montipora_capitata_HIv3___Scaffold_1:1128944-1147155	18211	133	0.2921	8.951193428039552	1	0.23396	noncoding

awk '$8 == "noncoding" {print $1}' cpc2output.txt.txt > Mcap_noncoding_transcripts_ids.txt
head Mcap_noncoding_transcripts_ids.txt
transcript::Montipora_capitata_HIv3___Scaffold_1:204190-204451
transcript::Montipora_capitata_HIv3___Scaffold_1:223365-223638
transcript::Montipora_capitata_HIv3___Scaffold_1:330241-330991
transcript::Montipora_capitata_HIv3___Scaffold_1:499117-506138
transcript::Montipora_capitata_HIv3___Scaffold_1:581275-583890
transcript::Montipora_capitata_HIv3___Scaffold_1:600468-626017
transcript::Montipora_capitata_HIv3___Scaffold_1:696118-715336
transcript::Montipora_capitata_HIv3___Scaffold_1:1058733-1069888
transcript::Montipora_capitata_HIv3___Scaffold_1:1128944-1147155
transcript::Montipora_capitata_HIv3___Scaffold_1:1152456-1156614
```

Extract putative lncRNAs from fasta 
```
sed 's/^transcript:://' Mcap_noncoding_transcripts_ids.txt > Mcap_noncoding_transcripts_ids_modified.txt

grep -A 1 -Ff Mcap_noncoding_transcripts_ids_modified.txt Mcap_lncRNA_candidates.fasta | sed '/^--$/d' > Mcap_lncRNAs_putative.fasta
```

I have putative lncRNAs. 

### 20241112

Download rfam sequences to blast against putative lncRNAs to remove any unwanted rRNAs, tRNAs, etc. First, make a refs folder for the sequences

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/mRNA
mkdir refs 
```

In the scripts folder: `nano wget_rfam.sh`

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/scripts            
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

echo "Download rfam sequences" $(date)

cd /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/refs

wget -r -np -nd -A .fa.gz ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/

echo "Download complete, cat all sequences" $(date)

cat *fa.gz > rfam_seqs.fa.gz
rm RF*

echo "Cat complete" $(date)
```

Submitted batch job 348568

After the sequences are done downloading, blast the putative lncRNAs against the rfam sequences. Make new folder for output. 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output
mkdir lnc_blast 
```

In the scripts folder: `nano blastn_rfam.sh`

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/scripts            
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load BLAST+/2.9.0-iimpi-2019b

echo "Blasting putative lncRNAs against rfam" $(date)
echo "Making blast db" $(date)

makeblastdb -in /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/refs/rfam_seqs.fa -out /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/refs/rfam_db -dbtype nucl

echo "Db creation complete, blast putative lncRNAs against rfam db" $(date)

blastn -query /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output/stringtie/Mcap_lncRNAs_putative.fasta -db /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/refs/rfam_db -out /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output/lnc_blast/Mcap_lnc_rfam_blastn.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

echo "Blast complete" $(date)

```

Submitted batch job 348600. Ran in about 5 mins. Look at output

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output/lnc_blast

head Mcap_lnc_rfam_blastn.tab 
Montipora_capitata_HIv3___Scaffold_10:50987091-51005681	MU825396.1/212773-212653	96.694	121	4	0	10869	10989	121	1	1.95e-48	202
Montipora_capitata_HIv3___Scaffold_1013:6609-13142	MU826834.1/1912895-1913082	95.531	179	6	2	3345	3522	188	11	6.68e-74	285
Montipora_capitata_HIv3___Scaffold_1013:13927-20177	MU827310.1/718302-718465	96.951	164	5	0	2981	3144	164	1	3.85e-71	276
Montipora_capitata_HIv3___Scaffold_11:3287015-3291500	MU827309.1/2350248-2350386	96.063	127	5	0	4035	4161	127	1	1.02e-50	207
Montipora_capitata_HIv3___Scaffold_13:24569367-24573028	MU827310.1/717380-717262	95.798	119	5	0	1580	1698	1	119	2.35e-46	193
Montipora_capitata_HIv3___Scaffold_13:24617531-24626060	LJWW01001071.1/6624-6820	95.690	116	3	2	2761	2876	142	29	9.11e-44	185
Montipora_capitata_HIv3___Scaffold_13:24738989-24746769	MU826834.1/1912895-1913082	95.506	178	5	2	3075	3251	188	13	1.03e-72	281
Montipora_capitata_HIv3___Scaffold_13:24798335-24808343	LSMT01000541.1/6279-6088	97.917	192	4	0	2639	2830	192	1	3.57e-88	333
Montipora_capitata_HIv3___Scaffold_4:36816332-36834523	MU827817.1/523255-523552	90.333	300	25	4	17249	17546	1	298	3.80e-105	390
Montipora_capitata_HIv3___Scaffold_4:31137676-31146429	MU826834.1/1905971-1906134	93.252	163	11	0	4963	5125	163	1	1.97e-60	241

wc -l Mcap_lnc_rfam_blastn.tab 
15 Mcap_lnc_rfam_blastn.tab
```

Only 15 sequences got hits as possible rfam contaminants. Make a list of these. 

```
awk '{print $1}' Mcap_lnc_rfam_blastn.tab > sequences_to_remove.txt
```

Remove sequences from fasta

```
cd ../stringtie
grep -v -F -f <(sed 's/^/>/' /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output/lnc_blast/sequences_to_remove.txt) Mcap_lncRNAs_putative.fasta | sed '/^$/d' > Mcap_lncRNAs_final.fasta

zgrep -c ">" Mcap_lncRNAs_putative.fasta
31519
zgrep -c ">" Mcap_lncRNAs_final.fasta
31504
```

### 20241113 

I need to quantify the lncRNAs with kallisto. Make kallisto output folder

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output
mkdir kallisto 
```

The lncRNA fasta (`Mcap_lncRNAs_final.fasta`) will be used to generate the kallisto index. Then, the RNAseq reads (located `/data/putnamlab/jillashey/DT_Mcap_2023/mRNA/data/trim`) will be aligned to the lncRNA index using kallisto. In the scripts folder: `nano kallisto_lnc.sh`

```
#!/bin/bash 
#SBATCH -t 100:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/scripts            
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load kallisto/0.48.0-gompi-2022a 

echo "Creating kallisto index from lncRNA fasta" $(date)

kallisto index -i /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output/kallisto/mcap_lncRNA_index /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output/stringtie/Mcap_lncRNAs_final.fasta

echo "lncRNA index creation complete, starting read alignment" $(date)

array=($(ls /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/data/trim/*R1_001.fastq.gz)) # call the clean sequences - make an array to align

for i in ${array[@]}; do
    # Extract just the filename from the input FASTQ file path
    filename=$(basename ${i})
    # Construct the output directory path
    output_dir="/data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output/kallisto/kallisto.${filename}"
    # Run kallisto quant
    kallisto quant -i /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output/kallisto/mcap_lncRNA_index -o ${output_dir} ${i} $(echo ${i} | sed 's/_R1/_R2/')
done

echo "lncRNA alignment complete!" $(date)
```

Submitted batch job 348664. Ran in about 3 hours. Run trinity to generate lncRNa count matrix. in the scripts folder: `nano trinity_gene_matrix.sh`

```
#!/bin/bash 
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=32GB --cpus-per-task=24
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Astrangia2021/lncRNA/scripts              
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load Trinity/2.15.1-foss-2022a

echo "Use trinity to generate lncRNA count matrix" $(date)

perl $EBROOTTRINITY/trinityrnaseq-v2.15.1/util/abundance_estimates_to_matrix.pl \
--est_method kallisto \
--gene_trans_map none \
--out_prefix /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output/kallisto/Mcap_lncRNA_count_matrix \
--name_sample_by_basedir \
/data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output/kallisto/*/abundance.tsv

echo "LncRNA count matrix created!" $(date)
```

Submitted batch job 348671. Ran super fast but no slurm error or output files...But we did get 4 output files: 

```
-rw-r--r--. 1 jillashey 8.1M Nov 13 12:20 Mcap_lncRNA_count_matrix.isoform.TPM.not_cross_norm
-rw-r--r--. 1 jillashey 6.6M Nov 13 12:20 Mcap_lncRNA_count_matrix.isoform.counts.matrix
-rw-r--r--. 1 jillashey    0 Nov 13 12:20 Mcap_lncRNA_count_matrix.isoform.TMM.EXPR.matrix
-rw-r--r--. 1 jillashey  684 Nov 13 12:20 Mcap_lncRNA_count_matrix.isoform.TPM.not_cross_norm.runTMM.R
```

According to the Trinity [wiki](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification#building-expression-matrices), these files mean the following: 

```
kallisto.isoform.counts.matrix  : the estimated RNA-Seq fragment counts (raw counts)
kallisto.isoform.TPM.not_cross_norm  : a matrix of TPM expression values (not cross-sample normalized)
kallisto.isoform.TMM.EXPR.matrix : a matrix of TMM-normalized expression values
```

Copy `Mcap_lncRNA_count_matrix.isoform.counts.matrix` onto local computer! 

### 20241118

Huang et al. 2019 used [RIblast](https://academic.oup.com/bioinformatics/article/33/17/2666/3778372) to compute interaction probabilities for each possible lncRNA-mRNA pair based on free energy required for hybridization. The github for RIblast is [here](https://github.com/fukunagatsu/RIblast). 

```
cd /data/putnamlab/jillashey
git clone https://github.com/fukunagatsu/RIblast.git
make # compile 
```

The first step is to generate a database. I'm going to generate a db of mRNAs and the lncRNAs will be my query. No idea how long this is supposed to take. 

```
interactive 

cd /data/putnamlab/jillashey/genome/Mcap/V3
gunzip Montipora_capitata_HIv3.genes.cds.fna.gz 

cd /data/putnamlab/jillashey/genome/RIblast
export PATH=$PATH:/data/putnamlab/jillashey/RIblast
RIblast db -i /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.genes.cds.fna -o mRNA_db
``` 

Taking a while so I am going to submit as a job. In the Mcap scripts folder (`/data/putnamlab/jillashey/DT_Mcap_2023/mRNA/scripts`): `nano RIblast_db.sh`

```
#!/bin/bash 
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=32GB --cpus-per-task=24
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/scripts             
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

echo "Make RI blast db for mRNAs" $(date)

cd /data/putnamlab/jillashey/genome/RIblast
export PATH=$PATH:/data/putnamlab/jillashey/RIblast
RIblast db -i /data/putnamlab/jillashey/genome/Mcap/V3/Montipora_capitata_HIv3.genes.cds.fna -o mRNA_db

echo "mRNA db creation complete" $(date)
```

Submitted batch job 349336

### 20241119

Took about 13 hours. Move output into refs folder. 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/scripts
mv mRNA_db.* ../refs/
```

Make output folder for RIblast 

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output
mkdir RIblast
```

Next step is to search the potential RNA-RNA interactions of the query sequence (lncRNAs) to the db sequences (mRNAs). In the scripts folder: `nano RIblast_ris.sh`

```
#!/bin/bash 
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=125GB --cpus-per-task=24
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/scripts             
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

echo "Start RNA-RNA interaction predictions" $(date)

export PATH=$PATH:/data/putnamlab/jillashey/RIblast
RIblast ris -i /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output/stringtie/Mcap_lncRNAs_final.fasta -o /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/output/RIblast/RIblast_lncRNA_output.txt -d /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/refs/mRNA_db -e -20.0

echo "RNA interaction prediction complete - lncRNA as query and mRNA as db" $(date)
```

Submitted batch job 349435. The `-e` flag refers to the hybridization energy threshold for seed search, default is -6.0 but I set it at -20.0 for increased stringency. 