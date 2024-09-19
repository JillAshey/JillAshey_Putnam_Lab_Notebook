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




