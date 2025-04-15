---
layout: post
title: Developmental 2023 Timeseries mRNA analysis - polyA
date: '2025-04-15'
categories: Analysis
tags: [Bioinformatics, mRNA, Montipora capitata]
projects: Developmental Timseries 2023
---

## Developmental 2023 Timeseries mRNA analysis - polyA

These data came from my developmental timeseries experiment in 2023 with *Montipora capitata* in Hawaii. In this experiment, *Montipora capitata* embryos and larvae were exposed to ambient and heat stress over 72 hours from embryo to swimming larvae. Samples were collected at 8 time points: 1, 4, 9, 14, 22, 28, 48 and 72 hpf. The github for this project is [here](https://github.com/JillAshey/Hawaii_Developmental_TimeSeries). 

We sequenced these samples using rRNA depletion library methods and we wanted to compare that with libraries prepared using polyA selection. The library prep was done by me using the [Zymo 3' Switchfree](https://www.zymoresearch.com/products/zymo-seq-switchfree-3-mrna-library-kit) library prep kit (see my notebook [posts](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook) for more information on library prep). The sequencing was done through Oklahoma Medical Research Foundation NGS Core, which I found on Genohub. Here is the sequencing information: 

- Instrument: Illumina NovaSeq X Plus - 10B - PE 150 Cycle
- Read length: 2 x 150bp (Paired End)
- Pricing unit: per lane
- Number of samples: 8 libraries
- Guaranteed number of pass filter PE reads/sample: 20M (10M in each direction)
- Deliverables: FastQ files uploaded to Genohub project bucket

Files were downloaded to this location on Unity: XXXX. Goodbye Andromeda, hello Unity. Time to analyze!

Check files to ensure transfer is complete and no files were corrupted. 

```
md5sum *fastq.gz > checkmd5.txt
md5sum -c checkmd5.txt 
M10_S136_R1_001.fastq.gz: OK
M10_S136_R2_001.fastq.gz: OK
M11_S137_R1_001.fastq.gz: OK
M11_S137_R2_001.fastq.gz: OK
M13_S138_R1_001.fastq.gz: OK
M13_S138_R2_001.fastq.gz: OK
M14_S139_R1_001.fastq.gz: OK
M14_S139_R2_001.fastq.gz: OK
M23_S140_R1_001.fastq.gz: OK
M23_S140_R2_001.fastq.gz: OK
M24_S141_R1_001.fastq.gz: OK
M24_S141_R2_001.fastq.gz: OK
M26_S142_R1_001.fastq.gz: OK
M26_S142_R2_001.fastq.gz: OK
M28_S143_R1_001.fastq.gz: OK
M28_S143_R2_001.fastq.gz: OK
M35_S144_R1_001.fastq.gz: OK
M35_S144_R2_001.fastq.gz: OK
M36_S145_R1_001.fastq.gz: OK
M36_S145_R2_001.fastq.gz: OK
M37_S146_R1_001.fastq.gz: OK
M37_S146_R2_001.fastq.gz: OK
M39_S147_R1_001.fastq.gz: OK
M39_S147_R2_001.fastq.gz: OK
M47_S148_R1_001.fastq.gz: OK
M47_S148_R2_001.fastq.gz: OK
M48_S149_R1_001.fastq.gz: OK
M48_S149_R2_001.fastq.gz: OK
M51_S150_R1_001.fastq.gz: OK
M51_S150_R2_001.fastq.gz: OK
M52_S151_R1_001.fastq.gz: OK
M52_S151_R2_001.fastq.gz: OK
M60_S152_R1_001.fastq.gz: OK
M60_S152_R2_001.fastq.gz: OK
M61_S153_R1_001.fastq.gz: OK
M61_S153_R2_001.fastq.gz: OK
M62_S154_R1_001.fastq.gz: OK
M62_S154_R2_001.fastq.gz: OK
M63_S155_R1_001.fastq.gz: OK
M63_S155_R2_001.fastq.gz: OK
M6_S132_R1_001.fastq.gz: OK
M6_S132_R2_001.fastq.gz: OK
M72_S156_R1_001.fastq.gz: OK
M72_S156_R2_001.fastq.gz: OK
M73_S157_R1_001.fastq.gz: OK
M73_S157_R2_001.fastq.gz: OK
M74_S158_R1_001.fastq.gz: OK
M74_S158_R2_001.fastq.gz: OK
M75_S159_R1_001.fastq.gz: OK
M75_S159_R2_001.fastq.gz: OK
M7_S133_R1_001.fastq.gz: OK
M7_S133_R2_001.fastq.gz: OK
M85_S160_R1_001.fastq.gz: OK
M85_S160_R2_001.fastq.gz: OK
M86_S161_R1_001.fastq.gz: OK
M86_S161_R2_001.fastq.gz: OK
M87_S162_R1_001.fastq.gz: OK
M87_S162_R2_001.fastq.gz: OK
M88_S163_R1_001.fastq.gz: OK
M88_S163_R2_001.fastq.gz: OK
M8_S134_R1_001.fastq.gz: OK
M8_S134_R2_001.fastq.gz: OK
M9_S135_R1_001.fastq.gz: OK
M9_S135_R2_001.fastq.gz: OK
```

Make directories on Unity. Sym link the raw data in the project directory to the work directory. 

```
cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/data/raw
ln -s /project/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/*fastq.gz .
```

Run fastqc on raw reads. `nano raw_fastqc.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/scripts 

# load modules needed
module load parallel/20240822
module load fastqc/0.12.1
module load uri/main
module load all/MultiQC/1.12-foss-2021b

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/data/raw

# Create an array of fastq files to process
files=($('ls' *001.fastq.gz)) 

# Run fastqc in parallel
echo "Starting fastqc..." $(date)
parallel -j 20 "fastqc {} -o /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/fastqc/raw && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)

cd /work/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/output/fastqc/raw

echo "Starting multiqc..." $(date)
multiqc *

echo "Initial QC of polyA data complete." $(date)
```

Submitted batch job 32423847

