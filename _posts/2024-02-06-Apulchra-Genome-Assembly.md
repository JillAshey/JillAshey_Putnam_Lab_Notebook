---
layout: post
title: Apulchra genome assembly 
date: '2024-02-06'
categories: Analysis, Genome Assembly
tags: [Bioinformatics, Genome, Assembly]
projects: Apulchra genome
---

## Apulchra genome assembly 

Sperm and tissue from adult *Acropora pulchra* colonies were collected from Moorea, French Polynesia and sequencing with PacBio (long reads) and Illumina (short reads). This post will detail the genome assembly notes. The github for this project is [here](https://github.com/hputnam/Apulchra_genome/tree/main). 

I'm going to write notes and code chronologically so that I can keep track of what I'm doing each day. When assembly is complete, I will compile the workflow in a separate post. 

### 20240206

Met w/ Ross and Hollie today re Apulchra genome assembly. We decided to move forward with the workflow from [Stephens et al. 2022](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac098/6815755) which assembled genomes for 4 Hawaiian corals. 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/Stephens_et_al_2022_workflow.png)

PacBio long reads were received in late Jan/early Feb 2024. According to reps from Genohub, the [PacBio raw output](https://github.com/hputnam/Apulchra_genome/blob/main/DNA_Seq_Info/20240129_Project_6693786_Acropora_pulchra_PacBio_Seq_Summary.pdf) looks good. 

We decided to move forward with [Canu](https://canu.readthedocs.io/en/latest/index.html) to assembly the genome. Canu is specialized to assemble PacBio sequences, operating in three phases: correction, trimming and assembly. According to the Canu website, "The correction phase will improve the accuracy of bases in reads. The trimming phase will trim reads to the portion that appears to be high-quality sequence, removing suspicious regions such as remaining SMRTbell adapter. The assembly phase will order the reads into contigs, generate consensus sequences and create graphs of alternate paths." 

The PacBio files that will be used for assembly are located here on Andromeda: `/data/putnamlab/KITT/hputnam/20240129_Apulchra_Genome_LongRead`. The files in the folder that we will use are: 

```
m84100_240128_024355_s2.hifi_reads.bc1029.bam
m84100_240128_024355_s2.hifi_reads.bc1029.bam.pbi
```

The bam file contains all of the read information in computer language and the pbi file is an index file of the bam. Both are needed in the same folder to run Canu. 

For Canu, input files must be fasta or fastq format. I'm going to use `bam2fastq`from the [PacBio github](https://github.com/pacificbiosciences/pbtk/). This module is not on Andromeda so I will need to install it via conda. 

The PacBio sequencing for the Apul genome were done with HiFi sequencing that are produced with circular consensus sequencing on PacoBio long read systems. Here's how HiFi reads are generated from the [PacBio website](https://www.pacb.com/technology/hifi-sequencing/):

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/hifi_read_generation.png)

Since Hifi sequencing was used, a specific HiCanu flag (`-pacbio-hifi`) must be used. Additionally, in the Canu [tutorial](https://canu.readthedocs.io/en/latest/quick-start.html#), it says that if this flag is used, the program is assuming that the reads are trimmed and corrected already. However, our reads are not. I'm going to try to run the first pass at Canu with the `-raw` flag. 

Before Canu, I will run `bam2fastq`. This is a part of the PacBio BAM toolkit package `pbtk`. I need to create the conda environment and install the package. Load miniconda module: `module load Miniconda3/4.9.2`. Create a conda environment. 

```
conda create --prefix /data/putnamlab/conda/pbtk
```

Install the package. Once the package is installed on the server, anyone in the Putnam lab group can use it. 

```
conda activate /data/putnamlab/conda/pbtk
conda install -c bioconda pbtk
```

In my own directory, make a new directory for genome assembly

```
cd /data/putnamlab/jillashey
mkdir Apul_Genome
cd Apul_Genome
mkdir assembly structural functional
cd assembly
mkdir scripts data output
```

Run code to make the PacBio bam file to a fastq file. In the scripts folder: `nano bam2fastq.sh`

```
#!/bin/bash -i
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Apul_Genome/assembly/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

conda activate /data/putnamlab/conda/pbtk

echo "Convert PacBio bam file to fastq file" $(date)

bam2fastq -o /data/putnamlab/jillashey/Apul_Genome/assembly/data/m84100_240128_024355_s2.hifi_reads.bc1029.fastq /data/putnamlab/KITT/hputnam/20240129_Apulchra_Genome_LongRead/m84100_240128_024355_s2.hifi_reads.bc1029.bam

echo "Bam to fastq complete!" $(date)

conda deactivate
```

Submitted batch job 294235

### 20240208

Job pended for about a day, then ran in 1.5 hours. I got this error message: `bash: cannot set terminal process group (-1): Function not implemented bash: no job control in this shell`, but not sure why. A fastq.gz file was produced! The file is pretty large (35G).

```
less m84100_240128_024355_s2.hifi_reads.bc1029.fastq.fastq.gz 

@m84100_240128_024355_s2/261887593/ccs
TATAAGTTTTACAGCTGCCTTTTGCTCAGCAAAGAAAGCAGCATTGTTATTGAACAGAAAAAGCCTTTTGGTGATATAAAGGTTTCTAAGGGACCAAAGTTTGATCTAGTATGCTAAGTGTGGTGGGTTAAAACTTTGTTTCACCTTTTTCCCGTGATGTTACAAAATTGGTGCAAATATTCATGACGTCGTACTCAAGTCTGACACTAAGAGCATGCAATTACTTAAACAAACAAGCCATACACCAATAAACTGAAGCTCTGTCAACTAGAAAACCTTTGAGTATTTTCATTGTAAAGTACAAGTGGTATATTGTCACTTGCTTTTACAACTTGAAGAACTCACAGTTAAGTTACTAGATTCACCATAGTGCTTGGCAATGAAGAAGCCAAATCACATAAAGTCGGAGCATGTGGTGTTTAGACCTAATCAAACAAGAACACAATATTTAGTACCTGCATCCTGTCTAAGGAGGAAATTTTAAGCTGCTTTCTTTTAAATTTTTTTTATTAGCATTTCAATGGTTGAGGTCGATTATAGGTGCTAGGCTTTAATTCCGTACTATGAAAGAAGAAAGGTCGTTGTTATTAACCATGTCAAACAGAGAAACACATGGTAAAAAATTGACTTCCTTTTCCTCTCGTTGCCACTTAAGCTTAATGATGGTGTTTGACCTGAAAGATGTTACAATTGTTTTAGATGAAAAGACTGTTCTGCGTAAAACAGTGAAGCCTCCCAACTTATTTTGTTATGTGGATTTTGTTGTCTTGTTAGTAACATGTATTGGACTATCTTTTGTGAGTACATAGCTTTTTTTCCATCAACTGACTATATACGTGGTGTAATTTGAGATCATGCCTCCAAGTGTTAGTCTTTTGTTTGGGGCTAACTCGTAAAAGACAAAGGGAGGGGGGTTGTCTAATTCCTAAGCAAAGCATTAAGTTTAACACAGGAAATTGTTTGCGTTGATATTGCTATCCTTTCAGCCCCAAACAAAAAATTTAATGGTTATTTTATTTTACATCTATTGTAAATATATTTTAACATTAATTTTTATTATTGCACTGTAAATACTTGTACTAATGTTCTGTTTGAATTAATTTTGATTCATTCCTTGTGCTTACAACAACAGGGATACAAAACCGATATGTATAATAATACTATTAGAGATGCTTATTTGCATTTTTAGCCCATACCATGAGTTTTAATAACGCCAGGCCATTGGAGATTTTATGGAGTGAGGATTCATTGTACAAACATGGTTGATTTAATATTAAAGTTGTATCCAAATAATTAATATCTGCTGTGATCAGTGAAAGATTGACCTTTCAGTTGTTTGGTTGCACCTTCATCTTATTGGAAACAACTGAATGGAGCATCTTTCCAGTTTAAAAATGTACCACTGCCCACTTTCATGAAGTTATGCCACATATTAATAATGACTATTAATTGTTGAAAACCCTTCTTCCAAAATGTTTCCATTTATTTGTAATAGCATATGTGGTCCATCAGAACAATAATTTAAATCATTACTATTAATAATTTTCCAATAACTGACTTTCAAACCTAGCCAACAAGCATAAGTCAGTAAGCCACAGAGTCCAGAGATACACTTACACTTTACTTTTCACTTCTGAAACATTTTATAATCTCAGTATGAGCATAGAACTTTTCAGTTGGGCAGCATGGAATAGAACCTTTGGACCCCTCTGTGAATATCAAAAATAGGCAACCACTTCCAGCATACACTCTAGCCTCCTTCATAAAGCAAGCCTTAGTGTTTTAGCTTCTACTAGTTAGATTCATTTTAAAAGAAGTTCAGTATACTTAATCTTATAGAAGCTGATTGTGATATAATTGCATAGGTGGATCTCAGAAAAGTGAGATGTAGCTGTCAAATTAAAAGAAGTCCTTTCCAAGCGTAGCTTCTGATAAACAATGCATTTTAGTTAACATTGGATTATGGTTTCAAAGGACTTGTAAAGCTAAATTCAAGTTTTTATGACAACTTGAAAGCCTTTGCCACAGTCTCCGCTGATTTAAGACTTCCATCAAAGTTAGAGTGGTGTGAATGCATCTCCACATGCAATTAATAAAGGTGAGGCAGACAACACAAAACACCCTGGTGCACCATCAACTCCACGGATCACTTGACTGTAACGCCATCTTATACAGCGACTGCCAATTGGAACTGGAAGATCAGGAATGATCTCTTTCACATGGGAAATGAGCATGGTCTTGATAATGCTTTCATCAGTATCCCAATGTTGAAGACAGAAGGACTTTGTTGAATGTACTAAAAGAGAGGGACCAACATCCACTGGATCTGGAACAGAATGCCAAAGCATATGCAAAATGTTATCATGTATAATAAAGAAAAAACCAAACTCATGAAGTTCAACGGTCACAGGCTTTGGTAAAAGACTGTACATGGTAAAGTACTCAGGGTGAAAGTTTGTATCCAAAAAACGAAAAGCTAACCTGTACCACATCTCTTCCGAGAGTCCACTGACACAAATGCAATATTAGGATTGCCAGTCGCATATTTACACGTCCATGGCACATTGATAACTGCTTCATGATCAAAAAAGTATCCTACTGCAAACCGTGATGAATATTCCACATTTTGTAAGGATTTGATTTCATTTTGTAAGAATGCTTGAATTGAACCTTGTAGCTGAAGAAGTTGTGGTACTGGAATAGTTACTATGACTGA
```

See how many `@m84100` are in the file. I'm not sure what these stand for, maybe contigs? 

```
zgrep -c "@m84100" m84100_240128_024355_s2.hifi_reads.bc1029.fastq.fastq.gz 
5898386
```

More than 5 million, so many not contigs? I guess it represents the number of HiFi reads generated. Now time to run Canu! Canu is already installed on the server, which is nice. 

In the scripts folder: `nano canu.sh`

```
#!/bin/bash 
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Apul_Genome/assembly/scripts
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load canu/2.2-GCCcore-11.2.0 

cd /data/putnamlab/jillashey/Apul_Genome/assembly/data

echo "Unzip paco-bio fastq file" $(date)

gunzip m84100_240128_024355_s2.hifi_reads.bc1029.fastq.fastq.gz

echo "Unzip complete, starting assembly" $(date)

canu -p apul -d /data/putnamlab/jillashey/Apul_Genome/assembly/data genomeSize=475m -raw -pacbio-hifi m84100_240128_024355_s2.hifi_reads.bc1029.fastq.fastq

echo "Canu assembly complete" $(date)
```

I'm not sure if the `raw` and `-pacbio-hifi` will be compatible, as the Canu tutorial says that the `-pacbio-hifi` assumes that the input is trimmed and corrected (still not sure what this means). Submitted batch job 294325


