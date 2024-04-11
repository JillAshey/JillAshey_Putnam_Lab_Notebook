---
layout: post
title: e5 deep dive ncRNA protein machinery 
date: '2024-04-11'
categories: Analysis
tags: [Bioinformatics, BLAST]
projects: e5 deep dive
---

## e5 deep dive - investigating ncRNA protein machinery in 3 species of corals 

The e5 deep dive project is examining ncRNA dynamics in 3 species of coral from Moorea, French Polynesia. The github for the project is [here](https://github.com/urol-e5/deep-dive). 

In this post, I will be assessing the presence of ncRNA-related proteins in the protein fasta files from the 3 species of interest (*Acropora pulchra*, *Porites evermanni*, and *Pocillopora tuahiniensis*). To assess presence, I gathered ncRNA-related protein sequences from related species on NCBI (*Pocillopora verrocosa*, *Orbicella faveolata*, *Acropora millepora*, *Stylophora pistillata*, *Nematostella vectensis*, and *Homo sapiens*) and will blast these sequences against the protein fasta files from the 3 coral species. 

The proteins I chose to assess are: 

- DNMT1
- DNMT3A
- Drosha 
- DGCR8
- XPO5
- Dicer
- AGO2
- Piwi
- PPK-1/PIP5K1A
- RNase P

A list of these proteins (from *Pocillopora verrocosa*, *Orbicella faveolata*, *Acropora millepora*, *Stylophora pistillata*, *Nematostella vectensis*, and *Homo sapiens*) with their NCBI accession number and links can be found in this [spreadsheet](https://docs.google.com/spreadsheets/d/1vyW-TaPRl1RgdYLJ4DIkl5RBm5e9WOKcBMWZ_sEe-r8/edit#gid=0). 

I already have an e5 folder on the HPC server but I am going to make new folders in it for this analysis. 

```
cd /data/putnamlab/jillashey/e5
mkdir ncRNA_prot scripts refs
cd ncRNA_prot 
```

In the `ncRNA_prot` folder, I am going to make a fasta file for each protein category with the sequence info from each species from NCBI included. For example, `dnmt1.fasta` would include all DNMT1 protein sequences from *Pocillopora verrocosa*, *Orbicella faveolata*, *Acropora millepora*, *Stylophora pistillata*, *Nematostella vectensis*, and *Homo sapiens*. The following fastas will be created: 

```
dnmt1.fasta
dnmt3a.fasta
drosha.fasta
dgcr8.fasta
xpo5.fasta
dicer.fasta
ago2.fasta
piwi.fasta
pip5k1a.fasta
rnase_p.fasta
```

The protein fasta files from *Acropora pulchra*, *Porites evermanni*, and *Pocillopora tuahiniensis* are already here XXXXX

In the scripts folder: `nano makeblastdb.sh`

```
#!/bin/bash
#SBATCH --job-name="makeblastdb"
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH -t 24:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/e5/scripts           
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load BLAST+/2.9.0-iimpi-2019b

```
