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
mkdir ncRNA_prot scripts refs output
cd ncRNA_prot 
```

In the `ncRNA_prot` folder, I am going to make a fasta file for each protein category with the sequence info from each species from NCBI included. For example, `dnmt1.fasta` would include all DNMT1 protein sequences from *Pocillopora verrocosa*, *Orbicella faveolata*, *Acropora millepora*, *Stylophora pistillata*, *Nematostella vectensis*, and *Homo sapiens*. The following fastas were created: 

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

zgrep -c ">" *
ago2.fasta:28
dgcr8.fasta:8
dicer.fasta:19
dnmt1.fasta:8
dnmt3a.fasta:7
drosha.fasta:7
pip5k1a.fasta:10
piwi.fasta:12
rnase_p.fasta:43
xpo5.fasta:17
```

The protein fasta files from *Acropora pulchra*, *Porites evermanni*, and *Pocillopora tuahiniensis* are already here `/data/putnamlab/jillashey/e5/ortho/protein_seqs`. These will be used to create the blast dbs. 

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

echo "Making blast dbs for e5 deep dive protein seqs" $(date)

#Amillepora for Apulchra
makeblastdb -in /data/putnamlab/jillashey/e5/ortho/protein_seqs/GCF_013753865.1_Amil_v2.1.protein.faa -out /data/putnamlab/jillashey/e5/refs/Amil_prot -dbtype prot

#Pmeandrina for Ptuhuensis
makeblastdb -in /data/putnamlab/jillashey/e5/ortho/protein_seqs/Pocillopora_meandrina_HIv1.genes.pep.faa -out /data/putnamlab/jillashey/e5/refs/Pmea_prot -dbtype prot

#Pevermanni
makeblastdb -in /data/putnamlab/jillashey/e5/ortho/protein_seqs/Porites_evermanni_v1.annot.pep.fa -out /data/putnamlab/jillashey/e5/refs/Peve_prot -dbtype prot

echo "Blast db creation complete" $(date)
```

Submitted batch job 311805. Dbs created! In the scripts folder: `nano ncRNA_prot_blastp.sh`

```
#!/bin/bash
#SBATCH --job-name="blastp"
#SBATCH --nodes=1 --ntasks-per-node=15
#SBATCH -t 72:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --mem=250GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/e5/scripts           
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load BLAST+/2.9.0-iimpi-2019b

echo "ncRNA protein machinery blast beginning" $(date)
echo "Starting first with AGO2" $(date)

echo "AGO2 for Apul blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/ago2.fasta -db /data/putnamlab/jillashey/e5/refs/Amil_prot -out /data/putnamlab/jillashey/e5/output/apul_ago2_blastp.tab -evalue 1E-40 -num_threads 15 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l apul_ago2_blastp.tab

echo "AGO2 for Ptuh blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/ago2.fasta -db /data/putnamlab/jillashey/e5/refs/Pmea_prot -out /data/putnamlab/jillashey/e5/output/ptuh_ago2_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l ptuh_ago2_blastp.tab

echo "AGO2 for Peve blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/ago2.fasta -db /data/putnamlab/jillashey/e5/refs/Peve_prot -out /data/putnamlab/jillashey/e5/output/peve_ago2_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l peve_ago2_blastp.tab

echo "Now doing DGCR8" $(date)

echo "DGCR8 for Apul blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/dgcr8.fasta -db /data/putnamlab/jillashey/e5/refs/Amil_prot -out /data/putnamlab/jillashey/e5/output/apul_dgcr8_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l apul_dgcr8_blastp.tab

echo "DGCR8 for Ptuh blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/dgcr8.fasta -db /data/putnamlab/jillashey/e5/refs/Pmea_prot -out /data/putnamlab/jillashey/e5/output/ptuh_dgcr8_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l ptuh_dgcr8_blastp.tab

echo "DGCR8 for Peve blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/dgcr8.fasta -db /data/putnamlab/jillashey/e5/refs/Peve_prot -out /data/putnamlab/jillashey/e5/output/peve_dgcr8_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l peve_dgcr8_blastp.tab

echo "Now doing Dicer" $(date)

echo "Dicer for Apul blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/dicer.fasta -db /data/putnamlab/jillashey/e5/refs/Amil_prot -out /data/putnamlab/jillashey/e5/output/apul_dicer_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l apul_dicer_blastp.tab

echo "Dicer for Ptuh blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/dicer.fasta -db /data/putnamlab/jillashey/e5/refs/Pmea_prot -out /data/putnamlab/jillashey/e5/output/ptuh_dicer_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l ptuh_dicer_blastp.tab

echo "Dicer for Peve blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/dicer.fasta -db /data/putnamlab/jillashey/e5/refs/Peve_prot -out /data/putnamlab/jillashey/e5/output/peve_dicer_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l peve_dicer_blastp.tab

echo "Now doing DNMT1" $(date)

echo "DNMT1 for Apul blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/dnmt1.fasta -db /data/putnamlab/jillashey/e5/refs/Amil_prot -out /data/putnamlab/jillashey/e5/output/apul_dnmt1_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l apul_dnmt1_blastp.tab

echo "DNMT1 for Ptuh blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/dnmt1.fasta -db /data/putnamlab/jillashey/e5/refs/Pmea_prot -out /data/putnamlab/jillashey/e5/output/ptuh_dnmt1_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l ptuh_dnmt1_blastp.tab

echo "DNMT1 for Peve blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/dnmt1.fasta -db /data/putnamlab/jillashey/e5/refs/Peve_prot -out /data/putnamlab/jillashey/e5/output/peve_dnmt1_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l peve_dnmt1_blastp.tab

echo "Now doing DNMT3A" $(date)

echo "DNMT3A for Apul blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/dnmt3a.fasta -db /data/putnamlab/jillashey/e5/refs/Amil_prot -out /data/putnamlab/jillashey/e5/output/apul_dnmt3a_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l apul_dnmt3a_blastp.tab

echo "DNMT3A for Ptuh blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/dnmt3a.fasta -db /data/putnamlab/jillashey/e5/refs/Pmea_prot -out /data/putnamlab/jillashey/e5/output/ptuh_dnmt3a_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l ptuh_dnmt3a_blastp.tab

echo "DNMT3A for Peve blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/dnmt3a.fasta -db /data/putnamlab/jillashey/e5/refs/Peve_prot -out /data/putnamlab/jillashey/e5/output/peve_dnmt3a_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l peve_dnmt3a_blastp.tab

echo "Now doing Drosha" $(date)

echo "Drosha for Apul blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/drosha.fasta -db /data/putnamlab/jillashey/e5/refs/Amil_prot -out /data/putnamlab/jillashey/e5/output/apul_drosha_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l apul_drosha_blastp.tab

echo "Drosha for Ptuh blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/drosha.fasta -db /data/putnamlab/jillashey/e5/refs/Pmea_prot -out /data/putnamlab/jillashey/e5/output/ptuh_drosha_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l ptuh_drosha_blastp.tab

echo "Drosha for Peve blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/drosha.fasta -db /data/putnamlab/jillashey/e5/refs/Peve_prot -out /data/putnamlab/jillashey/e5/output/peve_drosha_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l peve_drosha_blastp.tab

echo "Now doing pip5k1a" $(date)

echo "pip5k1a for Apul blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/pip5k1a.fasta -db /data/putnamlab/jillashey/e5/refs/Amil_prot -out /data/putnamlab/jillashey/e5/output/apul_pip5k1a_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l apul_pip5k1a_blastp.tab

echo "pip5k1a for Ptuh blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/pip5k1a.fasta -db /data/putnamlab/jillashey/e5/refs/Pmea_prot -out /data/putnamlab/jillashey/e5/output/ptuh_pip5k1a_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l ptuh_pip5k1a_blastp.tab

echo "pip5k1a for Peve blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/pip5k1a.fasta -db /data/putnamlab/jillashey/e5/refs/Peve_prot -out /data/putnamlab/jillashey/e5/output/peve_pip5k1a_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l peve_pip5k1a_blastp.tab

echo "Now doing piwi" $(date)

echo "piwi for Apul blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/piwi.fasta -db /data/putnamlab/jillashey/e5/refs/Amil_prot -out /data/putnamlab/jillashey/e5/output/apul_piwi_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l apul_piwi_blastp.tab

echo "piwi for Ptuh blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/piwi.fasta -db /data/putnamlab/jillashey/e5/refs/Pmea_prot -out /data/putnamlab/jillashey/e5/output/ptuh_piwi_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l ptuh_piwi_blastp.tab

echo "piwi for Peve blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/piwi.fasta -db /data/putnamlab/jillashey/e5/refs/Peve_prot -out /data/putnamlab/jillashey/e5/output/peve_piwi_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l peve_piwi_blastp.tab

echo "Now doing rnase P" $(date)

echo "rnase P for Apul blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/rnase_p.fasta -db /data/putnamlab/jillashey/e5/refs/Amil_prot -out /data/putnamlab/jillashey/e5/output/apul_rnase_p_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l apul_rnase_p_blastp.tab

echo "rnase P for Ptuh blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/rnase_p.fasta -db /data/putnamlab/jillashey/e5/refs/Pmea_prot -out /data/putnamlab/jillashey/e5/output/ptuh_rnase_p_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l ptuh_rnase_p_blastp.tab

echo "rnase P for Peve blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/rnase_p.fasta -db /data/putnamlab/jillashey/e5/refs/Peve_prot -out /data/putnamlab/jillashey/e5/output/peve_rnase_p_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l peve_rnase_p_blastp.tab

echo "Now doing xpo5" $(date)

echo "xpo5 for Apul blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/xpo5.fasta -db /data/putnamlab/jillashey/e5/refs/Amil_prot -out /data/putnamlab/jillashey/e5/output/apul_xpo5_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l apul_xpo5_blastp.tab

echo "xpo5 for Ptuh blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/xpo5.fasta -db /data/putnamlab/jillashey/e5/refs/Pmea_prot -out /data/putnamlab/jillashey/e5/output/ptuh_xpo5_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l ptuh_xpo5_blastp.tab

echo "xpo5 for Peve blastp" $(date)
blastp -query /data/putnamlab/jillashey/e5/ncRNA_prot/xpo5.fasta -db /data/putnamlab/jillashey/e5/refs/Peve_prot -out /data/putnamlab/jillashey/e5/output/peve_xpo5_blastp.tab -evalue 1E-40 -num_threads 10 -max_target_seqs 1 -max_hsps 1 -outfmt 6

wc -l peve_xpo5_blastp.tab

echo "Blasting complete!" $(date)
```

Submitted batch job 311809. Ran super fast but did not count lines. Let's look at the output.

```
cd /data/putnamlab/jillashey/e5/output

ls -othr
total 40K
-rw-r--r--. 1 jillashey 1.9K Apr 11 14:21 apul_ago2_blastp.tab
-rw-r--r--. 1 jillashey 2.8K Apr 11 14:21 ptuh_ago2_blastp.tab
-rw-r--r--. 1 jillashey 2.0K Apr 11 14:21 peve_ago2_blastp.tab
-rw-r--r--. 1 jillashey  493 Apr 11 14:21 apul_dgcr8_blastp.tab
-rw-r--r--. 1 jillashey  720 Apr 11 14:21 ptuh_dgcr8_blastp.tab
-rw-r--r--. 1 jillashey  479 Apr 11 14:21 peve_dgcr8_blastp.tab
-rw-r--r--. 1 jillashey 1.4K Apr 11 14:21 apul_dicer_blastp.tab
-rw-r--r--. 1 jillashey 1.9K Apr 11 14:21 ptuh_dicer_blastp.tab
-rw-r--r--. 1 jillashey 1.4K Apr 11 14:21 peve_dicer_blastp.tab
-rw-r--r--. 1 jillashey  482 Apr 11 14:21 apul_dnmt1_blastp.tab
-rw-r--r--. 1 jillashey  798 Apr 11 14:21 ptuh_dnmt1_blastp.tab
-rw-r--r--. 1 jillashey  563 Apr 11 14:21 peve_dnmt1_blastp.tab
-rw-r--r--. 1 jillashey  433 Apr 11 14:21 apul_dnmt3a_blastp.tab
-rw-r--r--. 1 jillashey  596 Apr 11 14:21 ptuh_dnmt3a_blastp.tab
-rw-r--r--. 1 jillashey  428 Apr 11 14:21 peve_dnmt3a_blastp.tab
-rw-r--r--. 1 jillashey  511 Apr 11 14:21 apul_drosha_blastp.tab
-rw-r--r--. 1 jillashey  722 Apr 11 14:21 ptuh_drosha_blastp.tab
-rw-r--r--. 1 jillashey  510 Apr 11 14:21 peve_drosha_blastp.tab
-rw-r--r--. 1 jillashey  680 Apr 11 14:21 apul_pip5k1a_blastp.tab
-rw-r--r--. 1 jillashey  988 Apr 11 14:21 ptuh_pip5k1a_blastp.tab
-rw-r--r--. 1 jillashey  681 Apr 11 14:21 peve_pip5k1a_blastp.tab
-rw-r--r--. 1 jillashey  830 Apr 11 14:21 apul_piwi_blastp.tab
-rw-r--r--. 1 jillashey 1.2K Apr 11 14:21 ptuh_piwi_blastp.tab
-rw-r--r--. 1 jillashey  842 Apr 11 14:21 peve_piwi_blastp.tab
-rw-r--r--. 1 jillashey 2.4K Apr 11 14:21 apul_rnase_p_blastp.tab
-rw-r--r--. 1 jillashey 3.7K Apr 11 14:21 ptuh_rnase_p_blastp.tab
-rw-r--r--. 1 jillashey 2.5K Apr 11 14:21 peve_rnase_p_blastp.tab
-rw-r--r--. 1 jillashey 1.1K Apr 11 14:21 apul_xpo5_blastp.tab
-rw-r--r--. 1 jillashey 1.4K Apr 11 14:21 ptuh_xpo5_blastp.tab
-rw-r--r--. 1 jillashey  933 Apr 11 14:22 peve_xpo5_blastp.tab

wc -l *
   28 apul_ago2_blastp.tab
    7 apul_dgcr8_blastp.tab
   19 apul_dicer_blastp.tab
    7 apul_dnmt1_blastp.tab
    6 apul_dnmt3a_blastp.tab
    7 apul_drosha_blastp.tab
   10 apul_pip5k1a_blastp.tab
   12 apul_piwi_blastp.tab
   34 apul_rnase_p_blastp.tab
   16 apul_xpo5_blastp.tab
   28 peve_ago2_blastp.tab
    7 peve_dgcr8_blastp.tab
   19 peve_dicer_blastp.tab
    8 peve_dnmt1_blastp.tab
    6 peve_dnmt3a_blastp.tab
    7 peve_drosha_blastp.tab
   10 peve_pip5k1a_blastp.tab
   12 peve_piwi_blastp.tab
   36 peve_rnase_p_blastp.tab
   13 peve_xpo5_blastp.tab
   28 ptuh_ago2_blastp.tab
    7 ptuh_dgcr8_blastp.tab
   19 ptuh_dicer_blastp.tab
    8 ptuh_dnmt1_blastp.tab
    6 ptuh_dnmt3a_blastp.tab
    7 ptuh_drosha_blastp.tab
   10 ptuh_pip5k1a_blastp.tab
   12 ptuh_piwi_blastp.tab
   38 ptuh_rnase_p_blastp.tab
   14 ptuh_xpo5_blastp.tab
  441 total
```

Is there a way to make some kind of phylogeny?



