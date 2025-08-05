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
  
head apul_ago2_blastp.tab
XP_022791662.1	XP_044181953.1	65.544	891	291	8	141	1030	133	1008	0.0	1206
PFX12762.1	XP_044181953.1	60.216	832	259	6	315	1079	182	1008	0.0	1031
XP_022809216.1	XP_044181220.1	71.823	362	97	4	2	361	1	359	0.0	511
XP_020617412.1	XP_044181953.1	72.414	493	130	3	1	493	522	1008	0.0	744
XP_020617413.1	XP_044181953.1	74.747	495	117	3	1	495	522	1008	0.0	766
XP_020617479.1	XP_044181952.1	54.524	431	175	8	109	538	73	483	1.73e-147	450
XP_020625223.1	XP_044171946.1	89.344	122	13	0	1	122	1	122	1.77e-78	235
XP_058954265.1	XP_044181953.1	63.545	993	326	11	48	1038	50	1008	0.0	1268
XP_058951070.1	XP_044181952.1	71.528	144	40	1	1	144	674	816	7.32e-62	210
XP_058951065.1	XP_044181068.1	88.235	425	50	0	1	425	430	854	0.0	758
```

The first column is the subject sequence IDs (the ones that I compiled) and the second column is the query IDs (the species of interest). Select the second column from each file and remove any duplicates 

```
awk '{print $2}' apul_ago2_blastp.tab | sort | uniq > apul_ago2_genelist.txt
```

Use the gene list to subset the protein fasta file for the sequences of interest

```
grep -A 1 -f apul_ago2_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/GCF_013753865.1_Amil_v2.1.protein.faa | grep -v "^--$" > apul_ago2.fasta
```

I could write a for loop for this...but I am lazy and will just do it manually. 

##### Ago2

Apul 

```
awk '{print $2}' apul_ago2_blastp.tab | sort | uniq > apul_ago2_genelist.txt

grep -A 1 -f apul_ago2_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/GCF_013753865.1_Amil_v2.1.protein.faa | grep -v "^--$" > apul_ago2.fasta
```

Peve

```
awk '{print $2}' peve_ago2_blastp.tab | sort | uniq > peve_ago2_genelist.txt

grep -A 1 -f peve_ago2_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/Porites_evermanni_v1.annot.pep.fa | grep -v "^--$" > peve_ago2.fasta
```

Ptuh 

```
awk '{print $2}' ptuh_ago2_blastp.tab | sort | uniq > ptuh_ago2_genelist.txt

grep -A 1 -f ptuh_ago2_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/Pocillopora_meandrina_HIv1.genes.pep.faa | grep -v "^--$" > ptuh_ago2.fasta
```

##### DGCR8

Apul 

```
awk '{print $2}' apul_dgcr8_blastp.tab | sort | uniq > apul_dgcr8_genelist.txt

grep -A 1 -f apul_dgcr8_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/GCF_013753865.1_Amil_v2.1.protein.faa | grep -v "^--$" > apul_dgcr8.fasta
```

Peve

```
awk '{print $2}' peve_dgcr8_blastp.tab | sort | uniq > peve_dgcr8_genelist.txt

grep -A 1 -f peve_dgcr8_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/Porites_evermanni_v1.annot.pep.fa | grep -v "^--$" > peve_dgcr8.fasta
```

Ptuh 

```
awk '{print $2}' ptuh_dgcr8_blastp.tab | sort | uniq > ptuh_dgcr8_genelist.txt

grep -A 1 -f ptuh_dgcr8_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/Pocillopora_meandrina_HIv1.genes.pep.faa | grep -v "^--$" > ptuh_dgcr8.fasta
```

##### Dicer

Apul 

```
awk '{print $2}' apul_dicer_blastp.tab | sort | uniq > apul_dicer_genelist.txt

grep -A 1 -f apul_dicer_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/GCF_013753865.1_Amil_v2.1.protein.faa | grep -v "^--$" > apul_dicer.fasta
```

Peve

```
awk '{print $2}' peve_dicer_blastp.tab | sort | uniq > peve_dicer_genelist.txt

grep -A 1 -f peve_dicer_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/Porites_evermanni_v1.annot.pep.fa | grep -v "^--$" > peve_dicer.fasta
```

Ptuh 

```
awk '{print $2}' ptuh_dicer_blastp.tab | sort | uniq > ptuh_dicer_genelist.txt

grep -A 1 -f ptuh_dicer_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/Pocillopora_meandrina_HIv1.genes.pep.faa | grep -v "^--$" > ptuh_dicer.fasta
```

##### DNMT1

Apul 

```
awk '{print $2}' apul_dnmt1_blastp.tab | sort | uniq > apul_dnmt1_genelist.txt

grep -A 1 -f apul_dnmt1_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/GCF_013753865.1_Amil_v2.1.protein.faa | grep -v "^--$" > apul_dnmt1.fasta
```

Peve

```
awk '{print $2}' peve_dnmt1_blastp.tab | sort | uniq > peve_dnmt1_genelist.txt

grep -A 1 -f peve_dnmt1_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/Porites_evermanni_v1.annot.pep.fa | grep -v "^--$" > peve_dnmt1.fasta
```

Ptuh 

```
awk '{print $2}' ptuh_dnmt1_blastp.tab | sort | uniq > ptuh_dnmt1_genelist.txt

grep -A 1 -f ptuh_dnmt1_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/Pocillopora_meandrina_HIv1.genes.pep.faa | grep -v "^--$" > ptuh_dnmt1.fasta
```

##### DNMT3A

Apul 

```
awk '{print $2}' apul_dnmt3a_blastp.tab | sort | uniq > apul_dnmt3a_genelist.txt

grep -A 1 -f apul_dnmt3a_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/GCF_013753865.1_Amil_v2.1.protein.faa | grep -v "^--$" > apul_dnmt3a.fasta
```

Peve

```
awk '{print $2}' peve_dnmt3a_blastp.tab | sort | uniq > peve_dnmt3a_genelist.txt

grep -A 1 -f peve_dnmt3a_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/Porites_evermanni_v1.annot.pep.fa | grep -v "^--$" > peve_dnmt3a.fasta
```

Ptuh 

```
awk '{print $2}' ptuh_dnmt3a_blastp.tab | sort | uniq > ptuh_dnmt3a_genelist.txt

grep -A 1 -f ptuh_dnmt3a_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/Pocillopora_meandrina_HIv1.genes.pep.faa | grep -v "^--$" > ptuh_dnmt3a.fasta
```

##### Drosha

Apul 

```
awk '{print $2}' apul_drosha_blastp.tab | sort | uniq > apul_drosha_genelist.txt

grep -A 1 -f apul_drosha_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/GCF_013753865.1_Amil_v2.1.protein.faa | grep -v "^--$" > apul_drosha.fasta
```

Peve

```
awk '{print $2}' peve_drosha_blastp.tab | sort | uniq > peve_drosha_genelist.txt

grep -A 1 -f peve_drosha_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/Porites_evermanni_v1.annot.pep.fa | grep -v "^--$" > peve_drosha.fasta
```

Ptuh 

```
awk '{print $2}' ptuh_drosha_blastp.tab | sort | uniq > ptuh_drosha_genelist.txt

grep -A 1 -f ptuh_drosha_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/Pocillopora_meandrina_HIv1.genes.pep.faa | grep -v "^--$" > ptuh_drosha.fasta
```

##### Pip5k1a

Apul 

```
awk '{print $2}' apul_pip5k1a_blastp.tab | sort | uniq > apul_pip5k1a_genelist.txt

grep -A 1 -f apul_pip5k1a_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/GCF_013753865.1_Amil_v2.1.protein.faa | grep -v "^--$" > apul_pip5k1a.fasta
```

Peve

```
awk '{print $2}' peve_pip5k1a_blastp.tab | sort | uniq > peve_pip5k1a_genelist.txt

grep -A 1 -f peve_pip5k1a_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/Porites_evermanni_v1.annot.pep.fa | grep -v "^--$" > peve_pip5k1a.fasta
```

Ptuh 

```
awk '{print $2}' ptuh_pip5k1a_blastp.tab | sort | uniq > ptuh_pip5k1a_genelist.txt

grep -A 1 -f ptuh_pip5k1a_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/Pocillopora_meandrina_HIv1.genes.pep.faa | grep -v "^--$" > ptuh_pip5k1a.fasta
```

##### Piwi

Apul 

```
awk '{print $2}' apul_piwi_blastp.tab | sort | uniq > apul_piwi_genelist.txt

grep -A 1 -f apul_piwi_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/GCF_013753865.1_Amil_v2.1.protein.faa | grep -v "^--$" > apul_piwi.fasta
```

Peve

```
awk '{print $2}' peve_piwi_blastp.tab | sort | uniq > peve_piwi_genelist.txt

grep -A 1 -f peve_piwi_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/Porites_evermanni_v1.annot.pep.fa | grep -v "^--$" > peve_piwi.fasta
```

Ptuh 

```
awk '{print $2}' ptuh_piwi_blastp.tab | sort | uniq > ptuh_piwi_genelist.txt

grep -A 1 -f ptuh_piwi_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/Pocillopora_meandrina_HIv1.genes.pep.faa | grep -v "^--$" > ptuh_piwi.fasta
```

##### RNase P

Apul 

```
awk '{print $2}' apul_rnase_p_blastp.tab | sort | uniq > apul_rnase_p_genelist.txt

grep -A 1 -f apul_rnase_p_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/GCF_013753865.1_Amil_v2.1.protein.faa | grep -v "^--$" > apul_rnase_p.fasta
```

Peve

```
awk '{print $2}' peve_rnase_p_blastp.tab | sort | uniq > peve_rnase_p_genelist.txt

grep -A 1 -f peve_rnase_p_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/Porites_evermanni_v1.annot.pep.fa | grep -v "^--$" > peve_rnase_p.fasta
```

Ptuh 

```
awk '{print $2}' ptuh_rnase_p_blastp.tab | sort | uniq > ptuh_rnase_p_genelist.txt

grep -A 1 -f ptuh_rnase_p_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/Pocillopora_meandrina_HIv1.genes.pep.faa | grep -v "^--$" > ptuh_rnase_p.fasta
```

##### Xpo5

Apul 

```
awk '{print $2}' apul_xpo5_blastp.tab | sort | uniq > apul_xpo5_genelist.txt

grep -A 1 -f apul_xpo5_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/GCF_013753865.1_Amil_v2.1.protein.faa | grep -v "^--$" > apul_xpo5.fasta
```

Peve

```
awk '{print $2}' peve_xpo5_blastp.tab | sort | uniq > peve_xpo5_genelist.txt

grep -A 1 -f peve_xpo5_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/Porites_evermanni_v1.annot.pep.fa | grep -v "^--$" > peve_xpo5.fasta
```

Ptuh 

```
awk '{print $2}' ptuh_xpo5_blastp.tab | sort | uniq > ptuh_xpo5_genelist.txt

grep -A 1 -f ptuh_xpo5_genelist.txt /data/putnamlab/jillashey/e5/ortho/protein_seqs/Pocillopora_meandrina_HIv1.genes.pep.faa | grep -v "^--$" > ptuh_xpo5.fasta
```

Copy the fasta files onto my local computer. 

I will use the following programs to analyze these data: 

- [Muscle](https://www.ebi.ac.uk/jdispatcher/msa/muscle)
- Jalview (application downloaded to my computer)

Links for the muscle alignment to use w/ Jalview. Copy these links into the Jalview application.

- Ago2: https://www.ebi.ac.uk/Tools/services/rest/muscle/result/muscle-I20240507-025347-0712-98646927-p1m/aln-clustalw
- DGCR8: https://www.ebi.ac.uk/Tools/services/rest/muscle/result/muscle-I20240507-024813-0571-38010587-p1m/aln-clustalw
- Dicer: https://www.ebi.ac.uk/Tools/services/rest/muscle/result/muscle-I20240507-025605-0874-10941366-p1m/aln-clustalw
- DNMT1: https://www.ebi.ac.uk/Tools/services/rest/muscle/result/muscle-I20240507-025843-0531-38950325-p1m/aln-clustalw
- DNMT3a: https://www.ebi.ac.uk/Tools/services/rest/muscle/result/muscle-I20240507-030001-0467-66490567-p1m/aln-clustalw
- Drosha: https://www.ebi.ac.uk/Tools/services/rest/muscle/result/muscle-I20240507-030250-0520-91955590-p1m/aln-clustalw
- Pip5k1a: https://www.ebi.ac.uk/Tools/services/rest/muscle/result/muscle-I20240507-030359-0442-62910879-p1m/aln-clustalw
- Piwi: https://www.ebi.ac.uk/Tools/services/rest/muscle/result/muscle-I20240507-030503-0247-70454786-p1m/aln-clustalw
- RNase P: https://www.ebi.ac.uk/Tools/services/rest/muscle/result/muscle-I20240507-030615-0015-90999376-p1m/aln-clustalw
- Xpo5: https://www.ebi.ac.uk/Tools/services/rest/muscle/result/muscle-I20240507-030715-0512-87536698-p1m/aln-clustalw

### Redo Apul using new Apul genome 

ON UNITY! 

Copied the fasta files of the query seqs onto Unity. Now build an Apul protein db and blast the fastas against the db. 

`nano apul_ncRNA_prot_blastp.sh`

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --partition=uri-cpu
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 100:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error
#SBATCH -D /work/pi_hputnam_uri_edu/jillashey/e5/ncRNA_prots

# Load modules 
module load uri/main
module load BLAST+/2.15.0-gompi-2023a

cd /work/pi_hputnam_uri_edu/jillashey/e5/ncRNA_prots

echo "Making Apul blast db for e5 deep dive expression protein seqs" $(date)

makeblastdb -in /work/pi_hputnam_uri_edu/jillashey/e5/Apulchra-genome.pep.faa -out /work/pi_hputnam_uri_edu/jillashey/e5/ncRNA_prots/Apul_prot -dbtype prot

echo "Apul Blast db creation complete" $(date)

echo "AGO2 for Apul blastp" $(date)
blastp -query ago2.fasta -db Apul_prot -out apul_ago2_blastp.tab -evalue 1E-40 -num_threads 15 -max_target_seqs 1 -max_hsps 1 -outfmt 6

echo "DGCR8 for Apul blastp" $(date)
blastp -query dgcr8.fasta -db Apul_prot -out apul_dgcr8_blastp.tab -evalue 1E-40 -num_threads 15 -max_target_seqs 1 -max_hsps 1 -outfmt 6

echo "Dicer for Apul blastp" $(date)
blastp -query dicer.fasta -db Apul_prot -out apul_dicer_blastp.tab -evalue 1E-40 -num_threads 15 -max_target_seqs 1 -max_hsps 1 -outfmt 6

echo "DNMT1 for Apul blastp" $(date)
blastp -query dnmt1.fasta -db Apul_prot -out apul_dnmt1_blastp.tab -evalue 1E-40 -num_threads 15 -max_target_seqs 1 -max_hsps 1 -outfmt 6

echo "DNMT3a for Apul blastp" $(date)
blastp -query dnmt3a.fasta -db Apul_prot -out apul_dnmt3a_blastp.tab -evalue 1E-40 -num_threads 15 -max_target_seqs 1 -max_hsps 1 -outfmt 6

echo "Drosha for Apul blastp" $(date)
blastp -query drosha.fasta -db Apul_prot -out apul_drosha_blastp.tab -evalue 1E-40 -num_threads 15 -max_target_seqs 1 -max_hsps 1 -outfmt 6

echo "pip5k1a for Apul blastp" $(date)
blastp -query pip5k1a.fasta -db Apul_prot -out apul_pip5k1a_blastp.tab -evalue 1E-40 -num_threads 15 -max_target_seqs 1 -max_hsps 1 -outfmt 6

echo "piwi for Apul blastp" $(date)
blastp -query piwi.fasta -db Apul_prot -out apul_piwi_blastp.tab -evalue 1E-40 -num_threads 15 -max_target_seqs 1 -max_hsps 1 -outfmt 6

echo "rnase p for Apul blastp" $(date)
blastp -query rnase_p.fasta -db Apul_prot -out apul_rnase_p_blastp.tab -evalue 1E-40 -num_threads 15 -max_target_seqs 1 -max_hsps 1 -outfmt 6

echo "xpo5 for Apul blastp" $(date)
blastp -query xpo5.fasta -db Apul_prot -out apul_xpo5_blastp.tab -evalue 1E-40 -num_threads 15 -max_target_seqs 1 -max_hsps 1 -outfmt 6

echo "Blast complete!" $(date)
```

Submitted batch job 40369015

Similar to above, select the second column from each file and remove any duplicates. Then use the gene list to subset the protein fasta for the seqs of interest 

```
# ago2 
awk '{print $2}' apul_ago2_blastp.tab | sort | uniq > apul_ago2_genelist.txt
grep -A 1 -f apul_ago2_genelist.txt /work/pi_hputnam_uri_edu/jillashey/e5/Apulchra-genome.pep.faa | grep -v "^--$" > apul_ago2.fasta

# dgcr8
awk '{print $2}' apul_dgcr8_blastp.tab | sort | uniq > apul_dgcr8_genelist.txt
grep -A 1 -f apul_dgcr8_genelist.txt /work/pi_hputnam_uri_edu/jillashey/e5/Apulchra-genome.pep.faa | grep -v "^--$" > apul_dgcr8.fasta

# dicer
awk '{print $2}' apul_dicer_blastp.tab | sort | uniq > apul_dicer_genelist.txt
grep -A 1 -f apul_dicer_genelist.txt /work/pi_hputnam_uri_edu/jillashey/e5/Apulchra-genome.pep.faa | grep -v "^--$" > apul_dicer.fasta

# dnmt1
awk '{print $2}' apul_dnmt1_blastp.tab | sort | uniq > apul_dnmt1_genelist.txt
grep -A 1 -f apul_dnmt1_genelist.txt /work/pi_hputnam_uri_edu/jillashey/e5/Apulchra-genome.pep.faa | grep -v "^--$" > apul_dnmt1.fasta

# dnmt3a
awk '{print $2}' apul_dnmt3a_blastp.tab | sort | uniq > apul_dnmt3a_genelist.txt
grep -A 1 -f apul_dnmt3a_genelist.txt /work/pi_hputnam_uri_edu/jillashey/e5/Apulchra-genome.pep.faa | grep -v "^--$" > apul_dnmt3a.fasta

# drosha
awk '{print $2}' apul_drosha_blastp.tab | sort | uniq > apul_drosha_genelist.txt
grep -A 1 -f apul_drosha_genelist.txt /work/pi_hputnam_uri_edu/jillashey/e5/Apulchra-genome.pep.faa | grep -v "^--$" > apul_drosha.fasta

# pip5k1a
awk '{print $2}' apul_pip5k1a_blastp.tab | sort | uniq > apul_pip5k1a_genelist.txt
grep -A 1 -f apul_pip5k1a_genelist.txt /work/pi_hputnam_uri_edu/jillashey/e5/Apulchra-genome.pep.faa | grep -v "^--$" > apul_pip5k1a.fasta

# piwi
awk '{print $2}' apul_piwi_blastp.tab | sort | uniq > apul_piwi_genelist.txt
grep -A 1 -f apul_piwi_genelist.txt /work/pi_hputnam_uri_edu/jillashey/e5/Apulchra-genome.pep.faa | grep -v "^--$" > apul_piwi.fasta

# rnase p
awk '{print $2}' apul_rnase_p_blastp.tab | sort | uniq > apul_rnase_p_genelist.txt
grep -A 1 -f apul_rnase_p_genelist.txt /work/pi_hputnam_uri_edu/jillashey/e5/Apulchra-genome.pep.faa | grep -v "^--$" > apul_rnase_p.fasta

# xpo5
awk '{print $2}' apul_xpo5_blastp.tab | sort | uniq > apul_xpo5_genelist.txt
grep -A 1 -f apul_xpo5_genelist.txt /work/pi_hputnam_uri_edu/jillashey/e5/Apulchra-genome.pep.faa | grep -v "^--$" > apul_xpo5.fasta
```