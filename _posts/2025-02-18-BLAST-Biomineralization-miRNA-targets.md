---
layout: post
title: BLAST biomin genes 
date: '2025-02-18'
categories: Analysis
tags: [Bioinformatics, BLAST]
projects: 
---

## Potential biomineralization targets of miRNAs

Investigating if any biomineralization genes are targeted by miRNAs in four species (AST, ACR, POC, and POR). A biomineralization protein list was created by FS (primarily from Stylophora pistillata). To see if any of the biomin genes are targeted by miRNAs, I am going to blast the biomin protein list against the protein fasta files for each species. 

```
cd /data/putnamlab/jillashey
mkdir biomin_blast
cd biomin_blast 
```

The sequences in fasta format are already on the server here: `/data/putnamlab/jillashey/Pacuta_HI_2022/data/blast`.

`nano biomin_blast.sh`

```
#!/bin/bash
#SBATCH --job-name="biomin_blast"
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --mem=100GB
#SBATCH --error="blast_out_error"
#SBATCH --output="blast_out"
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/jillashey/Pacuta_HI_2022/scripts            
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

module load BLAST+/2.9.0-iimpi-2019b

echo "Make blast dbs for ACR, AST, POC, POR" $(date)

makeblastdb -in /data/putnamlab/jillashey/Astrangia_Genome/apoculata_proteins_v2.0.fasta -out /data/putnamlab/jillashey/biomin_blast/AST_prot -dbtype prot

makeblastdb -in /data/putnamlab/jillashey/genome/Pmea/Pocillopora_meandrina_HIv1.genes.pep.faa -out /data/putnamlab/jillashey/biomin_blast/POC_prot -dbtype prot

makeblastdb -in /data/putnamlab/jillashey/genome/Peve/Porites_evermanni_v1.annot.pep.fa -out /data/putnamlab/jillashey/biomin_blast/POR_prot -dbtype prot

## not sure where ACR protein fasta file is

echo "Blast biomin seqs against protein dbs" $(date)

blastp -query /data/putnamlab/jillashey/Pacuta_HI_2022/data/blast/Biomineraliztion_Toolkit_FScucchia_ZDrefmt.fasta -db /data/putnamlab/jillashey/biomin_blast/AST_prot -out /data/putnamlab/jillashey/biomin_blast/AST_biomin_blast_results_tab.txt -outfmt 6 -max_target_seqs 2

blastp -query /data/putnamlab/jillashey/Pacuta_HI_2022/data/blast/Biomineraliztion_Toolkit_FScucchia_ZDrefmt.fasta -db /data/putnamlab/jillashey/biomin_blast/POC_prot -out /data/putnamlab/jillashey/biomin_blast/POC_biomin_blast_results_tab.txt -outfmt 6 -max_target_seqs 2

blastp -query /data/putnamlab/jillashey/Pacuta_HI_2022/data/blast/Biomineraliztion_Toolkit_FScucchia_ZDrefmt.fasta -db /data/putnamlab/jillashey/biomin_blast/POR_prot -out /data/putnamlab/jillashey/biomin_blast/POR_biomin_blast_results_tab.txt -outfmt 6 -max_target_seqs 2

echo "Blast complete!" $(date)
```

Submitted batch job 361991

Following the BLAST, the results in the .txt file will need to be compared to the miranda results for each species to see if the genes are targeted by any miRNAs. 

- [AST miranda data](https://github.com/JillAshey/Astrangia_repo/blob/main/output/Molecular/interactions/miranda_strict_all_1kb_apoc_shortstack_parsed.txt)
- [POC miranda data](https://github.com/urol-e5/deep-dive-expression/blob/main/F-Ptuh/output/11-Ptuh-mRNA-miRNA-interactions/three_prime_interaction/miranda_strict_all_1kb_parsed_ptuh_updated.txt)
- [POR miranda data](https://github.com/urol-e5/deep-dive-expression/blob/main/E-Peve/output/10-Peve-mRNA-miRNA-interactions/miranda_strict_all_1kb_parsed_peve_updated.txt)

The miranda data is comparing the 3'UTR sequence to the mature miRNA sequence, invoking strict binding. 