---
layout: post
title: Pacuta SRA uploads to NCBI
date: '2024-05-01'
categories: Protocol
tags: [Protocol, Bioinformatics]
projects: Pacuta HI 2022
---

## Genome submission to NCBI

This post details the NCBI Genome Submission upload for the assembled *Acropora pulchra* genome. The github for that project is [here](https://github.com/hputnam/Apulchra_genome). More information on genome submission on NCBI can be found [here](https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/). Below is the following information for this submission. 

### Overview 

- Genome submission: SUB14718394
- Submitting single genome 

### Submitter 

- Jill Ashey; 120 Flagg Road, Kingston RI 02881

### General info 

- Not associated with existing bioproject or biosamples 
- Release date: 2024-09-30
- Assembly date: 2024-07
- Assembly method: Hifiasm, run in JULY 2024
- Assembly name: Apul_v1.1
- Genome coverage: 100
- PacBio sequencing technology 
- Only submitting pacbio reads 
- Sample is the full genome 
- It is the final version - will not be doing re-assembly 
- De novo assembly 
- Do not automatically trim or remove sequences identified as contamination
- Submission category: original 
- Submission title: Acropora pulchra NCBI genome submission

### Bioproject general info 

Public description (4000 characters): This submission provides the genome assembly of Acropora pulchra, a scleractinian coral, from Moorea, French Polynesia. 

Links 
- OSF https://osf.io/y8963/
- Github https://github.com/hputnam/Apulchra_genome

### Biosample type 

- Invertebrate 
- Sample name: Acropora_pulchra_JA
- Organism: Acropora pulchra 
- Isolate: Coral host 
- Isolation source: Coral reef 
- collection date: 2022-10-23
- Location: French Polynesia: Moorea
- Tissue: sperm
- Developmental stage: adult 
- Broad scale env context: coral reef [ENVO:00000150]
- Sex: hermaphrodite 

### Files 

One or more chromosomes are still in multiple pieces and/or some sequences are not assembled into chromosomes

Command line upload 

Make new file with only sequence data 

```
cd /data/putnamlab/jillashey/Apul_Genome/ncbi
ln -s /data/putnamlab/tconn/repeats/apul_softmasked/apul.hifiasm.s55_pa.p_ctg.fa.k32.w100.z1000.ntLink.5rounds.fa.masked
```

Use FTP command line file upload to provide files. Activate FTP on the command line. Command prompt now changes to ftp>. Immediately login using the username and password from the NCBI submission portal.

```
ftp -i
open ftp-private.ncbi.nlm.nih.gov
USERNAME
PASSWORD
```

Go into the folder that NCBI provided from the NCBI submission portal. Make a new directory for the files and go into that folder. Put all the sequences from Andromeda into this new folder (they will still be on Andromeda after transfer is complete).

```
cd uploads/jillashey_uri.edu_XXXXX
mkdir apul_2024
cd apul_2024
mput *
local: apul.hifiasm.s55_pa.p_ctg.fa.k32.w100.z1000.ntLink.5rounds.fa.masked remote: apul.hifiasm.s55_pa.p_ctg.fa.k32.w100.z1000.ntLink.5rounds.fa.masked
227 Entering Passive Mode (130,14,250,5,196,151).
150 Opening BINARY mode data connection for apul.hifiasm.s55_pa.p_ctg.fa.k32.w100.z1000.ntLink.5rounds.fa.masked
226 Transfer complete
528682350 bytes sent in 7.06 secs (74889.03 Kbytes/sec)
```

Note: it takes at least 10 minutes for uploaded files to become available for selection within a submission.

Once transfers are complete, click select preload folder in the submission portal. Wait until all files have been uploaded and select the pacuta_2022 folder. Click submit!

Fasta contigs: `apul.hifiasm.s55_pa.p_ctg.fa.k32.w100.z1000.ntLink.5rounds.fa.masked`

### Assignment 

- Do any sequences belong to a chromosome? No
- Do any sequences belong to an organelle, eg mitochondrion or chloroplast? No
- Does any sequence belong to a plasmid? No 

### References 

- Sequence Authors: Jill Ashey 
- Reference
	- Unpublished 
	- Title: Genome assembly and annotation of Acropora pulchra from Moorea, French Polynesia 
	- Reference authors: Jill Ashey, Trinity Conn, Ross Cunning, Hollie M. Putnam

### Submission

BioProject: Processed
- PRJNA1162071 : Acropora pulchra genome sequencing (TaxID: 140239)
- Locus Tag Prefixes: ACE5DV

BioSample: Processed 
- Successfully loaded 
- SAMN43800006: Acropora_pulchra_JA (TaxID: 140239)
