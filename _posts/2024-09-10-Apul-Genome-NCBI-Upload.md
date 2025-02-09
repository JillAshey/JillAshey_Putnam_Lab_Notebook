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

On 10/25/24, I received an email from NCBI saying that they have received my submission and assigned it an accession number. The submission has passed initial QC checks and will now be manually reviewed by the indexing staff. Once that is complete, the genome will be released immediately. 

From NCBI:

```
We have assigned the following accession number to your submission:

SUBID           BioProject      BioSample       Accession       Organism
SUB14718394     PRJNA1162071    SAMN43800006    JBIQNU000000000 Acropora pulchra JA

Please cite the accession number JBIQNU000000000 like this:

This Whole Genome Shotgun project has been deposited at DDBJ/ENA/GenBank
under the accession JBIQNU000000000. The version described
in this paper is version JBIQNU010000000.

While you can cite the BioProject, we recommend you include the WGS
accession(s) to refer to the specific WGS assembly, especially for
BioProjects that include more than one genome assembly.  Note that we
prefer not to change the BioProject/BioSample links after the WGS assembly
is released.
```

## Raw PacBio sequencing data submission to NCBI

Gigabyte requires that our raw PacBio data also be stored on NCBI. I'm going to start a new SRA submission and link it to the existing BioProject (PRJNA1162071). 

### Overview 

- Submission: SUB15071615
- Submitting raw PacBio Hifi reads  

### Submitter 

- Jill Ashey; 120 Flagg Road, Kingston RI 02881

### General info 

- Associated with existing bioproject: PRJNA1162071
- Release date: Immediately 
- Not associated with existing biosample 

### Biosample type 

- Invertebrate 
- Sample name: Acropora_pulchra_raw_reads_JA
- BioProject accession: PRJNA1162071
- Organism: Acropora pulchra 
- Isolate: Coral host 
- Isolation source: Coral reef 
- collection date: 2022-10-23
- Location: French Polynesia: Moorea
- Tissue: sperm
- Developmental stage: adult 
- Broad scale env context: coral reef [ENVO:00000150]
- Sex: hermaphrodite 

### SRA metadata 

- Sample name: Acropora_pulchra_raw_reads_JA
- Library ID: m84100_240128_024355_s2
- Title: PacBio raw HiFi seqences of Acropora pulchra
- Library strategy (drop down options): WGS
- Library source (drop down options): Genomic 
- Library selection (drop down options): size fractionation 
- Library layout (drop down options): single 
- Platform (drop down options): PACBIO_SMRT
- Instrument model (drop down options): Revio 
- Design description: DNA was extracted by DNA Sequencing Center at Brigham Young University using the Qiagen Genomic Tip protocol and buffers (Qiagen Cat \# 10223). The samples were ethanol (2x) precipitated post column elution, put in the -20°C freezer overnight and then spun for 30 minutes at 14K rcf the following day. Ethanol was removed and the DNA pellets were suspended in low TE buffer. The resulting DNA was cleaned prior to library prep with the PacBio SRE Kit (PacBio Cat \# 102-208-300) to remove fragments under 25kb. Following extraction, DNA was sheared to ~17kb using a Diagenode Megaruptor (Diagenode Cat \# B06010003) and checked on an Agilent Femto Pulse system (Agilent Part \# M5330AA) to assess size. The DNA was then cleaned and concentrated post-shearing using a 1x AMPure bead cleaning (AMPure Cat \#20805800). The DNA was then put into a library using the PacBio SMRTbell prep kit 3.0 (PacBio Cat \# 102-141-700), following the instructions provided with the kit. The final sizing of the library was performed using the 35\% v/v dilution of AMPure PB beads (AMPure Part \#  100-265-900). The single SMRTbell library was then sequenced using one 8M SMRT Revio Cell, and run for 29 hours on a PacBio Revio sequencer. Consensus accuracy circular consensus sequencing (CCS) processing was used to generate HiFi reads. -- from methods section of paper 
- File type: bam
- Reference assembly: unaligned
- File name: m84100_240128_024355_s2.hifi_reads.bc1029.bam

### Files 

Command line upload 

Make new file with only sequence data 

```
cd /data/putnamlab/jillashey/Apul_Genome/ncbi
ln -s /data/putnamlab/KITT/hputnam/20240129_Apulchra_Genome_LongRead/m84100_240128_024355_s2.hifi_reads.bc1029.bam
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
mkdir apul_raw_bam_2024
cd apul_raw_bam_2024
mput m84100_240128_024355_s2.hifi_reads.bc1029.bam
local: m84100_240128_024355_s2.hifi_reads.bc1029.bam remote: m84100_240128_024355_s2.hifi_reads.bc1029.bam
227 Entering Passive Mode (130,14,250,6,195,171).
150 Opening BINARY mode data connection for m84100_240128_024355_s2.hifi_reads.bc1029.bam
226 Transfer complete
38666252812 bytes sent in 536 secs (72146.58 Kbytes/sec)
```

## Submit!

SAMN46708315 - assigned biosample number. SRA currently processing submission as of 2/6/24. Downloaded attributes file and will put on github repo. 
