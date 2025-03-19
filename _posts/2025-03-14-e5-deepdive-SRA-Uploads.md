---
layout: post
title: e5 SRA uploads to NCBI
date: '2025-03-14'
categories: Protocol
tags: [Protocol, Bioinformatics]
projects: e5 deep dive
---

## e5 deep dive SRA uploads to NCBI

This post details the NCBI Sequence Read Archive upload for the sequence data from the e5 deep dive experiment. The github for that project is [here](https://github.com/urol-e5/deep-dive). This protocol is based on A. Huffmyer's [SRA upload post](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/SRA-Uploads-10-November-2022/) and Putnam Lab [SRA upload protocol](https://github.com/Putnam-Lab/Lab_Management/blob/master/Bioinformatics_%26_Coding/Data_Mangament/SRA-Upload_Protocol.md). This post includes information on both the RNAseq and smRNAseq SRA uploads. 

### Overview for RNAseq uploads 

Below is the following information for this submission. 

Sequence Read Archive (SRA) submission: SUB15175231

#### Submitter

- Jill Ashey; 120 Flagg Road, Kingston RI 02881

#### General info 

- Not associated with current BioProject 
- Not associated with current BioSample
- Release date: immediately 

#### Project info

- Project title: RNAseq in Acropora pulchra, Pocillopora tuahiniensis, and Porites evermanni (Moorea, March 2020)
- Public description: Sample fragments of three coral species (Acropora pulchra, Pocillopora tuahiniensis, and Porites evermanni) were collected from the lagoon backreef in Moorea, French Polynesia on March 5, 2020 on the north shore (-17.476872, -149.80594). RNA-Sequencing was performed by Azenta Life Sciences. Strand-specific RNA sequencing libraries were prepared using NEBNext Ultra II Directional RNA Library Prep Kit for Illumina following the manufacturer’s instructions. These sequences were used to analyze mRNA and lncRNA expression in these species.
- Not associated with umbrella project

#### Biosample type: invertebrate 

[Link to submission spreadsheet](XXXXXXX)

Important note so that errors do not occur during upload: After filling out values for attributes provided in the template, your samples are not distinguishable by at least one, or a combination of attributes. **Make the combined value of all attributes unique for each sample while taking into account that "sample name", "sample title" and "description" are not included in this check for "uniqueness".**

#### SRA metadata 

LINK

#### Files 

Upload the fastq.gz files to NCBI with the following code. First, make a new folder that only contains sequence data.

```
cd /data/putnamlab/jillashey/e5
mkdir data 
cd data 

wget XXXX
```

Copy sequences from Roberts server (ACR [location](https://owl.fish.washington.edu/nightingales/A_pulchra/30-789513166/); POC [location](https://owl.fish.washington.edu/nightingales/P_meandrina/30-789513166/); POR [location](https://owl.fish.washington.edu/nightingales/P_evermanni/30-852430235/)). 

Use FTP command line file upload to provide files. Activate FTP on the command line. Command prompt now changes to ftp>. Immediately login using the username and password from the NCBI submission portal.

```
ftp -i
open ftp-private.ncbi.nlm.nih.gov
USERNAME
PASSWORD 
```

Go into the folder that NCBI provided from the NCBI submission portal. Make a new directory for the files and go into that folder. Put all the sequences from Andromeda into this new folder (they will still be on Andromeda after transfer is complete).

```
cd uploads/jillashey_uri.edu_XXXXXXXXXXX
mkdir e5_rnaseq
cd e5_rnaseq
mput *
```

BioProject ID: PRJNA1236658

### Overview for smRNAseq uploads 

Below is the following information for this submission. 

Sequence Read Archive (SRA) submission: SUB15177197

#### Submitter

- Jill Ashey; 120 Flagg Road, Kingston RI 02881

#### General info 

- Not associated with current BioProject 
- Not associated with current BioSample
- Release date: immediately 

#### Project info

- Project title: smRNAseq in Acropora pulchra, Pocillopora tuahiniensis, and Porites evermanni (Moorea, March 2020)
- Public description: Sample fragments of three coral species (Acropora pulchra, Pocillopora tuahiniensis, and Porites evermanni) were collected from the lagoon backreef in Moorea, French Polynesia on March 5, 2020 on the north shore (-17.476872, -149.80594). Small RNA sequencing libraries were prepared using NEB Small RNA Library Prep Kit (NEB CAT: E7560S) following the manufacturer’s instructions. These reads were used to analyze piRNAs and miRNAs. 
- Not associated with umbrella project

#### Biosample type: invertebrate 

[Link to submission spreadsheet](XXXXXXX)

Important note so that errors do not occur during upload: After filling out values for attributes provided in the template, your samples are not distinguishable by at least one, or a combination of attributes. **Make the combined value of all attributes unique for each sample while taking into account that "sample name", "sample title" and "description" are not included in this check for "uniqueness".**

#### SRA metadata 

LINK

#### Files 

Upload the fastq.gz files to NCBI with the following code. First, make a new folder that only contains sequence data.

```
cd /data/putnamlab/jillashey/e5
mkdir data 
cd data 

wget XXXX
```

Copy smRNA sequences from Roberts server (ACR [location](https://owl.fish.washington.edu/nightingales/A_pulchra/30-852430235/); POC [location](https://owl.fish.washington.edu/nightingales/P_meandrina/30-852430235/); POR [location](https://owl.fish.washington.edu/nightingales/P_evermanni/30-852430235/)). 

Use FTP command line file upload to provide files. Activate FTP on the command line. Command prompt now changes to ftp>. Immediately login using the username and password from the NCBI submission portal.

```
ftp -i
open ftp-private.ncbi.nlm.nih.gov
USERNAME
PASSWORD 
```

Go into the folder that NCBI provided from the NCBI submission portal. Make a new directory for the files and go into that folder. Put all the sequences from Andromeda into this new folder (they will still be on Andromeda after transfer is complete).

```
cd uploads/jillashey_uri.edu_XXXXXXXXXXX
mkdir e5_smrnaseq
cd e5_smrnaseq
mput *
```
