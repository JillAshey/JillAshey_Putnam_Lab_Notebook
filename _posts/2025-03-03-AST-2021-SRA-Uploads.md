---
layout: post
title: AST 2021 SRA uploads to NCBI
date: '2025-03-03'
categories: Protocol
tags: [Protocol, Bioinformatics]
projects: AST 2021
---

## AST 2021 SRA uploads to NCBI

This post details the NCBI Sequence Read Archive upload for the sequence data from the AST 2021 experiment. The github for that project is [here](https://github.com/JillAshey/Astrangia_repo). This protocol is based on A. Huffmyer's [SRA upload post](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/SRA-Uploads-10-November-2022/) and Putnam Lab [SRA upload protocol](https://github.com/Putnam-Lab/Lab_Management/blob/master/Bioinformatics_%26_Coding/Data_Mangament/SRA-Upload_Protocol.md). This post includes information on both the RNAseq and smRNAseq SRA uploads. 

Because sample AST-1105 did not sequence well, I am not including it in the data upload. The file itself seems to be corrupt. 

### Overview for RNAseq uploads 

Below is the following information for this submission. 

Sequence Read Archive (SRA) submission: SUB15150405

#### Submitter

- Jill Ashey; 120 Flagg Road, Kingston RI 02881

#### General info 

- Not associated with current BioProject 
- Not associated with current BioSample
- Release date: immediately 

#### Project info

- Project title: Astrangia poculata RNA-seq across seasons
- Public description: RNA-seq data from Astrangia poculata 2021 experiment exposing adults to chronic thermal stress across seasons (February - August)
- Not associated with umbrella project

#### Biosample type: invertebrate 

[Link to submission spreadsheet](https://github.com/JillAshey/Astrangia_repo/blob/main/data/NCBI/AST_invertebrate.xlsx)

Important note so that errors do not occur during upload: After filling out values for attributes provided in the template, your samples are not distinguishable by at least one, or a combination of attributes. **Make the combined value of all attributes unique for each sample while taking into account that "sample name", "sample title" and "description" are not included in this check for "uniqueness".**

#### SRA metadata

Complete 

#### Files 

Upload the fastq.gz files to NCBI with the following code. First, make a new folder that only contains sequence data.

```
cd /data/putnamlab/jillashey/Astrangia2021/mRNA
mkdir raw_file_rnaseq_sra
cd raw_file_rnaseq_sra

ln -s /data/putnamlab/jillashey/Astrangia2021/mRNA/data/raw/*fastq.gz .
rm *1105* # sample that did not sequence well
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
cd uploads/jillashey_uri.edu_XXXXXXXXXXX
mkdir astrangia_2021_rnaseq
cd astrangia_2021_rnaseq
mput *
```

This takes a few minutes. Once transfers are complete, click select preload folder in the submission portal. Wait until all files have been uploaded and select the `astrangia_2021_rnaseq` folder. Click submit!

Bioproject: PRJNA1231118

### Overview for smRNAseq uploads 

Below is the following information for this submission. 

Sequence Read Archive (SRA) submission: SUB15151122

#### Submitter

- Jill Ashey; 120 Flagg Road, Kingston RI 02881

#### General info 

- Not associated with current BioProject 
- Not associated with current BioSample
- Release date: immediately 

#### Project info

- Project title: Astrangia poculata smRNA-seq across seasons
- Public description: smRNA-seq data from Astrangia poculata 2021 experiment exposing adults to chronic thermal stress across seasons (February - August)
- Not associated with umbrella project

#### Biosample type: invertebrate 

[Link to submission spreadsheet](XXXX)

Important note so that errors do not occur during upload: After filling out values for attributes provided in the template, your samples are not distinguishable by at least one, or a combination of attributes. **Make the combined value of all attributes unique for each sample while taking into account that "sample name", "sample title" and "description" are not included in this check for "uniqueness".**

#### SRA metadata

Complete 

#### Files 

Upload the fastq.gz files to NCBI with the following code. First, make a new folder that only contains sequence data.

```
cd /data/putnamlab/jillashey/Astrangia2021/smRNA
mkdir raw_file_smrnaseq_sra
cd raw_file_smrnaseq_sra

ln -s /data/putnamlab/jillashey/Astrangia2021/smRNA/data/raw/*fastq.gz .
rm *1105* # sample that did not sequence well
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
cd uploads/jillashey_uri.edu_XXXXXXXXXXX
mkdir astrangia_2021_smrnaseq
cd astrangia_2021_smrnaseq
mput *
```

This takes a few minutes. Once transfers are complete, click select preload folder in the submission portal. Wait until all files have been uploaded and select the `astrangia_2021_smrnaseq` folder. Click submit!

Bioproject: PRJNA1231129