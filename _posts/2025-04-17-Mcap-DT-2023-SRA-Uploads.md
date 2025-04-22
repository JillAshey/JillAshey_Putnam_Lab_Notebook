---
layout: post
title: Mcap 2023 SRA uploads to NCBI
date: '2025-04-17'
categories: Protocol
tags: [Protocol, Bioinformatics]
projects: Mcap DT 2023
---

## Mcap Developmental Timeseries 2023 SRA uploads to NCBI

This post details the NCBI Sequence Read Archive upload for the sequence data from the Mcap developmental timeseries 2023 experiment. The github for that project is [here](https://github.com/JillAshey/DevelopmentalTimeseries). This protocol is based on A. Huffmyer's [SRA upload post](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/SRA-Uploads-10-November-2022/) and Putnam Lab [SRA upload protocol](https://github.com/Putnam-Lab/Lab_Management/blob/master/Bioinformatics_%26_Coding/Data_Mangament/SRA-Upload_Protocol.md). This post includes information on both the RNAseq and smRNAseq SRA uploads. 

### Overview for rRNA depleted RNAseq uploads 

Below is the following information for this submission. 

Sequence Read Archive (SRA) submission: SUB15265814

#### Submitter

- Jill Ashey; 120 Flagg Road, Kingston RI 02881

#### General info 

- Not associated with current BioProject 
- Not associated with current BioSample
- Release date: immediately 

#### General info 

- Not associated with current BioProject 
- Not associated with current BioSample
- Release date: immediately 

#### Project info

- Project title: Montipora capitata developmental timeseries 2023 
- Public description: RNA-seq data from rRNA depleted libraries from a 2023 experiment with Montipora capitata collecting across 8 developmental stages 
- Not associated with umbrella project

#### Biosample type: invertebrate 

[Link to submission spreadsheet](https://github.com/JillAshey/DevelopmentalTimeseries/blob/main/data/NCBI/rRNA_depleted_RNAseq_Invertebrate.1.0.xlsx)

Important note so that errors do not occur during upload: After filling out values for attributes provided in the template, your samples are not distinguishable by at least one, or a combination of attributes. **Make the combined value of all attributes unique for each sample while taking into account that "sample name", "sample title" and "description" are not included in this check for "uniqueness".**

#### SRA metadata

[Link to SRA metadata](https://github.com/JillAshey/DevelopmentalTimeseries/blob/main/data/NCBI/rRNA_depleted_RNAseq_SRA_data.tsv)

#### Files 

Upload the fastq.gz files to NCBI with the following code. First, make a new folder that only contains sequence data.

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/mRNA
mkdir raw_file_rnaseq_sra
cd raw_file_rnaseq_sra

ln -s /data/putnamlab/jillashey/DT_Mcap_2023/mRNA/data/raw/*fastq.gz .
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
mkdir mcap_2023_rRNA_depletion_rnaseq
cd mcap_2023_rRNA_depletion_rnaseq
mput *
```

This takes a few minutes (depending on how many files you have). Once transfers are complete, click select preload folder in the submission portal. Wait until all files have been uploaded and select the `mcap_2023_rRNA_depletion_rnaseq` folder. Click submit!

Bioproject: PRJNA1252744

### Overview for polyA RNAseq uploads 

Below is the following information for this submission. 

Sequence Read Archive (SRA) submission: SUB15267960

#### Submitter

- Jill Ashey; 120 Flagg Road, Kingston RI 02881

#### General info 

- Not associated with current BioProject 
- Not associated with current BioSample
- Release date: immediately 

#### General info 

- Not associated with current BioProject 
- Not associated with current BioSample
- Release date: immediately 

#### Project info

- Project title: Montipora capitata developmental timeseries 2023 - polyA libraries
- Public description: RNA-seq data from polyA libraries from a 2023 experiment with Montipora capitata collecting across 8 developmental stages 
- Not associated with umbrella project

#### Biosample type: invertebrate 

[Link to submission spreadsheet](https://github.com/JillAshey/DevelopmentalTimeseries/blob/main/data/NCBI/polyA_RNAseq_Invertebrate.1.0.xlsx)

Important note so that errors do not occur during upload: After filling out values for attributes provided in the template, your samples are not distinguishable by at least one, or a combination of attributes. **Make the combined value of all attributes unique for each sample while taking into account that "sample name", "sample title" and "description" are not included in this check for "uniqueness".**

#### SRA metadata

[Link to SRA metadata](https://github.com/JillAshey/DevelopmentalTimeseries/blob/main/data/NCBI/polyA_RNAseq_SRA_data.tsv)

#### Files 

Upload the fastq.gz files to NCBI with the following code. First, make a new folder that only contains sequence data.

```
cd /project/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA
mkdir raw_file_rnaseq_sra
cd raw_file_rnaseq_sra

ln -s /project/pi_hputnam_uri_edu/jillashey/Mcap_2023/polyA/*fastq.gz .
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
mkdir mcap_2023_rnaseq
cd mcap_2023_rnaseq
mput *
```

This takes a few minutes (depending on how many files you have). Once transfers are complete, click select preload folder in the submission portal. Wait until all files have been uploaded and select the `mcap_2023_zymo_smrnaseq ` folder. Click submit!

Bioproject: XXXXXXXX











### Overview for smRNAseq uploads - Zymo libraries 

**NOTE: These sequences should NOT be used!!! Library prep and sequencing was not successful; based on QC, it appears that adapter content was primarily sequenced. I am uploading these sequences to NCBI as a backup but they will likely not be released on NCBI. They will also not be used in analyses; we are resequencing the smRNA samples from this expeirment.**

Below is the following information for this submission. 

Sequence Read Archive (SRA) submission: SUB15267867

#### Submitter

- Jill Ashey; 120 Flagg Road, Kingston RI 02881

#### General info 

- Not associated with current BioProject 
- Not associated with current BioSample
- Release date: immediately 

#### General info 

- Not associated with current BioProject 
- Not associated with current BioSample
- Release date: 2026-05-01 

#### Project info

- Project title: Montipora capitata developmental timeseries 2023 - smRNAseq using Zymo libraries 
- Public description: smRNA-seq data from Zymo libraries from a 2023 experiment with Montipora capitata collecting across 8 developmental stages. DO NOT USE 
- Not associated with umbrella project

#### Biosample type: invertebrate 

[Link to submission spreadsheet](https://github.com/JillAshey/DevelopmentalTimeseries/blob/main/data/NCBI/miRNA_Zymo_smRNAseq_Invertebrate.1.0.xlsx)

Important note so that errors do not occur during upload: After filling out values for attributes provided in the template, your samples are not distinguishable by at least one, or a combination of attributes. **Make the combined value of all attributes unique for each sample while taking into account that "sample name", "sample title" and "description" are not included in this check for "uniqueness".**

#### SRA metadata

[Link to SRA metadata](https://github.com/JillAshey/DevelopmentalTimeseries/blob/main/data/NCBI/miRNA_Zymo_smRNAseq_SRA_data.tsv)

#### Files 

Upload the fastq.gz files to NCBI with the following code. First, make a new folder that only contains sequence data.

```
cd /data/putnamlab/jillashey/DT_Mcap_2023/smRNA
mkdir raw_file_smrnaseq_sra
cd raw_file_smrnaseq_sra

ln -s /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/*fastq.gz .
ln -s /data/putnamlab/jillashey/DT_Mcap_2023/smRNA/data/raw/first_batch/*fastq.gz .
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
mkdir mcap_2023_zymo_smrnaseq
cd mcap_2023_zymo_smrnaseq
mput *
```

This takes a few minutes (depending on how many files you have). Once transfers are complete, click select preload folder in the submission portal. Wait until all files have been uploaded and select the `mcap_2023_zymo_smrnaseq ` folder. Click submit!

Bioproject: PRJNA1252766


