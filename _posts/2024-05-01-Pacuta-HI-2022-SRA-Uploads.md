---
layout: post
title: Pacuta SRA uploads to NCBI
date: '2024-05-01'
categories: Protocol
tags: [Protocol, Bioinformatics]
projects: Pacuta HI 2022
---

## SRA uploads to NCBI

This post details the NCBI Sequence Read Archive upload for the sequence data from the Pacuta HI 2022 experiment. The github for that project is [here](https://github.com/fscucchia/Hawaii2022_pH_Temp_Mcap_Pacu/tree/main). This protocol is based on A. Huffmyer's [SRA upload post](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/SRA-Uploads-10-November-2022/) and Putnam Lab [SRA upload protocol](https://github.com/Putnam-Lab/Lab_Management/blob/master/Bioinformatics_%26_Coding/Data_Mangament/SRA-Upload_Protocol.md). 

### Overview 

Below is the following information for this submission. 

Sequence Read Archive (SRA) submission: SUB14419354

#### Submitter

- Jill Ashey; 120 Flagg Road, Kingston RI 02881

#### General info 

- Not associated with current BioProject 
- Not associated with current BioSample
- Release date: immediately 

#### Project info 
- project title: Pocillopora acuta juvenile RNA-seq
- public description: RNA-seq data from Pocillopora acuta juveniles to future warming and acidification conditions from 2022 Hawaii experiments
- Not associated with umbrella project 

#### Biosample type: invertebrate 

[Link to submission spreadsheet](https://github.com/fscucchia/Hawaii2022_pH_Temp_Mcap_Pacu/blob/main/Raw_data/NCBI/Pacuta_HI_2022_Invertebrate.1.0.xlsx)

Important note so that errors do not occur during upload: After filling out values for attributes provided in the template, your samples are not distinguishable by at least one, or a combination of attributes. **Make the combined value of all attributes unique for each sample while taking into account that "sample name", "sample title" and "description" are not included in this check for "uniqueness".**
 
#### SRA metadata 

[Link to SRA metadata spreadsheet](https://github.com/fscucchia/Hawaii2022_pH_Temp_Mcap_Pacu/blob/main/Raw_data/NCBI/Pacuta_HI_2022_SRA_metadata.tsv)

#### Files  

Upload the fastq.gz files to NCBI with the following code. First, make a new folder that only contains sequence data.
 
```
cd /data/putnamlab/jillashey/Pacuta_HI_2022
mkdir raw_file_rnaseq_sra
cd raw_file_rnaseq_sra

ln -s /data/putnamlab/KITT/hputnam/20231127_Scucchia_HI/*fastq.gz .
``` 

Use FTP command line file upload to provide files. Activate FTP on the command line. Command prompt now changes to `ftp>`. Immediately login using the username and password from the NCBI submission portal. 

```
ftp -i
open ftp-private.ncbi.nlm.nih.gov
USERNAME
PASSWORD 
```

Go into the folder that NCBI provided from the NCBI submission portal. Make a new directory for the files and go into that folder. Put all the sequences from Andromeda into this new folder (they will still be on Andromeda after transfer is complete). 

```
cd uploads/jillashey_uri.edu_XXXXXXXXXXX
mkdir pacuta_2022
cd pacuta_2022
mput *
```

This takes a few minutes. Once transfers are complete, click select preload folder in the submission portal. Wait until all files have been uploaded and select the `pacuta_2022` folder. Click submit!