---
layout: post
title: Developmental 2023 Timeseries mRNA analysis - egg and sperm analysis 
date: '2025-04-15'
categories: Analysis
tags: [Bioinformatics, mRNA, Montipora capitata]
projects: Developmental Timseries 2023
---

## Developmental 2023 Timeseries mRNA analysis - egg and sperm

I have sequenced and analyzed samples from the following time points: 1, 4, 9, 14, 22, 28, 48 and 72 hpf. The github for this project is [here](https://github.com/JillAshey/Hawaii_Developmental_TimeSeries). During the closed portion of my defense, we discussed incorporating unfertilized egg and sperm samples into the analyses to understand the contribution of the sperm v. egg to the mRNA complement. In my 2023 experiment, I did not collect sperm or unfertilized egg samples so I am going to use the sperm and egg samples that were used in [Van et Etten et al. 2020](https://peerj.com/articles/9739/#supp-1), which analyzed sperm and unfertilized egg samples in Mcap. In the paper, they state the following on accessing the sample fastq files: "The egg data are publicly available under NCBI BioProject [PRJNA616341](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA616341) (SAMN14486762, SAMN14486763, SAMN14486764) and the sperm data are publicly available under NCBI BioProject [PRJNA339779](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA339779)." They also state 

Here are the samples that I need to download from NCBI: 

- [SAMN05607941](https://www.ncbi.nlm.nih.gov/sra/SRX2039373[accn]) - sperm
- [SAMN14486762](https://www.ncbi.nlm.nih.gov/sra/?term=SAMN14486762) - egg
- [SAMN14486763](https://www.ncbi.nlm.nih.gov/sra/?term=SAMN14486763) - egg 
- [SAMN14486764](https://www.ncbi.nlm.nih.gov/sra/?term=SAMN14486764) - egg

The egg and sperm libraries were generated using [Illumina  TruSeq RNA Library Preparation Kit v2](https://www.illumina.com/products/by-type/sequencing-kits/library-prep-kits/truseq-rna-v2.html#tabs-fbcf765f10-item-fb7f2f9963-documentation), which uses polyA selection. They were sequenced on an Illumina MiSeq (150 cycles, PE). 