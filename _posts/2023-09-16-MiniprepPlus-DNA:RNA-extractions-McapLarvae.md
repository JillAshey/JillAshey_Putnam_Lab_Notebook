---
layout: post
title: MiniPrep Plus DNA/RNA extractions
date: '2023-09-16'
categories:
tags: [DNA, RNA, Extractions, Protocols]
projects: Larvae C/D 
---

# Extractions for Mcap larvae C/D Hawaii 2023 experiment 

This protocol is based on Maggie’s coral egg/larvae extraction [protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/Larvae-Ex-Protocol/) and Kevin’s coral larvae extraction [protocol](https://kevinhwong1.github.io/KevinHWong_Notebook/DNA-RNA-Extractions-on-P.-astreoides-larvae-BEAD-BEATING/). This post details the info about the extraction steps for the test samples that were collected in HI 2023 for A. Huffmyer's experiment on symbiont identity on *Montipora capitata* coral larval thermal performance. The github for that project is linked [here](https://github.com/AHuffmyer/larval_symbiont_TPC). 

### Samples 

The samples run today were R13, R21, R32, R42, R48, R54, R60, R64, R68, R74, R88, R96, R98 and R106. These samples had 50 larvae in them and were preserved with 700 uL of DNA/RNA Shield. This photo is post-bead beating: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/larvae_cd_mcap2023/samples_20230916.JPG)

### Materials 

- Zymo Quick-DNA/RNA Miniprep Plus Kit [HERE](https://files.zymoresearch.com/protocols/_d7003t_d7003_quick-dna-rna_miniprep_plus_kit.pdf) Protocol Booklet
- Tris-Ethylenediaminetetraacetic acid (EDTA) 1X buffer for DNA elution
- Heating block capable of heating to 70ºC
- Centrifuge and rotor capable of spinning at 15,000 rcf
- Plastics 
	- 3 1.5 mL microcentrifuge tubes per sample
	- 2 PCR tubes per sample
	- 2 Qubit tubes per sample 
	- 1 5 mL tube per sample 

### Protocol

Protocol was followed according to this [post](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2023-07-21-MiniprepPlus-DNA%3ARNA-extractions-McapLarvae.md). A few differences: 

- I added 300 uL of DNA/RNA shield to the samples after thawing them so that the volume in the tubes was 1 mL of shield. I aliquoted out 500 uL for extractions and saved the remaining fraction in case we need to re-extract. 
- In the RNA extraction, I added another wash step with 400 uL of wash buffer after the DNase incubation (i.e., after the incubation, I did 400 uL of prep buffer, 700 uL of wash buffer, 400 uL of wash buffer, and another 400 uL of wash buffer).
- I eluted in 80 uL of either Tris (DNA) or DNA/RNA free water (RNA) instead of 100 uL. 

### QC 

#### Qubit & Nanodrop results 

| Sample ID | Qubit DNA1 | Qubit DNA2 | Qubit DNA average | Nanodrop DNA concentration | Nanodrop DNA 260/280 | Nanodrop DNA 260/230 | Qubit RNA1 | Qubit RNA2 | Qubit RNA average | Nanodrop RNA concentration | Nanodrop RNA 260/280 | Nanodrop RNA 260/230 |
| --------- | ---------- | ---------- | ----------------- | -------------------------- | -------------------- | -------------------- | ---------- | ---------- | ----------------- | -------------------------- | -------------------- | -------------------- |
| R13       | 8.88       | 8.98       | 8.93              | 12.3                       | 1.92                 | 66.42                | 10.6       | 10.6       | 10.6              | 11.4                       | 2.06                 | 1.03                 |
| R21       | 7.58       | 7.64       | 7.61              | 12.1                       | 1.86                 | 11.75                | 11.2       | 11         | 11.1              | 11.1                       | 1.95                 | 0.53                 |
| R32       | 8.18       | 8.26       | 8.22              | 11.6                       | 1.92                 | 1.09                 | 11.8       | 11.8       | 11.8              | 12.7                       | 1.95                 | 0.86                 |
| R42       | 11         | 11.1       | 11.05             | 17.4                       | 1.93                 | 5.32                 | 10.6       | 10.6       | 10.6              | 9.4                        | 2.1                  | 1.08                 |
| R48       | 8.58       | 8.64       | 8.61              | 14.6                       | 1.85                 | 1.79                 | NA         | NA         | NA                | 7.2                        | 2.23                 | 1.13                 |
| R54       | 6.9        | 6.94       | 6.92              | 11.4                       | 1.96                 | 0.35                 | NA         | NA         | NA                | 7.8                        | 2.16                 | 1.06                 |
| R60       | 10.1       | 10.1       | 10.1              | 16.2                       | 1.95                 | 4.46                 | NA         | NA         | NA                | 8.9                        | 2.11                 | 1.1                  |
| R64       | 11         | 11.1       | 11.05             | 14.1                       | 1.98                 | 6.77                 | 10.4       | 10.4       | 10.4              | 11.5                       | 2.01                 | 1.07                 |
| R68       | 13.9       | 14.1       | 14                | 17.3                       | 1.96                 | 15.38                | 11         | 11         | 11                | 12.6                       | 1.93                 | 1                    |
| R74       | 9.36       | 9.38       | 9.37              | 10.2                       | 2.12                 | 1.46                 | 11.2       | 11.2       | 11.2              | 11.4                       | 2.04                 | 1.1                  |
| R88       | 9.2        | 9.28       | 9.24              | 13                         | 1.99                 | 8.79                 | NA         | NA         | NA                | 9.2                        | 2.07                 | 1.06                 |
| R96       | 8.18       | 8.3        | 8.24              | 11.1                       | 1.97                 | 24.67                | NA         | NA         | NA                | 8.5                        | 2.1                  | 1.25                 |
| R98       | 8.8        | 8.86       | 8.83              | 12.3                       | 1.96                 | 0.07                 | NA         | NA         | NA                | 14.8                       | 1.9                  | 0.77                 |
| R106      | 12.6       | 12.7       | 12.65             | 17.3                       | 1.94                 | 7.19                 | NA         | NA         | NA                | 10.7                       | 2.16                 | 1.18                 |

Some of these values are really weird... Nanodrop DNA 260/230 have weirdly high values for R13, R21, R42, R60, R64, R68, R88, R96 and R106 and weirdly low value for R98. Some of the Nanodrop RNA 260/230 values are lower than I would like (R21, R32, and R98). I got DNA concentrations on the Qubit for all samples, but I only got RNA concentrations on the Qubit for R13, R21, R32, R42, R64, R68 and R74 (7/14 of the samples). However, I got RNA concentrations when I ran the Nanodrop, though that is less reliable.

#### Gel 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/larvae_cd_mcap2023/gel_20230916.JPG)

Even though I didn't get RNA concentrations for all of the samples, I got distinct RNA bands for all samples. Need to discuss this further w/ Ariana. I'm thinking that I should run a tapestation on the samples that failed from this batch and the batch before (9/11/23). 