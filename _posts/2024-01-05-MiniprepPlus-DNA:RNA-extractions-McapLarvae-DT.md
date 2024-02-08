---
layout: post
title: MiniPrep Plus DNA/RNA extractions
date: '2023-10-27'
categories:
tags: [DNA, RNA, Extractions, Protocols]
projects: Mcap developmental time series 
---

# Extractions for Mcap developmental time series Hawaii 2023

This protocol is based on Maggie’s coral egg/larvae extraction [protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/Larvae-Ex-Protocol/) and Kevin’s coral larvae extraction [protocol](https://kevinhwong1.github.io/KevinHWong_Notebook/DNA-RNA-Extractions-on-P.-astreoides-larvae-BEAD-BEATING/). This post details the info about the extraction steps for the *M. capitata* developmental time series experiment in Hawaii 2023. The github for that project is linked [here](https://github.com/JillAshey/Hawaii_Developmental_TimeSeries).  

### Samples 

The samples run today were M10, M11, M14, M24, M26, M28, M49, M50, M51, M59, M62 and M63. These samples had between 700-1400 larvae per tube. M10, M11, and M14 were preserved in 1 mL of DNA/RNA Shield; M24, M26, M28, M49, M50, M51, M59, M62 and M63 were preserved in 700 uL of DNA/RNA Shield. This photo is post-bead beating: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/DT_mcap2023/samples_20240105.JPG)

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

- I added 300 uL of DNA/RNA shield to the samples that had <1mL of shield in them after thawing them so that the volume in the tubes was 1 mL of shield. I aliquoted out 500 uL for extractions and saved the remaining fraction in case we need to re-extract. 
- In the RNA extraction, I added another wash step with 400 uL of wash buffer after the DNase incubation (i.e., after the incubation, I did 400 uL of prep buffer, 700 uL of wash buffer, 400 uL of wash buffer, and another 400 uL of wash buffer).
- I eluted in 80 uL of either Tris (DNA) or DNA/RNA free water (RNA) instead of 100 uL. 

### QC 

#### Qubit results 

| Tube ID | Qubit DNA1 | Qubit DNA2 | Qubit DNA Average | Qubit RNA1 | Qubit RNA2 | Qubit RNA Average |
| ------- | ---------- | ---------- | ----------------- | ---------- | ---------- | ----------------- |
| M10     | NA         | NA         | NA                | 101        | 101        | 101               |
| M11     | NA         | NA         | NA                | 123        | 124        | 123.5             |
| M14     | NA         | NA         | NA                | 252        | 252        | 252               |
| M24     | 2.4        | 2.44       | 2.42              | 162        | 162        | 162               |
| M26     | NA         | NA         | NA                | 230        | 230        | 230               |
| M28     | NA         | NA         | NA                | 42.2       | 42.8       | 42.5              |
| M49     | 8.22       | 8.44       | 8.33              | 86.8       | 86.4       | 86.6              |
| M50     | 9.88       | 10.1       | 9.99              | 72.4       | 71.8       | 72.1              |
| M51     | 8.96       | 9.2        | 9.08              | 55.4       | 55.8       | 55.6              |
| M59     | 13.2       | 13.5       | 13.35             | 108        | 109        | 108.5             |
| M62     | 12.3       | 12.6       | 12.45             | 98.6       | 98.6       | 98.6              |
| M63     | 8.3        | 8.5        | 8.4               | 89         | 89.4       | 89.2              |

#### Gel 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/DT_mcap2023/gel_20240105.JPG)

Not much DNA in the samples, but a lot of RNA. The RNA bands look very blurry and some samples don't have distinct bands. I may nanodrop and/or tapestation these samples to check out the quality of RNA. 