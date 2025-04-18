---
layout: post
title: MiniPrep Plus DNA/RNA extractions
date: '2024-02-10'
categories:
tags: [DNA, RNA, Extractions, Protocols]
projects: Mcap developmental time series 
---

# Extractions for Mcap developmental time series Hawaii 2023

This protocol is based on Maggie’s coral egg/larvae extraction [protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/Larvae-Ex-Protocol/) and Kevin’s coral larvae extraction [protocol](https://kevinhwong1.github.io/KevinHWong_Notebook/DNA-RNA-Extractions-on-P.-astreoides-larvae-BEAD-BEATING/). This post details the info about the extraction steps for the *M. capitata* developmental time series experiment in Hawaii 2023. The github for that project is linked [here](https://github.com/JillAshey/Hawaii_Developmental_TimeSeries). 

**Important changes in this protocol**: I did not add beads to samples for beating step and I used 300uL as my input volume instead of 500uL. See the bottom of this [post](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-02-08-MiniprepPlus-DNA%3ARNA-extractions-McapLarvae-DT.md) for reasoning behind those choices.  

### Samples 

The samples run today were M13, M23, M35, M52, M60, M72, and M85. These samples had between 300-1600 eggs/larvae per tube. Sample were preserved in 500 ul of DNA/RNA Shield. This photo is post-beating: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/DT_mcap2023/samples_20240210.JPG)

I was originally going to extract M5 as well, but I spilled all of it on the ground post-beating. Good thing I collected n=6 samples per timepoint!

### Materials 

- Zymo Quick-DNA/RNA Miniprep Plus Kit [HERE](https://files.zymoresearch.com/protocols/_d7003t_d7003_quick-dna-rna_miniprep_plus_kit.pdf) Protocol Booklet
- Tris-Ethylenediaminetetraacetic acid (EDTA) 1X buffer for DNA elution
- Heating block capable of heating to 70ºC
- Centrifuge and rotor capable of spinning at 15,000 rcf
- Plastics 
	- 5 1.5 mL microcentrifuge tubes per sample
	- 2 PCR tubes per sample
	- 2 Qubit tubes per sample 

### Protocol

Protocol was followed according to this [post](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2023-07-21-MiniprepPlus-DNA%3ARNA-extractions-McapLarvae.md). A few differences: 

- I added between 200-500 uL of DNA/RNA shield to the samples that had <1mL of shield in them after thawing them so that the volume in the tubes was 1 mL of shield. I aliquoted out 300 uL for extractions and saved the remaining fraction in case we need to re-extract. 
- I did not add beads when mixing the samples after they thawed. 
- In the RNA extraction, I added another wash step with 400 uL of wash buffer after the DNase incubation (i.e., after the incubation, I did 400 uL of prep buffer, 700 uL of wash buffer, 400 uL of wash buffer, and another 400 uL of wash buffer).
- I eluted in 80 uL of either Tris (DNA) or DNA/RNA free water (RNA) instead of 100 uL. 

### QC 

#### Qubit results 

| Tube ID | Qubit DNA1 | Qubit DNA2 | Qubit DNA Average | Qubit RNA1 | Qubit RNA2 | Qubit RNA Average |
| ------- | ---------- | ---------- | ----------------- | ---------- | ---------- | ----------------- |
| M13     | NA         | NA         | NA                | 58.6       | 58.4       | 58.5              |
| M23     | NA         | NA         | NA                | 16.4       | 14.4       | 15.4              |
| M35     | NA         | NA         | NA                | 29         | 29         | 29                |
| M52     | NA         | NA         | NA                | 19.2       | 18.8       | 19                |
| M60     | 3.54       | 3.58       | 3.56              | 36.8       | 36.2       | 36.5              |
| M72     | 4.58       | 4.66       | 4.62              | 24.2       | 24         | 24.1              |
| M85     | 5.74       | 5.84       | 5.79              | 30.4       | 30.2       | 30.3              |

#### Gel 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/DT_mcap2023/gel_20240210.JPG)

This gel looks so much better!!! The 4hpf sample still doesn't look great, but I knew these samples might not be viable because I added water to shield when doing the sampling. Much happier with these extractions, so this is the method that I will be using moving forward with these samples for the developmental timeseries project. 