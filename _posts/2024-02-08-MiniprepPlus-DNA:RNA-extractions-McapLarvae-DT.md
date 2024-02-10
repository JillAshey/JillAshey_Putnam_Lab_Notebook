---
layout: post
title: MiniPrep Plus DNA/RNA extractions
date: '2024-02-08'
categories:
tags: [DNA, RNA, Extractions, Protocols]
projects: Mcap developmental time series 
---

# Extractions for Mcap developmental time series Hawaii 2023

This protocol is based on Maggie’s coral egg/larvae extraction [protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/Larvae-Ex-Protocol/) and Kevin’s coral larvae extraction [protocol](https://kevinhwong1.github.io/KevinHWong_Notebook/DNA-RNA-Extractions-on-P.-astreoides-larvae-BEAD-BEATING/). This post details the info about the extraction steps for the *M. capitata* developmental time series experiment in Hawaii 2023. The github for that project is linked [here](https://github.com/JillAshey/Hawaii_Developmental_TimeSeries).  

### Samples 

The samples run today were M4, M12, M15, M25, M38, M40, M47, M48, M61, M64, M76, and M83. These samples had between 300-1600 eggs/larvae per tube. Sample were preserved in 500 ul - 1 ml of DNA/RNA Shield. This photo is post-bead beating: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/DT_mcap2023/samples_20240208.JPG)

After bead beating, some of the tubes (M4, M47, M64, and M76) still looked like they had intact eggs/larvae so these samples were bead beat for another minute. Still some intact-looking larvae but the liquid looked pretty pigmented so I moved forward with the extraction as is. 

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

- I added between 200-500 uL of DNA/RNA shield to the samples that had <1mL of shield in them after thawing them so that the volume in the tubes was 1 mL of shield. I aliquoted out 500 uL for extractions and saved the remaining fraction in case we need to re-extract. 
- In the RNA extraction, I added another wash step with 400 uL of wash buffer after the DNase incubation (i.e., after the incubation, I did 400 uL of prep buffer, 700 uL of wash buffer, 400 uL of wash buffer, and another 400 uL of wash buffer).
- I eluted in 80 uL of either Tris (DNA) or DNA/RNA free water (RNA) instead of 100 uL. 
- I accidently added >1mL of lysis buffer to M40 during the lysis step (was only supposed to add 575uL). I am going to proceed with this sample but I will keep an eye on it for QC

### QC 

#### Qubit results 

| Tube ID | Qubit DNA1 | Qubit DNA2 | Qubit DNA Average | Qubit RNA1 | Qubit RNA2 | Qubit RNA Average |
| ------- | ---------- | ---------- | ----------------- | ---------- | ---------- | ----------------- |
| M12     | NA         | NA         | NA                | 127        | 128        | 127.5             |
| M15     | NA         | NA         | NA                | 151        | 151        | 151               |
| M25     | NA         | NA         | NA                | 82.6       | 82.8       | 82.7              |
| M38     | 3.04       | 3.12       | 3.08              | 194        | 195        | 194.5             |
| M4      | NA         | NA         | NA                | 42.4       | 42.2       | 42.3              |
| M40     | 5.3        | 5.32       | 5.31              | 100        | 101        | 100.5             |
| M47     | 8.52       | 8.44       | 8.48              | 100        | 99.4       | 99.7              |
| M48     | 7.82       | 7.82       | 7.82              | 106        | 106        | 106               |
| M61     | 13.5       | 13.4       | 13.45             | 91.8       | 91.6       | 91.7              |
| M64     | 4.38       | 4.4        | 4.39              | 71.8       | 71.8       | 71.8              |
| M76     | 4.68       | 4.64       | 4.66              | 71         | 71.2       | 71.1              |
| M83     | 17.1       | 17.3       | 17.2              | 48.4       | 48.4       | 48.4              |

#### Gel 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/DT_mcap2023/gel_20240208.JPG)

Gel is really all over the place. Some of the samples look good, some not so good. There also seems to be some residual DNA at the top of some samples.

I discussed with Zoe and she recommended to not add beads when vortexting the sample when it thaws. Additionally, she recommended that I use less input volume, as I may be overwhelming the column with RNA, which prevents the DNase from digesting all of the RNA. 
