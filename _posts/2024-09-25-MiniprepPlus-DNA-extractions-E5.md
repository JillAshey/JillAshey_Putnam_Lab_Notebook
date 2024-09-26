---
layout: post
title: MiniPrep Plus DNA/RNA extractions
date: '2024-09-25'
categories:
tags: [DNA, Extractions, Protocols]
projects: e5
---

# Extractions for e5 time series 

This protocol is based on the Putnam lab Zymo Miniprep protocol (see example from Zoe [here](https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/Protocols_Zymo_Quick_DNA_RNA_Miniprep_Plus/)). This post details the info about the extraction steps for samples from the e5 timeseries (Acropora, Pocillopora, Porites spp.) in Moorea 2020. The e5 github is linked [here](https://github.com/urol-e5). 

### Samples 

The samples run today were re-extractions from the e5 deep dive and timeseries samples. I am having library prep issues with 6 of the deep dive samples and ~29 of the timeseries DNA samples do not have enough ng of DNA for Azenta. I am going to re-extract to see if we can get higher concentrations. Since the samples had already been through one extraction, there was <1mL of DNA/RNA shield in the tube. 

| Sample               | tube # | Purpose    | Notes                            |
| -------------------- | ------ | ---------- | -------------------------------- |
| 7-20220208           | 7      | timeseries |                                  |
| 13                   | 13     | timeseries |                                  |
| 23                   | 23     | timeseries |                                  |
| 29                   | 29     | timeseries |                                  |
| 195                  | 195    | timeseries |                                  |
| 20240805_POC-52_TP1  | 213    | timeseries | Might be okay with original DNA? |
| 235                  | 235    | timeseries |                                  |
| 247-20220208         | 247    | timeseries |                                  |
| 251-20220208         | 251    | timeseries |                                  |
| 297                  | 297    | timeseries |                                  |
| 319                  | 319    | timeseries |                                  |
| 327                  | 327    | timeseries |                                  |
| 385                  | 385    | deepdive   |                                  |
| 393                  | 393    | deepdive   |                                  |
| 401                  | 401    | deepdive   |                                  |
| 421                  | 421    | deepdive   |                                  |
| 469-20211122         | 469    | timeseries |                                  |
| 477                  | 477    | timeseries |                                  |
| 487                  | 487    | deepdive   |                                  |
| 489                  | 489    | deepdive   |                                  |
| 491                  | 491    | timeseries |                                  |
| 493                  | 493    | timeseries |                                  |
| 515-20211104         | 515    | timeseries |                                  |
| 20240805_POC-219_TP3 | 521    | timeseries | Might be okay with original DNA? |
| 529                  | 529    | timeseries |                                  |
| 557                  | 557    | timeseries |                                  |
| 633-20220201         | 633    | timeseries |                                  |
| 637                  | 637    | timeseries |                                  |
| 641                  | 641    | timeseries |                                  |
| 669                  | 669    | timeseries |                                  |
| 711                  | 711    | timeseries |                                  |
| 713                  | 713    | timeseries |                                  |
| 741                  | 741    | timeseries |                                  |
| 749                  | 749    | timeseries |                                  |
| 877                  | 877    | timeseries |

There were 35 samples total, so I split them into two batches because the centrifuge can only fit 24 tubes at a time. 213, 385, 393, 401, and 521 were POC, the rest of the samples were POR. 

This photo is before I removed shield from the sample tubes: 

Set 1

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/refs/heads/master/images/samples_set1_20240926.JPG)

Set 2 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/refs/heads/master/images/samples_set2_20240926.JPG)

### Materials 

ONLY DID DNA EXTRACTIONS

- Zymo Quick-DNA/RNA Miniprep Plus Kit [HERE](https://files.zymoresearch.com/protocols/_d7003t_d7003_quick-dna-rna_miniprep_plus_kit.pdf) Protocol Booklet
- Tris-Ethylenediaminetetraacetic acid (EDTA) 1X buffer for DNA elution
- Heating block capable of heating to 70ºC
- Centrifuge and rotor capable of spinning at 15,000 rcf
- Plastics 
	- 3 1.5 mL microcentrifuge tubes per sample
	- 1 PCR tubes per sample
	- 1 Qubit tubes per sample 

### Protocol 

This [protocol](https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/Protocols_Zymo_Quick_DNA_RNA_Miniprep_Plus/) was followed. A few differences 

- For POC samples, 300 uL of sample was moved from the original shield tube into a new 1.5mL tube for the prot-k digestion step. 
- For POR samples, 100 uL of sample and 200 uL of fresh DNA/RNA shield was moved into a new 1.5mL tube for the prot-k digestion step. 
- Samples were eluted in 90 uL of Tris instead of 100 uL. 

### QC 

#### Qubit results 

For Qubit, I used the DNA HS kit on most samples; on 385, 393, and 401, I used the DNA BR kit because these were out of range with the HS kit. I used 5 uL of DNA as input for the HS kit and 2 uL of DNA as input for the BR kit. 

| tube # | Colony ID | Time point | Species                  | Purpose    | DNA1_ng_uL | DNA2_ng_uL | DNA_Avg | Volume_uL | DNA_ng  |
| ------ | --------- | ---------- | ------------------------ | ---------- | ---------- | ---------- | ------- | --------- | ------- |
| 7      | POR-262   | TP1        | Porites evermanni        | timeseries | 0.39       | 0.38       | 0.385   | 76        | 29.26   |
| 13     | POR-245   | TP1        | Porites evermanni        | timeseries | 0.386      | 0.381      | 0.3835  | 76        | 29.146  |
| 23     | POR-236   | TP1        | Porites evermanni        | timeseries | 0.432      | 0.436      | 0.434   | 76        | 32.984  |
| 29     | POR-216   | TP1        | Porites evermanni        | timeseries | 0.9        | 0.892      | 0.896   | 76        | 68.096  |
| 195    | POR-69    | TP1        | Porites evermanni        | timeseries | 0.317      | 0.319      | 0.318   | 76        | 24.168  |
| 213    | POC-52    | TP1        | Pocillopora tuahiniensis | timeseries | 5.713      | 5.68       | 5.6965  | 76        | 432.934 |
| 235    | POR-72    | TP1        | Porites evermanni        | timeseries | 0.584      | 0.23       | 0.407   | 76        | 30.932  |
| 247    | POR-74    | TP1        | Porites evermanni        | timeseries | 0.183      | 0.19       | 0.1865  | 76        | 14.174  |
| 251    | POR-73    | TP1        | Porites evermanni        | timeseries | 0.193      | 0.193      | 0.193   | 76        | 14.668  |
| 297    | POR-245   | TP2        | Porites evermanni        | timeseries | 1.28       | 1.29       | 1.285   | 76        | 97.66   |
| 319    | POR-236   | TP2        | Porites evermanni        | timeseries | 2          | 2          | 2       | 76        | 152     |
| 327    | POR-260   | TP2        | Porites evermanni        | timeseries | 0.566      | 0.566      | 0.566   | 76        | 43.016  |
| 385    | POC-47    | TP2        | Pocillopora tuahiniensis | deepdive   | 26.4       | 25.6       | 26      | 74        | 1924    |
| 393    | POC-50    | TP2        | Pocillopora tuahiniensis | deepdive   | 21.9       | 21.2       | 21.55   | 74        | 1594.7  |
| 401    | POC-53    | TP2        | Pocillopora tuahiniensis | deepdive   | 26         | 25.4       | 25.7    | 74        | 1901.8  |
| 421    | POR-82    | TP2        | Porites evermanni        | deepdive   | 0.924      | 0.908      | 0.916   | 76        | 69.616  |
| 469    | POR-83    | TP2        | Porites evermanni        | timeseries | 0.52       | 0.52       | 0.52    | 76        | 39.52   |
| 477    | POR-69    | TP2        | Porites evermanni        | timeseries | 0.44       | 0.448      | 0.444   | 76        | 33.744  |
| 487    | POR-71    | TP2        | Porites evermanni        | deepdive   | 1.71       | 1.7        | 1.705   | 76        | 129.58  |
| 489    | POR-79    | TP2        | Porites evermanni        | deepdive   | 1.4        | 1.39       | 1.395   | 76        | 106.02  |
| 491    | POR-73    | TP2        | Porites evermanni        | timeseries | 0.604      | 0.588      | 0.596   | 76        | 45.296  |
| 493    | POR-72    | TP2        | Porites evermanni        | timeseries | 0.712      | 0.712      | 0.712   | 76        | 54.112  |
| 515    | POR-216   | TP3        | Porites evermanni        | timeseries | 0.48       | 0.48       | 0.48    | 76        | 36.48   |
| 521    | POC-219   | TP3        | Pocillopora tuahiniensis | timeseries | 2.71       | 2.7        | 2.705   | 76        | 205.58  |
| 529    | POR-262   | TP3        | Porites evermanni        | timeseries | 0.54       | 0.521      | 0.5305  | 76        | 40.318  |
| 557    | POR-260   | TP3        | Porites evermanni        | timeseries | 0.464      | 0.444      | 0.454   | 76        | 34.504  |
| 633    | POR-73    | TP3        | Porites evermanni        | timeseries | 1.75       | 1.74       | 1.745   | 76        | 132.62  |
| 637    | POR-72    | TP3        | Porites evermanni        | timeseries | 0.52       | 0.516      | 0.518   | 76        | 39.368  |
| 641    | POR-74    | TP3        | Porites evermanni        | timeseries | 0.484      | 0.48       | 0.482   | 76        | 36.632  |
| 669    | POR-83    | TP3        | Porites evermanni        | timeseries | 0.296      | 0.297      | 0.2965  | 76        | 22.534  |
| 711    | POR-83    | TP4        | Porites evermanni        | timeseries | 0.298      | 0.301      | 0.2995  | 76        | 22.762  |
| 713    | POR-69    | TP4        | Porites evermanni        | timeseries | 0.792      | 0.784      | 0.788   | 76        | 59.888  |
| 741    | POR-72    | TP4        | Porites evermanni        | timeseries | 0.832      | 0.828      | 0.83    | 76        | 63.08   |
| 749    | POR-74    | TP4        | Porites evermanni        | timeseries | 0.66       | 0.648      | 0.654   | 76        | 49.704  |
| 877    | POR-262   | TP4        | Porites evermanni        | timeseries | 0.924      | 0.92       | 0.922   | 76        | 70.072  |

The POC all look pretty good, but POR is still very low. The highest DNA ng POR value is 152. 

#### Gel 

I ran a 1.5% gel at 100 volts for 1 hour. Samples are in the same order as the Qubit. 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/refs/heads/master/images/gel_20240926.jpg)

Gel looks terrible, can't see anything. There are some faint bands for 385, 393, and 401 (all POC). No bands for POR. 
