---
layout: post
title: MiniPrep Plus DNA/RNA extractions
date: '2024-03-27'
categories:
tags: [DNA, RNA, Extractions, Protocols]
projects: Moorea Poc spawning 2023 
---

# Extractions for Moorea Pocillopora spawning 2023

This protocol is based on Maggie’s coral egg/larvae extraction [protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/Larvae-Ex-Protocol/) and Kevin’s coral larvae extraction [protocol](https://kevinhwong1.github.io/KevinHWong_Notebook/DNA-RNA-Extractions-on-P.-astreoides-larvae-BEAD-BEATING/). This post details the info about the extraction steps for larval *Pocillopora* samples from Moorea 2023. The github for that project is linked [here](https://github.com/hputnam/Pocillopora_Spawning_Moorea). 

### Samples 

The samples run today were 127, 129, 131, and 133. These samples had 400 ul of larvae per tube and were preserved in 500 uL of DNA/RNA Shield. They were transported in a freezer box that also included formalin-fixed samples and were stored at 4°C for ~3 months. This photo is post-bead beating: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/samples_20240327.JPG)

I also did a second extraction run today with samples 6, 7, 107, 108, 110, and 113. These samples either had 100 uL of concentrated embryos (samples 6 & 7) or 200 uL of concentrated larvae (107, 108, 110, 113). 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/samples_20240327_2.JPG) 

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

- I added 500 uL of DNA/RNA shield to the samples after thawing them so that the volume in the tubes was ~1 mL of shield. I aliquoted out 500 uL for extractions and saved the remaining fraction in case we need to re-extract. 
- In the RNA extraction, I added another wash step with 400 uL of wash buffer after the DNase incubation (i.e., after the incubation, I did 400 uL of prep buffer, 700 uL of wash buffer, 400 uL of wash buffer, and another 400 uL of wash buffer).
- I eluted in 80 uL of either Tris (DNA) or DNA/RNA free water (RNA) instead of 100 uL. 
- In the second extraction of the day, I tipped the RNA samples onto the floor when trying to close lids to the caps. This was during the second incubation for RNA elution. 

### QC 

#### Qubit results 

These results are for both extractions: 

| ID  | Extraction batch | Qubit DNA1 (ng/uL) | Qubit DNA2 (ng/uL) | Qubit DNA Average (ng/uL) | Qubit RNA1 (ng/uL) | Qubit RNA2 (ng/uL) | Qubit RNA Average (ng/uL) |
| --- | ---------------- | ------------------ | ------------------ | ------------------------- | ------------------ | ------------------ | ------------------------- |
| 127 | 1                | 13.7               | 13.8               | 13.75                     | NA                 | NA                 | NA                        |
| 129 | 1                | 15.8               | 15.7               | 15.75                     | 37.4               | 37.6               | 37.5                      |
| 131 | 1                | 3.04               | 2.96               | 3                         | NA                 | NA                 | NA                        |
| 133 | 1                | 3.98               | 3.58               | 3.78                      | 19.6               | 18.8               | 19.2                      |
| 6   | 2                | NA                 | NA                 | NA                        | 17.4               | 17.2               | 17.3                      |
| 7   | 2                | NA                 | NA                 | NA                        | 21.8               | 21.6               | 21.7                      |
| 107 | 2                | 2.34               | 2.24               | 2.29                      | NA                 | NA                 | NA                        |
| 108 | 2                | NA                 | NA                 | NA                        | NA                 | NA                 | NA                        |
| 110 | 2                | 4.02               | 3.92               | 3.97                      | 11                 | 11                 | 11                        |
| 113 | 2                | 2.8                | 2.7                | 2.75                      | NA                 | NA                 | NA                        |


#### Gel 

Gel from first extraction

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/gel_20240327.JPG)

This gel does not look great and the RNA seems very degraded. These samples were larvae collected 3 days post-fertilization. 127 and 131 were exposed to high temperatures (31°C), and 129 and 133 were exposed to ambient temperatures (28°C). Interesting that we only got RNA qubit values for the ambient treatment. Hollie hypothesized that maybe the larvae in the high treatment died in larger numbers so there was no good RNA. 

Gel from second extraction 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/gel_20240327_2.JPG)

This gel looks much better! There are bands for all samples, even if some of the bands are faint. The two brightest bands are the embryo samples (6 & 7) in ambient temperatures, and the others are larvae 8hpf in high (107 & 113) and ambient (108 & 110) temperatures. We did not get RNA qubit values for all of the samples. Tomorrow, I'm going to run an RNA tapestation on the 6 samples from the second extraction batch. Based on how the tapestation looks, we will determine if we want to use these samples for the Zymo 3' switchfree library prep kit. In order to use the library prep kit, we need to have some idea of the concentration of the samples. 

On the following day (3/28/24), I also ran a RNA tapestation on the six extractions from the second extraction batch to see if the lower concentration samples could be picked up. See full tapestation report [here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/images/tapestation/RNA_POC_2023-03-28.pdf). 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/RNA_TS_overview_20240328.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/RNA_TS_6_20240328.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/RNA_TS_7_20240328.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/RNA_TS_107_20240328.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/RNA_TS_108_20240328.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/RNA_TS_110_20240328.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/RNA_TS_113_20240328.png)

There is RNA in the samples with no qubit values, which is good. It is at a pretty low concentration and the peaks are not very high. Based on these concentrations, we can still likely move forward with library prep. 