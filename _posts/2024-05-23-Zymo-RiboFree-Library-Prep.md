---
layout: post
title: RiboFree library prep
date: '2024-05-23'
categories:
tags: [DNA, RNA, Library prep, Protocols]
projects: Mcap developmental time series 
---

# RiboFree library prep test for Mcap developmental time series Hawaii 2023

This post details the info about the library prep steps for the *M. capitata* developmental time series experiment in Hawaii 2023. The github for that project is linked [here](https://github.com/JillAshey/Hawaii_Developmental_TimeSeries). I'm using the [Zymo-Seq RiboFree Total RNA Library Kit](https://www.zymoresearch.com/products/zymo-seq-ribofree-total-rna-library-kit) for total RNA library prep (I got the 12 prep kit). See Zymo's protocol [here](https://files.zymoresearch.com/protocols/r3000_zymo-seq_ribofree_total_rna_library_kit.pdf). 

After sequencing n=1 from each time point from the Mcap 2023 project, we decided to move forward with sequencing n=3 from each time point in the ambient. I need to RiboFree library prep 24 more samples. Today, I used samples M6, M10, M24, M36, M47, M62, M73, M86, M7, M11, M26, and M37. M6, M10, M24, M47, M7, M11, and M26 were extracted on [2/22/24](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-02-22-MiniprepPlus-DNA%3ARNA-extractions-McapLarvae-DT.md) and M36, M62, M73, M86 and M37 were extracted on [5/2/24](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-05-02-MiniprepPlus-DNA%3ARNA-extractions-McapLarvae-DT.md). 

The kit needs a minimum input of 10 of total RNA or a maximum input of 250 ng of total RNA. I'm used 190 ng as my input concentration; Maggie did a similar input for her [ribofree preps](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/zribo-lib-RNA-second/)). Here's the breakdown of input RNA volumes for each sample: 

| TubeID | Qubit RNA Average (ng/uL) | RNA (uL) | DNA/RNA free water (uL) | Total starting volume (uL) | Primer |
| ------ | ------------------------- | -------- | ----------------------- | -------------------------- | ------ |
| M6     | 60.9                      | 3.12     | 4.88                    | 8                          | 24     |
| M10    | 72                        | 2.64     | 5.36                    | 8                          | 25     |
| M24    | 94.2                      | 2.02     | 5.98                    | 8                          | 26     |
| M36    | 36.8                      | 5.16     | 2.84                    | 8                          | 27     |
| M47    | 81.2                      | 2.34     | 5.66                    | 8                          | 28     |
| M62    | 45                        | 4.22     | 3.78                    | 8                          | 29     |
| M73    | 30.7                      | 6.19     | 1.81                    | 8                          | 30     |
| M86    | 75.9                      | 2.50     | 5.50                    | 8                          | 21     |
| M7     | 27.2                      | 6.99     | 1.01                    | 8                          | 32     |
| M11    | 67.6                      | 2.81     | 5.19                    | 8                          | 33     |
| M26    | 197                       | 0.96     | 7.04                    | 8                          | 34     |
| M37    | 34.9                      | 5.44     | 2.56                    | 8                          | 35     |

Here's the RiboFree library prep workflow: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/DT_mcap2023/ribofree_workflow.png)

### Materials 

- [Zymo-Seq RiboFree library kit](https://www.zymoresearch.com/products/zymo-seq-ribofree-total-rna-library-kit)
- PCR tubes 
- Thermocycler 
- Heating block 
- Mini centrifuge
- Vortex 
- Aluminum beads (to keep things on ice)
- Magnetic stand for PCR tubes 
- 1.5 mL tubes 
- Kit contents 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/DT_mcap2023/ribofree_kit_contents.png)

### Protocol 

Protocol was followed according to this [post](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-01-12-Zymo-RiboFree-Library-Prep.md). The only difference is that the number of PCR library amplification cycles was changed to 15 from 11. 

Sections 1, 2, and 3 were completed on 5/23/24. Section 4 was completed on 5/24/24. QC was completed on 5/27/24. 

### QC 

Run [DNA Tapestation](https://github.com/meschedl/MESPutnam_Open_Lab_Notebook/blob/master/_posts/2019-07-30-DNA-Tapestation.md) for visualize libraries. Here's an example of what the library should look like on a Tapestation: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/DT_mcap2023/ribofree_library_visual_example.png)

Here's how my samples turned out on 5/27/24. See full report [here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/images/tapestation/DNA_Ribofree_2024-05-27.pdf). 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_overview_20240527.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M6_20240527.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M10_20240527.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M24_20240527.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M36_20240527.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M47_20240527.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M62_20240527.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M73_20240527.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M86_20240527.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M7_20240527.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M11_20240527.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M26_20240527.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M37_20240527.png)

Some of the peaks are pretty low, but most look good. 