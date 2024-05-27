---
layout: post
title: RiboFree library prep
date: '2024-05-27'
categories:
tags: [DNA, RNA, Library prep, Protocols]
projects: Mcap developmental time series 
---

# RiboFree library prep test for Mcap developmental time series Hawaii 2023

This post details the info about the library prep steps for the *M. capitata* developmental time series experiment in Hawaii 2023. The github for that project is linked [here](https://github.com/JillAshey/Hawaii_Developmental_TimeSeries). I'm using the [Zymo-Seq RiboFree Total RNA Library Kit](https://www.zymoresearch.com/products/zymo-seq-ribofree-total-rna-library-kit) for total RNA library prep (I got the 12 prep kit). See Zymo's protocol [here](https://files.zymoresearch.com/protocols/r3000_zymo-seq_ribofree_total_rna_library_kit.pdf). 

After sequencing n=1 from each time point from the Mcap 2023 project, we decided to move forward with sequencing n=3 from each time point in the ambient. Today, I used samples M8, M14, M28, M39, M48, M51, M61, M63, M74, M75, M87, and M88. M8, M14, M28, M48, M51 were extracted on [2/22/24](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-02-22-MiniprepPlus-DNA%3ARNA-extractions-McapLarvae-DT.md) and M39, M61, M63, M74, M75, M87, and M88 were extracted on [5/2/24](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-05-02-MiniprepPlus-DNA%3ARNA-extractions-McapLarvae-DT.md). 

The kit needs a minimum input of 10 of total RNA or a maximum input of 250 ng of total RNA. I'm used 190 ng as my input concentration; Maggie did a similar input for her [ribofree preps](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/zribo-lib-RNA-second/)). Here's the breakdown of input RNA volumes for each sample: 

| TubeID | Qubit RNA Average (ng/uL) | RNA (uL) | DNA/RNA free water (uL) | Total starting volume (uL) | Primer |
| ------ | ------------------------- | -------- | ----------------------- | -------------------------- | ------ |
| M8     | 39.2                      | 4.8      | 3.2                     | 8                          | 36     |
| M14    | 115                       | 1.7      | 6.3                     | 8                          | 37     |
| M28    | 31.8                      | 6.0      | 2.0                     | 8                          | 38     |
| M39    | 21.3                      | 8.9      | \-0.9                   | 8                          | 39     |
| M48    | 114                       | 1.7      | 6.3                     | 8                          | 40     |
| M51    | 31.5                      | 6.0      | 2.0                     | 8                          | 41     |
| M61    | 29.8                      | 6.4      | 1.6                     | 8                          | 42     |
| M63    | 29.8                      | 6.4      | 1.6                     | 8                          | 43     |
| M74    | 32.4                      | 5.9      | 2.1                     | 8                          | 44     |
| M75    | 27.4                      | 6.9      | 1.1                     | 8                          | 45     |
| M87    | 71.6                      | 2.7      | 5.3                     | 8                          | 46     |
| M88    | 97                        | 2.0      | 6.0                     | 8                          | 47     |

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

### QC 

Run [DNA Tapestation](https://github.com/meschedl/MESPutnam_Open_Lab_Notebook/blob/master/_posts/2019-07-30-DNA-Tapestation.md) for visualize libraries. Here's an example of what the library should look like on a Tapestation: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/DT_mcap2023/ribofree_library_visual_example.png)

Here's how my samples turned out on 5/27/24. See full report here XXXXXX