---
layout: post
title: RiboFree library prep
date: '2024-02-26'
categories:
tags: [DNA, RNA, Library prep, Protocols]
projects: Mcap developmental time series 
---

# RiboFree library prep test for Mcap developmental time series Hawaii 2023

This post details the info about the library prep steps for the *M. capitata* developmental time series experiment in Hawaii 2023. The github for that project is linked [here](https://github.com/JillAshey/Hawaii_Developmental_TimeSeries). I'm using the [Zymo-Seq RiboFree Total RNA Library Kit](https://www.zymoresearch.com/products/zymo-seq-ribofree-total-rna-library-kit) for total RNA library prep (I got the 12 prep kit). See Zymo's protocol [here](https://files.zymoresearch.com/protocols/r3000_zymo-seq_ribofree_total_rna_library_kit.pdf). 

I'll be using samples M9, M13, M23, M35, M52 and M85. M9 was extracted on [2/22/24](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-02-22-MiniprepPlus-DNA%3ARNA-extractions-McapLarvae-DT.md), and all of the other samples were extracted on [2/10/24](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-02-10-MiniprepPlus-DNA%3ARNA-extractions-McapLarvae-DT.md). These samples, along with M60 and M72 (prepped on [2/19/24](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-02-19-Zymo-RiboFree-Library-Prep.md)), will be sent for sequencing. The samples represent n=1 from each time point during *M. capitata* development. The kit needs a minimum input of 10 of total RNA or a maximum input of 250 ng of total RNA. I'm going to use 190 ng as my input concentration; Maggie did a similar input for her [ribofree preps](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/zribo-lib-RNA-second/)). Here's the breakdown of input RNA volumes for each sample: 

| TubeID | Qubit RNA Average (ng/uL) | Strip tube # | RNA (uL) | DNA/RNA free water (uL) | Total starting volume (uL) | Primer |
| ------ | ------------------------- | ------------ | -------- | ----------------------- | -------------------------- | ------ |
| M9     | 16.6                      | 1            | 8        | 0                       | 8                          | 5      |
| M13    | 58.5                      | 2            | 3.5      | 4.5                     | 8                          | 6      |
| M23    | 15.4                      | 3            | 8        | 0                       | 8                          | 7      |
| M35    | 29                        | 4            | 6.5      | 1.5                     | 8                          | 8      |
| M52    | 19                        | 5            | 8        | 0                       | 8                          | 9      |
| M85    | 30.3                      | 6            | 6.5      | 1.5                     | 8                          | 10     |


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

Protocol was followed according to this [post](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-01-12-Zymo-RiboFree-Library-Prep.md). The only difference is that the number of PCR library amplification cycles was changed to 15 from 11. M85 also got ~8 uL of Depletion Reagent 1 instead of 10 uL due to a pipetting error. See table above for which primer number was used for each sample.  

### QC 

Run [DNA Tapestation](https://github.com/meschedl/MESPutnam_Open_Lab_Notebook/blob/master/_posts/2019-07-30-DNA-Tapestation.md) for visualize libraries. Here's an example of what the library should look like on a Tapestation: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/DT_mcap2023/ribofree_library_visual_example.png)

Here's how my samples turned out on 2/26/24. See full report [here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/images/tapestation/DNA_RiboFree_2024-02-26.pdf). 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_overview_20240226.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M9_20240226.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M13_20240226.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M23_20240226.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M35_20240226.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M52_20240226.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M85_20240226.png)

Peaks are similar in height and width to my previous [prep](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-02-19-Zymo-RiboFree-Library-Prep.md) and similar to Maggie's [results](https://github.com/meschedl/MESPutnam_Open_Lab_Notebook/blob/master/tapestation_pdfs/2019-09-11%20-%2009.28.34.pdf) in terms of the bp distribution. Overall, I am happy with these libraries and they will be sent for sequencing. 