---
layout: post
title: miRNA library prep
date: '2024-03-03'
categories:
tags: [DNA, RNA, Library prep, Protocols]
projects: Mcap developmental time series 
---

# miRNA library prep for Mcap developmental time series Hawaii 2023

This post details the info about the miRNA library prep steps for the *M. capitata* developmental time series experiment in Hawaii 2023. The github for that project is linked [here](https://github.com/JillAshey/Hawaii_Developmental_TimeSeries). I'm using the [Zymo-Seq miRNA library kit](https://www.zymoresearch.com/products/zymo-seq-mirna-library-kit) for small RNA library prep (I got the 12 prep kit). See Zymo's protocol [here](https://files.zymoresearch.com/protocols/r3006_r3007-zymo-seq_mirna_library_kit.pdf). 

I'll be using samples M60 and M72 for this test. They were extracted on [2/10/24](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-02-10-MiniprepPlus-DNA%3ARNA-extractions-McapLarvae-DT.md) and had an average of 36.5 and 24.1 ng/ul of RNA, respectively. I already did a test with these samples on [2/29/24](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-01-09-Zymo-miRNA-Library-Prep.md). I used ~100 ng of input RNA, but it appears to have been too low of an amount. Today, I'm going to add 6.5 uL of RNA without diluting it with ultrapure water, which should be between 150-200 ng. 

Here's the miRNA library prep workflow: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/miRNA_lib_prep_workflow.png)

### Materials 

- [Zymo-Seq miRNA library kit](https://www.zymoresearch.com/products/zymo-seq-mirna-library-kit)
- PCR tubes 
- Thermocycler 
- Heating block 
- Mini centrifuge
- Bench top centrifuge  
- Aluminum beads (to keep things on ice)
- Magnetic stand for PCR tubes 
- Kit contents 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/miRNA_lib_prep_contents.png)

### Materials 

Protocol was followed according to this [post](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-01-09-Zymo-miRNA-Library-Prep.md). The only change was the input RNA amount. M60 got primer 3; M72 got primer 4. 

### QC 

Run [DNA Tapestation](https://github.com/meschedl/MESPutnam_Open_Lab_Notebook/blob/master/_posts/2019-07-30-DNA-Tapestation.md) for visualize libraries. Here's an example of what the library should look like on a Tapestation: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/DT_mcap2023/ribofree_library_visual_example.png)

Here's how my samples turned out on 3/3/24. See full report [here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/images/tapestation/DNA_ZymomiRNA_2024-03-03.pdf).

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_overview_20240303.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M60_20240303.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M72_20240303.png)

Increasing the input ng amount did not seem to improve the results. The peaks are still very low, as are the concentrations. My next step is to increase the number of PCR amplification cycles. I have to be mindful with the samples, as I only have 8 preps left. 