---
layout: post
title: miRNA library prep
date: '2024-03-04'
categories:
tags: [DNA, RNA, Library prep, Protocols]
projects: Mcap developmental time series 
---

# miRNA library prep for Mcap developmental time series Hawaii 2023

This post details the info about the miRNA library prep steps for the *M. capitata* developmental time series experiment in Hawaii 2023. The github for that project is linked [here](https://github.com/JillAshey/Hawaii_Developmental_TimeSeries). I'm using the [Zymo-Seq miRNA library kit](https://www.zymoresearch.com/products/zymo-seq-mirna-library-kit) for small RNA library prep (I got the 12 prep kit). See Zymo's protocol [here](https://files.zymoresearch.com/protocols/r3006_r3007-zymo-seq_mirna_library_kit.pdf). 

I'll be using samples M60 and M72 for this test. They were extracted on [2/10/24](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-02-10-MiniprepPlus-DNA%3ARNA-extractions-McapLarvae-DT.md) and had an average of 36.5 and 24.1 ng/ul of RNA, respectively. I already did a test with these samples on [2/29/24](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-01-09-Zymo-miRNA-Library-Prep.md) and on [3/3/24](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-03-03-Zymo-RiboFree-Library-Prep.md). In the 2/29/24 test, I used ~100 ng of input RNA, but it appears to have been too low of an amount. In the 3/3/24 test, I used ~150-200 ng of input RNA, but it still appears to be too low. Today, I'm going to keep the same input RNA amount and increase the number of amplification cycles from 16 to 22.

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

Protocol was followed according to this [post](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-01-09-Zymo-miRNA-Library-Prep.md). I changed the input RNA amount to ~150-200 ng and increased the number of amplification cycles from 16 to 22. M60 got primer 5; M72 got primer 6. 

### QC 

Run [DNA Tapestation](https://github.com/meschedl/MESPutnam_Open_Lab_Notebook/blob/master/_posts/2019-07-30-DNA-Tapestation.md) for visualize libraries. Here's an example of what the library should look like on a Tapestation: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/DT_mcap2023/miRNA_library_visual_example.png)

Here's how my samples turned out on 3/4/24. See full report [here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/images/tapestation/DNA_ZymomiRNA_2024-03-04.pdf).

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_overview_20240304.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M60_20240304.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M72_20240304.png)

These results look much better! The peaks are in the correct place (~150 bp) and the concentrations are higher than the other tests. There appears to be some contamination in both samples aroun 470 bp and 720 bp, perhaps due to overamplification. I wonder if I could remove these through another bead clean-up? 

Discussed results with Amy and Willow. Both of them seemed confident that these libraries would be adequate for sequencing. I could do another bead clean-up, but I would lose some cDNA in the process. I think I am going to stick with this protocol for the other 6 samples. I only have 6 more preps with this kit, so I am very limited in troubleshooting capacity. I sent the tapestation results to Genohub to have them confirm if these libraries are adequate for sequencing in their facility. 