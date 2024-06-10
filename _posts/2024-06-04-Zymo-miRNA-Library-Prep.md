---
layout: post
title: miRNA library prep
date: '2024-06-04'
categories:
tags: [DNA, RNA, Library prep, Protocols]
projects: Mcap developmental time series 
---

# miRNA library prep for Mcap developmental time series Hawaii 2023

This post details the info about the miRNA library prep steps for the *M. capitata* developmental time series experiment in Hawaii 2023. The github for that project is linked [here](https://github.com/JillAshey/Hawaii_Developmental_TimeSeries). I'm using the [Zymo-Seq miRNA library kit](https://www.zymoresearch.com/products/zymo-seq-mirna-library-kit) for small RNA library prep (I got the 12 prep kit). See Zymo's protocol [here](https://files.zymoresearch.com/protocols/r3006_r3007-zymo-seq_mirna_library_kit.pdf). 

After sequencing n=1 from each time point from the Mcap 2023 project, we decided to move forward with sequencing n=3 from each time point in the ambient. I need to RiboFree library prep 24 more samples. Today, I used samples M6, M10, M24, M36, M47, M62, M73, M86, M7, M11, M26, and M37. M6, M10, M24, M47, M7, M11, and M26 were extracted on [2/22/24](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-02-22-MiniprepPlus-DNA%3ARNA-extractions-McapLarvae-DT.md) and M36, M62, M73, M86 and M37 were extracted on [5/2/24](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-05-02-MiniprepPlus-DNA%3ARNA-extractions-McapLarvae-DT.md). I will use 150-200 ng of RNA from each sample. 

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

Protocol was followed according to this [post](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-01-09-Zymo-miRNA-Library-Prep.md). A few differences: 

- I changed the input RNA amount to ~150-200 ng
- I increased the number of amplification cycles from 16 to 22
- I did Section 1 on 6/3/24, then stored overnight at -20°C. On 6/4/24, I did Sections 2, 3, 4, 5, and QC 

### QC 

Run [DNA Tapestation](https://github.com/meschedl/MESPutnam_Open_Lab_Notebook/blob/master/_posts/2019-07-30-DNA-Tapestation.md) for visualize libraries. Here's an example of what the library should look like on a Tapestation: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/DT_mcap2023/miRNA_library_visual_example.png)

Here's how my samples turned out on 6/4/24. See full report [here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/images/tapestation/DNA_miRNA_2024-06-04.pdf).

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_overview_20240604.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M6_20240604.png) 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M10_20240604.png) 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M24_20240604.png) 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M36_20240604.png) 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M47_20240604.png) 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M62_20240604.png) 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M73_20240604.png) 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M86_20240604.png) 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M7_20240604.png) 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M11_20240604.png) 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M26_20240604.png) 

Ran out of tape for M37. QC for these samples didn't look good at all. Little to no peaks and some had no concentrations at all. Those that did have concentrations were <1ng/uL. I'm guessing that they didn't amplify in the PCR step. 

Today (6/5/24), I will redo the PCR amplification and cleaning steps on my libraries to see if I can amplify them at all. During the original PCR step, sample volume was 29 uL; 6 uL of primer and 35 uL of Taq PreMix was added for a total of 70 uL. 

- 6 uL primer / 29 uL sample = 0.2 (ie primer was 20% of original volume)
- 35 uL Taq PreMix / 29 uL sample = 1.2 (ie premix was 120% of orginal volume

To calculate how much primer and premix volume I will need to reamp: 

- 11 uL final library x 0.2 = 2.2 uL of primer. Changed to 2.5 uL of primer to each sample
- 11 uL final library x 1.2 = 13.2 uL of premix. Changed to 13.5 uL of premix to each sample

I am also going to use the Taq Premix that came with the first miRNA kit I got (I received this kit in January and used it to process samples in MArch). The premix has ~200 uL, which is enough for the reamp. I'm using the first premix because it definitely worked on my first round of samples. I also accidently stored the Taq PreMix that I used yesterday at -80C for a week instead of -20C. Maybe storage at -80C messed the premix up? 

Redid sections 4, 5, and QC on 6/5/24. Here's how my samples turned out on 6/5/24. See full report [here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/images/tapestation/DNA_miRNA_2024-06-05.pdf).

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_overview_20240605.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M6_20240605.png) 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M10_20240605.png) 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M24_20240605.png) 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M36_20240605.png) 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M47_20240605.png) 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M62_20240605.png) 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M73_20240605.png) 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M86_20240605.png) 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M7_20240605.png) 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M11_20240605.png) 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M26_20240605.png) 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M37_20240605.png) 

Much better!!!! Concentrations are all relatively higher and there are peaks in all samples. M6, M86, M7, M11, and M26 all have more than one peak, meaning I probably have to do a bead clean-up on the samples and then re-QC them. I also need to reach out to Zymo to ask if they can send me more Taq PreMix (if the one from the second kit was a dud). 