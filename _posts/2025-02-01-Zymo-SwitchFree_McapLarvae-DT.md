---
layout: post
title: SwitchFree library prep
date: '2025-02-01'
categories:
tags: [DNA, RNA, Library prep, Protocols]
projects: Kit test 
---

# SwitchFree library prep test

This post details the info about the SwitchFree library prep. See Zymo's protocol [here](https://files.zymoresearch.com/protocols/_r3008_r3009__zymo_seq_switchfree_3_mrna_library_kit.pdf). 

I will be prepping my 32 Mcap DT 2023 samples (see project github [repo](https://github.com/JillAshey/DevelopmentalTimeseries) for more information about the experiment); today, I am prepping 13 samples. I prepped another 13 samples a few days ago with moderate success (see [post](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2025-01-30-Zymo-SwitchFree_McapLarvae-DT.md) for more details). Some samples had double peaks which may mean I need lower RNA ng input, but I am going to continue with 100ng for now. 

I already sequenced these using rRNA depletion and smRNA sequencing but we are going to do polyA sequencing as well to better elucidate the MZT. The kit needs a minimum of 10 ng of total RNA or a maximum of 500 ng of total RNA, which is a large range. In my Ribofree and miRNA preps with the same samples, I used 190 and 100 ng per sample, respectively. I am going to use 100ng as input for this library prep. 

| TubeID | Starting volume (uL) | Volume for 100 ng RNA (uL) | Volume ultrapure water (uL) | Primer |
| ------ | -------------------- | -------------------------- | --------------------------- | ------ |
| M9     | 5                    | 5.0                        | 0.0                         | 14     |
| M13    | 5                    | 1.7                        | 3.3                         | 15     |
| M23    | 5                    | 5.0                        | 0.0                         | 16     |
| M28    | 5                    | 3.1                        | 1.9                         | 17     |
| M36    | 5                    | 2.7                        | 2.3                         | 18     |
| M37    | 5                    | 2.9                        | 2.1                         | 19     |
| M52    | 5                    | 5.0                        | 0.0                         | 20     |
| M62    | 5                    | 2.2                        | 2.8                         | 21     |
| M73    | 5                    | 3.3                        | 1.7                         | 22     |
| M74    | 5                    | 3.1                        | 1.9                         | 23     |
| M85    | 5                    | 3.3                        | 1.7                         | 24     |
| M87    | 5                    | 1.4                        | 3.6                         | 25     |
| M88    | 5                    | 1.0                        | 4.0                         | 26     |

I am using the Zymo UDI primers for these preps. See more information and barcode sequences [here](https://www.zymoresearch.com/products/zymo-seq-udi-primer-sets?srsltid=AfmBOoqmYVsF5dEMxwuwu7L6mn6Ot93O6ldOc9wwDwvXUhCkKhg0WK5Y). 

Here's the SwitchFree library prep workflow: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/switchfree_lib_prep_workflow.png)

### Materials 

- [Zymo-Seq SwitchFree 3' mRNA Library Kit](https://www.zymoresearch.com/products/zymo-seq-switchfree-3-mrna-library-kit)
- PCR tubes 
- Thermocycler 
- Mini centrifuge
- Aluminum beads (to keep things on ice)
- Magnetic stand for PCR tubes 
- Kit contents 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/switchfree_lib_prep_contents.png)

### Buffer preperation 

Once buffers are prepared for a kit, they do not need to be prepared again. 

- Add 300 uL of the Select-a-Size MagBead concentrate to 10 mL of Select-a-Size MagBead Buffer. 
	- These should be prepared at least 5 days ahead of library prep. 
- Add 24 mL of 100% ethanol to the DNA Wash Buffer concentrate. Store at room temperature. 

### Best practices 

- Avoid multiple freeze thaws of all components, and aliquot components as necessary
- Remove enzymes from cold storage just before use and return to cold storage immediately after use
- Thaw and maintain all components on ice unless noted otherwise 
- Flick to mix thawed components and briefly centrifuge prior to use 
- After adding each component, mix by pipetting up and down 15-20 times. Briefly centrifuge after 
- Pre-program thermal cycler with lid heating ON set to >100-105°C 
- Turn on thermocycler day-of to preheat 
- Pre-calculate how much to dilute RNA samples 

### Protocol 

Protocol was followed according to this [post](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-03-29-Zymo-SwitchFree.md). 

I eluted in 15 uL of DNA elution buffer ([previous prep](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2025-01-30-Zymo-SwitchFree_McapLarvae-DT.md) eluted in 20uL). 

#### QC

Run [DNA Tapestation](https://github.com/meschedl/MESPutnam_Open_Lab_Notebook/blob/master/_posts/2019-07-30-DNA-Tapestation.md) to visualize libraries. Here's an example of what the library should look like on a Tapestation: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/switchfree_lib_prep_library_example.png)

I ran the tapestation on the libraries on 2/1/25. See the full report XXXXX
