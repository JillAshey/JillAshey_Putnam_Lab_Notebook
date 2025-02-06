---
layout: post
title: SwitchFree library prep
date: '2025-02-05'
categories:
tags: [DNA, RNA, Library prep, Protocols]
projects: Kit test 
---

# SwitchFree library prep

This post details the info about the SwitchFree library prep. See Zymo's protocol [here](https://files.zymoresearch.com/protocols/_r3008_r3009__zymo_seq_switchfree_3_mrna_library_kit.pdf). 

I will be prepping my 32 Mcap DT 2023 samples (see project github [repo](https://github.com/JillAshey/DevelopmentalTimeseries) for more information about the experiment); today, I am prepping 10 samples, all redos (see information about preps on [1/30/25](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2025-01-30-Zymo-SwitchFree_McapLarvae-DT.md), [2/1/25](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2025-02-01-Zymo-SwitchFree_McapLarvae-DT.md), and [2/2/25](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2025-02-02-Zymo-SwitchFree_McapLarvae-DT.md). 

I already sequenced these using rRNA depletion and smRNA sequencing but we are going to do polyA sequencing as well to better elucidate the MZT. The kit needs a minimum of 10 ng of total RNA or a maximum of 500 ng of total RNA, which is a large range. In my Ribofree and miRNA preps with the same samples, I used 190 and 100 ng per sample, respectively. I already tried to prep these samples using the switchfree kit using 100ng of input, but they didn't do well. I'm now going to try with 50 ng RNA to see if that improves library quality. 

| TubeID | Starting volume (uL) | Volume for 50 ng RNA (uL) | Volume ultrapure water (uL) | Primer |
| ------ | -------------------- | ------------------------- | --------------------------- | ------ |
| M7     | 5                    | 1.8                       | 3.2                         | UDI_39 |
| M28    | 5                    | 1.6                       | 3.4                         | UDI_40 |
| M36    | 5                    | 1.4                       | 3.6                         | UDI_41 |
| M37    | 5                    | 1.4                       | 3.6                         | UDI_42 |
| M51    | 5                    | 1.6                       | 3.4                         | UDI_43 |
| M63    | 5                    | 1.7                       | 3.3                         | UDI_44 |
| M73    | 5                    | 1.6                       | 3.4                         | UDI_45 |
| M74    | 5                    | 1.5                       | 3.5                         | UDI_46 |
| M85    | 5                    | 1.7                       | 3.3                         | UDI_47 |
| M86    | 5                    | 0.7                       | 4.3                         | UDI_48 |

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

- Used 50ng RNA input instead of 100ng
- I eluted in 15 uL of DNA elution buffer.

#### QC

Run [DNA Tapestation](https://github.com/meschedl/MESPutnam_Open_Lab_Notebook/blob/master/_posts/2019-07-30-DNA-Tapestation.md) to visualize libraries. Here's an example of what the library should look like on a Tapestation: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/switchfree_lib_prep_library_example.png)

I ran the tapestation on the libraries on 2/5/25. See the full report [here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/images/tapestation/DNA_switchfree_mcap2023_2025-02-05.pdf). 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/refs/heads/master/images/tapestation/DNA_TS_overview_20250205.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/refs/heads/master/images/tapestation/DNA_TS_M7_20250205.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/refs/heads/master/images/tapestation/DNA_TS_M28_20250205.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/refs/heads/master/images/tapestation/DNA_TS_M36_20250205.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/refs/heads/master/images/tapestation/DNA_TS_M37_20250205.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/refs/heads/master/images/tapestation/DNA_TS_M51_20250205.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/refs/heads/master/images/tapestation/DNA_TS_M63_20250205.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/refs/heads/master/images/tapestation/DNA_TS_M73_20250205.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/refs/heads/master/images/tapestation/DNA_TS_M74_20250205.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/refs/heads/master/images/tapestation/DNA_TS_M85_20250205.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/refs/heads/master/images/tapestation/DNA_TS_M86_20250205.png)