---
layout: post
title: Pico Methyl-Seq library prep
date: '2024-06-13'
categories:
tags: [DNA, Library prep, Protocols, WGBS]
projects: e5
---

# Pico Methy-Seq Library Prep test for E5 samples

This post details the info about the WGBS library prep steps for the e5 deep dive samples collected in Moorea in 2020. The github for that project is linked [here](https://github.com/urol-e5/deep-dive). I'm using the [Zymo Pico Methyl Seq Library Prep](https://www.zymoresearch.com/products/pico-methyl-seq-library-prep-kit) for library prep (we have several kits with 25 preps). See Zymo's protocol [here](https://files.zymoresearch.com/protocols/_d5455_d5456_picomethylseq.pdf). 

I prepped the 15 deep dive E5 samples in June 2024 (see [notebook post](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-06-13-Zymo-Pico-Methyl-Seq-Library-Prep.md)) but 6 samples failed to amplify (3 *Pocillopora tuahiniensis* and 3 *Porites evermanni*). I re-tried the 6 samples that failed in August 2024 (see [notebook post](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-08-07-Zymo-Pico-Methyl-Seq-Library-Prep.md)) but these failed as well. I got slight peaks for the POC samples but concentrations were low (<1ng/ul). I decided to retry these samples again doing different things for the POC and POR samples. Because the POC samples did have tiny peaks, I'm going to increase the library amplification cycles (Section 4) from 9 to 12 cycles. The POR samples do not look like they amplified at all. Kevin did some WGBS prep with Porites astreoides and he found that 10ng input for Porites worked well (see his [post](https://kevinhwong1.github.io/KevinHWong_Notebook/Thermal-Transplant-WGBS-PicoMethyl-Protocol/)). For the POR samples, I'm going to decrease the DNA input from 25ng to 10ng. 

| Number | Species                  | colony_id | Extraction Date | Extraction notebook post                                                                                                                                                                                                                                     | DNA (ng/uL) | Volume for 25 ng (POC) or 10 ng (POR) DNA (uL) | Tris (uL) | Starting volume (uL) | Primer |
| ------ | ------------------------ | --------- | --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ----------- | ---------------------------------------------- | --------- | -------------------- | ------ |
| 385    | Pocillopora tuahiniensis | POC-47    | 20211012        | [https://kterpis.github.io/Putnam_Lab_Notebook/20211012-RNA-DNA-extractions-from-E5-project/](https://kterpis.github.io/Putnam_Lab_Notebook/20211012-RNA-DNA-extractions-from-E5-project/)                                                                   | 25.5        | 1.0                                            | 19.0      | 20                   | 22     |
| 393    | Pocillopora tuahiniensis | POC-50    | 20210903        | [https://github.com/Kterpis/Putnam_Lab_Notebook/blob/master/_posts/2021-09-03-20210903-RNA-DNA-extractions-from-E5-project.md](https://github.com/Kterpis/Putnam_Lab_Notebook/blob/master/_posts/2021-09-03-20210903-RNA-DNA-extractions-from-E5-project.md) | 29.4        | 0.9                                            | 19.1      | 20                   | 23     |
| 401    | Pocillopora tuahiniensis | POC-53    | 20211118        | [https://kterpis.github.io/Putnam_Lab_Notebook/20211118-RNA-DNA-extractions-from-E5-project/](https://kterpis.github.io/Putnam_Lab_Notebook/20211118-RNA-DNA-extractions-from-E5-project/)                                                                   | 30.8        | 0.8                                            | 19.2      | 20                   | 24     |
| 421    | Porites evermanni        | POR-82    | 20211008        | [https://kterpis.github.io/Putnam_Lab_Notebook/20211008-RNA-DNA-extractions-from-E5-project/](https://kterpis.github.io/Putnam_Lab_Notebook/20211008-RNA-DNA-extractions-from-E5-project/)                                                                   | 3.41        | 2.9                                            | 17.1      | 20                   | 28     |
| 487    | Porites evermanni        | POR-71    | 20211122        | [https://github.com/Kterpis/Putnam_Lab_Notebook/blob/master/_posts/2021-11-22-20211122-RNA-DNA-extractions-from-E5-project.md](https://github.com/Kterpis/Putnam_Lab_Notebook/blob/master/_posts/2021-11-22-20211122-RNA-DNA-extractions-from-E5-project.md) | 2.31        | 4.3                                            | 15.7      | 20                   | 34     |
| 489    | Porites evermanni        | POR-79    | 20211129        | [https://github.com/Kterpis/Putnam_Lab_Notebook/blob/master/_posts/2021-11-29-20211129-RNA-DNA-extractions-from-E5-project.md](https://github.com/Kterpis/Putnam_Lab_Notebook/blob/master/_posts/2021-11-29-20211129-RNA-DNA-extractions-from-E5-project.md) | 3.72        | 2.7                                            | 17.3      | 20                   | 35     |


Here's the Pico Methyl-Seq library prep workflow: 

![](https://raw.githubusercontent.com/meschedl/MESPutnam_Open_Lab_Notebook/master/images/PMS-workflow.png) 

### Protocol 

I followed this [protocol](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-06-13-Zymo-Pico-Methyl-Seq-Library-Prep.md) with a few modifications: 

- Section 1 - incubation with warmed DNA elubtion buffer shortened from 5 to 3 minutes. Emma reported that she typically got less out of the column if she did the full 5 minutes. 
- Section 3 - used 14 uL of warmed DNA elution buffer instead of 12 (to make sure I have enough volume for PCR). Incubation with warmed DNA elubtion buffer shortened from 5 to 2.5 minutes.
- Section 4 - increased number of cycles for POC samples from 9 to 12. 
- Section 5 - incubation with warmed DNA elubtion buffer shortened from 5 to 2.5 minutes.
- Section 6 - thought we initially ran out of Library Master Mix so initially Zoe and I were only going to run 2 and 3 samples, respectively. I was only going to run my POR samples but we found more Library Master Mix in Sam Gurr's RRBS box.

### QC

Run [DNA Tapestation](https://github.com/meschedl/MESPutnam_Open_Lab_Notebook/blob/master/_posts/2019-07-30-DNA-Tapestation.md) for visualize libraries. Here's an example of what the library should look like on a Tapestation: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/pico_lib_prep_library_example.png)

Here's what my libraries looked like on 8/28/24. See full results [here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/images/tapestation/DNA_Pico-2024-08-28.pdf). 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_overview_20240828.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_385_20240828.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_393_20240828.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_401_20240828.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_421_20240828.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_487_20240828.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_489_20240828.png)

Still no amplification. Briefly talked with Hollie about it. She suggested that I nanodrop the DNA and library cDNA to see if there are any impurities before moving forward with any more troubleshooting. But good news is that Zoe's LCM samples worked!! 
