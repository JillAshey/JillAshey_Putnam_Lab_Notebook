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

I prepped the 15 deep dive E5 samples in June 2024 (see [notebook post](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-06-13-Zymo-Pico-Methyl-Seq-Library-Prep.md)) but 6 samples failed to amplify (3 *Pocillopora tuahiniensis* and 3 *Porites evermanni*). I decided to retry those samples with new PrepAmp Polymerase and some small modifications to the protocol based on conversations with Emma on what worked for her at GMGI. 

The kit needs a minimum input of 10 ng DNA or a maximum input of 100 ng DNA. Emma suggested using 20 ng, as that amount has yielded better libraries than 10 ng, so I'll go with 25 ng. Both [Maggie](https://github.com/meschedl/MESPutnam_Open_Lab_Notebook/blob/master/_posts/2020-09-18-WGBS-PMS-protocol.md) and [Emma](https://github.com/emmastrand/GMGI_Notebook/blob/main/posts/2023-08-24%20Zymo%20Pico%20Methyl%20Seq%20Kit%20Protocol.md) have done this protocol before with success. 


| Number | Species                  | colony_id | Extraction Date | Extraction notebook post                                                                                                                                                                                                                                     | DNA (ng/uL) | Volume for 25 ng DNA (uL) | Tris (uL) | Starting volume (uL) | Primer |
| ------ | ------------------------ | --------- | --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ----------- | ------------------------- | --------- | -------------------- | ------ |
| 385    | Pocillopora tuahiniensis | POC-47    | 20211012        | [https://kterpis.github.io/Putnam_Lab_Notebook/20211012-RNA-DNA-extractions-from-E5-project/](https://kterpis.github.io/Putnam_Lab_Notebook/20211012-RNA-DNA-extractions-from-E5-project/)                                                                   | 25.5        | 1.0                       | 19.0      | 20                   | 22     |
| 393    | Pocillopora tuahiniensis | POC-50    | 20210903        | [https://github.com/Kterpis/Putnam_Lab_Notebook/blob/master/_posts/2021-09-03-20210903-RNA-DNA-extractions-from-E5-project.md](https://github.com/Kterpis/Putnam_Lab_Notebook/blob/master/_posts/2021-09-03-20210903-RNA-DNA-extractions-from-E5-project.md) | 29.4        | 0.9                       | 19.1      | 20                   | 23     |
| 401    | Pocillopora tuahiniensis | POC-53    | 20211118        | [https://kterpis.github.io/Putnam_Lab_Notebook/20211118-RNA-DNA-extractions-from-E5-project/](https://kterpis.github.io/Putnam_Lab_Notebook/20211118-RNA-DNA-extractions-from-E5-project/)                                                                   | 30.8        | 0.8                       | 19.2      | 20                   | 24     |
| 421    | Porites evermanni        | POR-82    | 20211008        | [https://kterpis.github.io/Putnam_Lab_Notebook/20211008-RNA-DNA-extractions-from-E5-project/](https://kterpis.github.io/Putnam_Lab_Notebook/20211008-RNA-DNA-extractions-from-E5-project/)                                                                   | 3.41        | 7.3                       | 12.7      | 20                   | 28     |
| 487    | Porites evermanni        | POR-71    | 20211122        | [https://github.com/Kterpis/Putnam_Lab_Notebook/blob/master/_posts/2021-11-22-20211122-RNA-DNA-extractions-from-E5-project.md](https://github.com/Kterpis/Putnam_Lab_Notebook/blob/master/_posts/2021-11-22-20211122-RNA-DNA-extractions-from-E5-project.md) | 2.31        | 10.8                      | 9.2       | 20                   | 34     |
| 489    | Porites evermanni        | POR-79    | 20211129        | [https://github.com/Kterpis/Putnam_Lab_Notebook/blob/master/_posts/2021-11-29-20211129-RNA-DNA-extractions-from-E5-project.md](https://github.com/Kterpis/Putnam_Lab_Notebook/blob/master/_posts/2021-11-29-20211129-RNA-DNA-extractions-from-E5-project.md) | 3.72        | 6.7                       | 13.3      | 20                   | 35     |

Here's the Pico Methyl-Seq library prep workflow: 

![](https://raw.githubusercontent.com/meschedl/MESPutnam_Open_Lab_Notebook/master/images/PMS-workflow.png) 

### Protocol 

I followed this [protocol](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-06-13-Zymo-Pico-Methyl-Seq-Library-Prep.md) with a few modifications: 

- Brand new PrepAmp Polymerase! 
- Section 1 - incubation with warmed DNA elubtion buffer shortened from 5 to 3 minutes. Emma reported that she typically got less out of the column if she did the full 5 minutes. 
- Section 3 - used 14 uL of warmed DNA elution buffer instead of 12 (to make sure I have enough volume for PCR). Incubation with warmed DNA elubtion buffer shortened from 5 to 2.5 minutes.
	- For sample 487, I added the DNA elution buffer and then knocked it over. Tried to tap the tube on the benchtop to get the liquid back to the bottom of the column but was unsure if it worked so I added ~7uL more DNA elution buffer to this sample. 
- Section 5 - incubation with warmed DNA elubtion buffer shortened from 5 to 2.5 minutes.

### QC

Run [DNA Tapestation](https://github.com/meschedl/MESPutnam_Open_Lab_Notebook/blob/master/_posts/2019-07-30-DNA-Tapestation.md) for visualize libraries. Here's an example of what the library should look like on a Tapestation: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/pico_lib_prep_library_example.png)

Here's what my samples looked like on 8/7/24. See full results [here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/images/tapestation/DNA_Pico-2024-08-07.pdf). 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_overview_20240807.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_385_20240807.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_393_20240807.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_401_20240807.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_421_20240807.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_487_20240807.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_489_20240807.png)

Got a little bit of amplification for most of the samples but concentrations are low. 393 did not amplify at all. 385 and 401 have concentrations of 1.69 and 3.82 ng/ul, respectively, which might be okay for sequencing. The Porites samples (421, 487, 489) all had concentrations >1 ng/ul, not sure if that is enough for sequencing. I could try to redo the samples but with a higher number of PCR cycles. I am hesitant to do that though, as that may bias these samples or create a large batch effect. 