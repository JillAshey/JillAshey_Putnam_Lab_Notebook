---
layout: post
title: Pico Methyl-Seq library prep
date: '2024-10-15'
categories:
tags: [DNA, Library prep, Protocols, WGBS]
projects: e5
---

# Pico Methy-Seq Library Prep test for E5 samples

This post details the info about the WGBS library prep steps for the e5 deep dive samples collected in Moorea in 2020. The github for that project is linked [here](https://github.com/urol-e5/deep-dive). I'm using the [Zymo Pico Methyl Seq Library Prep](https://www.zymoresearch.com/products/pico-methyl-seq-library-prep-kit) for library prep (we have several kits with 25 preps). See Zymo's protocol [here](https://files.zymoresearch.com/protocols/_d5455_d5456_picomethylseq.pdf). 

I prepped the 15 deep dive E5 samples in June 2024 (see [notebook post](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-06-13-Zymo-Pico-Methyl-Seq-Library-Prep.md)) but 6 samples failed to amplify (3 *Pocillopora tuahiniensis* and 3 *Porites evermanni*). I re-tried the 6 samples that failed in August 2024 (see [notebook post](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-08-07-Zymo-Pico-Methyl-Seq-Library-Prep.md)) but these failed as well. I got slight peaks for the POC samples but concentrations were low (<1ng/ul). I re-tried again with these samples in late August 2024 (see [notebook post](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-08-28-Zymo-Pico-Methyl-Seq-Library-Prep.md)) but these failed once again, depsite increasing the number of cycles for POC samples and decreasing the DNA input amount for POR samples. 

Going to retry one more time with these samples. For POC samples, I will increase the DNA input amount slightly (35-50ng) and increase the number of cycles to 15. For the POR samples, I will use the cleaned samples (see [notebook post](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-10-14-Zymo-Cleanup-DNA-E5.md)) and keep the number of cycles the same (10). I am also running one of Zoe's samples from her LCM experiment. 

| Number | Species                  | Timepoint | Collection Date | Site  | colony_id | Extraction Date | Extraction notebook post                                                                                                                                                                                                                                     | DNA (ng/uL)     | Volume eluted (uL) | Total DNA (ng)  | Starting volume (uL) | Volume for DNA (uL) | Tris (uL) | DNA input amount (ng) |
| ------ | ------------------------ | --------- | --------------- | ----- | --------- | --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | --------------- | ------------------ | --------------- | -------------------- | ------------------- | --------- | --------------------- |
| 385    | Pocillopora tuahiniensis | TP2       | 20200305        | Site1 | POC-47    | 20211012        | [https://kterpis.github.io/Putnam_Lab_Notebook/20211012-RNA-DNA-extractions-from-E5-project/](https://kterpis.github.io/Putnam_Lab_Notebook/20211012-RNA-DNA-extractions-from-E5-project/)                                                                   | 25.5            | 90                 | 2295.00         | 20                   | 2.0                 | 18.0      | 51                    |
| 393    | Pocillopora tuahiniensis | TP2       | 20200305        | Site1 | POC-50    | 20210903        | [https://github.com/Kterpis/Putnam_Lab_Notebook/blob/master/_posts/2021-09-03-20210903-RNA-DNA-extractions-from-E5-project.md](https://github.com/Kterpis/Putnam_Lab_Notebook/blob/master/_posts/2021-09-03-20210903-RNA-DNA-extractions-from-E5-project.md) | 29.4            | 90                 | 2646.00         | 20                   | 2.0                 | 18.0      | 58.8                  |
| 401    | Pocillopora tuahiniensis | TP2       | 20200305        | Site1 | POC-53    | 20211118        | [https://kterpis.github.io/Putnam_Lab_Notebook/20211118-RNA-DNA-extractions-from-E5-project/](https://kterpis.github.io/Putnam_Lab_Notebook/20211118-RNA-DNA-extractions-from-E5-project/)                                                                   | 30.8            | 90                 | 2772.00         | 20                   | 2.0                 | 18.0      | 61.6                  |
| 421    | Porites evermanni        | TP2       | 20200305        | Site1 | POR-82    | 20211008        | [https://kterpis.github.io/Putnam_Lab_Notebook/20211008-RNA-DNA-extractions-from-E5-project/](https://kterpis.github.io/Putnam_Lab_Notebook/20211008-RNA-DNA-extractions-from-E5-project/)                                                                   | 0.580           | 80                 | 46.40           | 20                   | 20                  | 0         | 11.6                  |
| 487    | Porites evermanni        | TP2       | 20200305        | Site1 | POR-71    | 20211122        | [https://github.com/Kterpis/Putnam_Lab_Notebook/blob/master/_posts/2021-11-22-20211122-RNA-DNA-extractions-from-E5-project.md](https://github.com/Kterpis/Putnam_Lab_Notebook/blob/master/_posts/2021-11-22-20211122-RNA-DNA-extractions-from-E5-project.md) | 0.954           | 80                 | 76.32           | 20                   | 20                  | 0         | 19.08                 |
| 489    | Porites evermanni        | TP2       | 20200305        | Site1 | POR-79    | 20211129        | [https://github.com/Kterpis/Putnam_Lab_Notebook/blob/master/_posts/2021-11-29-20211129-RNA-DNA-extractions-from-E5-project.md](https://github.com/Kterpis/Putnam_Lab_Notebook/blob/master/_posts/2021-11-29-20211129-RNA-DNA-extractions-from-E5-project.md) | 0.804           | 80                 | 64.32           | 20                   | 20                  | 0         | 16.08                 |
| 421    | Porites evermanni        | TP2       | 20200305        | Site1 | POR-82    | 20240925        | [https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-09-25-MiniprepPlus-DNA-extractions-E5.md](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-09-25-MiniprepPlus-DNA-extractions-E5.md)   | 0.116           | 75                 | 8.70            | 20                   | 20                  | 0         | 2.32                  |
| 487    | Porites evermanni        | TP2       | 20200305        | Site1 | POR-71    | 20240925        | [https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-09-25-MiniprepPlus-DNA-extractions-E5.md](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-09-25-MiniprepPlus-DNA-extractions-E5.md)   | 0.240           | 75                 | 18.00           | 20                   | 20                  | 0         | 4.8                   |
| 489    | Porites evermanni        | TP2       | 20200305        | Site1 | POR-79    | 20240925        | [https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-09-25-MiniprepPlus-DNA-extractions-E5.md](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-09-25-MiniprepPlus-DNA-extractions-E5.md)   | 0.141           | 75                 | 10.58           | 20                   | 20                  | 0         | 2.82                  |
| 13     | Pocillopora acuta ZD     | NA        | NA              | NA    | NA        | 20240908        | [https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/LCM-Exp-Test-Extraction/](https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/LCM-Exp-Test-Extraction/)                                                                                                   | See ZD notebook | See ZD notebook    | See ZD notebook | 20                   | 5                   | 15        | See ZD notebook       |

Here's the Pico Methyl-Seq library prep workflow: 

![](https://raw.githubusercontent.com/meschedl/MESPutnam_Open_Lab_Notebook/master/images/PMS-workflow.png) 

### Protocol 

I followed this [protocol](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-06-13-Zymo-Pico-Methyl-Seq-Library-Prep.md) with a few modifications: 

- Section XXXXX



### QC

Run [DNA Tapestation](https://github.com/meschedl/MESPutnam_Open_Lab_Notebook/blob/master/_posts/2019-07-30-DNA-Tapestation.md) for visualize libraries. Here's an example of what the library should look like on a Tapestation: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/pico_lib_prep_library_example.png)

Here's what my libraries looked like on 10/15/24. See full results XXXXXX


