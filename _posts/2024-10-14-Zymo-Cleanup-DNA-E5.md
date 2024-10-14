---
layout: post
title: PCR adapter clean-up e5 samples 
date: '2024-10-14'
categories:
tags: [DNA, Extractions, e5, Protocols]
projects: e5
---

# PCR adapter clean-up w/ e5 deep dive DNA samples  

There are 3 deep dive samples (all POR) that keep failing the Zymo Pico prep (see previous notebook [posts](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-08-28-Zymo-Pico-Methyl-Seq-Library-Prep.md)). There may be some adapter or inhibitor in the POR DNA that is preventing the preps from working. To combat this, we decided to use the Zymo PCR inhibitor remonval kit (see protocol [here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/protocols/d6031_onestep_pcr_inhibitor_removal_kit.pdf)). We received 10 free samples of this kit. 

I am going to clean the DNA extractions that I did (9/25/24) and that Kristina did (fall 2021) to see what will provide a better yield. Here is the sample info: 

| Number | Species           | Timepoint | Collection Date | Site  | colony_id | Extraction Date | Extraction notebook post                                                                                                                                                                                                                                     |
| ------ | ----------------- | --------- | --------------- | ----- | --------- | --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| 421    | Porites evermanni | TP2       | 20200305        | Site1 | POR-82    | 20211008        | [https://kterpis.github.io/Putnam_Lab_Notebook/20211008-RNA-DNA-extractions-from-E5-project/](https://kterpis.github.io/Putnam_Lab_Notebook/20211008-RNA-DNA-extractions-from-E5-project/)                                                                   |
| 487    | Porites evermanni | TP2       | 20200305        | Site1 | POR-71    | 20211122        | [https://github.com/Kterpis/Putnam_Lab_Notebook/blob/master/_posts/2021-11-22-20211122-RNA-DNA-extractions-from-E5-project.md](https://github.com/Kterpis/Putnam_Lab_Notebook/blob/master/_posts/2021-11-22-20211122-RNA-DNA-extractions-from-E5-project.md) |
| 489    | Porites evermanni | TP2       | 20200305        | Site1 | POR-79    | 20211129        | [https://github.com/Kterpis/Putnam_Lab_Notebook/blob/master/_posts/2021-11-29-20211129-RNA-DNA-extractions-from-E5-project.md](https://github.com/Kterpis/Putnam_Lab_Notebook/blob/master/_posts/2021-11-29-20211129-RNA-DNA-extractions-from-E5-project.md) |
| 421    | Porites evermanni | TP2       | 20200305        | Site1 | POR-82    | 20240925        | [https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-09-25-MiniprepPlus-DNA-extractions-E5.md](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-09-25-MiniprepPlus-DNA-extractions-E5.md)   |
| 487    | Porites evermanni | TP2       | 20200305        | Site1 | POR-71    | 20240925        | [https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-09-25-MiniprepPlus-DNA-extractions-E5.md](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-09-25-MiniprepPlus-DNA-extractions-E5.md)   |
| 489    | Porites evermanni | TP2       | 20200305        | Site1 | POR-79    | 20240925        | [https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-09-25-MiniprepPlus-DNA-extractions-E5.md](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-09-25-MiniprepPlus-DNA-extractions-E5.md)   |

### Protocol 

- Loosen filter cap and break bottom tip of filter off 
- Put filter in collection tube 
- Centrifuge for 3 mins at 8000g
- Move filter into 1.5 mL tube and add sample 
- Centrifuge for 3 mins at 8000g

### Qubit results 

Used Qubit HS DNA kit

| Number | Extraction Date | DNA1_ng/uL | DNA2_ng/uL | DNA average |
| ------ | --------------- | ---------- | ---------- | ----------- |
| 421    | 20211008        | 0.586      | 0.574      | 0.58        |
| 487    | 20211122        | 0.967      | 0.941      | 0.954       |
| 489    | 20211129        | 0.81       | 0.798      | 0.804       |
| 421    | 20240925        | 0.117      | 0.114      | 0.1155      |
| 487    | 20240925        | 0.24       | 0.24       | 0.24        |
| 489    | 20240925        | 0.143      | 0.139      | 0.141       |

