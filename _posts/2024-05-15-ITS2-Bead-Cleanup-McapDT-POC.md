---
layout: post
title: ITS2 clean-up 
date: '2024-05-15'
categories:
tags: [DNA, ITS2, Extractions, Protocols]
projects: Mcap DT 2023; POC 2023 spawning 
---

This post details information on ITS2 clean-up for the my ambient Mcap developmental timeseries samples from 2023 and for the POC spawning experimental samples from 2023. Mcap 2023 github repo is [here](https://github.com/JillAshey/Hawaii_Developmental_TimeSeries) and POC github is [here](https://github.com/hputnam/Poc_RAPID). See Ariana's [post](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/ITS2-amplicon-PCR-and-preparation-for-sequencing-20240326/) about the ITS2 protocol for this project. 

### Samples 

Samples (44 total samples) were amplified on [5/13/24](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-05-13-ITS2-LowConcentration-McapDT-POC.md) and [5/14/24](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-05-14-ITS2-LowConcentration-McapDT-POC.md), but there appeared to be some [primer dimer](https://www.minipcr.com/primer-dimer-pcr/) present, as evidenced by bands ~100bp. Therefore, we decided to do a bead clean-up on all samples. 

### Equipment and materials 

- [Kapa Pure Beads](https://elabdoc-prod.roche.com/eLD/web/pi/en/products/SEQ-KAPA-0161?searchTerm=07983271001&catalog=Researcher&orderBy=Relevance)
- Tris-HCl (10 mM)
- Molecular grade ethanol
- Magnetic stand 
- PCR tubes 
- Gel rig + materials 

### Protocol 

Protocol was followed according to this [post](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/_posts/2024-04-01-KAPA-bead-clean-up-protocol-for-removal-of-primer-dimers-from-PCR-product.md). A few notes: 

- After adding Kapa Beads, I incubated tubes for 15 minutes at room temperature.
- After the ethanol washes, I dried beads for ~4 minutes until beads were matte/ethanol had evaporated. 
- After adding the elution buffer (10 mM Tris-HCl), I incubated tubes for 8 minutes at room temperature. 

### QC 

Following the clean-up protocol, I ran a 2% gel for 90 minutes at 100 volts + 100 amps to check if the primer dimers were removed. I also used a 1kb ladder and a 100 bp ladder. The machine that the gel rig is hooked up to has only been reaching ~70-80 volts. 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/gel_left_20240515.JPG)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/gel_right_20240515.JPG)

Gel looks great! No primer dimer from what I can see. Hollie will send to Janet to confirm if these are good enough for sequencing. 