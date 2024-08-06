---
layout: post
title: MiniPrep Plus DNA/RNA extractions
date: '2024-08-05'
categories:
tags: [DNA, RNA, Extractions, Protocols]
projects: Mcap developmental time series 
---

# Extractions for e5 time series 

This protocol is based on the Putnam lab Zymo Miniprep protocol (see example from Zoe [here](https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/Protocols_Zymo_Quick_DNA_RNA_Miniprep_Plus/)). This post details the info about the extraction steps for 12 samples from the e5 timeseries (Acropora, Pocillopora, Porites spp.) in Moorea 2020. The e5 github is linked [here](https://github.com/urol-e5). 

### Samples 

The samples run today were clipped from snap frozen fragments by Zoe and I on [8/2/24](https://github.com/zdellaert/ZD_Putnam_Lab_Notebook/blob/master/_posts/2024-08-01-E5-Time-Series-Reextractions.md). These were all samples that either needed to be re-extracted because the initial extraction did not have enough RNA for sequencing (10 samples) or had never been extracted in the first place (2 samples). Kristen Terpis did the vast majority of the original (see her posts [here](https://github.com/Kterpis/Putnam_Lab_Notebook/tree/master/_posts)). Her extraction information is in a [spreadsheet](https://docs.google.com/spreadsheets/d/1A764av1a3VORX6m9aDUEcoY9Bx9l0fvGtV5ycm2J9Wo/edit?gid=0#gid=0) but we have consolidated the samples into those we want to [sequence](https://docs.google.com/spreadsheets/d/1iFsVfp1vix9IfNcqSajQfLgXXDg0CCPT/edit?gid=1681650319#gid=1681650319). Samples were in 1mL of DNA/RNA shield. 

- POC-53, TP4 (Nov) - brand new extraction
- POR-262, TP2 (Mar) - brand new extraction
- ACR-173, TP4 (Nov)
- ACR-186, TP1 (Jan)
- POC-219, TP3 (Sep)
- POR-245, TP4 (Nov)
- POR-74, TP2 (Mar)
- POR-83, TP1 (Jan)
- ACR-150, TP1 (Jan)
- ACR-145, TP1 (Jan)
- ACR-173, TP1 (Jan)
- POC-52, TP1 (Jan)

This photo is post bead-beating: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/samples_20240805.JPG)

Kristina typically labelled her tubes with just a number but I labelled tubes with the colony ID and timepoint.  

### Materials 

- Zymo Quick-DNA/RNA Miniprep Plus Kit [HERE](https://files.zymoresearch.com/protocols/_d7003t_d7003_quick-dna-rna_miniprep_plus_kit.pdf) Protocol Booklet
- Tris-Ethylenediaminetetraacetic acid (EDTA) 1X buffer for DNA elution
- Heating block capable of heating to 70ºC
- Centrifuge and rotor capable of spinning at 15,000 rcf
- Plastics 
	- 5 1.5 mL microcentrifuge tubes per sample
	- 2 PCR tubes per sample
	- 2 Qubit tubes per sample 

### Protocol 

- Samples were snap frozen (in 2020) and clipped into a tube with 1 mL of DNA/RNA shield and 0.5mm glass beads. 
- Samples were bead beat for 1.5 minutes on max speed. 

After bead beating, this [protocol](https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/Protocols_Zymo_Quick_DNA_RNA_Miniprep_Plus/) was followed exactly and eluted in 100 uL of either Tris (DNA) or DNA/RNA free water (RNA). 

### QC 

#### Qubit results 

For Qubit, I used the DNA BR kit and the RNA HS kit. 

| Colony ID | Timepoint | Volume | DNA1_ng_uL | DNA2_ng_uL | RNA1_ng_uL | RNA2_ng_uL | DNA_average | RNA_average |
| --------- | --------- | ------ | ---------- | ---------- | ---------- | ---------- | ----------- | ----------- |
| POC-53    | TP4       | 90     | 28.2       | 27.8       | 41         | 40.6       | 28          | 40.8        |
| POR-262   | TP2       | 90     | 62.6       | 61.4       | 36.6       | 36         | 62          | 36.3        |
| ACR-173   | TP4       | 90     | 16.1       | 15.7       | 12.1       | 12.1       | 15.9        | 12.1        |
| ACR-186   | TP1       | 90     | 14.1       | 13.9       | 35         | 35         | 14          | 35          |
| POC-219   | TP3       | 90     | 3.64       | 3.58       | 26.8       | 26.8       | 3.61        | 26.8        |
| POR-245   | TP4       | 90     | 35.4       | 35         | 28.6       | 28         | 35.2        | 28.3        |
| POR-74    | TP2       | 90     | 28.2       | 27.8       | 9.62       | 9.48       | 28          | 9.55        |
| POR-83    | TP1       | 90     | 39.8       | 39.4       | 11.7       | 11.4       | 39.6        | 11.55       |
| ACR-150   | TP1       | 90     | 23.4       | 23         | 12.9       | 13         | 23.2        | 12.95       |
| ACR-145   | TP1       | 90     | 46.8       | 46.4       | 39.6       | 39.4       | 46.6        | 39.5        |
| ACR-173   | TP1       | 90     | 52.8       | 52         | 15.3       | 15.3       | 52.4        | 15.3        |
| POC-52    | TP1       | 90     | NA         | NA         | 11         | 10.5       | NA          | 10.75       |

Qubit looks great! All samples now have enough RNA for both RNA-seq and small RNA-seq. 

#### Gel 

I ran a 1% gel at 80 volts for 1 hour. 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/gel_20240805.JPG)

The gel doesn't really look like anything...tried to look at it a few different ways. Maybe I didn't add the dye properly? Not sure, but Zoe is going to re-run the gel for these samples tomorrow. 

Zoe reran the gel and got the same thing (will add gel picture). 