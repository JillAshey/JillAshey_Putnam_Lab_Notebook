---
layout: post
title: RNA clean and concentrate protocol 
date: '2023-03-02'
categories: Protocols, RNA
tags: [Protocols, RNA]
projects: GSO Astrangia 
---

This protocol describes how to clean/concentrate RNA with the [Zymo RNA Clean & Concentrator-5 kit](https://www.zymoresearch.com/products/rna-clean-concentrator-5). 

- [Astrangia 2021 master extraction sheet](https://docs.google.com/spreadsheets/d/1inGyWAlwnnSJXvDR5anFCuhebnuQkrmizTKb37Thg_0/edit#gid=1846471493)
- [Clean + Concentrate test sheet](https://docs.google.com/spreadsheets/d/1kxe8YRLCuoOeKmdYfhXtNh23e5-N2kOQmyhU63M7pk0/edit#gid=785181220) (JA & ZD samples)

### Materials & Equipment 

- Zymo RNA Clean & Concentrator -5 kit
	- From kit:
		- RNA binding buffer 
		- RNA prep buffer 
		- RNA wash buffer 
		- DNase/RNase free water
		- Spin IC columns and collection tubes 
- Centrifuge and rotor capable of spinning at 15,000 rcf
- Plastics 
	- 2 1.5 mL microcentrifuge tubes per sample
	- 1 quibit tube per sample 
	- Several 1000 µL, 200 µL, and 20 µL filtered pipette tips per sample

### Samples 

Today, I tested out the RNA clean/concentrate kit on 4 of my samples with low RNA concentrations. I need 50 ng/µL to send to Azenta for RNA sequencing (small, long, mRNA). 

These were the samples I used: 

| Sample   | Avg RNA ng/µL | µg    | Extraction date |
| -------- | ------------- | ----- | --------------- |
| AST-2523 | 20.1          | 1.809 | 11/7/22         |
| AST-1065 | 16.3          | 1.467 | 11/7/22         |
| AST-1761 | 45.8          | 4.122 | 1/12/23         |
| AST-2403 | 30.6          | 2.754 | 1/16/23         |

Each sample had ~90 µL of RNA in the tube. 

### Protocol 

- Add 2 volumes of RNA binding buffer to each sample and pipette up and down to mix 
	- 40 µL of sample; 80 µL of buffer
- Add equal volume of 100% ethanol and pipette up and down to mix 
	- 120 µL of ethanol 
- Transfer sample to spin column and spin at 15,000g for 60 seconds. Discard flow-through 
- Add 400 µL of RNA prep bugger and spin at 15,000g for 60 seconds. Discard flow-through 
- Add 700 µL of RNA wash bugger and spin at 15,000g for 60 seconds. Discard flow-through 
- Add 400 µL of RNA wash bugger and spin at 15,000g for 60 seconds. Discard flow-through 
- Move spin column to 1.5 mL tube
- Add 15 µL of DNase/RNase free water to column and spin at 15,000g for 60 seconds - this is the cleaned & concentrated RNA!

### Quality control 

#### Quibit results 

| Tube     | Type     | RNA 1 (ng/uL) | RNA 2 (ng/uL) | Avg (ng/uL) | PRIOR Avg (ng/uL) |
| -------- | -------- | ------------- | ------------- | ----------- | ----------------- |
| S1       | Standard | 409.79        | NA            |             |                   |
| S2       | Standard | 8736.12       | NA            |             |                   |
| AST-2523 | Sample   | 32            | 31.8          | 31.9        | 20.1              |
| AST-1065 | Sample   | 28.4          | 28            | 28.2        | 16.3              |
| AST-1761 | Sample   | 65.4          | 65.8          | 65.6        | 45.8              |
| AST-2403 | Sample   | 44.2          | 44            | 44.1        | 30.6              |

Here's a table that Zoe & I put together that outlines our expected v. observed results. 

| A) SampleID | B) AvgConcentration (ng/uL) | C) ng  | D) Take 40 uL from each | E) Expected concentration in 15 uL | F) Avg concentration post-kit | G) ng post-kit | % RNA recovered |
| ----------- | --------------------------- | ------ | ----------------------- | ---------------------------------- | ----------------------------- | -------------- | --------------- |
| AST-2523    | 20.1                        | 2010   | 804                     | 53.6                               | 31.9                          | 478.5          | 0.5951493       |
| AST-1065    | 16.3                        | 1630   | 652                     | 43.4666667                         | 28.2                          | 423            | 0.648773        |
| AST-1761    | 45.8                        | 4580   | 1832                    | 122.133333                         | 65.6                          | 984            | 0.5371179       |
| AST-2403    | 30.6                        | 3060   | 1224                    | 81.6                               | 44.1                          | 661.5          | 0.5404412       |
|             |                             |        |                         |                                    |                               |                |                 |
| *Calculations*             |                             | =B\*100 | =C\*0.4                  | =D/15                               |                               | =H\*15          | =G/D             |