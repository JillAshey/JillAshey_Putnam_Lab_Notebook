---
layout: post
title: Test Carbohydrate Assay w/ GSO Astrangia samples
date: '2022-04-06'
categories: Protocols, Physiology
tags: [Physiology, Carbohydrates]
projects: GSO Astrangia 
---

## Goal: Determine amount of carbohydrates in Astrangia samples from 2021 experiment 

See full protocol [here](https://github.com/Putnam-Lab/Lab_Management/tree/master/Lab_Resources/Physiology_Protocols/Carbohydrates/Bove_Baumann_96well_Protocol) from KW. Modified from this [protocol](https://www.protocols.io/view/coral-carbohydrate-assay-for-96-well-plates-j8nlk4ro1g5r/v1) from CB

Ariana and I tested out the carbohydrate assay this week. She used Mcap larval samples and I used a few of my adult Astrangia samples. 

Samples used: FLD-0031, AST-1109, AST-1138, AST-1152

*Important note*: Symbionts were NOT removed for these samples, reading carbs for the entire holobiont 

### Equipment/materials

- 96-well plates
- Water bath (room temperature)
- Vortex
- Fume hood
- Pipettes and tips (P1000 and P200)
- Plate reader (can read absorbance at 485 nm)
- Glucose standards 
- Sulfuric acid 
- Phenol 

### Protocol 

#### Preparing Standards 

In this test run, we used standards prepared by KW. Info on how to prepare the standards is in his [protocol](https://github.com/Putnam-Lab/Lab_Management/tree/master/Lab_Resources/Physiology_Protocols/Carbohydrates/Bove_Baumann_96well_Protocol)

#### Assay 

1. Take samples out of -80°C and thaw at room temperature 
2. Label 5 mL tubes (10 tubes for standards + # of samples in the run)
	- add location in lab
3. When samples are thawing, make the standards and blanks as shown below: 

| Tube ID | Concentration (mg/mL) | Vol water (uL) | Vol 1 mM Glucose (ul) | Vol 10mM Glucose (uL) |
|:-------:|:---------------------:|:--------------:|:---------------------:|:---------------------:|
|    B    |         0.0000        |      1000      |           0           |           0           |
|    1    |        0.00901        |       950      |           50          |           0           |
|    2    |        0.01802        |       900      |          100          |           0           |
|    3    |        0.02703        |       850      |          150          |           0           |
|    4    |        0.03604        |       800      |          200          |           0           |
|    5    |        0.05406        |       700      |          300          |           0           |
|    6    |         0.0901        |       500      |          500          |           0           |
|    7    |         0.1802        |        0       |          1000         |           0           |
|    8    |         0.3604        |       800      |           10          |          200          |
|    9    |         0.901         |       500      |           0           |          500          |

4. Vortex samples after thawed to mix well
5. Add 100 uL of coral homogenate and 900 uL DI water to pre-labeled tubes for all samples. 
6. Set up a room temperature water bath in the fume hood with test tube rack (i.e. DI water in a plastic bin)
7. In the fume hood, add 25 uL of phenol to first sample 
8. Add 2.5 mL sulfuric acid to sample
	- *Caution*: When the sulfuric acid is added, the tube will get very hot from the chemical reaction occurring. Hold the tube near the top to avoid heat.
9. Put tube in rack in the water bath to incubate
10. Repeat steps 7-9 for all samples/standards
11. When last sample has been placed in water bath, incubate all samples for 30 minutes. During this time, set up the plate layout for the 96-well plate
12. After 30 minutes, pipette 200 uL of sample/standard into 96-well plate following plate layout. Samples and standards will be run in triplicate.

Plate layout: 

13. Read on plate reader at 485 nm

#### Calculations 

1. Create standard curve with known standard concentrations and absorbance values (y = mx + b)
2. Using the resulting equation, convert sample absorbance to concentrations (mg/mL)
3. Multiply sample concentration (mg/mL) by total slurry volume (mL) and dilution factor (1000/v of sample, usually 100 mL), then divide by surface area (cm2) for resulting units: mg/cm2

add code here 














