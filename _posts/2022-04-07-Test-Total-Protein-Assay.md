---
layout: post
title: Test Total Protein Assay w/ GSO Astrangia samples
date: '2022-04-07'
categories: Protocols, Physiology
tags: [Physiology, Total Protein]
projects: GSO Astrangia 
---

## Goal: Determine total protein in Astrangia samples from 2021 experiment 

Protocol from manufacturer [here](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/MAN0011430_Pierce_BCA_Protein_Asy_UG.pdf). Modified from this [protocol](https://github.com/urol-e5/protocols/blob/master/2020-01-01-Total-Protein-Protocol.md)

Ariana and I tested out the total protein assay this week. She used Mcap larval samples and I used my adult Astrangia samples.

*Important note*: Symbionts were NOT removed for my samples, reading total protein for the entire holobiont 

### Equipment/materials

- 96-well plates
- [Pierce™ BCA Protein Assay Kit](https://www.thermofisher.com/order/catalog/product/23225?SID=srch-srp-23225)
- Incubator with range up to 37°C
- Plate reader 
- Pipettes P10, P200, P1000 and tips
- P200 multi-channel pipette
- 1.5 and 50 mL tubes
- DI water

### Protocol

1. Take samples out of -80°C and thaw at room temperature 
2. Label 1.5 mL tubes for standards (9 tubes, A-I). If samples need to be diluted, label 1.5 mL tubes for samples as well.
	- add location in lab
3. Prepare protein standards in labeled 1.5 mL tubes, as described in the BCA Protein Assay protocol. The diluent is DI water 

| Vial | Volume of Diluent (μL) | Volume of Source of BSA (μL) | Final BSA Concentration (μg/mL) |
|------|------------------------|------------------------------|---------------------------------|
| A    | 0                      | 300 of Stock                 | 2000                            |
| B    | 125                    | 375 of Stock                 | 1500                            |
| C    | 325                    | 325 of Stock                 | 1000                            |
| D    | 175                    | 175 of vial B dilution       | 750                             |
| E    | 325                    | 325 of vial C dilution       | 500                             |
| F    | 325                    | 325 of vial E dilution       | 250                             |
| G    | 325                    | 325 of vial F dilution       | 125                             |
| H    | 400                    | 100 of vial G dilution       | 25                              |
| I    | 400                    | 0 (Blank)                    | 0                               |

4. Prepare working reagent 
	- Using the following formula, determine the total volume of WR reagent needed: (# standards + # samples) x (# replicates) x (volume of WR per sample) = total volume WR required
	- We will use 3 replicates and 200 uL of WR per sample 
	> Example calculation from a full plate: (9 standards + 22 samples) x (3 replicates) x (200 uL WR) = 18,600 uL or 18.6 mL. I like to round up a little for the WR so there is some extra when loading the plate. Here, 18,600 uL will be rounded up to 19000 uL.
	- Make WR by mixing Reagent A with Reagent B in a 50:1 ratio. 
	> Example calculation of volume needed for Reagent A and B: 19000 uL WR needed / 51 parts total of A and B = 372.549 uL. 372.549 uL is one part, therefore Reagent B is 372.549 uL. Using the ratio, calculate Reagent A by multiplying B by 50 parts: 372.549 uL x 50 = 18627.451. Double check the math by adding the volumes of Reagents A and B, which should equal ~19000. 372.549 uL Reagent B + 18627.451 Reagent A = 19000. Hooray! This is the amount needed to run a full plate. 
 5. Make a plate map in order to keep track of where the standards and samples go on the plate.  
 6. Load the plate by pipetting 25 uL of standard/sample into each well. 
 7. Pour the WR into a small trough. Using a multi-channel pipette, add 200 uL to each well. Pipette up and down 3 times to mix. 
 8. Cover the plate and incubate at 37°C for 30 minutes.
 9. Let covered plate cool on benchtop for 15 minutes. 
 10. Read on plate reader at 562 nm. 
 11. In R, subtract the average 562 nm absorbance measurement of the Blank standard replicates from the 562 nm measurements of all other individual standard and unknown sample replicates.
 12. Pepare a standard curve by plotting the average Blank-corrected 562nm measurement for each BSA standard vs. its concentration in μg/mL. Use the standard curve to determine the protein concentration of each unknown sample.

add code links















