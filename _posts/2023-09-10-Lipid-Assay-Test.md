---
layout: post
title: Lipid assay
date: '2023-09-10'
categories: Protocol
tags: [Protocol, Lipids, Larvae, Mcap]
projects: Developmental time series 
---

## Test of lipid assay protocol 

This protocol was modified from [Bove & Baumann 2021](https://www.protocols.io/view/coral-lipid-assay-for-96-well-plates-q26g789pqlwz/v1) and L. Fuess lipid protocol. We are testing the lipid assay with coral eggs, embryos, and larvae. 

The samples used in this assay are test physiology samples from my developmental time series project with *M. capitata* in Hawaii 2023. For experimental details, see my Hawaii notebook posts [here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2023-06-15-Hawaii2023-DailyPosts.md). 

### Equipment 

- SpeedVac
- Vortex
- Shaker 
- Centrifuge
- 96-well PCR plate 
- Water bath (up to 70°C)
- Thermocycler 
- 96-well plate 
- Plate reader 
- Pipettes + tips 
- 1.5 mL tubes 
- Fume hood

### Reagents 

- CH<sub>3</sub>OH (Methanol)
- CHCl<sub>3</sub> (Chloroform)
- 0.05 M NaCl in water 
- H<sub>3</sub>PO<sub>4</sub> (17% phosphoric acid)
- 0.2 mg/mL vanillin in 17% phosphoric acid 
- Concentrated (18M) sulfuric acid (H<sub>2</sub>SO4)
- 1.5 mg/mL corn oil 

### Protocol 

#### Make reagents 

Prepare reagents according to Putnam lab lipid [protocol](https://github.com/Putnam-Lab/Lab_Management/blob/master/Lab_Resources/Physiology_Protocols/Lipids/Bove_Baumann_96well_Protocol/Coral_Lipid_Protocol.md). 

##### 0.05 NaCl 

In a 50 mL labeled falcon tube, add 0.1461g of NaCl in 50 mL DI water. 

##### Stock 1.5 mg/mL corn oil 

In a 15 mL labeled falcon tube, add 245 uL of corn oil in 14.755 mL of CHCl<sub>3</sub> (chloroform)

Calculations:

- Density = 0.9188 g/ml (Noureddini et al., 1992)
- Total volume of ampule = 1 g* (1ml / 0.9188g) =1.08837614 ml
- Known concentration of ampule= 1000mg / 1.088ml = 918.8 mg/ml (same as known density)
- Therefore: (1.5 mg/mL * 15mL) / 918mg/mL = 0.0245 mL (or 245 μL) of 918mg/mL Corn oil standard concentrate

##### Stock 0.2 mg/mL vanillin in 17% phosphoric acid (H3PO4):

- 20 mL of 17% H<sub>3</sub>PO<sub>4</sub>
	- In a labeled 50 mL Falcon tube, add 4 mL 85% H<sub>3</sub>PO<sub>4</sub> to 16 mL of DI water
- Stock vanillin solution
	- In a labeled 50 mL Falcon tube, add 4 mg vanillin to 20 mL 17% H<sub>3</sub>PO<sub>4</sub>

##### Make subsets of the following reagents in labeled 50 mL falcon tubes daily

- CH<sub>3</sub>OH (Methanol)
- CHCl<sub>3</sub> (Chloroform)
- Concentrated (18M) sulfuric acid (H<sub>2</sub>SO4)

#### Lipid assay 

1. Prior to the assay, homogenize the eggs/embryos/larvae in 700 uL of PBS. 
2. Aliquot out 150 uL of sample extract and desiccate samples overnight in a SpeedVac. 
3. Add 500 uL of a 2:1 chloroform:methanol mixture and 100 uL of 0.05 M NaCl to each desiccated sample.
	- 333.33 uL of chloroform and 166.67 uL of methanol 
4. Vortex samples to dissolve lipids. 
5. Place samples on a shaker for 1 hour, vortexting the samples every 15 minutes. 
6. Centrifuge samples at 3000 rpm for 5 minutes. 
7. Move 100 uL of the bottom layer of the centrifuged sample to a 96-well PCR plate in triplicate
	- Decide if we will run things in duplicate or triplicate 
8. Add 50 uL of methanol to each well. 
9. Put plate in hot water bath (70°C) for 15 minutes. 
10. Following solvent evaporation, add 100 uL of 18 M sulfuric acid to each well. 
11. Put the PCR plate in a thermocycler at 90°C for 20 minutes, followed by 4°C for 20 minutes.
12. While the PCR plate is running on the thermocycler, make the lipid standards: 

	a. Make a stock serial dulution in 7 1.5 mL tubes for each plate.
	
	b. Add 300 uL of chloroform to standard tubes 2-7.
	
	c. Add 600 uL of 1.5 mg/mL stock corn oil to standard tube 1. 
	
	d. Transfer 300 uL from tube 1 into tube 2. Pipette up and down to mix. 
	
	e. Transfer 300 uL from tube 2 into tube 3. Pipette up and down to mix. 
	
	f. Repeat this process for tubes 3 through 6. Do not add corn oil to tube 7! This is the blank. 
	
	g. Discard 300 uL from tube 6 so total volume equals 300 uL (optional)

	| Standard  | Corn oil concentration (mg/mL) |
	| --------- | ------------------------------ |
	| 1         | 1.5                            |
	| 2         | 0.75                           |
	| 3         | 0.375                          |
	| 4         | 0.188                          |
	| 5         | 0.094                          |
	| 6         | 0.047                          |
	| 7 (blank) | 0                              |

13. Transfer 75 uL from the PCR well plate into a new 96 well plate. 
14. Run an initial absorbance reading on the plate reader at 540 nm. This is a baseline measurement that will be used for correcting the final plate absorbance. 
	- Cover plate before reading?
15. Remove the plate from the plate reader and add 34.5 uLof 0.2 mg/mL vanillin in 17% phosphoric acid to each well. 
16. Incubate the plate in the dark at room temperature for 5 minutes. 
17. Read the plate again at 540 nm. 
	- Cover plate before reading?

#### Calculations

NEED TO FIGURE THIS OUT AND HOW TO STANDARDIZE IT PER LARVAE

#### Important notes 

- Do as much of this protocol as possible in the fume hood. Most of the reagents are extremely toxic. 
- Dispose of waste in labeled bags as hazardous waste and [schedule a hazardous waste pick up](https://web.uri.edu/ehs/online-pickup/) with URI EHS. 
- Dr. Kevin Wong attempted this lipid assay several times at URI. See his posts [here](https://kevinhwong1.github.io/KevinHWong_Notebook/) and search for lipid.