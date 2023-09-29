---
layout: post
title: Lipid assay
date: '2023-09-26'
categories: Protocol
tags: [Protocol, Lipids, Larvae, Mcap]
projects: Developmental time series 
---

### 20230915 Lipid assay test

This protocol was modified from [Bove & Baumann 2021](https://www.protocols.io/view/coral-lipid-assay-for-96-well-plates-q26g789pqlwz/v1) and L. Fuess lipid protocol. We are testing the lipid assay with coral eggs, embryos, and larvae. 

The samples used in this assay are test physiology samples from my developmental time series project with *M. capitata* in Hawaii 2023. For experimental details, see my Hawaii notebook posts [here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2023-06-15-Hawaii2023-DailyPosts.md). 

Lauren and I followed [this protocol](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2023-09-10-Lipid-Assay-Test.md) with a few modifications. I am going to review the steps in this post and add more details/pictures. 

### Protocol 

**Remember** to do this protocol in fume hood whenever possible! Reagents are very toxic!

#### Make reagents according to the linked post above 

#### Lipid assay 

1. Prior to the assay, homogenize the eggs/embryos/larvae in 700 uL of PBS. 
	- Separate out 300 uL of holobiont fraction into a new tube. 
	- Centrifuge for 5 minutes at 3000 rpm to separate host and symbiont fraction.
	- Move the supernatent (host) into a new tube. 
	- Resuspend pellet (symbiont) with 300-400 uL of 1x PBS. 

2. In past attempts, we found that 150 uL of sample extract was too concentrated. We decided to do a dilution to see what volume of sample would be best to get a good signal. We did the following: 
- 10 uL of sample + 140 uL DI water 
- 25 uL of sample + 125 uL DI water 
- 50 uL of sample + 100 uL DI water 

3. Desiccate samples overnight in a SpeedVac. 

This is what the SpeedVac settings should be set at: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/desiccator_settings.JPG)

- Time - run infinitely or until it's turned off
- Brake - OFF or no spinning 
- Temp - room temp
- Mode vent - D-AQ or desiccate aqueous solution 

This is what the samples looked like after dessication: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/lipid_post_dessc_20230926.JPG)

4. Add 500 uL of a 2:1 chloroform:methanol mixture and 100 uL of 0.05 M NaCl to each desiccated sample.
	- 333.33 uL of chloroform and 166.67 uL of methanol 

5. Vortex samples for 15 seconds to dissolve lipids. 
6. Place samples on a shaker for 1 hour (set shaker to 500 rpm), vortexting the samples for 10 seconds every 15 minutes. 
7. While the samples are shaking, make the standards: 
	- Label 7 1.5 mL tubes 1-7. 
	- Add 300 μL of CHCl<sub>3</sub> to standard tubes 2 through 7
	- Add 600 μL of 1.5 mg/ml corn oil stock to standard tube 1
	- Transfer 300 μL of tube 1 and place in tube 2. Pipette mix
	- Pull 300 μL from tube 2 and place in tube 3. Pipette mix
	- Repeat this process for tubes 3 through 6. Do not add corn oil to tube 7! This is the blank
	- Discard 300 μL from tube 6 so total volume equals 300 μL (optional)

|  Standard | Concentration (mg/mL) |
|:---------:|:---------------------:|
|     1     |          1.5          |
|     2     |          0.75         |
|     3     |         0.375         |
|     4     |         0.188         |
|     5     |         0.094         |
|     6     |         0.047         |
| 7 (Blank) |           0           |


8. Centrifuge samples at 3000 rpm for 5 minutes. 
9. Move 100 uL of the standard or bottom layer of the centrifuged sample to a 96-well PCR plate in duplicate. 

This is the plate map: 

| x | 1     | 2     | 3             | 4             | 5            | 6            | 7             | 8             | 9            | 10           |
| - | ----- | ----- | ------------- | ------------- | ------------ | ------------ | ------------- | ------------- | ------------ | ------------ |
| A | Std 1 | Std 1 | LZ Holo 10 uL | LZ Holo 10 uL | LZ Sym 10 uL | LZ Sym 10 uL | JA Holo 10 uL | JA Holo 10 uL | JA Sym 10 uL | JA Sym 10 uL |
| B | Std 2 | Std 2 | LZ Holo 25 uL | LZ Holo 25 uL | LZ Sym 25 uL | LZ Sym 25 uL | JA Holo 25 uL | JA Holo 25 uL | JA Sym 25 uL | JA Sym 25 uL |
| C | Std 3 | Std 3 | LZ Holo 50 uL | LZ Holo 50 uL | LZ Sym 50 uL | LZ Sym 50 uL | JA Holo 50 uL | JA Holo 50 uL | JA Sym 50 uL | JA Sym 50 uL |
| D | Std 4 | Std 4 | LZ Host 10 uL | LZ Host 10 uL |              |              | JA Host 10 uL | JA Host 10 uL |              |              |
| E | Std 5 | Std 5 | LZ Host 25 uL | LZ Host 25 uL |              |              | JA Host 25 uL | JA Host 25 uL |              |              |
| F | Std 6 | Std 6 | LZ Host 50 uL | LZ Host 50 uL |              |              | JA Host 50 uL | JA Host 50 uL |              |              |
| G | Std 7 | Std 7 |               |               |              |              |               |               |              |

10. Add 50 uL of methanol to each well. 
11. Put plate in hot bead bath (70°C) for 15 minutes. 
	- Bead bath is essentially the same as a water batch except it's aluminum beads instead of water. 
12. Following solvent evaporation, add 100 uL of 18 M sulfuric acid to each well. Add the sulfuric acid slowly, as it tends to bubble and splash as it comes into contact with the sample. If the sample splashes into another well and contiminates another sample, toss that plate and re-start. 

This is what the PCR plate looked like post-acid addition: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/lipid_post_acid_20230926.JPG)

13. Put the PCR plate in a thermocycler at 90°C for 20 minutes, followed by 4°C for 20 minutes.
14. Transfer 75 uL from the PCR well plate into a new 96 well plate (not a PCR plate). 
15. Run an initial absorbance reading on the plate reader at 540 nm. This is a baseline measurement that will be used for correcting the final plate absorbance. 

This is what the lipid plate looked like before the initial reading: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/initial_lipid_20230926.JPG)

These are the values from the initial reading:

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/initial_reading_lipid_20230926.JPG)

16. Remove the plate from the plate reader and add 34.5 uLof 0.2 mg/mL vanillin in 17% phosphoric acid to each well. 
17. Incubate the plate in the dark at room temperature for 5 minutes. 
18. Read the plate again at 540 nm with the lid off. 

This is what the lipid plate looked like after the incubation and before the final reading: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/final_lipid_20230926.JPG)

These are the values from the final reading:

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/final_reading_lipid_20230926.JPG)

Similarly to previous [attempts](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2023-09-15-Lipid-Assay-Test.md), the values were really high, but only for the standards. 

We used a 2:1 chloroform:methanol mixture to extract lipids following Cheng et al. 2011 and Bove & Baumann. This is what Cheng et al. 2011 plate looked like following the vanillin/phosphoric acid addition: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/cheng_et_al_2011_fig2.png)

Ours looked darker, more similar to the C:M 1:1 from their paper. They also used a higher vanillin and phosphoric acid ratio (they used 0.5 mg vanillin per mL in 68% phosphoric acid, while we used 0.2 mg vanillin per mL in 17% phosphoric acid). 

### Next steps 
- Check expiration dates on all reagents to make sure nothing is expired
- Discuss at lab meeting
- Maybe remake the vanillin in phosphoric acid? 