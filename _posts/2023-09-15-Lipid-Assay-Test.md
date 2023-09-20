---
layout: post
title: Lipid assay
date: '2023-09-10'
categories: Protocol
tags: [Protocol, Lipids, Larvae, Mcap]
projects: Developmental time series 
---

### 20230915 Lipid assay test

This protocol was modified from [Bove & Baumann 2021](https://www.protocols.io/view/coral-lipid-assay-for-96-well-plates-q26g789pqlwz/v1) and L. Fuess lipid protocol. We are testing the lipid assay with coral eggs, embryos, and larvae. 

The samples used in this assay are test physiology samples from my developmental time series project with *M. capitata* in Hawaii 2023. For experimental details, see my Hawaii notebook posts [here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2023-06-15-Hawaii2023-DailyPosts.md). 

Lauren, Flo and I followed [this protocol](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2023-09-10-Lipid-Assay-Test.md). I am going to review the steps in this post and add more details/pictures. 

### Protocol 

**Remember** to do this protocol in fume hood whenever possible! Reagents are very toxic!

#### Make reagents according to the linked post above 

#### Lipid assay 

1. Prior to the assay, homogenize the eggs/embryos/larvae in 700 uL of PBS. 
	- Separate out 300 uL of holobiont fraction into a new tube. 
	- Centrifuge for 5 minutes at 3000 rpm to separate host and symbiont fraction.
	- Move the supernatent (host) into a new tube. 
	- Resuspend pellet (symbiont) with 700 uL of 1x PBS. 
2. Aliquot out 150 uL of sample extract and desiccate samples overnight in a SpeedVac. 

This is what the SpeedVac settings should be set at: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/desiccator_settings.JPG)

- Time - run infinitely or until it's turned off
- Brake - OFF or no spinning 
- Temp - room temp
- Mode vent - D-AQ or desiccate aqueous solution 

This is what the samples looked like after desicatting overnight: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/desiccated_samples.JPG)

3. Add 500 uL of a 2:1 chloroform:methanol mixture and 100 uL of 0.05 M NaCl to each desiccated sample.
	- 333.33 uL of chloroform and 166.67 uL of methanol 

This is what the samples looked like after the addition of those regents: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/lipid_samples_post_CHCl3.JPG)

4. Vortex samples for 15 seconds to dissolve lipids. 
5. Place samples on a shaker for 1 hour (set shaker to 500 rpm), vortexting the samples for 10 seconds every 15 minutes. 
6. Centrifuge samples at 3000 rpm for 5 minutes. 
7. Move 100 uL of the bottom layer of the centrifuged sample to a 96-well PCR plate in triplicate. 
	- We did this round in triplicate, as they were just test samples but I will probably run my actual samples in duplicate. 
8. Add 50 uL of methanol to each well. 
9. Put plate in hot bead bath (70°C) for 15 minutes. 
	- Bead bath is essentially the same as a water batch except it's aluminum beads instead of water. 
10. Following solvent evaporation, add 100 uL of 18 M sulfuric acid to each well. Add the sulfuric acid slowly, as it tends to bubble and splash as it comes into contact with the sample. If the sample splashes into another well and contiminates another sample, toss that plate and re-start. 

This is what the PCR plate looked like post-acid addition: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/lipid_samples_post_acid.JPG)

11. Put the PCR plate in a thermocycler at 90°C for 20 minutes, followed by 4°C for 20 minutes.
12. While the PCR plate is running on the thermocycler, make the lipid standards. We didn't have corn oil today, so we proceeded w/o the standards. This was fine, as we were just doing practice samples. 
13. Transfer 75 uL from the PCR well plate into a new 96 well plate (not a PCR plate). 
	- Holobiont samples were in A1, B1, C1, D1, E1, F1, G1, H1, and A2. 
	- Host samples were in A6, B6, C6, D6, E6, F6, G6, H6, and A7. 
	- Symbiont samples were in A11, B11, C11, D11, E11, F11, G11, H11, and A12. 

14. Run an initial absorbance reading on the plate reader at 540 nm. This is a baseline measurement that will be used for correcting the final plate absorbance. 

This was the initial readings: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/initial_lipid_reading_test.JPG)

15. Remove the plate from the plate reader and add 34.5 uLof 0.2 mg/mL vanillin in 17% phosphoric acid to each well. 
16. Incubate the plate in the dark at room temperature for 5 minutes. 

This is what the plate looked like after incubation: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/lipid_test_plate.JPG)

17. Read the plate again at 540 nm with the lid off. 

This was the final reading: 

![](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/images/final_lipid_reading_test.JPG?raw=true)

Values were extremely high! Some were even too high for the plate reader to read. This means that we need to use a lot less material next time. 

#### Next steps 
- Run a plate with just the corn oil as the standard curve to get an idea of what absorbance the standard curve is at. This will help me to figure out how much material to use for the next trial 
- Use a few different volumes (ie, 10, 15, 25, 50 uL, etc) to see which would be best for this assay. 
- We didn't melt the plate, hooray!


