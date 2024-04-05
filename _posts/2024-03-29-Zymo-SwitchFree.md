---
layout: post
title: SwitchFree library prep
date: '2024-03-29'
categories:
tags: [DNA, RNA, Library prep, Protocols]
projects: Kit test 
---

# SwitchFree library prep test

This post details the info about the SwitchFree library prep steps. Hollie got this kit as a test to see if we can produce adequate libraries for RNA sequencing with different species. I'm using the [Zymo-Seq SwitchFree 3' mRNA Library Kit](https://www.zymoresearch.com/products/zymo-seq-switchfree-3-mrna-library-kit). See Zymo's protocol [here](https://files.zymoresearch.com/protocols/_r3008_r3009__zymo_seq_switchfree_3_mrna_library_kit.pdf). 

We will be using samples 6, 7, 107, 108, 110, and 113 from the POC 2023 spawning project. These samples were extracted on [3/27/24](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-03-27-MiniprepPlus-DNA%3ARNA-extractions-Poc2023.md). The kit needs a minimum of 10 ng of total RNA or a maximum of 500 ng of total RNA, which is a large range. Since some of the concentrations in these samples were low, we will be using 11 ng of total RNA as input. Here's a breakdown of input RNA volumes for each sample: 


| TubeID | TS RNA (ng/uL) | Strip tube # | RNA (uL) | Ultrapure water (uL) | Total starting volume (ul) | Primer |
| ------ | -------------- | ------------ | -------- | -------------------- | -------------------------- | ------ |
| 6      | 9.52           | 1            | 1.2      | 3.8                  | 5.0                        | 1      |
| 7      | 11.5           | 2            | 1.0      | 4.0                  | 5.0                        | 2      |
| 107    | 4.65           | 3            | 2.4      | 2.6                  | 5.0                        | 3      |
| 108    | 2.76           | 4            | 4.0      | 1.0                  | 5.0                        | 4      |
| 110    | 5.23           | 5            | 2.1      | 2.9                  | 5.0                        | 5      |
| 113    | 3.5            | 6            | 3.1      | 1.9                  | 5.0                        | 6      |

Here's the SwitchFree library prep workflow: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/switchfree_lib_prep_workflow.png)

### Materials 

- [Zymo-Seq SwitchFree 3' mRNA Library Kit](https://www.zymoresearch.com/products/zymo-seq-switchfree-3-mrna-library-kit)
- PCR tubes 
- Thermocycler 
- Mini centrifuge
- Aluminum beads (to keep things on ice)
- Magnetic stand for PCR tubes 
- Kit contents 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/switchfree_lib_prep_contents.png)

### Buffer preperation 

Once buffers are prepared for a kit, they do not need to be prepared again. 

- Add 300 uL of the Select-a-Size MagBead concentrate to 10 mL of Select-a-Size MagBead Buffer. 
	- These should be prepared at least 5 days ahead of library prep. 
- Add 24 mL of 100% ethanol to the DNA Wash Buffer concentrate. Store at room temperature. 

### Best practices 

- Avoid multiple freeze thaws of all components, and aliquot components as necessary
- Remove enzymes from cold storage just before use and return to cold storage immediately after use
- Thaw and maintain all components on ice unless noted otherwise 
- Flick to mix thawed components and briefly centrifuge prior to use 
- After adding each component, mix by pipetting up and down 15-20 times. Briefly centrifuge after 
- Pre-program thermal cycler with lid heating ON set to >100-105°C 
- Turn on thermocycler day-of to preheat 
- Pre-calculate how much to dilute RNA samples 

### Protocol 

#### Section 1: cDNA Synthesis 

Thermocycler settings: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/switchfree_lib_prep_section1_thermocycler.png)

- Thaw PolyA R1 Reagent and R2 Reagent on ice 
- Prepare samples to correct volume in strip tubes 
- Add 5 uL of PolyA R1 Reagent to each sample 
- Mix by pipetting and centrifuge briefly 
- Run Steps 1-2 in thermocycler 
	- SKIP this step if RNA is degraded/low quality 
- While in the thermocycler on the 4°C hold, add 10 uL of R2 Reagent to each sample 
- Mix by pipetting 
- Run Steps 3-4
- While Steps 3-4 are running, allow Reaction Clean-up Solution to equilibrate to room temperature 
- While in the thermocycler on the 42°C hold, add 4 uL of Reaction Clean-up Solution to each sample 
- Mix by pipetting and centrifuge briefly 
- Immediately return to thermocycler 
- Run Steps 5-7
- Remove tubes from thermocycler 
- Add 26 uL of 95% ethanol to each sample 
- Do bead clean-up with Select-a-Size MagBeads 
	- Allow beads to equilibrate to room temperature >30 minutes before use 
	- Resuspend beads immediately before use by shaking or vortexing until homogenous 
	- Add 110 uL of Select-a-Size MagBeads to each sample 
	- Pipette to mix until homogenous
	- Incubate for 5 minutes at room temperature 
	- Put samples on magnetic stand for 5 mins or until beads have fully separated from solution 
	- Aspirate slowly and discard supernatent 
	- While the sample is still on the magnetic stand, add 200 uL of DNA Wash Buffer without disturbing the pellet 
	- Aspirate slowly and discard supernatent 
	- Repeat this wash step again 
	- Keep the caps open to air-dry bead for 1 minute 
	-  Aspirate any residual wash buffer
	-  Continue to air dry pellet until it appears matte without cracking (see picture below)
		- ![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/DT_mcap2023/ribofree_beads.png)
	- Remove sample from magnetic stand 
	- Add 10 uL of DNA Elution Buffer to each sample 
	- Pipette to mix until homogenous
	- Put samples on magnetic stand for 1-2 minutes or until elute is clear
	- Move elute to new PCR tubes 
	
#### Section 2: Adapter Ligation 

Thermocycler settings: 

Step 1: 95°C for 3 minutes, 4°C for 2 minutes, 4°C hold 

Step 2: 25°C for 45 minutes (lid temperature at 40°C), 4°C hold 

- Thaw 3'mRNA L Reagent on ice 
- Put samples in thermocycler and run Step 1 
- Remove samples from thermocycler and put on ice 
- Add 15 uL of 3'mRNA L Reagent to each sample 
- Mix by pipetting and centrifuge briefly 
- Run Step 2 
- Remove samples from thermocycler 
- Add 25 uL of DNase/RNase free water to each sample 
- Do bead clean-up with Select-a-Size MagBeads 
	- Allow beads to equilibrate to room temperature >30 minutes before use 
	- Resuspend beads immediately before use by shaking or vortexing until homogenous 
	- Add 100 uL of Select-a-Size MagBeads to each sample 
	- Pipette to mix until homogenous
	- Incubate for 5 minutes at room temperature 
	- Put samples on magnetic stand for 5 mins or until beads have fully separated from solution 
	- Aspirate slowly and discard supernatent 
	- While the sample is still on the magnetic stand, add 200 uL of DNA Wash Buffer without disturbing the pellet 
	- Aspirate slowly and discard supernatent 
	- Repeat this wash step again 
	- Keep the caps open to air-dry bead for 1 minute 
	-  Aspirate any residual wash buffer
	-  Continue to air dry pellet until it appears matte without cracking (see picture above)
	- Remove sample from magnetic stand 
	- Add 15 uL of DNA Elution Buffer to each sample 
	- Pipette to mix until homogenous
	- Put samples on magnetic stand for 1-2 minutes or until elute is clear
	- Move elute to new PCR tubes 

#### Section 3: Library amplification 

Thermocycler settings: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/switchfree_lib_prep_section2_thermocycler.png) 

Given the low input for these preps and my previous experience with other Zymo kits, we are going to do 21 PCR amplification cycles. 

- Thaw Index Primers and ZymoTaq Premix on ice 
- Allow Library Binding Solution to equilibrate to room temperature >30 minutes before use 
- Add 10 uL of unique Index Primer to each sample
- Mix by pipetting 
- Add 25 uL of Zymo TaqPremix to each sample 
- Mix by pipetting and centrifuge briefly 
- Run thermocycler program above 
- Remove samples from thermocycler 
- Add 50 uL of DNase/RNase free water to each sample 
- Add 80 uL of Select-a-Size MagBeads to each sample 
- Pipette until homogenous 
- Incubate for 5 minutes at room temperature 
- Place samples on magnetic stand for 5 minutes or until beads have separated from solution 
- Aspirate slowly and discard supernatent 
- Remove samples from magnetic stand 
- Add 100 uL of DNA Elution Buffer to beads 
- Pipette until homogenous 
- Add 80 uL of Library Binding Solution to each sample 
- Pipette until homogenous 
- Incubate for 5 minutes at room temperature 
- Do bead clean-up with Select-a-Size MagBeads (beads already added to samples)
	- Put samples on magnetic stand for 5 mins or until beads have fully separated from solution 
	- Aspirate slowly and discard supernatent 
	- While the sample is still on the magnetic stand, add 200 uL of DNA Wash Buffer without disturbing the pellet 
	- Aspirate slowly and discard supernatent 
	- Repeat this wash step again 
	- Keep the caps open to air-dry bead for 1 minute 
	-  Aspirate any residual wash buffer
	-  Continue to air dry pellet until it appears matte without cracking (see picture above)
	- Remove sample from magnetic stand 
	- Add 15-20 uL of DNA Elution Buffer to each sample 
	- Pipette to mix until homogenous
	- Put samples on magnetic stand for 1-2 minutes or until elute is clear
	- Move elute to new PCR tubes 

THIS IS THE FINAL 3' mRNA-seq LIBRARY. STORE AT -20°C. 

#### QC

Run [DNA Tapestation](https://github.com/meschedl/MESPutnam_Open_Lab_Notebook/blob/master/_posts/2019-07-30-DNA-Tapestation.md) to visualize libraries. Here's an example of what the library should look like on a Tapestation: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/switchfree_lib_prep_library_example.png)

I ran the tapestation on the libraries on 3/31/24. See the full report [here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/images/tapestation/DNA_switchfree_POC_2024-03-31.pdf). 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_overview_20240331.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_6_20240331.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_7_20240331.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_107_20240331.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_108_20240331.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_110_20240331.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_113_20240331.png)

Tapestation looks great!! Sample 108 looks a little wonky, perhaps due to a PCR bubble product (read more [here](https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000001918)). Libraries with a bubble product can still be sequenced, which is good. Need to discuss next steps with Hollie + Ariana. 