---
layout: post
title: SwitchFree library prep
date: '2024-02-29'
categories:
tags: [DNA, RNA, Library prep, Protocols]
projects: Kit test 
---

# SwitchFree library prep test

This post details the info about the SwitchFree library prep steps. Hollie got this kit as a test to see if we can produce adequate libraries for RNA sequencing with different species. I'm using the [Zymo-Seq SwitchFree 3' mRNA Library Kit](https://www.zymoresearch.com/products/zymo-seq-switchfree-3-mrna-library-kit). See Zymo's protocol [here](https://files.zymoresearch.com/protocols/_r3008_r3009__zymo_seq_switchfree_3_mrna_library_kit.pdf). 

We will be using XXXXXX samples 

The kit needs a minimum of 10 ng of total RNA or a maximum of 500 ng of total RNA, which is a large range. Add info about starting input amt XXXXXXXXX

5 uL starting volume 

Here's the SwitchFree library prep workflow: 

XXXXXX add 

### Materials 

- [Zymo-Seq SwitchFree 3' mRNA Library Kit](https://www.zymoresearch.com/products/zymo-seq-switchfree-3-mrna-library-kit)
- PCR tubes 
- Thermocycler 
- Mini centrifuge
- Aluminum beads (to keep things on ice)
- Magnetic stand for PCR tubes 
- Kit contents 

XXXXXX

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

XXXXXXXXX

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

- 