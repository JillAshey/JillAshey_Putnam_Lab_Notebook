---
layout: post
title: miRNA library prep
date: '2023-10-27'
categories:
tags: [DNA, RNA, Library prep, Protocols]
projects: Mcap developmental time series 
---

# miRNA library prep for Mcap developmental time series Hawaii 2023

This post details the info about the library prep steps for the *M. capitata* developmental time series experiment in Hawaii 2023. The github for that project is linked [here](https://github.com/JillAshey/Hawaii_Developmental_TimeSeries).  

I'm using the [Zymo-Seq miRNA library kit](https://www.zymoresearch.com/products/zymo-seq-mirna-library-kit) for small RNA library prep (I got the 12 prep kit). Because there are only 12 preps in the kit, I want to stretch it to see if I can get more preps from the kit. **Therefore, in my protocol, I will be halving the volumes for each step, including the input.** See Zymo's protocol [here](https://files.zymoresearch.com/protocols/r3006_r3007-zymo-seq_mirna_library_kit.pdf). 

The kit needs a minimum input of 1 ng of total RNA or a maximum input of 200 ng of total RNA. 

Here's the miRNA library prep workflow: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/miRNA_lib_prep_workflow.png)

### Materials 

- Zymo-Seq miRNA library kit, which contains: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/miRNA_lib_prep_contents.png)

- PCR tubes 
- Thermocycler 
- Heating block 
- Mini centrifuge
- Bench top centrifuge  
- Aluminum beads (to keep things on ice)
- Magnetic stand for PCR tubes 

### Best practices 

- Avoid multiple freeze thaws of all components, and aliquot components as necessary
- Remove enzymes from cold storage just before use and return to cold storage immediately after use
- Thaw and maintain all components on ice unless noted otherwise 
- Flick to mix thawed components and briefly centrifuge prior to use 
- After adding each component, mix by pipetting up and down 15-20 times. Briefly centrifuge after 
- Pre-program thermal cycler with lid heating ON set to >95°C 
- Turn on thermocycler day-of to preheat 
- Pre-calculate how much to dilute RNA samples 

### Protocol 

#### Summary of thermocycler programs 

Create these programs beforehand with the lid heading ON and set to >95°C

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/miRNA_lib_prep_thermocycler.png)

#### Section 1: Adapter ligation and blocking 

- Turn on thermocycler and make sure programs are set / cycler is preheated
- Thaw RNA samples on ice 
- Thaw components with orange caps and red caps on ice
	- Orange caps (adapter ligation):
		- RNase Inhibitor
		- miRNA Adapter 
		- Ligation Enhancer -- THAW THIS AT ROOM TEMPERATURE 
		- miRNA Ligation Buffer 
		- miRNA Ligase 
		- DNase/RNase free water 
	- Red caps (adapter blocking):
		- Blocking Agent 
		- Blocking Enzyme 
- Once components have thawed, flick to mix and briefly centrifuge 
- Prior to use, pre-heat Ligation Enhancer at 37°C for 5 minutes to improve ease of pipetting. Place at room temperature after pre-heating and while in use. 
- Adjust each RNA sample volume to 3.25 with DNase/RNase free water in a 0.2 mL PCR tube. 
- Add 2 uL of Ligation to each sample slowly. 
- After adding Ligation Enhancer to each sample, set pipette to 10 uL and mix Ligation Enhancer by pipetting up and down slowly 15-20 times. Centrifuge and put on ice 
- Place samples in pre-heated thermocycler and run Steps 1-2 (Adapter Prep)
- Remove samples and put on ice
- Add 0.5 uL of RNase Inhibitor to each sample on ice 
- Add 0.75 uL of miRNA ligation buffer to each sample on ice 
- Add 0.5 uL of miRNA Ligase to each sample on ice 
- Pipette up and down 15-20 times and briefly centrifuge 
- Place samples in pre-heated thermocycler and run Steps 3-5 (Adapter Ligation)
- Once thermocycler reaches the step 5 hold, **do not remove** samples from thermocycler
- While still in thermocycler, add 0.5 uL of Blocking Agent to each sample. Mix with pipette (centrifuge if needed)
- Run Steps 6-8 (Blocking Agent) in thermocycler 
- Once thermocycler reaches the step 8 hold, **do not remove** samples from thermocycler
- While still in thermocycler, add 0.5 uL of Blocking Enzyme to each sample. Mix with pipette (centrifuge if needed)
- Run Steps 9-11 (Blocking Enzyme) in thermocycler 
- Remove samples from thermocycler and put on ice 

**Safe stopping point: samples can be stored at < -20°C overnight**

#### Section 2: Circularization and Dimer Removal 

- Turn on thermocycler and make sure programs are set / cycler is preheated
- Equilibrate Dimer Removal Beads (clear cap) to room temp for 30 mins (don't place on ice)
- Thaw components with yellow caps and the Dimer Removal Agent (clear cap) on ice. 
	- Yellow caps 
		- Circularization Enzyme 
		- Circularization Buffer 
- Once components have thawed, flick to mix and briefly centrifuge 
- For each sample, prepare a master mix of the Circlarization components as described below:
	
	XXXXXXX ADD 
	
- Add 1 uL of the master mix to each sample on ice. Pipette to mix and briefly centrifuge 
- Place samples in pre-heated thermocycler and run Steps 1-2 (Circularization)
- Once thermocycler reaches the step 2 hold, **do not remove** samples from thermocycler
- While still in thermocycler, add 0.5 uL of Dimer Removal Agent to each sample. Mix with pipette (centrifuge if needed)
- Run Steps 3-4 (Dimer Anneal)
- Prepare the Dimer Removal Beads 
	- Vortex the Dimer Remval Beads for 30 seconds or until homogenous 
	- For each sample, aliquot 20 uL of Dimer Remval Beads into new PCR tube 
	- Just before Step 3 is done on the thermocycler, place PCR tubes containing Dimer Remval Beads on a magnetic stand at room temp for 1 min or until supernatent is clear 
	- Remove and discard cleared supernatent 
	- Move the PCR tubes containing prepared Dimer Remval Beads from the magnetic stand to a PCR tube rack at room temp
- Once Step 4 has been reached, remove samples from thermocycler and briefly centrifuge 
- Put samples on PCR rack at room temperature 
- Transfer the entire sample reaction to a tube of prepared Dimer Remval Beads and pipette until homogenous 
- Immediately return samples to thermocycler 
- Run Steps 5-6 (Dimer Removal)
- Once Step 6 has been reached, remove samples from thermocycler and centrifuge for 10 seconds 
- Place samples on magnetic rack for 1 min or until supernatent is clear 
- Transfer cleared supernatent to a new PCR tube and place on ice

#### Section 3: Reverse Transcription

- Turn on thermocycler and make sure programs are set / cycler is preheated
- Thaw the RNase Inhibitor (orange) and all components with brown caps on ice. 
	- Brown caps
		- RT Primer Mix 
		- RT Buffer 
		- RT Enzyme 
- Once components have thawed, flick to mix and briefly centrifuge 
- Add 2 uL of RT Primer Mix to each sample on ice 
- Pipette to mix and briefly centrifuge 
- Put samples in pre-heated thermocycler and run Steps 1-2 (RT Primer Prep). 
- Once Step 2 has been reached, remove samples from thermocycler and put on ice 
- For each sample, prepare a master mix of the Reverese Transcription components as described below

	XXXXXXX ADD 

- Add 2.25 uL of the master mix to each sample 
- Pipette to mix and briefly centrifuge 
- Place samples in thermo cycler and run Steps 3-5 (Reverse Transcription)

**Safe stopping point: samples can be stored at < -20°C overnight**

#### Section 4: Index PCR 

- Turn on thermocycler and make sure programs are set / cycler is preheated
- Thaw the following components with green caps on ice: 
	- Zymo Taq Premix 
	- miRNA UDI Primers (individual tubes) 
- Once components have thawed, flick to mix and briefly centrifuge 
- Add 3 uL of a unique miRNA UDI Primer to each sample 

NEED TO LOOK MORE INTO THIS TO ENSURE THAT I AM NOT DUPLICATING ON SAMPLES. May need to sequence on different lanes 

- Add 17.5 of Zymo Taq Premix to each sample 
- Pipette to mix and briefly centrifuge 
- Place samples in pre-heated thermocycler and run the Index PCR Program 

NEED TO FIGURE OUT HOW MANY CYCLES 

- Remove samples from thermocycler and place on ice

**Safe stopping point: samples can be stored at < -20°C overnight**

#### Section 5: Library Purification 

NEED TO FIGURE OUT IF I SHOULD HALVE THESE VOLUMES AS WELL PRIOR TO FINAL ELUTION

- Add 300 uL of Select-a-Size Magbead Concentrate to 10 mL Select-a-Size Magbead Buffer for first time usage 
- Allow Select-a-Size Magbeads to equilibrate to room temperature for at least 30 mins prior to use 
- Add 24 mL of 100% EtOH to the DNA Wash Buffer Concentrate for first time usage 
- Just before use, vortex Select-a-Size Magbeads until homogenous 
- Add 62.5 uL of Select-a-Size Magbeads to each sample 
- Pipette to mix until homogenous 
- Incubate for 5 mins at room temperature 
- Place the samples on a magnetic stand until the beads have fully separated from the solution, leaving the solution clear 
- Without disturbing the beads, remove and discard the cleared supernatent 
- Leave the samples on the magnetic stand
- Add 100 uL of DNA Wash Buffer to each sample
- Without disturbing the beads, remove and discard the supernatent
- REPEAT this step for a total of 2 washes 
- Keep the samples on the magnetic stand with the tube caps open to air dry the beads 
- After 1 min, remove any residual DNA Wash Buffer 
- Air dry until the bead pellet is no longer glossy (5-10 mins)
	- The time will depend on the humidity and temperature in the room. Optimally dried beads appear matte without cracking 
	- Start with 3 mins of air-dry time and increase the time as needed 
- Remove samples from magnetic stand and resuspend beads with 12.5 uL of DNase/RNase-Free Water. Pipette to mix 
	- 12.5 uL is NOT halved. Need to decide if I should halve it 
- Put the samples on the magnetic stand for 3 mins or until the solution is clear 
- Move 10 uL of the eluate to a new tube, avoiding bead carryover 
- This is the final miRNA library!!!!! Store at < -20°C

#### Library Validation and Quantification 

Libraries can be visualized using a gel or Tapestation and quantified using Qubit or Tapestation. 