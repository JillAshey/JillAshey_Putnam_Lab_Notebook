---
layout: post
title: miRNA library prep
date: '2023-10-27'
categories:
tags: [DNA, RNA, Library prep, Protocols]
projects: Mcap developmental time series 
---

# miRNA library prep for Mcap developmental time series Hawaii 2023

This post details the info about the miRNA library prep steps for the *M. capitata* developmental time series experiment in Hawaii 2023. The github for that project is linked [here](https://github.com/JillAshey/Hawaii_Developmental_TimeSeries). I'm using the [Zymo-Seq miRNA library kit](https://www.zymoresearch.com/products/zymo-seq-mirna-library-kit) for small RNA library prep (I got the 12 prep kit). See Zymo's protocol [here](https://files.zymoresearch.com/protocols/r3006_r3007-zymo-seq_mirna_library_kit.pdf). 

The kit needs a minimum input of 1 ng of total RNA or a maximum input of 200 ng of total RNA in 6.5 uL. In my Ribofree library preps, I used ~190 ng of RNA from each sample. For this kit, I have no starting point with which to base my input amount. I will likely start with 100 ng as a happy middle ground. 

Here's the miRNA library prep workflow: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/miRNA_lib_prep_workflow.png)

### Materials 

- [Zymo-Seq miRNA library kit](https://www.zymoresearch.com/products/zymo-seq-mirna-library-kit)
- PCR tubes 
- Thermocycler 
- Heating block 
- Mini centrifuge
- Bench top centrifuge  
- Aluminum beads (to keep things on ice)
- Magnetic stand for PCR tubes 
- Kit contents 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/miRNA_lib_prep_contents.png)

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
- Pre-program thermal cycler with lid heating ON set to >95°C 
- Turn on thermocycler day-of to preheat 
- Pre-calculate how much to dilute RNA samples 

### Protocol 

#### Summary of thermocycler programs 

Create these programs beforehand with the lid heading ON and set to >95°C

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/miRNA_lib_prep_thermocycler.png)

#### Section 1: Adapter ligation and blocking 

- Pre-heat thermocycler 
- Thaw components with orange and red caps on ice 
- Thaw Ligation Enhancer at room temperature
	- Flick to mix and briefly centrifuge 
	- Pre-heat Ligation Enhancer at 37°C for 5 minutes to improve ease of pipetting, as its very viscous
- Prepare samples to correct volume in strip tubes 

| TubeID | Volume for 100 ng RNA (uL) | Volume ultrapure water (uL) | Total volume (uL) |
| ------ | -------------------------- | --------------------------- | ----------------- |
| M60    | 2.74                       | 3.76                        | 6.5               |
| M72    | 4.15                       | 2.35                        | 6.5               |

- Add 1 uL of miRNA Adapter to each sample on ice  
- Add 4 uL of Ligation Enhancer to each sample slowly (enhancer is viscous)
- Set pipette to 10 uL and mix by pipetting up and down slowly 15-20 times 
- Briefly centrifuge 
- Put samples in pre-heated thermocycler and run Steps 1-2 (Adapter Prep)
- After incubation is complete, move samples to ice 
- Add 1 uL of RNase Inhibitor to each sample 
- Add 1.5 uL of miRNA Ligation Buffer to each sample 
- Add 1 uL of miRNA Ligase to to each sample on ice 
- Set pipette to 10 uL and mix by pipetting up and down slowly 15-20 times 
- Briefly centrifuge 
- Run Steps 3-5 (Adapter Ligation)
- While in the thermocycler at the 68°C hold, add 1 uL of Blocking Agent to each sample 
- Pipette to mix 
- Run Steps 6-8 (Blocking Agent)
- While in the thermocycler at the 37°C hold, add 1 uL of Blocking Enzyme to each sample 
- Pipette to mix 
- Run Steps 9-11 (Blocking Enzyme)

#### Section 2: Circularization and Dimer Removal 

- Equilibrate Dimer Removal Beads to room temperature for >30 mins before use
	- Don't put Dimer Removal Beads on ice any time 
- Thaw components with yellow caps and Dimer Removal Agent on ice 
	- Flick to mix and briefly centrifuge 
- For each sample, prepare a master mix of the Circularization components in a PCR tube 
	- Circularization Buffer: 1.1 uL x # of samples 
	- Circularization Enzyme: 1.1 uL x # of samples 
- Add 2 uL of master mix to each sample on ice 
- Pipette to mix and briefly centrifuge 
- Run Steps 1-2 (Circularization)
- While in the thermocycler at the 37°C hold, add 1 uL of Dimer Removal Agent to each sample 
- Pipette to mix 
- Run Steps 3-4 (Dimer Anneal)
- Once Steps 3-4 are running, prepare Dimer Removal Beads 
	- Vortex Dimer Removal Beads for 30 seconds or until homogenous 
	- Aliquot 20 uL of Dimer Removal Beads into a new PCR tube for each sample 
	- Just before Step 3 (Dimer Anneal) finishes, put tubes containing beands on magnetic stand at room temperature for 1 minute or until supernatent is clear 
	- Remove and discard supernatent 
	- Move tubes with beads to a tube rack at room temperature 
- After Step 4 hold, remove samples from thermocycler and briefly centrifuge 
- Move entire sample into tube with beads 
- Pipette to mix until homogenous 
- Run Steps 5-6 (Dimer Removal)
- Once incubation is over, remove samples from thermocycler and centrifuge for 10 seconds 
- Put samples on magnetic stand at room temperature for 1 minute or until supernatent is clear 
- Move supernatent to new PCR tube and place on ice 

#### Section 3: Reverse Transcription 

- Thaw components with brown caps and the RNase Inhibitor (orange cap) on ice 
	- Flick to mix and briefly centrifuge 
- Add 4 uL of RT Primer Mix to each sample on ice 
- Pipette to mix and briefly centrifuge 
- Run Steps 1-2 (RT Primer Prep)
- Move samples to ice afterwards 
- For each sample, prepare a master mix of the Reverse Transcription components in a PCR tube 
	- RT Buffer: 3.3 uL x # of samples 
	- RNase Inhibitor: 1.1 uL x # of samples 
	- RT Enzyme: 1.65 uL x # of samples 
- Add 5.5 of master mix to each sample 
- Pipette to mix and briefly centrifuge 
- Run Steps 3-5 (Reverse Transcription)

#### Section 4: Index PCR

- Thaw ZymoTaq Premix and miRNA UDI Primers 
	- Flick to mix and briefly centrifuge 
- Add 6 uL of a unique miRNA UDI Primer to each sample 
- Add 35 uL of Zymo Taq Premix to each sample 
- Pipette to mix and briefly centrifuge 
- Run Index PCR Program 

#### Section 5: Library Purification 

- Equilibrate Select-a-Size MagBeads to room temperature at least 30 mins prior to use 
- Just before use, vortex Select-a-Size MagBeads until homogenous 
- Add 125 uL of Select-a-Size MagBeads to each sample 
- Mix thoroughly by pipetting until homogenous 
- Incubate for 5 minutes at room temperature 
- Put samples on a magnetic stand until elute is clear (~5 mins)
- Remove and discard supernatent 
- Add 200 uL of DNA Wash Buffer to each sample 
- Remove and discard supernatent 
- Repeat wash step again
- Keep tubes on magnetic stand and air dry beads for 1 minute 
- Remove any residual buffer 
- Continue to air dry the beads until matte (5-10 mins)
- Remove samples from magnetic stand and resuspend beads with 12.5 uL of DNase/RNase-free Water
- Pipette to mix 
- Incubate for 5 minutes at room temperature 
- Put samples on a magnetic stand until elute is clear (~3-5 mins)
- Transfer 10 uL of elute into a new tube 

THIS IS THE FINAL miRNA-SEQ LIBRARY. STORE AT -20°C

#### QC

Run [DNA Tapestation](https://github.com/meschedl/MESPutnam_Open_Lab_Notebook/blob/master/_posts/2019-07-30-DNA-Tapestation.md) for visualize libraries. Here's an example of what the library should look like on a Tapestation: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/DT_mcap2023/miRNA_library_visual_example.png)

