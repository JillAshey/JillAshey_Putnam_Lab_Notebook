---
layout: post
title: RiboFree library prep
date: '2023-10-27'
categories:
tags: [DNA, RNA, Library prep, Protocols]
projects: Mcap developmental time series 
---

# RiboFree library prep test for Mcap developmental time series Hawaii 2023

This post details the info about the library prep steps for the *M. capitata* developmental time series experiment in Hawaii 2023. The github for that project is linked [here](https://github.com/JillAshey/Hawaii_Developmental_TimeSeries). I'm using the [Zymo-Seq RiboFree Total RNA Library Kit](https://www.zymoresearch.com/products/zymo-seq-ribofree-total-rna-library-kit) for total RNA library prep (I got the 12 prep kit). See Zymo's protocol [here](https://files.zymoresearch.com/protocols/r3000_zymo-seq_ribofree_total_rna_library_kit.pdf). 

I'll be using samples M60 and M72 for this test. They were extracted on [2/10/24](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-02-10-MiniprepPlus-DNA%3ARNA-extractions-McapLarvae-DT.md) and had an average of 36.5 and 24.1 ng/ul of RNA, respectively. The kit needs a minimum input of 10 of total RNA or a maximum input of 250 ng of total RNA. I'm going to use 190 ng as my input concentration. When Maggie did this library prep with adult *M. capitata*, she used 195 ng of RNA as input for each sample (see her library prep post [here](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/zribo-lib-RNA-second/). She was working with similar RNA concentrations as me (see extraction post [here](https://emmastrand.github.io/EmmaStrand_Notebook/Holobiont-Integration-August-DNA-RNA-Extractions/), so I feel like it would be best to follow her input since she was able to produce libraries. 

Here's the RiboFree library prep workflow: 

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


### Materials 

- [Zymo-Seq RiboFree library kit](https://www.zymoresearch.com/products/zymo-seq-ribofree-total-rna-library-kit)
- PCR tubes 
- Thermocycler 
- Heating block 
- Mini centrifuge
- Vortex 
- Aluminum beads (to keep things on ice)
- Magnetic stand for PCR tubes 
- 1.5 mL tubes 

### Buffer preperation 

Once buffers are prepared for a kit, they do not need to be prepared again. 

- Add 300 uL of the Select-a-Size MagBead concentrate to 10 mL of Select-a-Size MagBead Buffer. 
	- It was difficult to pipette the concentrate, and I'm sure I missed some beads. In Maggie's protocol, she added 1mL of MagBead buffer to the bead concentrate, it was very hard to pipette the concentrated beads, which is why she did not add the beads to the buffer as it says in the protocol
	- Shake or vortex to mix. Store at 4°C. 
	- These should be prepared at least 5 days ahead of library prep. 
- Add 24 mL of 100% ethanol to the Zymo-Seq Wash Buffer concentrate. Store at room temperature. 


### Best practices 

- Preset the thermocycler programs 
- Set thermocycler lid to 105°C unless instructed otherwise. 
- Thaw and keep -80°C and -20°C components on ice unless instructed otherwise. Flick to mix and centrifuge before use
	- Avoid multiple freeze-thaws, make aliquots if needed
	- The Adapter Ligation Buffer 2, Adapter Ligation Buffer 3, and the Adapter Ligation Master Mix should only be thawed **4 times max**. Make aliquots.
- Allow Select-a-Size MagBeads to equilibrate to room temperature >30 mins before use 
- Resuspend magnetic particles immediately before each use by vigorously inverting and vortexing the Select-a-Size MagBeads until homogenous

### Protocol 

#### Section 1: cDNA Synthesis (yellow caps)

- Set thermocycler program (XXXXXX) 
- Thaw samples and reagents on ice (cDNA Synthesis Reagent 1 and cDNA Synthesis Reagent 2)
- Prepare samples to correct volume in strip tubes 

| TubeID | Volume for 190 ng RNA (uL) | Volume ultrapure water (uL) | Total volume (uL) |
| ------ | -------------------------- | --------------------------- | ----------------- |
| M60    | 5.21                       | 2.79                        | 8                 |
| M72    | 7.88                       | 0.12                        | 8                 |

- Add 2 uL of cDNA synthesis reagent 1 to each sample tube for a total of 10 uL. 
- Mix by flicking or pipetting. Centrifuge briefly
- Put samples in thermocycler and run Steps 1-2 (Primer Annealing)
- Add 10 uL of cDNA synthesis reagent 2 to each sample for a total of 20 uL while in thermocycler on the 4°C hold. 
- Mix by pipetting. Centrifuge briefly 
- Run Steps 3-5 (Reverse Transcription)

#### Section 2: RiboFree Universal Depletion (red caps)

- While in thermocycler on the 4°C hold, add 10 uL of Depletion Reagent 1 to each sample for a total of 30 uL. 
- Mix by pipetting. Centrifuge briefly 
- Put samples in thermocycler and run Steps 1-3 (Pre-Depletion Incubation)
- While in thermocycler on the 68°C hold, add 10 uL of Depletion Reagent 2 to each sample for a total of 40 uL. 
- Mix by pipetting 
- Run Step 4 (Depletion Reaction)
- While in thermocycler on the 68°C hold, add 10 uL of Depletion Reagent 3 to each sample for a total of 50 uL. 
- Mix by pipetting 
- Run Steps 6-7 (Stop Depletion)
- While in thermocycler on the 25°C hold, add 2 uL of Depletion Reagent 4 to to each sample for a total of 52 uL. 
- Flick to mix and centrifuge briefly. Put samples back in thermocycler fast 
- Run Steps 8-10 (Depletion Cleanup)
- Remove tubes from thermocycler and centrifuge briefly
- At room temperature, add 26 uL of 95% ethanol to each sample for a total of 78 uL 
- Mix by pipetting 
- Do clean-up with MagBeads 
	- Add 156 uL of Select-a-Size MagBeads to each sample 
	- Mix by pipetting until homogenized
	- Incubate for 5 mins at room temperature 
	- Put samples on magnetic stand for 5 mins or until beads have fully separated from solution 
	- Aspirate slowly and discard and supernatent 
	- While the sample is still on the magnetic stand, add 200 uL of Zymo-Seq Wash Buffer without disturbing the pellet 
	- Aspirate slowly and discard supernatent 
	- Repeat this wash step again 
	- Keep the caps open to air-dry bead
	-  





