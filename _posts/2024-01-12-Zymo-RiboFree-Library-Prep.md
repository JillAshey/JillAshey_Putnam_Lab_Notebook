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

I'll be using samples M60 and M72 for this test. They were extracted on [2/10/24](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-02-10-MiniprepPlus-DNA%3ARNA-extractions-McapLarvae-DT.md) and had an average of 36.5 and 24.1 ng/ul of RNA, respectively. The kit needs a minimum input of 10 of total RNA or a maximum input of 250 ng of total RNA. I'm going to use 190 ng as my input concentration. When Maggie did this library prep with adult *M. capitata*, she used 195 ng of RNA as input for each sample (see her library prep post [here](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/zribo-lib-RNA-second/). She was working with similar RNA concentrations as me (see extraction post [here](https://emmastrand.github.io/EmmaStrand_Notebook/Holobiont-Integration-August-DNA-RNA-Extractions/), so I feel like it would be best to follow her input since she was able to produce libraries). 

Here's the RiboFree library prep workflow: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/DT_mcap2023/ribofree_workflow.png)

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
- Kit contents 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/DT_mcap2023/ribofree_kit_contents.png)

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

Thermocycler settings: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/DT_mcap2023/ribofree_cDNA_synthesis_thermocycler.png)

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

Thermocycler settings: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/DT_mcap2023/ribofree_depletion_thermocycler.png)

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
	- Aspirate slowly and discard supernatent 
	- While the sample is still on the magnetic stand, add 200 uL of Zymo-Seq Wash Buffer without disturbing the pellet 
	- Aspirate slowly and discard supernatent 
	- Repeat this wash step again 
	- Keep the caps open to air-dry bead for 1 minute 
	-  Aspirate any residual wash buffer
	-  Continue to air dry pellet until it appears matte (see picture below)
		- ![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/DT_mcap2023/ribofree_beads.png)
	- Remove sample from magnetic stand 
	- Resuspend beads in 11 uL of DNA Elution buffer 
	- Put samples on magnetic stand for 1-2 minutes or until elute is clear
	- Move elute to new PCR tubes 
- Put samples in thermocycler and run Steps 11-12 (DNA Elution)
- Remove samples from thermocycler and centrifuge briefly 
- Put samples on magnetic stand and wait until elute is clear 
- Transfer 10 uL of elute to new PCR tubes 

#### Section 3: Adapter Ligation (green caps)

- Thaw Adapter Ligation Master Mix to room temperature. Once thawed, vortex for at least 30 seconds and invert to mix 
- Preheat thermocycler to 98°C for a 3 minute incubation (lid temperature at 105°C)
- Combine the following on ice in the PCR tube containing the 10 uL elute 

| Component                 | Volume |
| ------------------------- | ------ |
| Elute                     | 10     |
| Adapter Ligation Buffer 1 | 2      |
| DNase/RNase-Free Water    | 8      |
| **TOTAL**                     | **20**     |

- Mix by pipetting and centrifuge briefly 
- Incubate the tube on ice for 2 minutes 
- Heat shock by immediately placing tube in 98°C thermocycler (105°C lid) for 3 minutes 
- Return tube to ice and incubate for at least 2 minutes
	- Reset thermocycler to 37°C (lid temperature at 45°C) for a 1 hour incubation
- Vortex the Adapter Ligation Master Mix for at least 30 seconds and invert to mix
- Add the following on ice in the order defined below to the sample tube

| Component                   | Volume |
| --------------------------- | ------ |
| Reaction from steps above   | 20     |
| Adapter Ligation Buffer 2   | 2      |
| Adapter Ligation Buffer 3   | 2      |
| Adapter Ligation Master Mix | 26     |
| **TOTAL**                       | **50**     |

- Mix entire reaction thoroughly by vortexting for 1 minute to ensure complete homogenization 
- Centrifuge briefly 
- Incubate samples at 37°C (lid temperature at 45°C) for a 1 hour incubation
- Remove samples from thermocycler 
- At room temperature, add 85 uL of DNA Elution Buffer to each sample 
- Mix by pipetting 
- Do clean-up with MagBeads 
	- Add 48.6 uL of Select-a-Size MagBeads to each sample 
	- Mix by pipetting until homogenized
	- Incubate for 5 mins at room temperature 
	- Put samples on magnetic stand for 5 mins or until beads have fully separated from solution 
	- Aspirate slowly and discard supernatent 
	- While the sample is still on the magnetic stand, add 200 uL of Zymo-Seq Wash Buffer without disturbing the pellet 
	- Aspirate slowly and discard supernatent 
	- Repeat this wash step again 
	- Keep the caps open to air-dry bead for 1 minute 
	-  Aspirate any residual wash buffer
	-  Continue to air dry pellet until it appears matte (see picture below)
		- See photo above
	- Remove sample from magnetic stand 
	- Resuspend beads in 15 uL of DNA Elution buffer 
	- Put samples on magnetic stand for 1-2 minutes or until elute is clear
	- Move elute to new PCR tubes 

#### Section 4: Library Amplification (clear caps)

Thermocycler settings: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/DT_mcap2023/ribofree_library_amp_thermocycler.png)

- Allow Library Binding Solution to equilibrate to room temperature for >30 minutes before use
- Add 10 uL of specific Zymo-Seq UDI Primers to each sample 
- Mix by pipetting 
- Add 25 uL of Amplification PreMix to each sample 
- Mix by gently pipetting and centrifuge briefly 
- Put samples in thermocycler and run program above 
- Remove samples from thermocycler and centrifuge briefly 
- At room temperature, add 50 uL of DNA Elution Buffer 
- Add 80 uL of Select-a-Size MagBeads to each sample 
- Mix thoroughly by pipetting until homogenous 
- Incubate for 5 minutes at room temperature 
- Place samples on magnetic stand for 5 minutes or until solution is clear 
- Aspirate and discard supernatent 
- Remove samples from magnetic stand 
- Add 100 uL of DNA Elution Buffer to beads 
- Mix thoroughly by pipetting until homogenous 
- Add 80 uL of Library Binding Solution to each sample 
- Mix thoroughly by pipetting until homogenous 
- Incubate for 5 minutes at room temperature 
- Do clean-up with MagBeads 
	- Put samples on magnetic stand for 5 mins or until beads have fully separated from solution 
	- Aspirate slowly and discard supernatent 
	- While the sample is still on the magnetic stand, add 200 uL of Zymo-Seq Wash Buffer without disturbing the pellet 
	- Aspirate slowly and discard supernatent 
	- Repeat this wash step again 
	- Keep the caps open to air-dry bead for 1 minute 
	-  Aspirate any residual wash buffer
	-  Continue to air dry pellet until it appears matte (see picture below)
		- See photo above
	- Remove sample from magnetic stand 
	- Resuspend beads in 20 uL of DNA Elution buffer 
	- Put samples on magnetic stand for 1-2 minutes or until elute is clear
	- Move elute to new PCR tubes 

THIS IS THE FINAL RNA-SEQ LIBRARY. STORE AT -20°C

#### QC

Run [DNA Tapestation](https://github.com/meschedl/MESPutnam_Open_Lab_Notebook/blob/master/_posts/2019-07-30-DNA-Tapestation.md) for visualize libraries. Here's an example of what the library should look like on a Tapestation: 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/DT_mcap2023/ribofree_library_visual_example.png)

Here's how my samples turned out on 2/15/24. Decent concentration in M72, but not really in M60.

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_overview_20240215.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M60_20240215.png)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/tapestation/DNA_TS_M72_20240215.png)

The peaks are in the correct places but not as pronounced as they should be. They are also spread out a decent amount on the x-axis. Maybe I should increase the number of PCR amplification cycles? Discussed with Cassie and she recommended to increase the number of PCR cycles. Maggie (she used an older version of the ribofree kit) had a similar problem to me where she was getting low concentrations and the peaks were more like hills when there were 12 PCR amplification cycles (see her protocol [here](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/zribo-lib-RNA-second/)). Her tapestation results are [here](https://github.com/meschedl/MESPutnam_Open_Lab_Notebook/blob/master/tapestation_pdfs/2019-09-10%20-%2013.12.25.pdf). For the amplification step, she ended up increasing the number of cycles to 15 and tapestation [results](https://github.com/meschedl/MESPutnam_Open_Lab_Notebook/blob/master/tapestation_pdfs/2019-09-11%20-%2009.28.34.pdf) looked much better (ie higher peaks, not as wide). Using the same samples (M60 and M72), I will increase the number of PCR cycles to 15 and keep the rest of the protocol the same. 

Overall, I feel pretty good that I got peaks/hills in the correct locations! I think that increasing the number of amplification cycles will help remedy this issue. 


