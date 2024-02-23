---
layout: post
title: Pico Methyl-Seq library prep
date: '2024-02-26'
categories:
tags: [DNA, Library prep, Protocols, WGBS]
projects: e5
---

# Pico Methy-Seq Library Prep test for E5 samples

This post details the info about the WGBS library prep steps for the e5 deep dive samples collected in Moorea in 2020. The github for that project is linked [here](https://github.com/urol-e5/deep-dive). I'm using the [Zymo Pico Methyl Seq Library Prep](https://www.zymoresearch.com/products/pico-methyl-seq-library-prep-kit) for library prep (we have several kits with 25 preps). See Zymo's protocol [here](https://files.zymoresearch.com/protocols/_d5455_d5456_picomethylseq.pdf). 

Since there are 3 species for this project (*Acropora pulchra*, *Pocillopora spp.*, and *Porites evermanni*), I'll be doing a test of one of each species today. 

ADD INFO ABOUT EXTRACTIONS AND WHAT AMT OF DNA TO ADD (Emma suggested starting with 20 ng of DNA). Both Maggie and Emma have done this protocol before (add links to their github posts)

Here's the Pico Methyl-Seq library prep workflow: 

ADD 

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

ADD 

### Buffer preperation 

Once buffers are prepared for a kit, they do not need to be prepared again. 

- Add 26 mL of 95% EtOH to the M-Wash Buffer concentrate
- Add 26 mL of 95% EtOH to the DNA Wash Buffer concentrate

### Best practices 

- Preset the thermocycler programs 
- Thaw and keep -80°C and -20°C components on ice unless instructed otherwise. Flick to mix and centrifuge before use
	- Avoid multiple freeze-thaws, make aliquots if needed

#### Section 1: Bisulfite Conversion of DNA 

Thermocycler settings: 

XXXXX

- Thaw samples on ice 
- Prepare samples to correct volume in strip tubes 

XXXXXX add info 

- Add 130 uL of Lightning Conversion Reagent to each tube 
- Run thermocycler steps 
- Store at 4°C for <20 hours (overnight)

NEXT DAY!

- Label Zymo-spin IC columns and collection tubes for each sample 
- Add 600 uL of M-Binding Buffer to each column
- Add bisulfite-converted sample to column and invert several times to mix 
- Centrifuge at 10,000 g for 30 seconds 
- Discard flow through 
- Add 100 uL M-Wash Buffer 
- Centrifuge at 10,000 g for 30 seconds + discard flow through 
- Add 200 uL L-Desulphonation Buffer to column and incubate at room temperature for 15 minutes 
	- Near end of incubation, heat DNA Elution buffer to 56°C on thermoblock 
- After incubation, centrifuge at 10,000 g for 30 seconds + discard flow through
- Add 200 uL of M-Wash Buffer 
- Centrifuge at 10,000 g for 30 seconds + discard flow through 
- Repeat the wash step above 
- Add 8 uL of warmed DNA Elution buffer to column and let sit for 1 minute
- Centrifuge at 10,000 g for 30 seconds to elute bisulfite-converted DNA

#### Section 2: Amplification with PrepAmp Primer 

Thermocycler settings: 

XXXXX

Lid should be set at 25°C

- Thaw the following reagents on ice 
	- PreAmpBuffer (5x)
	- PreAmp Primer (40 uM)
	- PreAmp pre-mix 
	- PreAmp polymerase 
- Make Priming Reaction on ice for each sample in new strip tubes 

| Component                 | Volume (uL) |
| ------------------------- | ------ |
| PreAmpBuffer (5x)                     | 2     |
| PreAmp Primer (40 uM) | 1      |
| Bisulfite-converted DNA    | 7      |
| **TOTAL**                     | **10**     |

- Pipette up and down to mix, spin down
- In a separate tube, combine the following for the PreAmp Master Mix
	- Calculate how much of each component will be needed for the number of samples 
		- 1 uL PreAmpBuffer (5x) x number of samples 
		- 3.75 uL PreAmp Pre-mix x number of samples 
		- 0.3 uL PreAmp Polymerase x number of samples 
	- Each sample will be getting 5.05 uL of master mix 
- Make 'diluted' PreAmp Polymerase mix on ice to avoid adding less than 0.5 uL 
	- The original protocol asks you to add 0.3ul to the tubes in the thermocycler, sometimes that small of an amount does not come out of the tip so you can add DNA elution buffer to the enzyme to pipette 0.5 ul
	- 0.3 uL PreAmp Polymerase x number of samples 
	- 0.2 uL DNA Elution Buffer x number of samples 
- Run the thermocycler program with the lid temperature at 25°C with 2 cycles 
	- During the Step 2 8°C hold:
		- Cycle 1: Add 5.05 ul of PreAmp Master Mix. Pipette up and down to mix, centrifuge briefly, put back in thermocycler and proceed past hold 
		- Cycle 2: Add 0.5 uL of diluted PreAmp polymerase to each tube.Pipette up and down to mix, centrifuge briefly, put back in thermocycler and proceed past hold 

#### Section 3: Purification with DNA Clean-up and Concentrator (DCC)

- Label 1.5 mL tubes + IC columns for each sample
- Heat DNA Elution buffer to 56°C on thermoblock
- Move the sample from the PCR tube into the 1.5 mL tube 
- Add 108.85 uL of DNA binding buffer (7:1 ratio of DNA binding buffer:sample)
- Vortex and spin down to mix
- Move sample into a Zymo-Spin IC column 
- Centrifuge at 10,000 g for 30 seconds + discard flow through 
- Add 200 uL of DNA Wash Buffer
- Centrifuge at 10,000 g for 30 seconds + discard flow through 
- Repeat this wash step 
- Move spin columns to 1.5 mL tubes 
- Add 12 uL of warmed DNA Elution Buffer to the column matrix 
- Incubate for 1 minute at room temperature 
- Centrifuge at 10,000 g for 30 seconds to elute 

#### Section 4: Amplification 

Thermocycler settings: 

XXXXX

- In a new PCR tube, mix the following: 

| Component                 | Volume (uL) |
| ------------------------- | ------ |
| LibraryAmp Master Mix (2x)                     | 12.5     |
| LibraryAmp Primer (10 uM) | 1      |
| Sample elute    | 11.5      |
| **TOTAL**                     | **25**     |

- Pipette up and down to mix
- Centrifuge briefly 
- Run thermocycler program for a total of XXX cycles 
- Label 1.5 mL tubes + IC columns for each sample
- Heat DNA Elution buffer to 56°C on thermoblock
- Move the sample from the PCR tube into the 1.5 mL tube 
- Add 175 uL of DNA binding buffer (7:1 ratio of DNA binding buffer:sample)
- Vortex and spin down to mix
- Move sample into a Zymo-Spin IC column 
- Centrifuge at 10,000 g for 30 seconds + discard flow through 
- Add 200 uL of DNA Wash Buffer
- Centrifuge at 10,000 g for 30 seconds + discard flow through 
- Repeat this wash step 
- Move spin columns to 1.5 mL tubes 
- Add 12.5 uL of warmed DNA Elution Buffer to the column matrix 
- Incubate for 5 minutes at room temperature 
- Centrifuge at 10,000 g for 30 seconds to elute 

#### Section 5: Amplification with Index Primers 

Thermocycler settings: 

XXXXX

- In a new PCR tube, mix the following: 

| Component                 | Volume (uL) |
| ------------------------- | ------ |
| LibraryAmp Master Mix (2x)                     | 12.5     |
| Index Primer (10 uM) | 0.5      |
| Sample elute    | 12      |
| **TOTAL**                     | **25**     |

- Run thermocycler program for a total of 12 cycles 
	- The zymo protocol says 10 cycles but in Maggie's protocol, she did 12 cycles 
- Label 1.5 mL tubes + IC columns for each sample
- Heat DNA Elution buffer to 56°C on thermoblock
- Move the sample from the PCR tube into the 1.5 mL tube 
- Add 175 uL of DNA binding buffer (7:1 ratio of DNA binding buffer:sample)
- Vortex and spin down to mix
- Move sample into a Zymo-Spin IC column 
- Centrifuge at 10,000 g for 30 seconds + discard flow through 
- Add 200 uL of DNA Wash Buffer
- Centrifuge at 10,000 g for 30 seconds + discard flow through 
- Repeat this wash step 
- Move spin columns to 1.5 mL tubes 
- Add 12 uL of warmed DNA Elution Buffer to the column matrix 
- Incubate for 1 minute at room temperature 
- Centrifuge at 10,000 g for 30 seconds to elute 

NEED TO COMPLETE PROTOCOL AND COMPARE WITH EMMA + MAGGIE'S 

Emma suggested using 20ng of DNA, double check that 
