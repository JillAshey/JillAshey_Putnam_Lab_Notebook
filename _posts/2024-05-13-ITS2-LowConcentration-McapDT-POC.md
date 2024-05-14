---
layout: post
title: ITS2 
date: '2024-05-13'
categories:
tags: [DNA, ITS2, Extractions, Protocols]
projects: Mcap DT 2023; POC 2023 spawning 
---

# ITS2 for Mcap developmental timeseries Hawaii 2023 and POC spawning 2023 sample 

This post details information on ITS2 amplification for the my ambient Mcap developmental timeseries samples from 2023 and for the POC spawning experimental samples from 2023. Mcap 2023 github repo is [here](https://github.com/JillAshey/Hawaii_Developmental_TimeSeries) and POC github is [here](https://github.com/hputnam/Poc_RAPID). 

These samples had pretty low DNA concentrations (see DNA QC post [here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-05-08-DNA-QC.md) for sample concentration info). I did a test run with 6 samples (not in triplicate) diluted to 0.104 ng/uL to evaluate if the ITS2 region would amplify with low DNA input (see post [here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-05-09-ITS2-Test-LowConcentration.md)). The test was successful and all samples amplified. I can now move forward with processing all samples in triplicate diluted to <1 ng/uL. 

I followed Ariana's ITS2 [protocol](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/_posts/2024-03-30-ITS2-amplification-protocol.md) with some modifications for testing some samples with low concentrations. 

### Equipment and materials 

- 2X Phusion Mastermix
- Forward primer (ITSintfor2) at 100uM: TCG TCG GCA GCG TCA GAT GTG TAT AAG AGA CAG GAA TTG CAG AAC TCC GTG
- Reverse primer (ITS2_Reverse) at 100uM: GTC TCG TGG GCT CGG AGA TGT GTA TAA GAG ACA GGG GAT CCA TAT GCT TAA GTT CAG CGG GT
- Loading Dye NEB 6X Purple Loading Dye NEB Cat # B7024S
- Gel Stain Biotium GelGreen Nucleic Acid Gel Stain, 10,000X in Water Fisher Cat NC9728313
- DNA Ladder 1kb Gel ladder and 100bp gel ladder
- Agarose and 1XTAE buffer for gel electrophoresis
- Ultrapure water
- PCR strip tubes and 1.5 mL tubes

### Protocol 

Before starting the protocol, prepare the following: 

- Sample list. If you are performing preparations in multiple batches/PCRs, randomize the sample order for processing to remove batch effects confounding with treatment groups.

- Prepare list of appropriate controls. Controls will include a negative control (input = nuclease-free water) and a positive control (gDNA from a previously successfully amplified and sequenced sample). These controls should be included in each batch. If you perform preparations over multiple days, include a positive and negative control each day or with each PCR.

- Label 3 PCR strip tubes per sample (including negative and positive controls) with the sample ID
	- 96 well PCR plates can also be used for lots of samples 
	- I used a PCR plate for my samples since I had so many and they were in triplicate. 

- Labeled 1 PCR strip tube per sample for pooled PCR product

#### Sample Dilution

First, dilute gDNA from each of your samples to a consistent concentration. ITS2 amplifies readily, so input DNA into the first PCR can be at a low concentration. We dilute samples so that there is reduced bias towards high DNA concentration samples during PCR and subsequent sequencing. My lowest DNA concentration was 0.104 ng/uL so I will attempt to dilute all of my samples to that concentration. 

- Calculate the volume required to obtain ~0.1 ng/uL gDNA in a 11.5 uL total volume.
	- This step is specific to my samples. Usually, DNA is diluted in 10 uL
- Thaw gDNA on ice 
- Aliquot the required volume of water for each sample into the labeled gDNA strip tubes.
- For each sample, aliquot the required volume of gDNA into the labeled gDNA strip tubes
- After adding gDNA and water to each tube, vortex for ~5 sec and then spin down.
- Samples can be stored at -20°C until proceeding to the next step. 

Here are the samples that I did today: 

| TubeID              | hpf     | Treatment            | Project                       | Qubit HS (ng/uL) | DNA for dilution (uL) | Water (uL) | Total volume | Expected final concentration (ng/uL) | Actual DNA input | Actual water input | Actual final concentration (ng/uL) | Batch |
| ------------------- | ------- | -------------------- | ----------------------------- | ---------------- | --------------------- | ---------- | ------------ | ------------------------------------ | ---------------- | ------------------ | ---------------------------------- | ----- |
| M48                 | 22 hpf  | Ambient              | Developmental Timeseries 2023 | 0.949            | 1.10                  | 10.40      | 11.50        | 0.09043478261                        | 1.10             | 10.40              | 0.09043478261                      | 1     |
| M39                 | 14 hpf  | Ambient              | Developmental Timeseries 2023 | 0.46             | 2.26                  | 9.24       | 11.50        | 0.09043478261                        | 2.26             | 9.24               | 0.09043478261                      | 1     |
| M73                 | 48 hpf  | Ambient              | Developmental Timeseries 2023 | 8.62             | 0.12                  | 11.38      | 11.50        | 0.09043478261                        | 0.5              | 11                 | 0.3747826087                       | 1     |
| M9                  | 1 hpf   | Ambient              | Developmental Timeseries 2023 | 0.104            | 10.00                 | 1.50       | 11.50        | 0.09043478261                        | 10.00            | 1.50               | 0.09043478261                      | 1     |
| 107                 | Larvae  | High                 | POC                           | 0.582            | 1.79                  | 9.71       | 11.50        | 0.09043478261                        | 1.79             | 9.71               | 0.09043478261                      | 1     |
| 108                 | Larvae  | Ambient              | POC                           | 0.549            | 1.89                  | 9.61       | 11.50        | 0.09043478261                        | 1.89             | 9.61               | 0.09043478261                      | 1     |
| M51                 | 22 hpf  | Ambient              | Developmental Timeseries 2023 | 3.53             | 0.29                  | 11.21      | 11.50        | 0.09043478261                        | 0.5              | 11                 | 0.1534782609                       | 1     |
| M14                 | 4 hpf   | Ambient              | Developmental Timeseries 2023 | 0.263            | 3.95                  | 7.55       | 11.50        | 0.09043478261                        | 3.95             | 7.55               | 0.09043478261                      | 1     |
| M60                 | 28 hpf  | Ambient              | Developmental Timeseries 2023 | 0.292            | 3.56                  | 7.94       | 11.50        | 0.09043478261                        | 3.56             | 7.94               | 0.09043478261                      | 1     |
| 8                   | Embryos | Ambient              | POC                           | 0.439            | 2.37                  | 9.13       | 11.50        | 0.09043478261                        | 2.37             | 9.13               | 0.09043478261                      | 1     |
| M23                 | 9 hpf   | Ambient              | Developmental Timeseries 2023 | 0.116            | 8.97                  | 2.53       | 11.50        | 0.09043478261                        | 8.97             | 2.53               | 0.09043478261                      | 1     |
| M72                 | 48 hpf  | Ambient              | Developmental Timeseries 2023 | 0.848            | 1.23                  | 10.27      | 11.50        | 0.09043478261                        | 1.23             | 10.27              | 0.09043478261                      | 1     |
| M52                 | 22 hpf  | Ambient              | Developmental Timeseries 2023 | 0.412            | 2.52                  | 8.98       | 11.50        | 0.09043478261                        | 2.52             | 8.98               | 0.09043478261                      | 1     |
| M8                  | 1 hpf   | Ambient              | Developmental Timeseries 2023 | 0.155            | 6.71                  | 4.79       | 11.50        | 0.09043478261                        | 6.71             | 4.79               | 0.09043478261                      | 1     |
| 110                 | Larvae  | Ambient              | POC                           | 0.944            | 1.10                  | 10.40      | 11.50        | 0.09043478261                        | 1.10             | 10.40              | 0.09043478261                      | 1     |
| M35                 | 14 hpf  | Ambient              | Developmental Timeseries 2023 | 0.275            | 3.78                  | 7.72       | 11.50        | 0.09043478261                        | 3.78             | 7.72               | 0.09043478261                      | 1     |
| M63                 | 28 hpf  | Ambient              | Developmental Timeseries 2023 | 4.47             | 0.23                  | 11.27      | 11.50        | 0.09043478261                        | 0.5              | 11                 | 0.1943478261                       | 1     |
| M28                 | 9 hpf   | Ambient              | Developmental Timeseries 2023 | NA               | NA                    | NA         | NA           | NA                                   | 11.5             | 0                  | 0                                  | 1     |
| 106                 | Larvae  | Ambient              | POC                           | 0.788            | 1.32                  | 10.18      | 11.50        | 0.09043478261                        | 1.32             | 10.18              | 0.09043478261                      | 1     |
| 111                 | Larvae  | High                 | POC                           | 1.81             | 0.57                  | 10.93      | 11.50        | 0.09043478261                        | 0.57             | 10.93              | 0.09043478261                      | 1     |
| M6                  | 1 hpf   | Ambient              | Developmental Timeseries 2023 | 0.333            | 3.12                  | 8.38       | 11.50        | 0.09043478261                        | 3.12             | 8.38               | 0.09043478261                      | 1     |
| M13                 | 4 hpf   | Ambient              | Developmental Timeseries 2023 | 0.989            | 1.05                  | 10.45      | 11.50        | 0.09043478261                        | 1.05             | 10.45              | 0.09043478261                      | 1     |
| 112                 | Larvae  | Ambient              | POC                           | 1.71             | 0.61                  | 10.89      | 11.50        | 0.09043478261                        | 0.61             | 10.89              | 0.09043478261                      | 1     |
| M83                 | 72 hpf  | Ambient              | Developmental Timeseries 2023 | 19.45            | 0.05                  | 11.45      | 11.50        | 0.09043478261                        | 0.5              | 11                 | 0.8456521739                       | 1     |
| M24                 | 9 hpf   | Ambient              | Developmental Timeseries 2023 | 0.114            | 9.12                  | 2.38       | 11.50        | 0.09043478261                        | 9.12             | 2.38               | 0.09043478261                      | 1     |
| M36                 | 14 hpf  | Ambient              | Developmental Timeseries 2023 | 1.085            | 0.96                  | 10.54      | 11.50        | 0.09043478261                        | 0.96             | 10.54              | 0.09043478261                      | 1     |
| R55 (POS)           | Larvae  | Cladocopium; Ambient | Mcap 2023 AH                  | 14.3             | 0.07                  | 11.43      | 11.50        | 0.09043478261                        | 0.5              | 11                 | 0.6217391304                       | 1     |
| Ultrapure H2O (NEG) | NA      | NA                   | NA                            | NA               | NA                    | 11.5       | 11.5         | 0                                    | 0                | 11.5               | 0                                  | 1     |

The smallest pipette we have is 0.5, so I had to round up for some samples (expected v. actual), which made the concentrations a little all over the place (but all still <1 ng/uL). I could do a serial dilution in the future, but as long as the concentration is <1ng/uL, I'm happy. 

#### Make master mix 

Master mix will be prepared to allow for triplicate reactions for each sample. Reactions are performed in triplicate to 1) increase total PCR product for each sample, 2) provide individual replicates in the case that some reactions do not amplify or do not amplify well and 3) increase diversity of sequences amplified by chance. 

- Prepare master mix components. Thaw the Phusion Master Mix on ice. Thaw the primers on ice. If required, hydrate lyophilized primers if not done so already. Thaw nuclease-free water on ice. Thaw diluted gDNA samples on ice.
- Calculate the volume of each regent needed by first calculating the total number of reactions requires. The total reactions will be the number of samples (including positive and negative control) x 3. Include an extra 10% volume for pipetting error.
	- Once again, only doing 1 replicate today so my total number of reactions + error is 9.
- Label a 1.5 mL tube for master mix 
- Label two 1.5 mL tubes for dilution of primers (if required). This protocol calls for 10uM primers. Primers are often received or hydrated at 100uM. To dilute, add 1 part primer to 9 parts nuclease free water. For example, to make 10 total uL of primer stock, add 1 uL of the respective primer to 9 uL of H20. See the table below for the calculation for total primer stock needed to make the master mix. Do this for both the F and R primer separately. Vortex and spin down after diluting.
- Add in volume of Ultra pure/ nuclease free H20 required to the master mix tube. I did NOT do this today, as I added 11.5 uL of diluted DNA to the tube. The protocol typically requires 10.5 uL of ultrapure water per reaction and 1 ng of sample DNA. 
- Add the necessary volume of Phusion master mix to the 1.5 mL tube
- Add the necessary 10 uM F and R primers to the 1.5 mL tube

Master mix calculations for today: 

|                       | Per rxn (uL) | Rxns | Total volume (uL) |
| --------------------- | ------------ | ---- | ------------ |
| 2x Phusion Master Mix | 12.5         | 94    | 1175        |
| F Primer (10um)       | 0.5          | 94    | 47          |
| R Primer (10um)       | 0.5          | 94    | 47          |
| Ultrapure H2O         | 0            | 94    | 0            |
| DNA                   | 11.5         |      |              |
| TOTAL VOLUME PER RXN  | 25           |      |

- Add 13.5 uL of master mix to the 11.5 uL of DNA sample

**IMPORTANT NOTE FOR TODAY'S AMP**: Instead of adding 13.5 uL of master mix to the diluted sample, I added 24 uL uL of master mix.  I got mixed up with the traditional ITS2 [protocol](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/_posts/2024-03-30-ITS2-amplification-protocol.md) where 24 uL of master mix is added to 1 uL of dilute DNA. I realized this about halfway through loading the plate with master mix. At that point, I decided to just load them all with 24 uL and see what happens. The ratios of the master mix, primers and DNA will be off but the ITS2 region still should amplify. I anticipate there will be more primer dimer or other sort of contamination because there is an excess of primers and master mix. 

- Vortex and spin down. 
- Keep tubes on ice. Samples are now ready for PCR

#### PCR 

Conduct PCR to amplify the ITS2 region. Load all reaction tubes (including negative and positive controls) to the thermocycler. Run the following program (pre-programmed "ITS2" PCR protocol under PUTNAM):

- 1 cycle at 95°C for 3 min
- 35 cycles at 95°C for 30 sec, 52°C for 30 sec, and 72°C for 30 sec
- 1 cycle at 72°C for 2 min
- Hold at 4°C

PCR will run for approx. 1.5 hours. Samples will stay at 4°C on the thermocycler until running a QC gel or store samples at 4°C in fridge overnight.

This is a good stopping point if needed.

#### Gel QC for PCR product 

Run a 2% gel at 80-100V and 100 amps for 1.5 hours. Use the 1 kb DNA ladder and 100 bp ladder to each row. Today, I used the large gel box with 100 spaces. For some reason, the voltage will not go past ~70 volts. Since the gel was so large, I had to break it into separate images: 

Top left 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/gel_top_left_20240513.png)

Top right 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/gel_top_right_20240513.png)

Bottom left 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/gel_bottom_left_20240513.png)

Bottom right 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/gel_bottom_right_20240513.png)

Everything amplified successfully! I poked holes in the bottom of the wells that had samples 107, 108, M52 and M8 so that just looks like a blob. I do have a decent amount of primer dimer in some of my samples. To remedy this issue, I will do a bead cleanup (example [here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-04-18-ITS2-Bead-Cleanup-McapLarvae.md)). 