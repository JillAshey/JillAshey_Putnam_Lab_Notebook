---
layout: post
title: ITS2 
date: '2024-05-14'
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

| TubeID              | hpf     | Treatment            | Project                       | Qubit HS (ng/uL) | DNA for dilution (uL) | Water (uL) | Total volume | Expected final concentration (ng/uL) | Actual DNA input | Actual water input | Actual final concentration (ng/uL) | Batch |
| ------------------- | ------- | -------------------- | ----------------------------- | ---------------- | --------------------- | ---------- | ------------ | ------------------------------------ | ---------------- | ------------------ | ---------------------------------- | ----- |
| M61                 | 28 hpf  | Ambient              | Developmental Timeseries 2023 | 8.66             | 0.12                  | 11.38      | 11.50        | 0.09043478261                        | 0.5              | 11                 | 0.3765217391                       | 2     |
| M88                 | 72 hpf  | Ambient              | Developmental Timeseries 2023 | 67.2             | 0.02                  | 11.48      | 11.50        | 0.09043478261                        | 0.62             | 10.88              | 3.622956522                        | 2     |
| M74                 | 48 hpf  | Ambient              | Developmental Timeseries 2023 | 6.28             | 0.17                  | 11.33      | 11.50        | 0.09043478261                        | 0.5              | 11                 | 0.2730434783                       | 2     |
| M75                 | 48 hpf  | Ambient              | Developmental Timeseries 2023 | 7.82             | 0.13                  | 11.37      | 11.50        | 0.09043478261                        | 0.5              | 11                 | 0.34                               | 2     |
| M47                 | 22 hpf  | Ambient              | Developmental Timeseries 2023 | 2.22             | 0.47                  | 11.03      | 11.50        | 0.09043478261                        | 0.5              | 11                 | 0.09652173913                      | 2     |
| M62                 | 28 hpf  | Ambient              | Developmental Timeseries 2023 | 9.52             | 0.11                  | 11.39      | 11.50        | 0.09043478261                        | 0.5              | 11                 | 0.4139130435                       | 2     |
| M86                 | 72 hpf  | Ambient              | Developmental Timeseries 2023 | 41.1             | 0.03                  | 11.47      | 11.50        | 0.09043478261                        | 1.01             | 10.49              | 3.609652174                        | 2     |
| M7                  | 1 hpf   | Ambient              | Developmental Timeseries 2023 | 0.129            | 8.06                  | 3.44       | 11.50        | 0.09043478261                        | 8.06             | 3.44               | 0.09043478261                      | 2     |
| M37                 | 14 hpf  | Ambient              | Developmental Timeseries 2023 | 0.977            | 1.06                  | 10.44      | 11.50        | 0.09043478261                        | 1.06             | 10.44              | 0.09043478261                      | 2     |
| M10                 | 4 hpf   | Ambient              | Developmental Timeseries 2023 | 0.132            | 7.88                  | 3.62       | 11.50        | 0.09043478261                        | 7.88             | 3.62               | 0.09043478261                      | 2     |
| M26                 | 9 hpf   | Ambient              | Developmental Timeseries 2023 | 0.11             | 9.45                  | 2.05       | 11.50        | 0.09043478261                        | 9.45             | 2.05               | 0.09043478261                      | 2     |
| 9                   | Embryos | Ambient              | POC                           | 0.27             | 3.85                  | 7.65       | 11.50        | 0.09043478261                        | 3.85             | 7.65               | 0.09043478261                      | 2     |
| 6                   | Embryos | Ambient              | POC                           | 0.131            | 7.94                  | 3.56       | 11.50        | 0.09043478261                        | 7.94             | 3.56               | 0.09043478261                      | 2     |
| M11                 | 4 hpf   | Ambient              | Developmental Timeseries 2023 | 0.136            | 7.65                  | 3.85       | 11.50        | 0.09043478261                        | 7.65             | 3.85               | 0.09043478261                      | 2     |
| 113                 | Larvae  | High                 | POC                           | 1.245            | 0.84                  | 10.66      | 11.50        | 0.09043478261                        | 0.84             | 10.66              | 0.09043478261                      | 2     |
| 7                   | Embryos | Ambient              | POC                           | 0.327            | 3.18                  | 8.32       | 11.50        | 0.09043478261                        | 3.18             | 8.32               | 0.09043478261                      | 2     |
| 109                 | Larvae  | High                 | POC                           | 2.77             | 0.38                  | 11.12      | 11.50        | 0.09043478261                        | 0.5              | 11                 | 0.1204347826                       | 2     |
| M87                 | 72 hpf  | Ambient              | Developmental Timeseries 2023 | 29.7             | 0.04                  | 11.46      | 11.50        | 0.09043478261                        | 1.41             | 10.09              | 3.641478261                        | 2     |
| R55 (POS)           | Larvae  | Cladocopium; Ambient | Mcap 2023 AH                  | 14.3             | 0.07                  | 11.43      | 11.50        | 0.09043478261                        | 0.5              | 11                 | 0.6217391304                       | 2     |
| Ultrapure H2O (NEG) | NA      | NA                   | NA                            | NA               | NA                    | 11.5       | 11.5         | 0                                    | 0                | 11.5               | 0                                  | 2     |

The smallest pipette we have is 0.5, so I had to round up for some samples (expected v. actual), which made the concentrations a little all over the place (but all still <1 ng/uL). I did a initial dilution for M86, M87, and M88 because they all had high concentrations. 

- M86: 1 uL of DNA at 29.7 ng/uL + 39 uL H20 = 40 uL at 0.74 ng/uL. Use this as the input concentration
- M87: 1 uL of DNA at 41.1 ng/uL + 39 uL H20 = 40 uL at 1.03 ng/uL. Use this as the input concentration
- M88: 1 uL of DNA at 67.2 ng/uL + 39 uL H20 = 40 uL at 1.68 ng/uL. Use this as the input concentration

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
| 2x Phusion Master Mix | 12.5         | 66    | 750        |
| F Primer (10um)       | 0.5          | 66    | 33          |
| R Primer (10um)       | 0.5          | 66    | 33          |
| Ultrapure H2O         | 0            | 66    | 0            |
| DNA                   | 11.5         |      |              |
| TOTAL VOLUME PER RXN  | 25           |      |

- Add 13.5 uL of master mix to the 11.5 uL of DNA sample
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

Run a 2% gel at 80-100V and 100 amps for 2 hours. Use the 1 kb DNA ladder and 100 bp ladder to each row. Today, I used one small gel and one medium gel. For some reason, the voltage will not go past ~80 volts. 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/gel_20240514.JPG)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/gel_smol_20240514.JPG)

Both gels are really faint, along with the ladders...I added ~3-4 uL of sample/ladder to each well. Maybe I didn't add enough volume to each well? 

I reran the gels on 5/15 and they looked much better. 

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/gel_med_20240515.JPG)

![](https://raw.githubusercontent.com/JillAshey/JillAshey_Putnam_Lab_Notebook/master/images/gel_small_20240515.JPG)

The bands to the left and right of the numbers are the same sample as the number, as each sample was run in triplicate. 
