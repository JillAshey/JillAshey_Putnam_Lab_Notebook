---
layout: post
title: Hawaii 2023 daily notebook posts 
date: '2023-06-15'
categories: Field work, Mcap
tags: [Field work, Hawaii, Montipora capitata, Mcap]
projects: Developmental timeseries - HI 2023
---

This post includes the full notebook for the June-July 2023 *Montipora capitata* coral spawning and field expedition in Hawaii at the Hawaii Institute of Marine Biology. My project will focus on development of *M capitata* larvae through their use of energetic stores and expression of coding and non-coding transcripts. Dr. Ariana Huffmyer is also conducting experiments on *M. capitata* larvae; her experiments are detailed in her Github notebook [here](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/Hawaii-2023-Coral-Spawning-and-Field-Expedition-Daily-Entry/). Github for the developmental timeseries is [here](https://github.com/JillAshey/Hawaii_Larval_TimeSeries/tree/main). 

## 20230615

Arrived! Slowly but surely. 

## 20230616

### Set-up 

Supplies and equipment were unpacked and lab space was set up. AH filled the dry shipper up with liquid nitrogen at UH Manoa to act as a dewar in the field. Incubators were unpacked and set up. Incubators were set to ambient (27°C) and high (30°C) temperature treatments.

### Tube prep & labeling 

JA prepped 50 mL falcon tubes and 1.5 mL tubes for sampling. The 50 mL falcon tubes will be used to hold the larvae in the incubators. The 1.5 mL tubes will be used for sampling for physiology (carbs, lipids, protein), molecular (mRNA and ncRNA) and size. 1.5 mL tubes were labeled in the following way: 

| Color  | Category                     | Number Range    | Buffer/preservation | Storage                                   | Destination |
| ------ | ---------------------------- | --------------- | ------------------- | ----------------------------------------- | ----------- |
| Blue   | Size       | S4-S94 (no S16) | 20% Z-fix           | 4C; transport room temp (or on ice packs) | URI         |
| Orange | Molecular  | M1-M94 (no M16) | RNA/DNA shield      | \-80C; transport ln2                | URI         |
| Green  | Physiology | P1-P94 (no P16) | Snap freeze         | \-80C; transport lN2                      | URI         |

### Logger deployment 

2 loggers were deployed in each incubator. We found that there was variation in the temperatures on the different shelves of the incubators, with the top shelf being the warmest and the lowest shelf the coolest. The incubators were filled with 50 mL falcon tubes at this time to test airflow. We will need to reduce the number of falcon tubes in the incubators for the developmental time series experiment. Below is the metadata for each logger: 

| SerialNumber | Type               | DateDeployed | Incubator |
| ------------ | ------------------ | ------------ | --------- |
| 21723476     | HOBO Tidbit MX2203 | 20230616     | 1         |
| 21723475     | HOBO Tidbit MX2203 | 20230616     | 1         |
| 21723473     | HOBO Tidbit MX2203 | 20230616     | 2         |
| 21723474     | HOBO Tidbit MX2203 | 20230616     | 2         |
| 21723472     | HOBO Tidbit MX2203 | 20230616     | 3         |
| 21723477     | HOBO Tidbit MX2203 | 20230616     | 4         |
| 20987156     | HOBO U22           | 20230617     | 3         |
| 20444042     | HOBO U22           | 20230617     | 4         |

### Spawning 

The CRL team went to Reef 13 in Kaneohe Bay to collect gametes from wildtype spawning slicks. They reported that spawning was light and they were able to collect some gametes for preliminary testing. Two falcon tubes of gametes (5 mL gametes per tube in 40 mL FSW) were used for testing. Tubes incubated on the bench top until all bundles were broken up around 23:00. Tubes were rocked every 10-15 minutes to facilitate bundle breakage. After the majority of bundles were broken up, the fertilized eggs were rinsed 3x with FSW to clean the sperm. 

### Developmental time series testing

Since spawning was light, AH and JA decided to run a test for the developmental time series experiment to see if the embryos could survive in the 50 mL falcon tubes with minimal water changes. Stocking density was also assessed. 

To calculate the concentration of eggs, 5 uL was taken from the concentrated egg stock and counted under the micrscope (average of 76 eggs in 5 uL concentrated stock). Based on this estimate, tubes were stocked with either 10, 50, 100, 200 or 500 uL of the concentrated egg stock (n=3 tubes per density per temperature treatment) in 40-45 mL FSW. Stocking took place from 22:50-23:15. Below is the metadata from the prelim 50 mL falcon tubes: 

| FalconID | DateLoaded | Volume_uL | Treatment |
| -------- | ---------- | --------- | --------- |
| 30       | 20230616   | 10        | Ambient   |
| 43       | 20230616   | 10        | Ambient   |
| 46       | 20230616   | 10        | Ambient   |
| 50       | 20230616   | 50        | Ambient   |
| 41       | 20230616   | 50        | Ambient   |
| 55       | 20230616   | 50        | Ambient   |
| 65       | 20230616   | 100       | Ambient   |
| 154      | 20230616   | 100       | Ambient   |
| 94       | 20230616   | 100       | Ambient   |
| 21       | 20230616   | 200       | Ambient   |
| 176      | 20230616   | 200       | Ambient   |
| 88       | 20230616   | 200       | Ambient   |
| 104      | 20230616   | 500       | Ambient   |
| 191      | 20230616   | 500       | Ambient   |
| 162      | 20230616   | 10        | High      |
| 143      | 20230616   | 10        | High      |
| 109      | 20230616   | 10        | High      |
| 139      | 20230616   | 50        | High      |
| 140      | 20230616   | 50        | High      |
| 142      | 20230616   | 50        | High      |
| 148      | 20230616   | 100       | High      |
| 130      | 20230616   | 100       | High      |
| 170      | 20230616   | 100       | High      |
| 181      | 20230616   | 200       | High      |
| 164      | 20230616   | 200       | High      |
| 105      | 20230616   | 200       | High      |
| 63       | 20230616   | 500       | High      |
| 224      | 20230616   | 500       | High      |

## 20230617

### Prelim tube checks 

at 03:00: 

- Ambient has some clumping in all tubes, but more so in the higher density tubes.
- All in high treatment look like they have some clumping and death. Water is more opaque in these tubes than the ambient tubes. Also less evidence of fertilization in high treatments. Less debris and water is cloudier. 
- Rocking tubes every 30-60 minutes will help with clumping. As for high treatment, maybe JA can use a glass pipette to remove debris periodically? Don't want to do frequent water changes, as I don't want to disturb them too much 

At 08:00:

- Water changes performed on all tubes. Water changes performed by pouring contents through a cell strainer, refilling the tube with fresh FSW and then adding the embryos back in. 
- Put filled tubes in middle rack in incubator.
- Changed temps to be set at 26°C for ambient and 29°C for high. 
- Calculated out what I'll need for my actual sampling: 6 time points x 2 treatments x 7 replicates = 84 tubes

### Sampling prelim test

JA did a sampling run-through to prep for actual sampling. She decided to use one 50 mL tube to sample for all 3 metrics (physiology, molecular and size). Here's the sampling protocol: 

- Remove a 50 mL tube from the incubator and record the tube ID and time it was removed from the incubator. 
- Pour the contents through a cell strainer and set the cell strainer in a 6 well plate with a little bit of FSW. 
- Sample 20 uL of the concentrated embryos for size. Preserve in 20% Zfix. Note the time that this sample was preserved. Store at 4°C.
- Use a P1000 to remove the rest of the embryos. Record how much volume total it took to collect the embryos. Put these embryos in a 20 mL glass vial. 
- Subset 20 uL of the embryos and count how many embryos are in that volume. Record this number. Using the volume of the embryos, the total number of embryos in the 50 mL falcon tube can be calculated. 
- Move half of the total volume to the molecular tube and the other half to the physiology tube. This should split the embryos evenly between the physiology and molecular sample tubes. 
- Preserve the molecular tube by adding 500 uL of DNA/RNA shield. Record the time preserved. Store at -80°C.
- Preserve the physiology tube by snap freezing in liquid nitrogen. Record the time preserved. Store at -80°C.
- Sampling complete for one falcon tube! Repeat for all other falcon tubes for that sampling time point. 

### Spawning (part 2)

Much more spawning tonight! CRL team brought back 6 50 mL falcon tubes (5 mL gametes in 40 mL FSW) at 21:45. Tubes were left on the benchtop to fertilize and rocked gently every 5 minutes. After an hour, the fertilized eggs were rinsed 2x with FSW. 

### Loading tubes 

Decided to load the 50 mL falcon tubes with 200 uL stocking density. The embryos seemed to survive well at this density and it will provide a sufficient number of larvae to sample from each tube. 200 uL of concentrated eggs were loaded into 50 mL tubes with 40-45 FSW.  Tubes were inverted gently and loaded into an incubator with the tubes laying on their sides. Incubators were set to either 26°C or 29°C. 24 tubes with 200 uL of concentrated eggs were loaded into each incubator (4 incubators total; 2 ambient, 2 high). 1 tube per incubator was loaded with no eggs to serve as a "blank". 

| Incubator | Treatment | Time Tubes Loaded | Tubes in Incubator |
| --------- | --------- | ----------------- | ------------------ |
| 1         | Ambient   | 22:50             | 23:05              |
| 2         | High      | 23:39             | 23:53              |
| 3         | Ambient   | 23:05             | 23:20              |
| 4         | High      | 23:25             | 23:39              |

### Sampling 

Sampling of the fertilized egg stock was done at 23:30 for physiology, molecular and size (all at ambient). Sampling was done according to the protocol above.  

## 20230618

### Tube upkeep 

After being loaded into their respective incubators, tubes were gently rocked every hour to prevent embryos from clumping. This helped reduce die-off and kept the embryos developing at a regular pace. One 1L Nalgene of FSW was kept in each incubator for water changes so that the water was already at the treatment temperatures. Water changes on all tubes were done at 06:30 and 16:45. Water quality measurements (temperature, pH, salinity) was taken at 00:31, 05:45, and 16:20 on 3 tubes per incubator. Water quality measurements for today are below: 

| Date     | Time  | FalconID | pH (mV) | pH_nbs | temp (C) | sal (psu) | Tris date | Notes                                               |
| -------- | ----- | -------- | ------- | ------ | -------- | --------- | --------- | --------------------------------------------------- |
| 20230618 | 0:31  | Nalgene1 | \-52.8  | 8.08   | 23.87    | 33        | 20230617  | Nalgene1 was sitting on bench prior to measurements |
| 20230618 | 0:34  | Nalgene2 | \-51.1  | 8.05   | 28.3     | 33        | 20230617  |                                                     |
| 20230618 | 0:37  | Nalgene3 | \-51.2  | 8.05   | 28.51    | 33        | 20230617  |                                                     |
| 20230618 | 0:39  | Nalgene4 | \-52.3  | 8.07   | 25.5     | 33        | 20230617  |                                                     |
| 20230618 | 5:46  | 22       | \-45.2  | 7.94   | 24.6     | 34        | 20230617  |                                                     |
| 20230618 | 5:47  | 184      | \-45.5  | 7.95   | 24.51    | 34        | 20230617  |                                                     |
| 20230618 | 5:49  | 69       | \-47.7  | 7.99   | 25       | 34        | 20230617  |                                                     |
| 20230618 | 5:51  | 81       | \-46.6  | 7.97   | 29.15    | 34        | 20230617  |                                                     |
| 20230618 | 5:53  | 53       | \-45    | 7.94   | 29.92    | 34        | 20230617  |                                                     |
| 20230618 | 5:54  | 18       | \-41.8  | 7.89   | 29.53    | 34        | 20230617  |                                                     |
| 20230618 | 16:19 | 87       | \-30.4  | 7.69   | 28.37    | 34        | 20230617  |                                                     |
| 20230618 | 16:20 | 75       | \-6.4   | 7.27   | 29.34    | 34        | 20230617  |                                                     |
| 20230618 | 16:22 | 27       | \-23.5  | 7.57   | 29.25    | 34        | 20230617  |                                                     |
| 20230618 | 16:24 | 76       | \-40.2  | 7.86   | 28.03    | 34        | 20230617  |                                                     |
| 20230618 | 16:26 | 51       | \-38.9  | 7.84   | 28.51    | 34        | 20230617  |                                                     |
| 20230618 | 16:28 | 202      | \-40.2  | 7.86   | 27.51    | 34        | 20230617  |                                                     |
| 20230618 | 16:30 | 52       | \-36.8  | 7.8    | 25.23    | 34        | 20230617  |                                                     |
| 20230618 | 16:32 | 70       | \-28.1  | 7.65   | 25.18    | 34        | 20230617  |                                                     |
| 20230618 | 16:33 | 71       | \-41.2  | 7.88   | 25.33    | 34        | 20230617  |                                                     |
| 20230618 | 16:35 | 36       | \-32.3  | 7.72   | 25.35    | 34        | 20230617  |                                                     |
| 20230618 | 16:37 | Nalgene1 | \-50.1  | 8.03   | 24.61    | 34        | 20230617  |                                                     |
| 20230618 | 16:38 | Nalgene2 | \-47.1  | 7.98   | 28       | 34        | 20230617  |                                                     |
| 20230618 | 16:39 | Nalgene3 | \-45.6  | 7.95   | 28.8     | 34        | 20230617  |                                                     |
| 20230618 | 16:40 | Nalgene4 | \-47.9  | 7.99   | 25.14    | 34        | 20230617  |

### Sampling

Sampling took place at several time points today: 4, 9, 14, and 22 hours post fertilization (hpf). Sampling was done as detailed above. At each time point, n=6 falcon tubes were sampled from each treatment (n=3 falcon tubes per incubator). Total embryos were calculated for each tube, but this number is approximate. Some molecular samples (tubes M10, M11, M12, M17 and M18) may not be viable. In these tubes, there was already shield in the tubes when adding the sample. 

## 20230619

### Tube upkeep 

Tubes were checked every 3-6 hours and rocked gently to prevent build up. Embryos are looking good! JA periodically took an aliquot out of the tubes to check development under the microscope. Water changes occurred at 08:00. Water quality data was also collected: 

| Date     | Time  | FalconID | pH (mV) | pH_nbs | temp (C) | sal (psu) | Tris date | Notes                                               |
| -------- | ----- | -------- | ------- | ------ | -------- | --------- | --------- | --------------------------------------------------- |                                                    |
| 20230619 | 22:40 | 14       | \-45    | 7.94   | 26.08    | 34        | 20230617  |                                                     |
| 20230619 | 22:41 | 39       | \-35.3  | 7.78   | 25.64    | 34        | 20230617  |                                                     |
| 20230619 | 22:42 | 47       | \-45.5  | 7.95   | 28.59    | 34        | 20230617  |                                                     |
| 20230619 | 22:43 | 89       | \-36.2  | 7.79   | 28.68    | 34        | 20230617  |                                                     |
| 20230619 | 22:45 | 44       | \-43.1  | 7.91   | 29.63    | 34        | 20230617  |                                                     |
| 20230619 | 22:46 | 75       | \-38.2  | 7.83   | 29.54    | 34        | 20230617  |                                                     |
| 20230619 | 22:48 | 72       | 47.8    | 7.99   | 24.92    | 34        | 20230617  |                                                     |
| 20230619 | 22:49 | 52       | \-35.9  | 7.79   | 25.52    | 34        | 20230617  |

### Sampling

Sampling took place at several time points today: 28 and 48 hpf. Sampling was done as detailed above. At each time point, n=6 falcon tubes were sampled from each treatment (n=3 falcon tubes per incubator). Total embryos were calculated for each tube, but this number is approximate. 

## 20230620

### Tube upkeep 

Tubes were checked every 3-6 hours and rocked gently to prevent build up. Embryos are looking good! JA periodically took an aliquot out of the tubes to check development under the microscope. 

### Sampling

Final sampling today!! Occurred at 72 hpf. Sampling was done as detailed above. At each time point, n=6 falcon tubes were sampled from each treatment (n=3 falcon tubes per incubator). Total embryos were calculated for each tube, but this number is approximate. 

The following table has the hpf, as well as the start and end times for each sampling time point. 

| Date     | hpf            | TimeStart | TimeStop | NumberOfSamples |
| -------- | -------------- | --------- | -------- | --------------- |
| 20230617 | Fertilized egg | 23:30     | 23:39    | 18              |
| 20230618 | 4              | 2:45      | 4:20     | 36              |
| 20230618 | 9              | 8:41      | 10:38    | 36              |
| 20230618 | 14             | 13:53     | 15:51    | 36              |
| 20230618 | 22             | 20:57     | 22:54    | 36              |
| 20230619 | 28             | 2:15      | 4:05     | 36              |
| 20230619 | 48             | 21:36     | 23:21    | 36              |
| 20230620 | 72             | 21:26     | 23:25    | 36              |

Each time point took me about ~2 hours to sample, so for each sample, its not exactly 4 or 9 hpf. I alternated sampling between high and ambient treatments, so there should be no effect of time on treatment if that makes sense. 

Hooray! Sampling for my project complete! See A. Huffmyer's notebook [posts](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/Hawaii-2023-Coral-Spawning-and-Field-Expedition-Daily-Entry/) and project [github](https://github.com/AHuffmyer/larval_symbiont_TPC) for the other experiments that were completed during this trip. 
