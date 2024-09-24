---
layout: post
title: Qubit e5 samples 
date: '2024-09-23'
categories:
tags: [DNA, Extractions, e5, Protocols]
projects: e5
---

# Qubit w/ e5 timeseries DNA samples  

We are sending out the e5 timeseries DNA samples for WBGS but some samples still have low concentrations. There are ~20 samples that have <200 ng of DNA (see spreadsheet [here](https://docs.google.com/spreadsheets/d/1iFsVfp1vix9IfNcqSajQfLgXXDg0CCPT/edit?gid=1140341238#gid=1140341238); mostly POR). There are 10 samples that had too low concentration for Qubit. I'm going to re-Qubit these 20 samples with the HS Qubit kit to get a better idea of concentrations. 

**Note**: All samples except for 20240805_POC-52_TP1 are POR. 
 
| SampleName          | QC_Number | DNA1_ng_uL | DNA2_ng_uL | DNA_Average | DNA_uL | DNA_ng  |
| ------------------- | --------- | ---------- | ---------- | ----------- | ------ | ------- |
| 7-20220208          | 1         | 0.0604     | 0.0632     | 0.0618      | 82     | 5.0676  |
| 195                 | 2         | 0.0632     | 0.0624     | 0.0628      | 85     | 5.338   |
| 749                 | 3         | 0.448      | 0.287      | 0.3675      | 85     | 31.2375 |
| 557                 | 4         | 0.116      | 0.0892     | 0.1026      | 85     | 8.721   |
| 637                 | 5         | 0.305      | 0.293      | 0.299       | 85     | 25.415  |
| 247-20220208        | 6         | 0.102      | 0.0912     | 0.0966      | 85     | 8.211   |
| 13                  | 7         | 0.121      | 0.106      | 0.1135      | 82     | 9.307   |
| 469-20211122        | 8         | 0.334      | 0.334      | 0.334       | 82     | 27.388  |
| 251-20220208        | 9         | 0.118      | 0.104      | 0.111       | 85     | 9.435   |
| 20240805_POC-52_TP1 | 10        | 0.11       | 0.101      | 0.1055      | 85     | 8.9675  |
| 235                 | 11        | 1.4        | 1.4        | 1.4         | 82     | 114.8   |
| 23                  | 12        | 0.222      | 0.204      | 0.213       | 85     | 18.105  |
| 711                 | 13        | 0.848      | 0.84       | 0.844       | 85     | 71.74   |
| 327                 | 14        | 0.391      | 0.375      | 0.383       | 85     | 32.555  |
| 491                 | 15        | 0.432      | 0.448      | 0.44        | 85     | 37.4    |
| 669                 | 16        | 0.808      | 0.8        | 0.804       | 85     | 68.34   |
| 477                 | 17        | 0.272      | 0.272      | 0.272       | 85     | 23.12   |
| 493                 | 18        | 0.508      | 0.5        | 0.504       | 85     | 42.84   |
| 515-20211104        | 19        | 0.324      | 0.323      | 0.3235      | 85     | 27.4975 |
| 529                 | 20        | 0.308      | 0.304      | 0.306       | 85     | 26.01   |

Not great. 7 samples have <10 ng of DNA, and 19 samples have <100 ng of DNA. Azenta said: "Our lab prefers to have at least ~300ng of input DNA given our internal workflow optimizations. For samples under 200ng or so, we can proceed at best-efforts. Like I said, we have a lot of experience working with less-than-optimal samples and often see good results with starting materials as low as 10ng. However, if you were able to re-extract and/or pool samples to have more starting material, that would be a good option."

Options moving forward

- Send as-is and hope Azenta can successfully library prep 
- Re-extract samples that have <200ng to see if we can get higher concentrations
- Pool samples to get to higher concentrations 

