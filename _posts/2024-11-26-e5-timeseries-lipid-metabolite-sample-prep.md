---
layout: post
title: Metabolite and lipid sample prep
date: '2024-11-26'
categories:
tags: [Protocols, Lipids, Metabolites]
projects: e5
---

# e5 sample prep for lipid and metabolite analysis

To characterize lipids and metabolites in the e5 timeseries data, we have to homogenize and separate the samples before sending them to the [UW facility](https://northwestmetabolomicsorg.wpcomstaging.com/). We first sent them 4 test samples to see if these samples would work. The lipids worked great but the metabolite data was lacking some key metabolites --  we may not have given them enough sample as input and the instruments may have gotten disrupted by salt since we used 10x PBS to airbrush the samples. See the [former protocol](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-08-29-e5-Airbrushing-Metabolomics-Test-Samples.md) e5 timeseries molecular [github](https://github.com/urol-e5/timeseries_molecular/tree/main/M-multi-species) for information about test samples. 

To remedy these issues, Ariana and I developed the following [protocol](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/E5-metabolomics-lipidomics-tissue-preps/), which I followed for all samples unless otherwise noted. 

The major changes of the protocol: 

- Use LCMS grade water for airbrushing 
- Airbrush fragment primarily with air and remove any residual tissue with water + airbrush, which should lead to a more concentrated sample 
- Collect homogenate from each sample for protein assay
- Run protein assay to quantify the ug/mL of protein in tissue 

### Sample prep 

Below are the details for the samples I prepped. There are 96 samples total. I will continue to add to this table as I complete each sample. 

| Sample  | Timepoint | Species                  | Site            | Frozen Frag | Fragment bag (M1 or M2) | Bag notes                                                                      | Jill finding note | Date of airbrushing | Airbrush notes                                                          | Serial number | Metabolite volume (uL) | Lipid volume (uL) | Backup volume (uL) | Protein volume (uL) | Total volume (mL) |
| ------- | --------- | ------------------------ | --------------- | ----------- | ----------------------- | ------------------------------------------------------------------------------ | ----------------- | ------------------- | ----------------------------------------------------------------------- | ------------- | ---------------------- | ----------------- | ------------------ | ------------------- | ----------------- |
| POR-73  | TP2       | Porites evermanni        | Manava (site 1) | yes         | green                   | No date on whirlpak, but from March - based on handwriting. Zoe has picture,   |                   | 20241126 ?          | Could be from TP1 or TP2 - need to check handwriting and photos of bags | 1             | 1000                   | 1000              | 1000               | 1000                | 6                 |
| POR-73  | TP1       | Porites evermanni        | Manava (site 1) | yes         | orange                  | From "undated" January bag - Zoe has picture                                   |                   | 20241126 ?          | Could be from TP1 or TP2 - need to check handwriting and photos of bags | 1             | 1000                   | 1000              | 1000               | 1000                | 6                 |
| POC-42  | TP3       | Pocillopora tuahiniensis | Manava (site 1) | yes         | pink                    |                                                                                |                   | 20241126            |                                                                         | 2             | 1000                   | 1000              | 1000               | 1000                | 5                 |
| ACR-145 | TP4       | Acropora pulchra         | Manava (site 1) | yes         | grey                    |                                                                                |                   | 20241126            |                                                                         | 3             | 1000                   | 1000              | 1000               | 1000                | 5.5               |
| POC-53  | TP3       | Pocillopora tuahiniensis | Manava (site 1) | yes         | grey                    |                                                                                |                   | 20241126            |                                                                         | 4             | 1000                   | 1000              | 1000               | 1000                | 4                 |
| POC-52  | TP1       | Pocillopora tuahiniensis | Manava (site 1) | yes         | orange                  | From "undated" January bag - Zoe has picture                                   |                   | 20241126            |                                                                         | 5             | 1000                   | 1000              | 1000               | 1000                | 4                 |
| POC-259 | TP4       | Pocillopora tuahiniensis | Mahana (site 2) | yes         | pink                    |                                                                                |                   | 20241126            |                                                                         | 6             | 1000                   | 1000              | 1000               | 1000                | 5                 |
| POC-42  | TP4       | Pocillopora tuahiniensis | Manava (site 1) | yes         | blue                    |                                                                                |                   | 20241126            |                                                                         | 7             | 1000                   | 1000              | 1000               | 1000                | 9                 |
| POC-53  | TP1       | Pocillopora tuahiniensis | Manava (site 1) | yes         | orange                  | No date on whirlpak, but from January - based on handwriting. Zoe has picture, |                   | 20241126            |                                                                         | 8             | 1000                   | 1000              | 1000               | 1000                | 4                 |
| POR-216 | TP1       | Porites evermanni        | Mahana (site 2) | yes         | green                   |                                                                                |                   | 20241126            |                                                                         | 9             | 1000                   | 1000              | 1000               | 1000                | 9                 |
| POR-72  | TP3       | Porites evermanni        | Manava (site 1) | yes         | grey                    |                                                                                |                   | 20241126            |                                                                         | 10            | 1000                   | 1000              | 1000               | 1000                | 4                 |
| ACR-237 | TP3       | Acropora pulchra         | Mahana (site 2) | yes         | green                   |                                                                                |                   | 20241126            |                                                                         | 11            | 1000                   | 1000              | 1000               | 500                 | 3.5               |
| POR-216 | TP3       | Porites evermanni        | Mahana (site 2) | yes         | green                   |                                                                                |                   | 20241126            |                                                                         | 12            | 1000                   | 1000              | 1000               | 500                 | 3.5               |

In the table, the serial number is the number on the tube that will be sent to UW for analysis. 