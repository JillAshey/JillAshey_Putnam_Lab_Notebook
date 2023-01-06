---
layout: post
title: SRA uploads to NCBI
date: '2022-12-07'
categories: Protocol
tags: [Protocol, Bioinformatics]
projects: Sediment Stress
---

## SRA uploads to NCBI

This post details the NCBI Sequence Read Archive upload for J. Ashey's sediment stress project. This protocol is based on A. Huffmyer's [SRA upload post](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/SRA-Uploads-10-November-2022/) and Putnam Lab [SRA upload protocol](https://github.com/Putnam-Lab/Lab_Management/blob/master/Bioinformatics_%26_Coding/Data_Mangament/SRA-Upload_Protocol.md). 

### Overview 

These sequences are from a project assessing sedimentation stress in Caribbean and Pacific corals. The project's github can be found [here](https://github.com/JillAshey/SedimentStress).

Pacific corals (*Montipora capitata*, *Pocillopora acuta*, and *Porites lobata*) from Kāneʻohe Bay, Oʻahu, Hawaiʻi, were exposed to unsterilized terrigenous red soil for up to 7 days. Caribbean corals (*Acropora cervicornis*, *Montestraea cavernosa*, and *Orbicella faveolata*) from Key West, Florida, were exposed to sterilized white carbonate sediment for 18 days. After each experiment, fragments from each species were frozen and stored at -80°C until extraction.

RNA extractions were done with the Direct-zol™ RNA MiniPrep  (Zymo Research; Cat# R2070) kit. For library prep, samples were diluted and processed following the TruSeq stranded mRNA Library Prep for NeoPrep kit (Document # 15049725 v03, Illumina) protocol. Quality controlled libraries were sequenced on HiSeq 50 cycle single read sequencing v4 by the High Throughput Genomics Core Facility at the University of Utah.

#### 1. BioProject

I created a new submission on [NCBI Submission Portal](https://submit.ncbi.nlm.nih.gov/subs/bioproject/) for a new BioProject.

- Provided all my info as the submitter 
- Project type 
	- Project data type is raw sequence reads 
	- Sample scope is multispecies
- Target 
	- Not putting organism name because I am submitting reads for multiple species
	- Under multispecies description, add species names from the study (*Montipora capitata*, *Pocillopora acuta*, *Porites lobata*, *Acropora cervicornis*, *Montestraea cavernosa*, and *Orbicella faveolata*)
- General info
	- Release date set as Dec. 20, 2022 to allow time for edits 
	- Project title is "Characterizing transcriptomic responses to sedimentation across location and morphology in reef-building corals"
	- Project description is "Gene expression in response to sedimentation across location (Florida, Hawai'i) and morphology (branching, intermediate, massive). Data includes RNAseq (gene expression) sequences from 6 reef-building coral species (*Montipora capitata*, *Pocillopora acuta*, *Porites lobata*, *Acropora cervicornis*, *Montestraea cavernosa*, and *Orbicella faveolata*)"
	- Project is not part of a larger initiative already registered with NCBI
	- Grants associated with this project: 
		- # 1939795 and # 1939263 Harnessing the Data Revolution; National Science Foundation 

Submitted at 2:18PM 20221213 to NCBI. Submission ID is SUB12414011; BioProject ID is PRJNA911752.

2. BioSamples 

Using the Invertebrate attribute table b/c I'm uploading adult coral samples. I deleted some of the columns that I wasn't using. 

Trying to submit Invertebrate attribute file, but I keep getting errors. This is an example of the rows in my attribute table: 

| \*sample\_name    | sample\_title | bioproject\_accession | \*organism         | isolate    | breed          | host           | isolation\_source | \*collection\_date | \*geo\_loc\_name                        | \*tissue       | collected\_by   | dev\_stage | env\_broad\_scale            | host\_tissue\_sampled | identified\_by     | lat\_lon               | description                                                                                                                                                                                                                                                                                                                              | Treatment              | SequencingID | SedimentType                     | AnalyzedBy | Permit No.      | Grants                               |
| ----------------- | ------------- | --------------------- | ------------------ | ---------- | -------------- | -------------- | ----------------- | ------------------ | --------------------------------------- | -------------- | --------------- | ---------- | ---------------------------- | --------------------- | ------------------ | ---------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------- | ------------ | -------------------------------- | ---------- | --------------- | ------------------------------------ |
| 17\_ctl2\_Of\_ZTH | JA1           | PRJNA911752           | Obicella faveolata | Coral host | Not applicable | Not applicable | Not applicable    | 9-Jun-16           | Florida Keys; USA: Florida: Fort Pierce | Whole organism | Francois Seneca | Adult      | coral reef \[ENVO:00000150\] | Whole host organism   | Ellis and Solander | 27.460350N, -80.311073 | Adult life stage sample of Orbicella faveolata collected from Key West National Marine Sanctuary, Florida coral reef environment. Individuals were transported to Fort Pierce, Florida, for exposure to sterilized white carbonate sediment for 18 days. Samples were preserved in liquid nitrogen and stored at -80°C until processing. | Control                | 14004RX1     | Sterilized coral rubble sediment | Jill Ashey | #FKNMS-2016-017 | NSF HDR Awards #1939795 and #1939263 |
| 18\_T3.3\_Of\_VLL | JA2           | PRJNA911752           | Obicella faveolata | Coral host | Not applicable | Not applicable | Not applicable    | 9-Jun-16           | Florida Keys; USA: Florida: Fort Pierce | Whole organism | Francois Seneca | Adult      | coral reef \[ENVO:00000150\] | Whole host organism   | Ellis and Solander | 27.460350N, -80.311073 | Adult life stage sample of Orbicella faveolata collected from Key West National Marine Sanctuary, Florida coral reef environment. Individuals were transported to Fort Pierce, Florida, for exposure to sterilized white carbonate sediment for 18 days. Samples were preserved in liquid nitrogen and stored at -80°C until processing. | Treatment 3 (300 mg/L) | 14005RX1     | Sterilized coral rubble sediment | Jill Ashey | #FKNMS-2016-017 | NSF HDR Awards #1939795 and #1939263 |


I kept the green fields (mandatory), filled out some blue fields or put "Not applicable", and filled in some of the yellow (optional) fields. I deleted the yellow columns that I wasn't using. I also added some of my own attribute columns that are specific to my experiment: Treatment, SequencingID, SedimentType, AnalyzedBy, Permit No., Grants

My attribute file is [here](https://github.com/JillAshey/SedimentStress/blob/master/Data/NCBI_upload/Invertebrate.1.0-2.xlsx). When I try to upload it to the Attributes page on the BioSample submission page, it gives me this **error**: 

"Your table upload failed because multiple BioSamples cannot have identical attributes. You should have one BioSample for each specimen, and each of your BioSamples must have differentiating information (excluding sample name, title, bioproject accession and description). This check was implemented to encourage submitters to include distinguishing information in their samples. If the distinguishing information is in the sample name, title or description, please recode it into an appropriate attribute, either one of the predefined attributes or a custom attribute you define. If it is necessary to represent true biological replicates as separate BioSamples, you might add an 'aliquot' or 'replicate' attribute, e.g., 'replicate = biological replicate 1', as appropriate. Note that multiple assay types, e.g., RNA-seq and ChIP-seq data may reference the same BioSample if appropriate." And then it gives me a list of sample names.

I'm not sure what this means?? I deleted the unused columns, tried relabeling the sample names, added unique letters to the end of the sample name, etc. This is the confusing part: "**You should have one BioSample for each specimen, and each of your BioSamples must have differentiating information (excluding sample name, title, bioproject accession and description)**". ??????????????

20230106
Seemed to have fixed this problem while zooming with Ariana in December. Issue was mostly with formatting, as Excel sometimes applies its own formatting to the data. I had to redo some formatting things (dates, locations) to make sure that it could be properly be uploaded to NCBI. This [link](https://www.ncbi.nlm.nih.gov/biosample/docs/attributes/) is a good overview of the possible BioSample Attributes and how to format them in the Excel sheet. 

This is an updated example of the rows in my attribute table:

| \*sample_name  | sample_title | bioproject_accession | \*organism         | isolate    | breed          | host           | isolation_source | \*collection_date | \*geo_loc_name                        | \*tissue       | collected_by    | dev_stage | env_broad_scale            | host_tissue_sampled | identified_by      | lat_lon                   | description                                                                                                                                                                                                                                                                                                                              | Treatment              | SequencingID | medium                           | AnalyzedBy | Permit No.      | Grants                               | Identifier     |
| -------------- | ------------ | -------------------- | ------------------ | ---------- | -------------- | -------------- | ---------------- | ----------------- | ------------------------------------- | -------------- | --------------- | --------- | -------------------------- | ------------------- | ------------------ | ------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------- | ------------ | -------------------------------- | ---------- | --------------- | ------------------------------------ | -------------- |
| 17_ctl2_Of_ZTH | JA1          | PRJNA911752          | Obicella faveolata | Coral host | Not applicable | Not applicable | Not applicable   | 9-Jun-16          | USA: Florida Keys Fort Pierce Florida | Whole organism | Francois Seneca | Adult     | coral reef [ENVO:00000150] | Whole host organism | Ellis and Solander | 27.460350 N, -80.311073 W | Adult life stage sample of Orbicella faveolata collected from Key West National Marine Sanctuary, Florida coral reef environment. Individuals were transported to Fort Pierce, Florida, for exposure to sterilized white carbonate sediment for 18 days. Samples were preserved in liquid nitrogen and stored at -80°C until processing. | Control                | 14004RX1     | Sterilized coral rubble sediment | Jill Ashey | #FKNMS-2016-017 | NSF HDR Awards #1939795 and #1939263 | 17_ctl2_Of_ZTH |
| 18_T3.3_Of_VLL | JA2          | PRJNA911752          | Obicella faveolata | Coral host | Not applicable | Not applicable | Not applicable   | 9-Jun-16          | USA: Florida Keys Fort Pierce Florida | Whole organism | Francois Seneca | Adult     | coral reef [ENVO:00000150] | Whole host organism | Ellis and Solander | 27.460350 N, -80.311073 W | Adult life stage sample of Orbicella faveolata collected from Key West National Marine Sanctuary, Florida coral reef environment. Individuals were transported to Fort Pierce, Florida, for exposure to sterilized white carbonate sediment for 18 days. Samples were preserved in liquid nitrogen and stored at -80°C until processing. | Treatment 3 (300 mg/L) | 14005RX1     | Sterilized coral rubble sediment | Jill Ashey | #FKNMS-2016-017 | NSF HDR Awards #1939795 and #1939263 | 18_T3.3_Of_VLL |

Submitted BioSamples at 12:05 on 1/6/23 under submission number SUB12414115.
 

