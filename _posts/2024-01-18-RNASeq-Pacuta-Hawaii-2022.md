---
layout: post
title: Pacuta HI 2022
date: '2024-01-18'
categories: Analysis
tags: [Bioinformatics, Pacuta, mRNA]
projects: Pacuta HI 2022
---

## Pacuta 2022 mRNA analysis

These data came from the Pacuta 2022 experiment in Hawaii, done by Federica and myself. In this experiment, larval and spat *Pocillopora acuta* were subjected to a combination of high pH and temperature treatments. The github for that project is [here](https://github.com/fscucchia/Hawaii2022_pH_Temp_Mcap_Pacu). 

Files were downloaded to this location: `/data/putnamlab/KITT/hputnam/20231127_Scucchia_HI`

### 20240118

Make a new directory in my own folder on Andromeda for this experiment

```
cd /data/putnamlab/jillashey
mkdir Pacuta_HI_2022
cd Pacuta_HI_2022
mkdir data scripts output
cd data 
mkdir raw trim 
cd raw
```

Copy the files into the raw data folder 

```
cp /data/putnamlab/KITT/hputnam/20231127_Scucchia_HI/* .
```

Check md has already been done by Hollie. Let's count how many reads each file has. 