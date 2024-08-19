---
layout: post
title: e5 closest
date: '2024-08-19'
categories: Analysis
tags: [Bioinformatics]
projects: e5 deep dive
---

## e5 deep dive - investigating ncRNA protein machinery in 3 species of corals 

The e5 deep dive project is examining ncRNA dynamics in 3 species of coral from Moorea, French Polynesia. The github for the project is [here](https://github.com/urol-e5/deep-dive). 

In this post, I will be assessing what genomic features are closest to the ncRNAs that we have identified through the deep dive. I will need gffs and fasta files for this task. 

##### Important miRNA files

- [Apul miRNA fasta](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/mir.fasta); [Apul miRNA gff](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/output/13.2.1-Apul-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Results.gff3)
- [Ptuh miRNA fasta](https://github.com/urol-e5/deep-dive/blob/main/F-Pmea/output/13.2.1-Pmea-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/mir.fasta); [Ptuh miRNA gff](https://github.com/urol-e5/deep-dive/blob/main/F-Pmea/output/13.2.1-Pmea-sRNAseq-ShortStack-31bp-fastp-merged-cnidarian_miRBase/ShortStack_out/Results.gff3)
- [Peve miRNA fasta](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/output/08.2-Peve-sRNAseq-ShortStack-31bp-fastp-merged/ShortStack_out/mir.fasta); [Peve miRNA gff](https://github.com/urol-e5/deep-dive/blob/main/E-Peve/output/08.2-Peve-sRNAseq-ShortStack-31bp-fastp-merged/ShortStack_out/Results.gff3)

Because I don't own the deep dive repo, I copied all needed files onto my local computer and then put them on Andromeda. I also added the species prefix to each file. 

```
cd /data/putnamlab/jillashey/e5
mkdir ncRNA_gff
cd ncRNA_gff
```

Let's try to convert the gff to bed files for Apul. 

```
awk 'BEGIN {OFS="\t"} {if ($1 !~ /^#/) print $1, $4-1, $5, $3, $6, $7}' Apul_Results.gff3 | sort -k1,1 -k2,2n > Apul_Results.bed
```

Success! Here's what the bed file looks like: 

```
head Apul_Results.bed
NC_058066.1	152482	152910	Unknown_sRNA_locus	140	-
NC_058066.1	161063	161674	Unknown_sRNA_locus	549	.
NC_058066.1	172072	172496	Unknown_sRNA_locus	105	-
NC_058066.1	203241	203651	Unknown_sRNA_locus	100	.
NC_058066.1	204534	205150	Unknown_sRNA_locus	313	.
NC_058066.1	205744	206966	Unknown_sRNA_locus	1930	.
NC_058066.1	210840	211344	Unknown_sRNA_locus	1247	.
NC_058066.1	349654	351297	Unknown_sRNA_locus	3279	+
NC_058066.1	351490	353439	Unknown_sRNA_locus	8889	.
NC_058066.1	598650	599068	Unknown_sRNA_locus	114	+
```

I'm going to remove the unknowns. 

```
awk 'BEGIN {OFS="\t"} $4 != "Unknown_sRNA_locus"' Apul_Results.bed > filtered_Apul_Results.bed

head filtered_Apul_Results.bed 
NC_058066.1	5224791	5225215	siRNA24_locus	1264	+
NC_058066.1	7563376	7563797	siRNA21_locus	118	+
NC_058066.1	8905067	8905484	siRNA22_locus	102	-
NC_058066.1	12757124	12757218	MIRNA_hairpin	8293	-
NC_058066.1	12757146	12757168	mature_miRNA	1413	-
NC_058066.1	12757176	12757198	miRNA-star	85	-
NC_058066.1	20088629	20088720	MIRNA_hairpin	148	+
NC_058066.1	20088649	20088671	miRNA-star	7	+
NC_058066.1	20088678	20088700	mature_miRNA	102	+
NC_058066.1	20346226	20346321	MIRNA_hairpin	249212	-
```

I want to find the genomic features that are closest to the miRNAs. Looks like [bedtools closest](https://bedtools.readthedocs.io/en/latest/content/tools/closest.html) doesn't care if the input is bed or gff. Oh well. 

```
interactive 
module load BEDTools/2.30.0-GCC-11.3.0
bedtools closest -a filtered_Apul_Results.bed -b /data/putnamlab/jillashey/genome/Amil_v2.01/amil_GFFannotation.CDS_sorted.gff > output.tab

***** WARNING: File /data/putnamlab/jillashey/genome/Amil_v2.01/amil_GFFannotation.CDS_sorted.gff has inconsistent naming convention for record:
chr1	maker	CDS	15503	15530	.	+	0	ID=Amillepora00001-RA:cds;Parent=Amillepora00001-RA;

ERROR: Database file /data/putnamlab/jillashey/genome/Amil_v2.01/amil_GFFannotation.CDS_sorted.gff contains chromosome Sc0000015, but the query file does not.
       Please rerun with the -g option for a genome file.
       See documentation for details.

# Retry 
bedtools closest -a filtered_Apul_Results.bed -b /data/putnamlab/jillashey/genome/Amil_v2.01/amil_GFFannotation.CDS_sorted.gff -g /data/putnamlab/jillashey/genome/Amil.v2.01.chrs.fasta -sorted > output.tab
Error: Can't open genome file/data/putnamlab/jillashey/genome/Amil.v2.01.chrs.fastaExiting...

# original code from above 
bedtools closest -a filtered_Apul_Results.bed -b /data/putnamlab/jillashey/genome/Amil_v2.01/amil_GFFannotation.CDS_sorted.gff > output.tab
***** WARNING: File /data/putnamlab/jillashey/genome/Amil_v2.01/amil_GFFannotation.CDS_sorted.gff has inconsistent naming convention for record:
chr1	maker	CDS	15503	15530	.	+	0	ID=Amillepora00001-RA:cds;Parent=Amillepora00001-RA;

ERROR: Database file /data/putnamlab/jillashey/genome/Amil_v2.01/amil_GFFannotation.CDS_sorted.gff contains chromosome Sc0000015, but the query file does not.
       Please rerun with the -g option for a genome file.
       See documentation for details.

wc -l output.tab 
203 output.tab

head output.tab 
NC_058066.1	5224791	5225215	siRNA24_locus	1264	+	.	.	.	-1	-1	.	.	.	.
NC_058066.1	7563376	7563797	siRNA21_locus	118	+	.	.	.	-1	-1	.	.	.	.
NC_058066.1	8905067	8905484	siRNA22_locus	102	-	.	.	.	-1	-1	.	.	.	.
NC_058066.1	12757124	12757218	MIRNA_hairpin	8293	-	.	.	.	-1	-1	.	.	.	.
NC_058066.1	12757146	12757168	mature_miRNA	1413	-	.	.	.	-1	-1	.	.	.	.
NC_058066.1	12757176	12757198	miRNA-star	85	-	.	.	.	-1	-1	.	.	.	.
NC_058066.1	20088629	20088720	MIRNA_hairpin	148	+	.	.	.	-1	-1	.	.	.	.
NC_058066.1	20088649	20088671	miRNA-star	7	+	.	.	.	-1	-1	.	.	.	.
NC_058066.1	20088678	20088700	mature_miRNA	102	+	.	.	.	-1	-1	.	.	.	.
NC_058066.1	20346226	20346321	MIRNA_hairpin	249212	-	.	.	.	-1	-1	.	.	.	.
```

Unsure what this means...its giving me information from the miRNA stuff but no information about nearby features in the genome. Not sure why its not compatible. Trying with Ptuh instead which used the Pmea genome in this project. Need to sort the gff files. Let's do the genome gff file first. 

```
cd /data/putnamlab/jillashey/genome/Pmea
sort -k1,1 -k4,4n Pocillopora_meandrina_HIv1.genes.gff3 > Pocillopora_meandrina_HIv1.genes_sorted.gff3
```

Sort miRNA gff file and run `bed closest`. 

```
cd /data/putnamlab/jillashey/e5/ncRNA_gff

sort -k1,1 -k4,4n Ptuh_Results.gff3 > Ptuh_Results_sorted.gff3

interactive 
module load BEDTools/2.30.0-GCC-11.3.0
bedtools closest -a Ptuh_Results_sorted.gff3 -b /data/putnamlab/jillashey/genome/Pmea/Pocillopora_meandrina_HIv1.genes_sorted.gff3 > Ptuh_output.bed
```

Success!!!

```
wc -l Ptuh_output.bed 
21604 Ptuh_output.bed

head Ptuh_output.bed 
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	Unknown_sRNA_locus	9092	9521	10813	+	.	ID=Cluster_1;DicerCall=N;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	CDS	10771	11117	.	+	0	Parent=Pocillopora_meandrina_HIv1___RNAseq.g20902.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	Unknown_sRNA_locus	9092	9521	10813	+	.	ID=Cluster_1;DicerCall=N;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	exon	10771	11117	.	+	0	Parent=Pocillopora_meandrina_HIv1___RNAseq.g20902.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	Unknown_sRNA_locus	9092	9521	10813	+	.	ID=Cluster_1;DicerCall=N;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	10771	23652	.	+	.	ID=Pocillopora_meandrina_HIv1___RNAseq.g20902.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	Unknown_sRNA_locus	53578	53997	287	+	.	ID=Cluster_2;DicerCall=N;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	52050	53624	.	-	.	ID=Pocillopora_meandrina_HIv1___RNAseq.g20906.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	Unknown_sRNA_locus	53578	53997	287	+	.	ID=Cluster_2;DicerCall=N;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	CDS	53573	53624	.	-	0	Parent=Pocillopora_meandrina_HIv1___RNAseq.g20906.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	Unknown_sRNA_locus	53578	53997	287	+	.	ID=Cluster_2;DicerCall=N;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	exon	53573	53624	.	-	0	Parent=Pocillopora_meandrina_HIv1___RNAseq.g20906.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	Unknown_sRNA_locus	150243	150718	2549	-	.	ID=Cluster_3;DicerCall=N;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	143552	155669	.	-	.	ID=Pocillopora_meandrina_HIv1___RNAseq.g20914.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	Unknown_sRNA_locus	150243	150718	2549	-	.	ID=Cluster_3;DicerCall=N;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	CDS	150290	150371	.	-	1	Parent=Pocillopora_meandrina_HIv1___RNAseq.g20914.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	Unknown_sRNA_locus	150243	150718	2549	-	.	ID=Cluster_3;DicerCall=N;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	exon	150290	150371	.	-	1	Parent=Pocillopora_meandrina_HIv1___RNAseq.g20914.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	Unknown_sRNA_locus	150243	150718	2549	-	.	ID=Cluster_3;DicerCall=N;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	CDS	150573	150661	.	-	0	Parent=Pocillopora_meandrina_HIv1___RNAseq.g20914.t1
```

Remove everything that is unknown. 

```
awk 'BEGIN {OFS="\t"} $3 != "Unknown_sRNA_locus"' Ptuh_output.bed > filtered_Ptuh_output.bed
wc -l filtered_Ptuh_output.bed 
405 filtered_Ptuh_output.bed

head filtered_Ptuh_output.bed
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	siRNA22_locus	173728	174150	1257	+	.	ID=Cluster_4;DicerCall=22;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	CDS	174509	175333	.	-	0	Parent=Pocillopora_meandrina_HIv1___RNAseq.g20918.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	siRNA22_locus	173728	174150	1257	+	.	ID=Cluster_4;DicerCall=22;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	exon	174509	175333	.	-	0	Parent=Pocillopora_meandrina_HIv1___RNAseq.g20918.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	siRNA22_locus	173728	174150	1257	+	.	ID=Cluster_4;DicerCall=22;MIRNA=N	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	174509	176444	.	-	.	ID=Pocillopora_meandrina_HIv1___RNAseq.g20918.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	MIRNA_hairpin	818027	818120	12096	+	.	ID=Cluster_19;DicerCall=23;MIRNA=Y	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	816355	820160	.	+	.	ID=Pocillopora_meandrina_HIv1___RNAseq.g21001.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	mature_miRNA	818049	818070	3240	+	.	ID=Cluster_19.mature;Parent=Cluster_19	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	816355	820160	.	+	.	ID=Pocillopora_meandrina_HIv1___RNAseq.g21001.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	miRNA-star	818079	818100	9	+	.	ID=Cluster_19.star;Parent=Cluster_19	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	816355	820160	.	+	.	ID=Pocillopora_meandrina_HIv1___RNAseq.g21001.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	MIRNA_hairpin	2872019	2872110	177	+	.	ID=Cluster_34;DicerCall=21;MIRNA=Y	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	2868586	2871318	.	+	.	ID=Pocillopora_meandrina_HIv1___TS.g25957.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	MIRNA_hairpin	2872019	2872110	177	+	.	ID=Cluster_34;DicerCall=21;MIRNA=Y	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	CDS	2870522	2871318	.	+	2	Parent=Pocillopora_meandrina_HIv1___TS.g25957.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	MIRNA_hairpin	2872019	2872110	177	+	.	ID=Cluster_34;DicerCall=21;MIRNA=Y	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	exon	2870522	2871318	.	+	2	Parent=Pocillopora_meandrina_HIv1___TS.g25957.t1
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	mature_miRNA	2872041	2872061	110	+	.	ID=Cluster_34.mature;Parent=Cluster_34	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	2868586	2871318	.	+	.	ID=Pocillopora_meandrina_HIv1___TS.g25957.t1
```

Do any features overlap?

```
bedtools intersect -a Ptuh_Results_sorted.gff3 -b /data/putnamlab/jillashey/genome/Pmea/Pocillopora_meandrina_HIv1.genes_sorted.gff3 > Ptuh_overlaps.bed

awk 'BEGIN {OFS="\t"} $3 != "Unknown_sRNA_locus"' Ptuh_overlaps.bed > filtered_Ptuh_overlaps.bed

wc -l filtered_Ptuh_overlaps.bed
174 filtered_Ptuh_overlaps.bed

head filtered_Ptuh_overlaps.bed
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	MIRNA_hairpin	818027	818120	12096	+	.	ID=Cluster_19;DicerCall=23;MIRNA=Y
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	mature_miRNA	818049	818070	3240	+	.	ID=Cluster_19.mature;Parent=Cluster_19
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	miRNA-star	818079	818100	9	+	.	ID=Cluster_19.star;Parent=Cluster_19
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	siRNA24_locus	5624826	5625249	125	+	.	ID=Cluster_73;DicerCall=24;MIRNA=N
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	siRNA24_locus	7205683	7206107	1021	+	.	ID=Cluster_100;DicerCall=24;MIRNA=N
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	siRNA22_locus	13015507	13015928	82	+	.	ID=Cluster_275;DicerCall=22;MIRNA=N
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	siRNA22_locus	16197870	16198683	17179	+	.	ID=Cluster_306;DicerCall=22;MIRNA=N
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	MIRNA_hairpin	20372416	20372510	522228	-	.	ID=Cluster_356;DicerCall=22;MIRNA=Y
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	miRNA-star	20372436	20372456	19166	-	.	ID=Cluster_356.star;Parent=Cluster_356
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	mature_miRNA	20372469	20372490	769	-	.	ID=Cluster_356.mature;Parent=Cluster_356
```

This doesn't tell me if the files are overlapping...it looks exactly like the results gff file. Does this indicate that most of the features have some kind of overlap??

I think if I subtract the 5th column (end of miRNA feature) by the 13th column (beginning of genomic feature), this will tell me if there is overlap in the two features. If the number is positive, it means there is overlap. If it is negative, it means there is no overlap. 

```
awk '{ $(NF+1) = $13 - $5 }1' OFS='\t' filtered_Ptuh_output.bed > filtered_Ptuh_output_overlap.bed

head filtered_Ptuh_output_overlap.bed
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	siRNA22_locus	173728	174150	1257	+	.	ID=Cluster_4;DicerCall=22;MIRNA=Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	CDS	174509	175333	.	-	0	Parent=Pocillopora_meandrina_HIv1___RNAseq.g20918.t1	359
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	siRNA22_locus	173728	174150	1257	+	.	ID=Cluster_4;DicerCall=22;MIRNA=Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	exon	174509	175333	.	-	0	Parent=Pocillopora_meandrina_HIv1___RNAseq.g20918.t1	359
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	siRNA22_locus	173728	174150	1257	+	.	ID=Cluster_4;DicerCall=22;MIRNA=Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	174509	176444	.	-	.	ID=Pocillopora_meandrina_HIv1___RNAseq.g20918.t1	359
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	MIRNA_hairpin	818027	818120	12096	+	.	ID=Cluster_19;DicerCall=23;MIRNA=Y	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	816355	820160	.	+	.	ID=Pocillopora_meandrina_HIv1___RNAseq.g21001.t1	-1765
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	mature_miRNA	818049	818070	3240	+	.	ID=Cluster_19.mature;Parent=Cluster_19	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	816355	820160	.	+	.	ID=Pocillopora_meandrina_HIv1___RNAseq.g21001.t1	-1715
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	miRNA-star	818079	818100	9	+	.	ID=Cluster_19.star;Parent=Cluster_19	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	816355	820160	.	+	.	ID=Pocillopora_meandrina_HIv1___RNAseq.g21001.t1	-1745
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	MIRNA_hairpin	2872019	2872110	177	+	.	ID=Cluster_34;DicerCall=21;MIRNA=Y	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	2868586	2871318	.	+	.	ID=Pocillopora_meandrina_HIv1___TS.g25957.t1	-3524
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	MIRNA_hairpin	2872019	2872110	177	+	.	ID=Cluster_34;DicerCall=21;MIRNA=Y	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	CDS	2870522	2871318	.	+	2	Parent=Pocillopora_meandrina_HIv1___TS.g25957.t1	-1588
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	MIRNA_hairpin	2872019	2872110	177	+	.	ID=Cluster_34;DicerCall=21;MIRNA=Y	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	exon	2870522	2871318	.	+	2	Parent=Pocillopora_meandrina_HIv1___TS.g25957.t1	-1588
Pocillopora_meandrina_HIv1___Sc0000000	ShortStack	mature_miRNA	2872041	2872061	110	+	.	ID=Cluster_34.mature;Parent=Cluster_34	Pocillopora_meandrina_HIv1___Sc0000000	AUGUSTUS	transcript	2868586	2871318	.	+	.	ID=Pocillopora_meandrina_HIv1___TS.g25957.t1	-3475
```

Count number of positive and negative numbers

```
awk '{if ($NF < 0) neg++; else if ($NF > 0) pos++} END {print "Positive count:", pos; print "Negative count:", neg}' filtered_Ptuh_output_overlap.bed

Positive count: 99
Negative count: 306
```

Negative = no overlap though I should check my math with someone to make sure I'm not dumb. 

Questions to answer now that I have this data: 

- Which miRNA features are overlapping with their closest genomic feature? 
- On average, how far away are miRNAs from their closest genomic feature? 