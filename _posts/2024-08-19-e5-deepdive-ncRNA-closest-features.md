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

Let's try to convert the gff to bed files. 

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

Unsure what this means...its giving me information from the miRNA stuff but no information about nearby features in the genome. Not sure why its not compatible. 


